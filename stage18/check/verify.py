"""Rebuild the three rank-transfer tests from the original count tables."""
import csv
import gzip
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path
import re

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
HK = WORK / "stage13/external"


def read(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(name, rows):
    with (HERE / name).open("w", encoding="utf8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def percentiles(vector):
    ordered = np.sort(vector)
    left = np.searchsorted(ordered, vector, side="left")
    right = np.searchsorted(ordered, vector, side="right")
    return (left + right - 1) * (50.0 / (len(vector) - 1))


def hong_kong():
    by_id = defaultdict(set)
    for row in read(WORK / "stage5/paired/hgnc_ids.tsv.gz"):
        by_id[row["ensembl_gene_id"]].add(row["symbol"])
    fixed = read(WORK / "stage5/reference/features.tsv")
    by_symbol = defaultdict(list)
    source = HK / "raw/GSE173891_raw_gene_count.txt.gz"
    expect_hash = json.loads((HK / "source_manifest.json").read_bytes())[source.name]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == expect_hash
    metadata = sorted(read(HK / "sample_map.tsv"), key=lambda x:x["column"])
    wanted = [row["column"] for row in metadata]
    totals = np.zeros(10, dtype=np.int64)
    with gzip.open(source, "rt", encoding="utf8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        assert next(reader) == wanted
        for ens, *c in reader:
            vals = np.array(c, dtype=np.int64)
            totals += vals
            names = by_id.get(ens.partition(".")[0], ())
            if len(names) == 1:
                by_symbol[next(iter(names))].append(vals)
    np.testing.assert_equal(totals, [int(r["all_feature_total"]) for r in metadata])
    rows = [r for r in fixed if len(by_symbol[r["gene"]]) == 1]
    genes = [r["gene"] for r in rows]
    target = np.array([r["ecm_component"] == "other" for r in rows])
    assert len(genes) == 10805 and target.sum() == 101
    raw = np.stack([by_symbol[g][0] for g in genes])
    with np.load(HK / "paired_counts.npz", allow_pickle=False) as f:
        np.testing.assert_equal(genes, f["genes"])
        np.testing.assert_equal(raw, f["counts"])
        np.testing.assert_equal(target, f["target"])
    shifts = []
    for donor in sorted({r["donor"] for r in metadata}):
        ix = {r["arm"]:i for i,r in enumerate(metadata) if r["donor"] == donor}
        shifts.append((percentiles(raw[:,ix["tension"]]) - percentiles(raw[:,ix["control"]]))[target])
    return genes,target,np.stack(shifts)


def atlas(symbols):
    with np.load(WORK / "stage10/reference_counts.npz", allow_pickle=False) as f:
        lookup = {g:i for i,g in enumerate(f["genes"])}
        idx = [lookup[g] for g in symbols]
        baseline = f["CP"][idx] * 1e6 / f["CP_total"]
        loaded = f["TP3D"][idx] * 1e6 / f["TP3D_total"]
        residual = list(f["classes"]).index("residual")
        loaded[:,residual] = baseline[:,residual]
        weights = f["CP_total"] / sum(f["CP_total"])
        return (1+(loaded*weights).sum(axis=1)) / (1+(baseline*weights).sum(axis=1))


def load_profile(ids, rundir, runs, checksums):
    counts = np.zeros(len(ids))
    total = 0.
    for run in sorted(runs):
        raw = gzip.decompress((rundir / (run+".tsv.gz")).read_bytes())
        assert hashlib.md5(raw).hexdigest() == checksums[run]
        lines = csv.reader(raw.decode().splitlines(), delimiter="\t")
        assert next(lines) == ["Geneid",run+"_count"]
        data = {k:float(v) for k,v in lines}
        counts += np.array([data[i] for i in ids])
        total += sum(data.values())
    return 1e6 * counts / total


def destinations(genes, target):
    ids_by_gene = {r["gene"]:r["entrez_id"] for r in read(WORK / "stage5/zhao/feature_map.tsv")}
    ids = [ids_by_gene[g] for g in genes]
    factors = atlas(genes)
    wh = WORK / "stage17/wuhan_counts"
    zj = WORK / "stage17/zhejiang"
    levels = {}
    for location in ("Wuhan","Zhejiang"):
        groups = defaultdict(list)
        hours = {}
        if location == "Wuhan":
            directory = wh
            ena = json.loads((WORK/"stage17/culture/ena_runs.json").read_bytes())
            for r in ena:
                m = re.search(r"(GSM\d+): sh-Control, tension (\d+) h, rep", r["experiment_title"])
                if m:
                    groups[m[1]].append(r["run_accession"])
                    hours[m[1]] = int(m[2])
            assert len(groups) == 9 and sum(map(len,groups.values())) == 12
        else:
            directory = zj
            for r in read(zj / "samples.tsv"):
                label = {"normal PDLSCs":0,"stretch-induced PDLSCs":12}.get(r["sample_title"])
                if label is not None:
                    groups[r["sample_alias"]].append(r["run_accession"])
                    hours[r["sample_alias"]] = label
            assert len(groups) == 2
        sums = {r["run"]:r["md5"] for r in read(directory/"count_manifest.tsv")}
        base, pred = [], []
        observed = defaultdict(list)
        for sample in sorted(groups):
            x = load_profile(ids, directory/"counts",groups[sample],sums)
            r = percentiles(x)[target]
            if hours[sample] == 0:
                base.append(r)
                pred.append(percentiles(x*factors)[target]-r)
            else:
                observed[hours[sample]].append(r)
        for h,values in sorted(observed.items()):
            levels[location+str(h)] = np.mean(values,axis=0)-np.mean(base,axis=0), np.mean(pred,axis=0)
    return levels


def loss(y,p):
    delta = y-p
    return float(abs(np.mean(delta))), float(np.mean(abs(delta)))


def main():
    genes, target, paired = hong_kong()
    levels = destinations(genes,target)
    result, orient, summary = [], [], []
    expected = HERE.parent / "transfer"
    destination_key = {"Wuhan24":"Wuhan_24h","Wuhan48":"Wuhan_48h","Zhejiang12":"GSE106887_12h"}
    model_key = {"zero":"zero","Hong_Kong":"hong_kong_gene","Hong_Kong_common":"hong_kong_common","atlas":"atlas_tissue"}
    gene_rows = read(expected / "genes.tsv")
    assert len(gene_rows) == 3*4*101
    published = {(r["destination"],r["method"],r["gene"]):[float(r["observed"]),float(r["forecast"])]
                 for r in gene_rows}
    source_rows = read(expected / "hong_kong_donor_changes.tsv")
    source_pairs = {(r["gene"],int(r["donor"])):float(r["observed_rank_change"]) for r in source_rows}
    symbols = np.array(genes)[target]
    assert len(source_pairs) == 505
    for i,pair in enumerate(paired,start=1):
        np.testing.assert_allclose(pair,[source_pairs[g,i] for g in symbols],rtol=0,atol=2e-12)
    for name,(y,atlas_pred) in levels.items():
        hk = paired.mean(axis=0)
        metrics = {}
        for method,p in [("zero",np.zeros_like(hk)),("Hong_Kong",hk),
                         ("Hong_Kong_common",np.full_like(hk,hk.mean())),("atlas",atlas_pred)]:
            pe,ge = loss(y,p)
            result.append(dict(destination=name,model=method,panel_error=pe,gene_mae=ge,
                               observed_mean=y.mean(),predicted_mean=p.mean()))
            metrics[method]=(pe,ge)
            values = [published[destination_key[name],model_key[method],g] for g in symbols]
            np.testing.assert_allclose(np.array(values).T,[y,p],rtol=0,atol=2e-12)
        errors = []
        scaled_errors = []
        spread = np.abs(hk-hk.mean()).mean()
        for signs in itertools.product((-1,1),repeat=len(paired)):
            p = np.mean(np.array(signs)[:,None]*paired,axis=0)
            a,b = loss(y,p)
            centred = p-p.mean()
            width = abs(centred).mean()
            assert width > 0
            matched = hk.mean()+centred*(spread/width)
            ma,mb = loss(y,matched)
            np.testing.assert_allclose([ma,abs(matched-matched.mean()).mean()],
                                       [metrics["Hong_Kong"][0],spread],rtol=0,atol=2e-12)
            errors.append([a,b])
            scaled_errors.append(mb)
            orient.append(dict(destination=name,donors_1_to_5="".join("+" if s>0 else "-" for s in signs),
                               predicted_mean=p.mean(), gene_spread=width, panel_error=a, gene_mae=b,
                               matched_mean_and_spread_gene_mae=mb))
        errors = np.array(errors)
        actual = np.array(metrics["Hong_Kong"])
        np.testing.assert_allclose(actual,errors[-1],atol=1e-12,rtol=0)
        z = np.array(metrics["zero"])
        count=(errors <= actual+1e-12).sum(axis=0)
        summary.append(dict(destination=name, total_assignments=len(errors),
                            panel_no_worse_than_real=int(count[0]), gene_no_worse_than_real=int(count[1]),
                            assignments_with_panel_gain=int(sum(errors[:,0]<z[0]-1e-12)),
                            assignments_with_gene_gain=int(sum(errors[:,1]<z[1]-1e-12)),
                            matched_mean_spread_no_worse_than_real=int(sum(v<=actual[1]+1e-12 for v in scaled_errors)),
                            matched_mean_spread_beating_zero=int(sum(v<z[1]-1e-12 for v in scaled_errors)),
                            reversed_panel_error=errors[0,0],reversed_gene_mae=errors[0,1],
                            correct_panel_error=actual[0],correct_gene_mae=actual[1]))
    save("scores_independent.tsv",result)
    save("source_orientations.tsv",orient)
    save("orientation_summary.tsv",summary)
    print("Verified the original Hong Kong table, 14 destination run checksums and the integer atlas. Scored",len(orient),"orientation/destination combinations.")


if __name__ == "__main__":
    main()
