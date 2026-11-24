"""Score the frozen tissue and fibroblast models in the Wuhan tensile series."""
import csv
import gzip
import hashlib
import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def table(path):
    with path.open(encoding="utf8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(path, rows):
    with path.open("w", encoding="utf8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_counts(genes):
    mappings = {r["gene"]:r["entrez_id"] for r in table(WORK/"stage5/zhao/feature_map.tsv")}
    identifiers = [mappings[g] for g in genes]
    assert len(set(identifiers)) == len(genes) == 10806
    manifest = {r["run"]:r["md5"] for r in table(HERE/"count_manifest.tsv")}
    runs = table(HERE/"libraries.tsv")
    assert set(manifest) == {r["run"] for r in runs}
    columns, totals, metadata = {}, {}, {}
    for row in runs:
        run, sample, time = row["run"], row["gsm"], int(row["hours"])
        payload = gzip.decompress((HERE/"counts"/(run+".tsv.gz")).read_bytes())
        assert hashlib.md5(payload).hexdigest() == manifest[run]
        it = iter(payload.decode().splitlines())
        assert next(it).split("\t") == ["Geneid",run+"_count"]
        data = {a:float(b) for a,b in (line.split("\t") for line in it)}
        vector = np.array([data[e] for e in identifiers])
        columns[sample] = columns.get(sample,0) + vector
        totals[sample] = totals.get(sample,0) + sum(data.values())
        if sample in metadata:
            assert metadata[sample] == time
        metadata[sample] = time
    samples = sorted(columns)
    times = np.array([metadata[s] for s in samples])
    assert list(times) == [0]*3+[24]*3+[48]*3
    return np.array(samples), times, np.stack([columns[s]*1e6/totals[s] for s in samples])


def losses(p,y):
    diff = p-y
    return dict(observed_mean=y.mean(), predicted_mean=p.mean(),
                panel_abs_error=abs(diff.mean()), zero_panel_abs_error=abs(y.mean()),
                panel_gain=abs(y.mean())-abs(diff.mean()), gene_mae=np.abs(diff).mean(),
                zero_gene_mae=np.abs(y).mean(), gene_gain=np.abs(y).mean()-np.abs(diff).mean())


def main():
    with np.load(WORK/"stage11/specificity/rank_changes.npz",allow_pickle=False) as data:
        genes, target = data["gene"], data["component"] == "other"
    assert target.sum() == 101
    spec = importlib.util.spec_from_file_location("frozen", WORK/"stage13/external/score.py")
    frozen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(frozen)
    samples, times, cpm = load_counts(genes)
    ranks = np.stack([frozen.ranks(a)[target] for a in cpm])
    obs = {"rank":ranks, "log2":np.log2(1+cpm[:,target])}
    summary, deletions, targets = [], [], []
    for model, coeff in frozen.atlas(genes).items():
        if model == "zero":
            continue
        rank_preds = np.stack([frozen.ranks(b*np.exp2(coeff))[target] - ranks[j] for j,b in enumerate(cpm[:3])])
        for time in (24,48):
            selected = np.flatnonzero((times == 0) | (times == time))
            for omit in [None, *selected]:
                keep = selected if omit is None else selected[selected != omit]
                baseline, exposed = keep[times[keep] == 0], keep[times[keep] == time]
                for scale in ("rank","log2"):
                    y = obs[scale][exposed].mean(axis=0) - obs[scale][baseline].mean(axis=0)
                    p = rank_preds[baseline].mean(axis=0) if scale == "rank" else coeff[target]
                    forecasts = {model:p}
                    if model == "tissue":
                        forecasts["tissue_common_shift"] = np.full_like(p,p.mean())
                    for method,prediction in forecasts.items():
                        row = dict(hours=time, model=method, scale=scale, genes=len(y), background=len(genes),
                                   **losses(prediction,y))
                        if omit is None:
                            summary.append(row)
                            targets.extend(dict(hours=time, model=method, scale=scale, gene=g,
                                                observed=a, predicted=b) for g,a,b in zip(genes[target],y,prediction))
                        else:
                            deletions.append(dict(omitted=samples[omit], **row))
    save(HERE/"scores.tsv", summary)
    save(HERE/"deletions.tsv", deletions)
    save(HERE/"targets.tsv", targets)
    for r in summary:
        if r["model"] == "tissue" and r["scale"] == "rank":
            print(r)


if __name__ == "__main__":
    main()
