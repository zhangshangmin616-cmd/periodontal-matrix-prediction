"""Transfer measured donor-paired tensile responses to independent human cultures."""
import csv
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def module(name, file):
    spec = spec_from_file_location(name, WORK / file)
    loaded = module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def write(name, records):
    with (HERE / name).open("w", encoding="utf8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def destination_inputs(shared_genes, target):
    # The archived count readers validate the official file checksums. Wuhan's
    # reader sums multiple SRR objects within a library before normalization.
    wuhan = module("wuhan_counts", "stage17/wuhan_counts/score.py")
    zhejiang = module("zhejiang_counts", "stage17/zhejiang/score.py")
    with np.load(WORK / "stage11/specificity/rank_changes.npz", allow_pickle=False) as frozen:
        all_genes = frozen["gene"]
        original_panel = all_genes[frozen["component"] == "other"]
    assert set(original_panel) == set(shared_genes[target]) and len(all_genes) == 10806
    lookup = {g:i for i,g in enumerate(all_genes)}
    indices = [lookup[g] for g in shared_genes]
    samples, times, cpm = wuhan.load_counts(all_genes)
    common = cpm[:,indices]
    result = []
    for time in (24, 48):
        mask = (times == 0) | (times == time)
        result.append((f"Wuhan_{time}h", samples[mask], times[mask] == time, common[mask]))
    result.append(("GSE106887_12h", np.array(["GSM2856840","GSM2856842"]), np.array([False,True]),
                   zhejiang.counts(shared_genes)))
    return result


def losses(p, y):
    residual = p-y
    sign_eligible = (np.abs(p) > 1e-10) & (np.abs(y) > 1e-10)
    concordant = ((p > 0) == (y > 0)) & sign_eligible
    return dict(observed_panel_mean=y.mean(), forecast_panel_mean=p.mean(),
                zero_panel_error=abs(y.mean()), panel_error=abs(residual.mean()),
                panel_gain=abs(y.mean())-abs(residual.mean()),
                zero_gene_mae=np.abs(y).mean(), gene_mae=np.abs(residual).mean(),
                gene_gain=np.abs(y).mean()-np.abs(residual).mean(),
                same_gene_direction=int(concordant.sum()), nonzero_directions=int(sign_eligible.sum()))


def main():
    previous = module("frozen_score", "stage13/external/score.py")
    with np.load(WORK / "stage13/external/paired_counts.npz", allow_pickle=False) as source:
        genes, target = source["genes"], source["target"]
        cpm = (source["counts"] * 1e6 / source["total"]).T
    assert cpm.shape == (10,10805) and target.sum() == 101
    rank_target = lambda x: previous.ranks(x)[target]
    hk_rank = np.stack([rank_target(lib) for lib in cpm])
    hk_deltas = hk_rank[5:] - hk_rank[:5]
    multipliers = np.exp2(previous.atlas(genes)["tissue"])
    all_scores, all_genes, deletions = [], [], []
    source_rows = [dict(gene=g, donor=d+1, observed_rank_change=v)
                   for d, delta in enumerate(hk_deltas) for g, v in zip(genes[target], delta, strict=True)]
    for label, samples, loaded, dest in destination_inputs(genes, target):
        actual = np.stack([rank_target(lib) for lib in dest])
        atlas = np.stack([rank_target(lib*multipliers)-base for lib,base in zip(dest,actual,strict=True)])
        units = [("none",None), *(("source_donor",d) for d in range(5))]
        if loaded.sum() > 1 and (~loaded).sum() > 1:
            units.extend([("destination_library", j) for j in range(len(dest))])
        for category,index in units:
            donors = np.ones(5,dtype=bool)
            libraries = np.ones(len(dest),dtype=bool)
            deleted_id = ""
            if category == "source_donor":
                donors[index] = False
                deleted_id = f"HK_donor_{index+1}"
            if category == "destination_library":
                libraries[index] = False
                deleted_id = samples[index]
            base, after = libraries & ~loaded, libraries & loaded
            y = actual[after].mean(axis=0)-actual[base].mean(axis=0)
            hk = hk_deltas[donors].mean(axis=0)
            predictions = {"zero":np.zeros_like(hk), "hong_kong_gene":hk,
                           "hong_kong_common":np.full_like(hk,hk.mean()),
                           "atlas_tissue":atlas[base].mean(axis=0)}
            for method,p in predictions.items():
                row = dict(destination=label, method=method, background=len(genes), targets=len(y),
                           source_donors=int(donors.sum()), baseline_libraries=int(base.sum()),
                           exposed_libraries=int(after.sum()), **losses(p,y))
                if category != "none":
                    deletions.append(dict(deletion=category, omitted=deleted_id, **row))
                else:
                    all_scores.append(row)
                    all_genes.extend(dict(destination=label, method=method, gene=g, observed=yi, forecast=pi)
                                     for g,yi,pi in zip(genes[target],y,p,strict=True))
    assert len(all_scores) == 12 and len(all_genes) == 1212 and len(deletions) == 108
    write("scores.tsv", all_scores)
    write("genes.tsv", all_genes)
    write("deletions.tsv", deletions)
    write("hong_kong_donor_changes.tsv", source_rows)


if __name__ == "__main__":
    main()
