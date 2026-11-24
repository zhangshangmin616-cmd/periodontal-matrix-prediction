"""Evaluate frozen atlas forecasts on five independent paired culture contrasts."""
import csv
import gzip
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def ranks(x):
    unique, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    midrank = np.cumsum(counts) - (counts + 1) / 2
    return 100 * midrank[inverse] / (x.size - 1)


def atlas(genes):
    with np.load(WORK / "stage10/reference_counts.npz", allow_pickle=False) as ref:
        ix = {g: i for i, g in enumerate(ref["genes"])}
        sel = np.array([ix[g] for g in genes])
        k = list(ref["classes"]).index("fibroblast_like")
        cp = ref["CP"][sel, k] * 1e6 / ref["CP_total"][k]
        tp = ref["TP3D"][sel, k] * 1e6 / ref["TP3D_total"][k]
    with gzip.open(WORK / "stage8/axial_log2_multipliers.tsv.gz", "rt", encoding="utf8") as stream:
        tissue = {r["gene"]: float(r["E3"]) for r in csv.DictReader(stream, delimiter="\t")}
    return {"zero": np.zeros(len(genes)), "tissue": np.array([tissue[g] for g in genes]),
            "fibroblast": np.log2((tp + 1) / (cp + 1))}


def save(name, rows):
    with (HERE / name).open("w", encoding="utf8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    with np.load(HERE / "paired_counts.npz", allow_pickle=False) as data:
        genes, target = data["genes"], data["target"]
        abundance = 1e6 * data["counts"] / data["total"]
    models = atlas(genes)
    rows = []
    effects = []
    for donor in range(5):
        control, loaded = abundance[:, donor], abundance[:, donor + 5]
        log_obs = np.log2(1 + loaded) - np.log2(1 + control)
        rank_base = ranks(control)
        rank_obs = ranks(loaded) - rank_base
        for method, shift in models.items():
            rank_est = ranks(control * np.exp2(shift)) - rank_base
            for unit, p, y in [("log2", shift, log_obs), ("rank", rank_est, rank_obs)]:
                for name, sel in [("target", target), ("background", np.ones(len(genes), dtype=bool))]:
                    err = (p-y)[sel]
                    rows.append(dict(donor=donor+1, scale=unit, scope=name, method=method, genes=int(sel.sum()),
                                gene_mae=np.mean(abs(err)), panel_abs_error=abs(err.mean()),
                                signed_error=err.mean(), observed_mean=y[sel].mean(), predicted_mean=p[sel].mean()))
            for i in np.flatnonzero(target):
                effects.append(dict(donor=donor+1, gene=genes[i], method=method, log2_observed=log_obs[i],
                                    log2_predicted=shift[i], rank_observed=rank_obs[i], rank_predicted=rank_est[i]))
    summary = []
    for unit in ["log2", "rank"]:
        for scope in ["target", "background"]:
            for method in models:
                sample = [r for r in rows if (r["scale"],r["scope"],r["method"]) == (unit,scope,method)]
                baseline = [r for r in rows if (r["scale"],r["scope"],r["method"]) == (unit,scope,"zero")]
                summary.append(dict(scale=unit, scope=scope, method=method, pairs=len(sample), genes=sample[0]["genes"],
                                    gene_mae=np.mean([r["gene_mae"] for r in sample]),
                                    panel_abs_error=np.mean([r["panel_abs_error"] for r in sample]),
                                    signed_error=np.mean([r["signed_error"] for r in sample]),
                                    donors_better_than_zero=sum(r["gene_mae"] < z["gene_mae"] for r,z in zip(sample, baseline))))
    save("hongkong_scores.tsv", summary)
    save("hongkong_donor_scores.tsv", rows)
    save("hongkong_target_effects.tsv", effects)


if __name__ == "__main__":
    main()
