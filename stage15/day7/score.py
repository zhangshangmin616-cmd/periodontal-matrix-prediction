"""Evaluate published human progenitor qPCR bar changes at the seventh day."""
import csv
import gzip
from itertools import product
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def write(name, rows):
    with (HERE / name).open("w", newline="", encoding="utf8") as file:
        table = csv.DictWriter(file, list(rows[0]), delimiter="\t", lineterminator="\n")
        table.writeheader()
        table.writerows(rows)


def read_measurements():
    with (HERE / "coordinates.tsv").open(encoding="utf8") as table:
        rows = list(csv.DictReader(table, delimiter="\t"))
    genes = [row["gene"] for row in rows]
    effects = []
    for row in rows:
        b, c, t, tol = (float(row[k]) for k in ("baseline_y", "control_top", "loaded_top", "tolerance"))
        obs = np.log2((b-t) / (b-c))
        varied = [np.log2((b+db-t-dt) / (b+db-c-dc)) for db, dc, dt in product((-tol, tol), repeat=3)]
        effects.append((obs, min(varied), max(varied)))
    return genes, np.array(effects).T


def side_factors(genes):
    symbols = genes + ["GAPDH"]
    with gzip.open(WORK / "stage10/state_multipliers.tsv.gz", "rt", encoding="utf8") as src:
        rows = {r["gene"]: r for r in csv.DictReader(src, delimiter="\t") if r["gene"] in symbols}
    sides = ("TP7D", "PP7D")
    tissue = np.exp2([[float(rows[g][s]) for g in symbols] for s in sides])
    with np.load(WORK / "stage10/reference_counts.npz") as atlas:
        f = list(atlas["classes"]).index("fibroblast_like")
        select = [list(atlas["genes"]).index(g) for g in symbols]
        def fibro_cpm(arm):
            return atlas[arm][select, f] * 1e6 / atlas[arm + "_total"][f]
        fibro = np.array([(1 + fibro_cpm(s)) / (1 + fibro_cpm("CP")) for s in sides])
    return {"fixed_shares": tissue, "fibroblast_like": fibro}


def reference_adjust(sides, weights):
    linear = (1-weights[:, None])*sides[0, :] + weights[:, None]*sides[1, :]
    return np.log2(linear[:, :-1] / linear[:, -1, None])


def bound_gains(curve, obs_low, obs_high):
    # A ratio of two positive affine factors is monotone or constant on [0, 1].
    # Each pair of adjacent pressure grid values thus encloses the full subinterval.
    p_low = np.minimum(curve[1:], curve[:-1])
    p_high = np.maximum(curve[1:], curve[:-1])
    zero_low = np.maximum.reduce([np.zeros_like(obs_low), obs_low, -obs_high])
    zero_high = np.maximum(np.abs(obs_low), np.abs(obs_high))
    largest_error = np.maximum(np.abs(obs_low - p_high), np.abs(obs_high - p_low))
    smallest_error = np.maximum(np.maximum(obs_low - p_high, p_low - obs_high), 0.)
    return np.min((zero_low - largest_error).mean(axis=1)), np.max((zero_high - smallest_error).mean(axis=1))


def main():
    genes, (obs, lo, hi) = read_measurements()
    with (WORK / "stage5/reference/features.tsv").open(encoding="utf8") as src:
        labels = {r["gene"]: r["ecm_component"] for r in csv.DictReader(src, delimiter="\t")}
    assert all(labels.get(g) != "other" for g in genes)
    models = side_factors(genes)
    predictions = []
    scores = []
    ranges = []
    weights = np.arange(1001) / 1000.
    zero_gene = np.abs(obs).mean()
    zero_panel = abs(obs.mean())
    for name, factors in models.items():
        curve = reference_adjust(factors, weights)
        lower, upper = bound_gains(curve, lo, hi)
        ranges.append(dict(model=name, pressure_interval="0 to 1", full_gene_gain_lower=lower,
                           full_gene_gain_upper=upper, note="conservative including graph-reading uncertainty"))
        for i in (0, 500, 1000):
            gene_pred = curve[i]
            for gene_index, gene in enumerate(genes):
                predictions.append(dict(gene=gene, model=name, pressure=weights[i],
                                        measured=obs[gene_index], measured_lo=lo[gene_index],
                                        measured_hi=hi[gene_index], predicted=gene_pred[gene_index],
                                        in_101=labels.get(gene) == "other"))
            for rule, pred in [("zero", np.zeros_like(gene_pred)), ("full_gene", gene_pred),
                               ("uniform_model_mean", np.full_like(gene_pred, gene_pred.mean()))]:
                e = pred - obs
                scores.append(dict(model=name, pressure=weights[i], rule=rule, n=len(genes),
                                   gene_mae=np.abs(e).mean(), gene_gain=zero_gene-np.abs(e).mean(),
                                   panel_mae=abs(e.mean()), panel_gain=zero_panel-abs(e.mean()),
                                   source_mean=pred.mean()))
    write("targets.tsv", predictions)
    write("scores.tsv", scores)
    write("ranges.tsv", ranges)


if __name__ == "__main__":
    main()
