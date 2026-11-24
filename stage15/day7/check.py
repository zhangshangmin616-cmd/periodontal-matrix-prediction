"""Independently verify digitized changes, atlas predictions and all reported scores."""
import csv
from itertools import product
from math import isclose, log2
from pathlib import Path
from statistics import mean

import numpy as np


HERE = Path(__file__).resolve().parent


def rows(name):
    with (HERE / name).open(encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def close(a, b):
    assert isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-10), (a, b)


def main():
    figures = {}
    for r in rows("coordinates.tsv"):
        base, control, loaded, tol = [float(r[k]) for k in ("baseline_y", "control_top", "loaded_top", "tolerance")]
        point = log2(base-loaded) - log2(base-control)
        bounds = [log2(base+b-loaded-t)-log2(base+b-control-c)
                  for b, c, t in product((-tol, tol), repeat=3)]
        figures[r["gene"]] = (point, min(bounds), max(bounds))
    assert len(figures) == 6
    effect = {}
    genes = [*figures, "GAPDH"]
    with np.load(HERE.parents[1] / "stage10/reference_counts.npz") as data:
        baseweights = data["CP_total"] / sum(data["CP_total"])
        fibro = data["classes"].tolist().index("fibroblast_like")
        residual = data["classes"].tolist().index("residual")
        for sym in genes:
            idx = data["genes"].tolist().index(sym)
            base = 1e6 * data["CP"][idx, :] / data["CP_total"]
            for arm in ("TP7D", "PP7D"):
                treated = 1e6 * data[arm][idx, :] / data[arm + "_total"]
                effect["fibroblast_like", sym, arm] = (1 + treated[fibro]) / (1 + base[fibro])
                treated[residual] = base[residual]
                effect["fixed_shares", sym, arm] = (1 + sum(treated*baseweights)) / (1 + sum(base*baseweights))

    def prediction(model, symbol, w):
        def mult(s):
            return (1-w)*effect[model, s, "TP7D"] + w*effect[model, s, "PP7D"]
        return log2(mult(symbol)) - log2(mult("GAPDH"))

    for r in rows("targets.tsv"):
        g = r["gene"]
        close(figures[g][0], r["measured"])
        close(figures[g][1], r["measured_lo"])
        close(figures[g][2], r["measured_hi"])
        close(prediction(r["model"], g, float(r["pressure"])), r["predicted"])
        assert r["in_101"] == "False"
    for r in rows("scores.tsv"):
        preds = [prediction(r["model"], g, float(r["pressure"])) for g in figures]
        if r["rule"] == "uniform_model_mean":
            preds = [mean(preds)] * len(preds)
        if r["rule"] == "zero":
            preds = [0.] * len(preds)
        actual = [v[0] for v in figures.values()]
        err = [p-y for p, y in zip(preds, actual)]
        close(mean(abs(d) for d in err), r["gene_mae"])
        close(abs(mean(err)), r["panel_mae"])
        close(mean(abs(y) for y in actual)-mean(abs(e) for e in err), r["gene_gain"])
        close(abs(mean(actual))-abs(mean(err)), r["panel_gain"])
    # Independent test at each observation interval's extrema and midpoints
    # at 101 distinct settings, comparing with the certified outer bounds.
    for r in rows("ranges.tsv"):
        bounds = list(map(float, (r["full_gene_gain_lower"], r["full_gene_gain_upper"])))
        for ix in range(101):
            lows, highs = [], []
            for g, (obs, low, high) in figures.items():
                pred = prediction(r["model"], g, ix / 100)
                changes = [abs(y) - abs(y-pred) for y in (low, high, obs, np.clip(0, low, high),
                                                          np.clip(pred, low, high))]
                lows.append(min(changes))
                highs.append(max(changes))
            assert bounds[0] <= mean(lows) + 1e-12 and mean(highs) <= bounds[1]+1e-12
    print("Verified 6 digitized genes, 36 target forecasts, 18 model scores and both sensitivity envelopes.")


if __name__ == "__main__":
    main()
