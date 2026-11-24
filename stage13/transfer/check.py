"""Check cross-study forecasts, loss aggregation and all recorded single deletions."""
import csv
from collections import defaultdict
from math import fsum, isclose
from pathlib import Path
from statistics import mean

import numpy as np


HERE = Path(__file__).resolve().parent


def rows(name):
    with (HERE / name).open(encoding="utf8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def near(x, y):
    assert isclose(float(x), float(y), abs_tol=2e-10, rel_tol=1e-12), (x,y)


def main():
    work = HERE.parents[1]
    with np.load(work / "stage11/specificity/rank_changes.npz") as d:
        col_genes = d["gene"].tolist()
        col = {g:k for k,g in enumerate(col_genes)}
        units = d["units"].tolist()
        obs = dict(zip(units, d["observed"].tolist()))
        atlas = dict(zip(units, d["fitted"].tolist()))
        calibration = np.flatnonzero(d["component"]=="outside_ecm").tolist()
    labels = {"Seoul_week3":[u for u in units if u.startswith("Kim_")],
              "Leuven_day7":["Zhao_day7"]}
    lookup = {(r["source"],r["destination"]):r for r in rows("coefficients.tsv")}
    for (source,dest),record in lookup.items():
        source_units = labels[source]
        numerator = fsum(atlas[u][j]*obs[u][j] for j in calibration for u in source_units)
        denominator = fsum(atlas[u][j]**2 for j in calibration for u in source_units)
        near(record["unconstrained_slope"],numerator/denominator)
        near(record["bounded_slope"],max(0,min(1,numerator/denominator)))

    archive = rows("forecasts.tsv")
    all_models = list(archive[0])[6:]
    assert len(all_models) == 7
    by_unit = defaultdict(list)
    for r in archive:
        s, dest, u, gene = r["source"],r["destination"],r["unit"],r["gene"]
        j = col[gene]
        empirical = mean(obs[t][j] for t in labels[s])
        near(empirical,r["clinical_full"])
        near(.25*empirical,r["clinical_quarter"])
        near(obs[u][j],r["observed_change"])
        near(atlas[u][j],r["atlas_full"])
        near(.25*atlas[u][j],r["atlas_quarter"])
        near(float(lookup[s,dest]["bounded_slope"])*atlas[u][j],r["atlas_source_calibrated"])
        near(0,r["zero"])
        by_unit[s,dest,u].append(r)
    unit_scores = defaultdict(list)
    for (source,destination,unit), group in by_unit.items():
        assert len(group) == 101
        uniform = mean(float(r["clinical_full"]) for r in group)
        for r in group:
            near(uniform, r["clinical_uniform"])
        for model in all_models:
            for bounds in ["original","bounded"]:
                errors = []
                for r in group:
                    p, b, y = map(float,(r[model],r["untreated_rank"],r["observed_change"]))
                    if bounds == "bounded":
                        p = max(0.,min(100.,b+p))-b
                    errors.append(p-y)
                unit_scores[source,destination,model,bounds].append(
                    (mean(map(abs,errors)),abs(mean(errors))))
    scores = rows("scores.tsv")
    assert len(scores) == 42
    for r in scores:
        key = tuple(r[k] for k in ["source","destination","model","bounds"])
        values = unit_scores[key]
        near(mean(v[0] for v in values),r["gene_mae"])
        near(mean(v[1] for v in values),r["panel_mae"])

    group_metrics = rows("group_metrics.tsv")
    assert len(group_metrics) == 63
    for r in group_metrics:
        key = r["source"],r["destination"]
        matching = [g for k,g in by_unit.items() if k[:2] == key]
        b = np.array([[float(q["untreated_rank"]) for q in g] for g in matching])
        y = np.array([[float(q["observed_change"]) for q in g] for g in matching])
        p = np.array([[float(q[r["model"]]) for q in g] for g in matching])
        if r["bounds"] == "bound_unit_then_mean":
            estimate = np.clip(p,-b,100.-b).mean(axis=0)
        elif r["bounds"] == "mean_then_bound_group":
            estimate = np.clip(p.mean(axis=0),-b.mean(axis=0),100.-b.mean(axis=0))
        else:
            estimate = p.mean(axis=0)
        e = estimate - y.mean(axis=0)
        near(np.abs(e).mean(),r["group_gene_mae"])
        near(abs(e.mean()),r["group_panel_abs_error"])

    deletions = defaultdict(list)
    for r in rows("deletions_by_unit.tsv"):
        key = tuple(r[k] for k in ["source","destination","model","exclusion_type"])
        deletions[key].append(float(r["gene_gain_vs_zero"]))
    summary = rows("deletions_summary.tsv")
    assert len(summary) == 18
    for r in summary:
        key = tuple(r[k] for k in ["source","destination","model","exclusion_type"])
        v = deletions[key]
        near(len(v),r["deletions"])
        near(sum(gain>0 for gain in v),r["improved"])
        near(min(v),r["worst_gain"])
        near(max(v),r["best_gain"])
    print(f"Recalculated {len(archive)} forecast rows, {len(scores)} unit-first and {len(group_metrics)} group-first scores.")
    print(f"Verified {sum(map(len,deletions.values()))} single-deletion gain entries against the summaries.")


if __name__ == "__main__":
    main()
