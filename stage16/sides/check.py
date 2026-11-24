"""Reconstruct clinical ranks and gain bounds independently of the two scoring scripts."""
import csv
import gzip
from pathlib import Path
import numpy as np


HERE = Path(__file__).resolve().parent
STAGES = HERE.parents[1]


def table(path):
    with path.open(encoding="utf8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def percentiles(values):
    """Sort once, then assign the midpoint of each equal-value block."""
    values = np.asarray(values)
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    starts = np.r_[0, 1 + np.flatnonzero(sorted_values[:-1] != sorted_values[1:])]
    ends = np.r_[starts[1:], values.size]
    mid = (starts + ends - 1) * (50.0 / (values.size - 1))
    result = np.empty(values.size)
    result[order] = np.repeat(mid, ends-starts)
    return result


def close(a, b):
    np.testing.assert_allclose(a, b, rtol=0, atol=3e-12)


def from_counts(genes, target, controls):
    zhao = STAGES/"stage5/zhao"
    ids = {r["gene"]:r["entrez_id"] for r in table(zhao/"feature_map.tsv")}
    mapped = [ids[g] for g in genes]
    assert len(set(mapped)) == 10806
    groups = {"control": [], "day7": []}
    for sample in table(zhao/"joined_manifest.tsv"):
        arm, run = sample["group"], sample["run"]
        if sample["batch"] != "4" or arm not in groups:
            continue
        with gzip.open(zhao/"counts"/(run + ".tsv.gz"), "rt", encoding="utf8") as stream:
            r = csv.reader(stream, delimiter="\t")
            assert next(r) == ["Geneid", run+"_count"]
            counts = {key:float(value) for key,value in r}
        ranks = percentiles([counts[key] for key in mapped])
        if arm == "control":
            close(ranks, percentiles(controls[run]))
        groups[arm].append(ranks[target])
    assert len(groups["day7"]) == 19 and len(groups["control"]) == 5
    return np.mean(groups["day7"],axis=0) - np.mean(groups["control"],axis=0)


def clinical_forecasts():
    with np.load(STAGES/"stage8/mixture_inputs.npz",allow_pickle=False) as v:
        genes, target = v["symbols"], v["idx"]
        b, pressure, tension = v["control"], v["pressure"], v["tension"]
        controls = dict(zip(v["samples"], b, strict=True))
    assert len(set(genes[target])) == 101
    observed = from_counts(genes, target, controls)
    with np.load(STAGES/"stage11/specificity/rank_changes.npz",allow_pickle=False) as v:
        source_by_gene = {g:i for i,g in enumerate(v["gene"])}
        indices = [source_by_gene[g] for g in genes[target]]
        assert np.all(v["component"][indices]=="other")
        unit = np.flatnonzero(v["units"]=="Zhao_day7").item()
        close(observed, v["observed"][unit, indices])
        previous = v["fitted"][unit, indices]
    mean_target_rank = lambda rows:np.mean([percentiles(row)[target] for row in rows],axis=0)
    baseline = mean_target_rank(b)
    pressure_effect = mean_target_rank(pressure) - baseline
    tension_effect = mean_target_rank(tension) - baseline
    pooled_effect = mean_target_rank(0.5*pressure + 0.5*tension) - baseline
    close(pooled_effect,previous)
    models = {"original":(pressure_effect,tension_effect),
              "side_labels_exchanged":(tension_effect,pressure_effect),
              "same_on_both_sides":(pooled_effect,pooled_effect)}
    by_gene = {(r["assignment"],r["gene"]):r for r in table(HERE/"forecasts.tsv")}
    metrics = {r["assignment"]:r for r in table(HERE/"metrics.tsv")}
    assert len(by_gene) == 303 and set(models) == set(metrics)
    for label,(left,right) in models.items():
        diff = left-right
        for i,g in enumerate(genes[target]):
            row = by_gene[label,g]
            expected = (left[i],right[i],diff[i],pooled_effect[i],observed[i])
            stored = [float(row[k]) for k in ("predicted_pressure","predicted_tension",
                "predicted_pressure_minus_tension","predicted_pooled","observed_unsplit_change")]
            close(stored,expected)
        m = metrics[label]
        assert int(m["targets"])==len(target) and int(m["targets_with_side_contrast"])==np.sum(abs(diff)>1e-10)
        delta = pooled_effect-observed
        expected = (observed.mean(), pooled_effect.mean(), abs(delta.mean()),abs(delta).mean(),
                    abs(observed.mean()),abs(observed).mean(),abs(diff).mean())
        close([float(m[k]) for k in ("observed_panel_change","predicted_panel_change","panel_error",
              "gene_mae","zero_panel_error","zero_gene_mae","predicted_side_mean_absolute_contrast")], expected)
    return genes[target]


def log_effects(panel):
    with np.load(STAGES/"stage10/reference_counts.npz",allow_pickle=False) as v:
        genes = {g:i for i,g in enumerate(v["genes"])}
        take = [genes[g] for g in panel]
        total = v["CP_total"]
        weights = total / total.sum()
        cp = 1e6 * v["CP"][take] / total
        residual = np.flatnonzero(v["classes"]=="residual").item()
        ref = np.flatnonzero(panel=="ACTB").item()
        effects = {}
        for lib in ("PP7D","TP7D"):
            post = 1e6*v[lib][take] / v[lib+"_total"]
            post[:, residual] = cp[:, residual]
            fold = (1+(post*weights).sum(axis=1)) / (1+(cp*weights).sum(axis=1))
            effects[lib] = np.log2(fold/fold[ref])
        return effects


def gain(x,y):
    return abs(y)-abs(x-y)


def check_bounds(panel):
    atlas = log_effects(panel)
    contrasts = {"pressure_vs_untreated":atlas["PP7D"],"tension_vs_untreated":atlas["TP7D"],
                 "pressure_vs_tension":atlas["PP7D"]-atlas["TP7D"]}
    tg = np.flatnonzero(panel=="TGFB1").item()
    hk = np.flatnonzero(panel=="ACTB").item()
    result = table(HERE/"bounds.tsv")
    assert len(result)==len(contrasts)==3
    for row in result:
        name = row["contrast"]
        p = contrasts[name]
        assert p[hk]==0 and int(row["targets"])==len(panel) and row["target_predictions_relative_to"]=="ACTB"
        assert int(row["biological_targets_with_reported_assay"]) == 1
        allowed = (name!="pressure_vs_tension")
        limit = 2+max(abs(p))
        lo, hi = [], []
        # Absolute loss is linear between its knots in y and constant in either tail.
        # Evaluating the knots and one point in each allowed tail gives sharp closure bounds.
        for i,x in enumerate(p):
            knots = [-limit,0.,x,limit]
            if i==hk:
                knots=[0.]
            elif i==tg and allowed:
                knots=[y for y in knots if y>=0.]
            values = [gain(x,y) for y in knots]
            lo.append(min(values))
            hi.append(max(values))
        close([float(row[key]) for key in ("TGFB1_prediction","TGFB1_min_gain","TGFB1_max_gain",
              "full101_min_gain","full101_max_gain")],[p[tg],lo[tg],hi[tg],np.mean(lo),np.mean(hi)])
        if allowed and p[tg]>0:
            # With a strict positive observed response, the lower endpoint is a limit.
            eps = 1e-8
            close(gain(p[tg],eps) - lo[tg],2*eps)
        if not allowed:
            for delta in (-limit,0.,limit):
                left,right = 1+max(0.,delta),1+max(0.,-delta)
                assert left>0 and right>0
                close(left-right, delta)


def main():
    close(percentiles([2,0,1,1,0]), [100,12.5,62.5,62.5,12.5])
    panel=clinical_forecasts()
    check_bounds(panel)
    print("Verified 24 observed clinical libraries, 303 gene-by-assignment forecasts, three aggregate scores, and three gain intervals from the original atlas counts.")


if __name__=="__main__":
    main()
