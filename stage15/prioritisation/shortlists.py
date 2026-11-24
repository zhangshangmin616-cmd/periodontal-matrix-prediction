"""Retrieve measured large responses with a fixed-size signed gene shortlist."""
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def rows(path):
    with path.open(encoding="utf8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, data):
    with path.open("w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=data[0], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(data)


def column(records, name):
    return np.array([float(r[name]) for r in records])


def top_weights(x, k):
    a = np.abs(x)
    cut = np.sort(a)[-k]
    greater, tied = a > cut, a == cut
    return greater.astype(float) + tied * ((k - greater.sum()) / tied.sum())


def baseline_ranks(raw):
    ordered = np.sort(raw)
    ranks = (ordered.searchsorted(raw) + ordered.searchsorted(raw, side="right") - 1) / 2
    return ranks * 100 / (len(raw) - 1)


def clinical():
    blocks = defaultdict(list)
    for row in rows(WORK / "stage13/transfer/forecasts.tsv"):
        blocks[row["destination"], row["unit"]].append(row)
    seoul = []
    for (dst, unit), values in sorted(blocks.items()):
        values = sorted(values, key=lambda x: x["gene"])
        g = np.array([v["gene"] for v in values])
        b, y, a, m = [column(values, key) for key in ("untreated_rank", "observed_change", "atlas_full", "clinical_full")]
        for model, p in [("atlas_tissue", a), ("other_clinical", m)]:
            yield dst, unit, model, g, b, y, p
        if dst == "Seoul_week3":
            seoul.append((g, b, y, a, m))
    assert len(seoul) == 11
    g = seoul[0][0]
    assert all(np.array_equal(r[0], g) for r in seoul)
    b, y, a, m = [np.mean([r[i] for r in seoul], axis=0) for i in range(1, 5)]
    yield "Seoul_week3", "cohort_mean", "atlas_tissue", g, b, y, a
    yield "Seoul_week3", "cohort_mean", "other_clinical", g, b, y, m


def cultures():
    by_hk = defaultdict(list)
    for row in rows(WORK / "stage13/external/hongkong_target_effects.tsv"):
        by_hk[int(row["donor"]), row["method"]].append(row)
    by_xa = defaultdict(list)
    for row in rows(WORK / "stage14/atlas_xian/targets.tsv"):
        if row["scale"] == "rank":
            by_xa[row["frequency_Hz"], row["method"]].append(row)
    with np.load(WORK / "stage13/external/paired_counts.npz") as raw:
        hk = raw["counts"]
        target = raw["target"]
        genes = raw["genes"][target]
    with np.load(WORK / "stage14/xian/xian_counts.npz") as raw:
        xa = raw["counts"][:, raw["arm"] == "0"]
        assert np.array_equal(genes, raw["genes"][raw["target"]])
    ix = np.argsort(genes)
    g = genes[ix]
    untreated_hk = [baseline_ranks(arr)[target][ix] for arr in hk[:, :5].T]
    untreated_xa = np.mean([baseline_ranks(arr)[target][ix] for arr in xa.T], axis=0)
    hk_observations = []
    for donor in range(1, 6):
        for model in ("tissue", "fibroblast"):
            r = sorted(by_hk[donor,model], key=lambda z:z["gene"])
            assert np.array_equal(g, [z["gene"] for z in r])
            y, p = column(r, "rank_observed"), column(r, "rank_predicted")
            yield "HongKong", str(donor), "atlas_"+model, g, untreated_hk[donor-1], y, p
        hk_observations.append(y)
    for i, y in enumerate(hk_observations):
        p = np.mean(np.delete(hk_observations, i, axis=0), axis=0)
        yield "HongKong", str(i+1), "other_four_donors", g, untreated_hk[i], y, p
    for freq in ("0.1", "0.5"):
        for model in ("tissue", "fibroblast"):
            r = sorted(by_xa[freq,model], key=lambda z:z["gene"])
            assert np.array_equal(g, [z["gene"] for z in r])
            y, p = column(r, "observed"), column(r, "prediction")
            yield "Xian_"+freq, "group", "atlas_"+model, g, untreated_xa, y, p
        yield "Xian_"+freq, "group", "HongKong_mean", g, untreated_xa, y, np.mean(hk_observations, axis=0)


def evaluate(records):
    results, detail = [], []
    for dst, unit, model, genes, b, y, p in records:
        assert len(set(genes)) == 101 and np.isfinite([b,y,p]).all()
        strata = np.searchsorted(np.quantile(b, [.2,.4,.6,.8]), b, side="left")
        for k in (5, 10, 20):
            z, w = top_weights(p,k), top_weights(y,k)
            hits = z*w
            same = (np.sign(p) == np.sign(y)) & (y != 0)
            ex = ex_signed = 0.
            for s in np.unique(strata):
                subset = strata == s
                zz, ww, n = z[subset], w[subset], subset.sum()
                ex += zz.sum()*ww.sum()/n
                for sign in (-1, 1):
                    ex_signed += zz[np.sign(p[subset])==sign].sum()*ww[np.sign(y[subset])==sign].sum()/n
            key = dict(destination=dst,unit=unit,forecast=model,k=k)
            results.append(dict(**key,overlap=hits.sum(),signed_overlap=hits[same].sum(),
                                expected=ex,expected_signed=ex_signed,unmatched_random=k*k/101))
            for i in np.flatnonzero((z+w)>0):
                detail.append(dict(**key,gene=genes[i],baseline_quintile=int(strata[i])+1,
                                   predicted=p[i],observed=y[i],forecast_weight=z[i],observed_weight=w[i]))
    return results, detail


def main():
    results, detail = evaluate([*clinical(), *cultures()])
    write(HERE/"scores.tsv",results)
    write(HERE/"shortlists.tsv",detail)
    groups = defaultdict(list)
    for r in results:
        unit = "cohort" if r["unit"]=="cohort_mean" else "participant" if r["destination"] in ("HongKong","Seoul_week3") else "group"
        groups[r["destination"],unit,r["forecast"],r["k"]].append(r)
    summary = []
    for (dst, unit, model, k), members in groups.items():
        summary.append(dict(destination=dst,unit_type=unit,units=len(members),forecast=model,k=k,
                            **{field:np.mean([r[field] for r in members])
                               for field in ("overlap","signed_overlap","expected","expected_signed","unmatched_random")}))
    write(HERE/"summary.tsv",summary)


if __name__ == "__main__":
    main()
