"""Independent check using archived forecasts, original counts, and sorted-tie ranks."""
import csv
import gzip
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def rows(path):
    with path.open(encoding="utf8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def ranks(values):
    order = np.sort(values)
    return 50*(order.searchsorted(values, side="left")+order.searchsorted(values, side="right")-1)/(len(order)-1)


def close(x, y):
    np.testing.assert_allclose(x, y, atol=1e-10, rtol=0)


def main():
    with np.load(HERE.parent / "xian/xian_counts.npz") as d:
        genes, mask, arm, run = [d[k] for k in ("genes","target","arm","run")]
        cpm = (d["counts"].T / d["total"][:,None]) * 1e6
    with np.load(WORK / "stage10/reference_counts.npz") as d:
        indices = {g:i for i,g in enumerate(d["genes"])}
        ix = [indices[g] for g in genes]
        cell = tuple(d["classes"]).index("fibroblast_like")
        models = {"fibroblast": np.log2(1+d["TP3D"][ix,cell]/d["TP3D_total"][cell]*1e6)
                               - np.log2(1+d["CP"][ix,cell]/d["CP_total"][cell]*1e6)}
        cp = d["CP"][ix,cell]*1e6/d["CP_total"][cell]
        tp = d["TP3D"][ix,cell]*1e6/d["TP3D_total"][cell]
        fibroblast_ratio = np.exp2(np.log2((tp+1)/(cp+1)))
    with gzip.open(WORK / "stage8/axial_log2_multipliers.tsv.gz","rt") as stream:
        tissue = {row["gene"]:float(row["E3"]) for row in csv.DictReader(stream,delimiter="\t")}
    models["tissue"] = np.array([tissue[g] for g in genes])
    models["zero"] = np.zeros(len(genes))
    archived = {r["gene"]:r for r in rows(WORK / "stage12/culture/predictions.tsv")}
    assert set(archived) == set(genes[mask]) and len(genes[mask])==101
    for name, prior in [("tissue","mixture"),("fibroblast","fibroblast"),("zero","zero")]:
        close(models[name][mask], [float(archived[g][prior]) for g in genes[mask]])
    observed_before = {(r["freq"],r["scale"],r["gene"]):float(r["xian"]) for r in rows(HERE.parent / "xian/target_changes.tsv")}
    effects = {(r["frequency_Hz"],r["scale"],r["method"],r["gene"]):r for r in rows(HERE / "targets.tsv")}
    assert len(effects)==2*2*3*101
    complete, deletions = rows(HERE / "scores.tsv"), rows(HERE / "deletions.tsv")
    assert len(complete)==12 and len(deletions)==66
    for r in complete+deletions:
        scale, mode, freq = r["scale"],r["method"],r["frequency_Hz"]
        include = run != r.get("omitted_run","")
        ctrl = cpm[include & (arm=="0")]
        load = cpm[include & (arm==freq)]
        assert (len(ctrl),len(load)) == (int(r["control_runs"]),int(r["loaded_runs"]))
        if scale=="log2":
            y = np.log2(1+load[:,mask]).mean(0)-np.log2(1+ctrl[:,mask]).mean(0)
            p = models[mode][mask]
        else:
            base = np.mean([ranks(c)[mask] for c in ctrl],axis=0)
            y = np.mean([ranks(c)[mask] for c in load], axis=0)-base
            ratios = fibroblast_ratio if mode=="fibroblast" else np.exp2(models[mode])
            p = np.mean([ranks(c*ratios)[mask] for c in ctrl], axis=0)-base
        d = p-y
        calculated = [np.abs(d).mean(),np.abs(y).mean()-np.abs(d).mean(),abs(d.mean()),abs(y.mean())-abs(d.mean())]
        close(calculated, [float(r[k]) for k in ("gene_mae","gain_vs_zero","panel_abs_error","panel_gain_vs_zero")])
        if "omitted_run" not in r:
            close(y, [observed_before[freq,scale,g] for g in genes[mask]])
            close(y, [float(effects[freq,scale,mode,g]["observed"]) for g in genes[mask]])
            close(p, [float(effects[freq,scale,mode,g]["prediction"]) for g in genes[mask]])
    print("Verified 12 primary scores, 66 single-run deletion scores and 1,212 target rows.")


if __name__ == "__main__":
    main()
