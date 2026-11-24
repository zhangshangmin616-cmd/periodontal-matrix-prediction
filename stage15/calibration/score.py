"""Train one amplitude on the other laboratory's non-matrix genes."""
import csv
import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def load_atlas(genes):
    spec = importlib.util.spec_from_file_location("atlas_models", WORK / "stage13/external/score.py")
    frozen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(frozen)
    models = frozen.atlas(genes)
    return {key: models[key] for key in ("tissue", "fibroblast")}


def load_data():
    with np.load(WORK / "stage13/external/paired_counts.npz") as pair:
        genes, target = pair["genes"], pair["target"]
        h = pair["counts"] * 1e6 / pair["total"]
    with np.load(WORK / "stage14/xian/xian_counts.npz") as src:
        assert np.array_equal(genes, src["genes"]) and np.array_equal(target, src["target"])
        x = src["counts"] * 1e6 / src["total"]
        arms, runs = src["arm"], src["run"]
    with (WORK / "stage5/reference/features.tsv").open(encoding="utf8") as fh:
        pathway = {r["gene"]: r["ecm_component"] for r in csv.DictReader(fh, delimiter="\t")}
    outside = np.array([pathway[g] == "outside_ecm" for g in genes])
    assert np.array_equal(target, [pathway[g] == "other" for g in genes])
    assert (len(genes), int(target.sum())) == (10805, 101)
    return genes, target, outside, h, x, arms, runs


def views(h, x, arms, drop_h=None, drop_x=None):
    hkeep = np.arange(5) != (-1 if drop_h is None else drop_h)
    xkeep = np.arange(len(arms)) != (-1 if drop_x is None else drop_x)
    hl, xl = np.log2(1+h), np.log2(1+x)
    paired = (hl[:, 5:]-hl[:, :5])[:, hkeep]
    baseline = xl[:, (arms == "0") & xkeep].mean(1)
    low, high = [xl[:, (arms == hz) & xkeep].mean(1)-baseline for hz in ("0.1","0.5")]
    return ({"HK": paired.mean(1), "XA": (low+high)/2},
            {"HK": paired.mean(1), "XA_0.1": low, "XA_0.5": high}, paired)


def lad_scale(a, y):
    nz = a != 0
    if not nz.any():
        return 0.0
    ratio = y[nz]/a[nz]
    order = np.argsort(ratio, kind="stable")
    weights = abs(a[nz][order])
    idx = np.searchsorted(np.cumsum(weights), weights.sum()/2)
    return float(max(0,ratio[order[idx]]))


def coefficients(a, y):
    scale = max(0, (a*y).sum()/(a*a).sum())
    ac, yc = a-a.mean(), y-y.mean()
    centred = max(0,(ac*yc).sum()/(ac*ac).sum())
    return {
        "zero": (0.,0.), "full": (0.,1.), "quarter": (0.,.25),
        "LS_scale": (0.,scale), "LAD_scale": (0.,lad_scale(a,y)),
        "constant_LS": (y.mean(),0.), "affine_LS": (y.mean()-centred*a.mean(),centred)}


def evaluate(y, p):
    signed = p-y
    mae = abs(signed).mean()
    const = p.mean()
    flat = abs(const-y).mean()
    active = (p!=0) & (y!=0)
    return dict(gene_mae=mae, gain=abs(y).mean()-mae, panel_error=abs(signed.mean()),
                panel_gain=abs(y.mean())-abs(signed.mean()), uniform_same_mean_mae=flat,
                gain_over_uniform=flat-mae, predicted_mean=const, observed_mean=y.mean(),
                nonzero_pairs=int(active.sum()), same_sign=int(np.sum((p[active]*y[active])>0)))


def write(name, rows):
    with (HERE / name).open("w", encoding="utf-8", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def main():
    genes,target,outside,h,x,arms,runs = load_data()
    models = load_atlas(genes)
    source,destination,hpairs = views(h,x,arms)
    scopes = {"target":target, "all":np.ones_like(target), "outside_ecm":outside}
    contrasts = [("HK","XA_0.1"),("HK","XA_0.5"),("XA","HK")]
    scores, params, donors, shuffle, oracle, effects = [], [], [], [], [], []
    controls = {"HK":h[:,:5].mean(1), "XA":x[:,arms=="0"].mean(1)}
    for model, atlas in models.items():
        for src, measured in source.items():
            fits = coefficients(atlas[outside], measured[outside])
            for rule,(intercept,slope) in fits.items():
                params.append(dict(source=src,model=model,rule=rule,training_genes=int(outside.sum()),
                                   intercept=intercept,slope=slope))
            ctrl = controls[src]
            bins = np.searchsorted(np.quantile(ctrl,np.arange(1,10)/10),ctrl[target],side="right")
            for _,dst in [pair for pair in contrasts if pair[0]==src]:
                obs = destination[dst]
                for i in np.flatnonzero(target):
                    effects.append(dict(source=src,destination=dst,model=model,gene=genes[i],
                                        source_observed=measured[i],destination_observed=obs[i],atlas=atlas[i],
                                        source_control_cpm=ctrl[i],target_decile=int(bins[np.sum(target[:i])])))
                fit_oracle = lad_scale(atlas[target],obs[target])
                oracle.append(dict(destination=dst,model=model,oracle_LAD_scale=fit_oracle,
                                   **evaluate(obs[target],fit_oracle*atlas[target])))
                for rule, (intercept,slope) in fits.items():
                    p = intercept+slope*atlas
                    names = dict(source=src,destination=dst,model=model,rule=rule)
                    for scope,mask in scopes.items():
                        scores.append(dict(**names,scope=scope,genes=int(mask.sum()),
                                           **evaluate(obs[mask],p[mask])))
                    if dst=="HK":
                        for j,yp in enumerate(hpairs.T,1):
                            donors.append(dict(**names,donor=j,**evaluate(yp[target],p[target])))
                    original = abs(p[target]-obs[target]).mean()
                    total = 0.
                    for band in np.unique(bins):
                        keep = bins==band
                        costs = abs(p[target][keep,None]-obs[target][None,keep])
                        total += costs.sum()/int(keep.sum())
                    expected = total/int(target.sum())
                    shuffle.append(dict(**names,genes=int(target.sum()),expected_shuffled_mae=expected,
                                        observed_mae=original,gain_over_shuffled=expected-original))
    deletions = []
    for lab,n in (("HK",5),("XA",len(arms))):
        for j in range(n):
            sv,dv,_ = views(h,x,arms,drop_h=j if lab=="HK" else None,drop_x=j if lab=="XA" else None)
            for src,dst in contrasts:
                if lab=="XA" and src=="HK" and arms[j] not in ("0",dst.split("_")[1]):
                    continue
                for model,atlas in models.items():
                    fitted = coefficients(atlas[outside],sv[src][outside])
                    for rule in ("full","quarter","LS_scale","LAD_scale","affine_LS"):
                        intercept,slope=fitted[rule]
                        preds = intercept+slope*atlas[target]
                        deletions.append(dict(source=src,destination=dst,omitted_lab=lab,
                                              omitted=f"donor_{j+1}" if lab=="HK" else runs[j],model=model,rule=rule,
                                              intercept=intercept,slope=slope,**evaluate(dv[dst][target],preds)))
    for path,rows in [("parameters.tsv",params),("scores.tsv",scores),("donors.tsv",donors),
                      ("assignment.tsv",shuffle),("oracle.tsv",oracle),("targets.tsv",effects),
                      ("deletions.tsv",deletions)]:
        write(path,rows)


if __name__=="__main__":
    main()
