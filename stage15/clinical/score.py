"""Score cross-study magnitude correction and MAE lower bounds for a fixed atlas."""
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOL = 1e-10


def write(name, rows):
    with (HERE / name).open("w", encoding="utf8", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def identify(p, y):
    """Lower weighted median for |y-b*p|; direction statistics retain zero y."""
    p, y = p.ravel(), y.ravel()
    use = p != 0
    weights = np.abs(p[use])
    if not weights.size:
        return 0.0
    ratio = y[use] / p[use]
    order = np.argsort(ratio, kind="stable")
    cum = np.cumsum(weights[order])
    return float(ratio[order[np.searchsorted(cum, .5 * cum[-1], side="left")]])


def diagnostics(x, y):
    magnitude = np.abs(x)
    valid = (magnitude > TOL) & (np.abs(y) > TOL)
    concord = valid & (np.sign(x) == np.sign(y))
    above = np.where(y == 0, magnitude, -x * np.sign(y))
    return dict(right_derivative=float(above.mean()),
                paired_nonzero=int(valid.sum()), concordant=int(concord.sum()),
                concordant_share=float(magnitude[concord].sum() / magnitude[valid].sum()))


def rowmetrics(forecast, outcome):
    diff = forecast - outcome
    return dict(gene_mae=float(np.abs(diff).mean()),
                panel_mae=float(np.abs(diff.mean(axis=1)).mean()))


def boundaries(base, forecast):
    clipped = np.clip(forecast, -base, 100 - base)
    return clipped, int((np.abs(clipped-forecast) > TOL).sum())


def target_views(study, labels, p, y, base, units=True):
    if study == "Seoul_week3":
        yield "participant_first", "all", p, y, base
        if units:
            for i, name in enumerate(labels):
                yield "one_participant", name, p[i:i+1], y[i:i+1], base[i:i+1]
    yield "group_mean", "all", p.mean(axis=0, keepdims=True), y.mean(axis=0, keepdims=True), base.mean(axis=0, keepdims=True)


def main():
    with np.load(ROOT / "stage11/specificity/rank_changes.npz") as z:
        labels, genes, kinds = z["units"], z["gene"], z["component"]
        fitted, observed = z["fitted"], z["observed"]
    target, train = kinds == "other", kinds == "outside_ecm"
    assert train.sum() == 10589 and target.sum() == 101
    groups = {"Seoul_week3": np.flatnonzero(np.char.startswith(labels,"Kim_")),
              "Leuven_day7": np.flatnonzero(labels == "Zhao_day7"),
              "Leuven_day28": np.flatnonzero(labels == "Zhao_day28")}
    with (ROOT / "stage12/review/control_ranks.tsv").open(encoding="utf8") as stream:
        lookup = {(r["unit"],r["gene"]):float(r["untreated_rank"]) for r in csv.DictReader(stream,delimiter="\t")}
    base = np.array([[lookup[lab,g] for g in genes[target]] for lab in labels])
    coefficients, scores, oracle, curves = [], [], [], []
    for source,dest in [("Seoul_week3","Leuven_day7"),("Seoul_week3","Leuven_day28"),("Leuven_day7","Seoul_week3")]:
        si, di = groups[source], groups[dest]
        sx, sy = fitted[si][:,train], observed[si][:,train]
        tp, ty = fitted[di][:,target], observed[di][:,target]
        scales = {"zero":0., "quarter":.25, "full":1.}
        for estimand,x,y in [("source_unit_first",sx,sy),
                             ("source_group_mean",sx.mean(axis=0),sy.mean(axis=0))]:
            signed = identify(x,y)
            scale = max(0.,signed)
            scales[estimand] = scale
            coefficients.append(dict(source=source,destination=dest,estimand=estimand,genes=int(train.sum()),
                source_units=len(si),signed_slope=signed,nonnegative_slope=scale,
                training_zero_mae=float(np.abs(y).mean()),training_gene_mae=float(np.abs(y-scale*x).mean()),
                **diagnostics(x,y)))
        for view,lab,p,y,baseline in target_views(dest,labels[di],tp,ty,base[di]):
            background = dict(destination=dest,estimand=view,unit=lab)
            zero = rowmetrics(np.zeros_like(y),y)
            m = np.broadcast_to(p.mean(axis=1,keepdims=True),p.shape)
            for rule,scale in scales.items():
                for identity,term in [("assigned",p),("common_panel_shift",m)]:
                    guess = scale * term
                    fit = rowmetrics(guess,y)
                    bounded,n = boundaries(baseline,guess)
                    scores.append(dict(source=source,**background,rule=rule,assignment=identity,
                        slope=scale,**fit,gene_gain=zero["gene_mae"]-fit["gene_mae"],
                        panel_gain=zero["panel_mae"]-fit["panel_mae"],outside_bounds=n,
                        clipped_gene_mae=rowmetrics(bounded,y)["gene_mae"]))
            for family,offset,x in [("amplitude",np.zeros_like(p),p),("preserve_mean",m,p-m)]:
                signed = identify(x,y-offset)
                b = max(0.,signed)
                optimum = offset + b*x
                fit = rowmetrics(optimum,y)
                clip,n = boundaries(baseline,optimum)
                anchor = rowmetrics(offset,y)["gene_mae"]
                bounded = offset + min(1.,b)*x
                oracle.append(dict(**background,family=family,n_observations=y.shape[0],targets=y.shape[1],
                    signed_optimum=signed,nonnegative_optimum=b,**fit,
                    zero_gene_mae=zero["gene_mae"],gain_vs_zero=zero["gene_mae"]-fit["gene_mae"],
                    anchor_gene_mae=anchor,gain_vs_anchor=anchor-fit["gene_mae"],
                    bounded_01_optimum=min(1.,b),bounded_01_gene_mae=rowmetrics(bounded,y)["gene_mae"],
                    outside_rank_bounds=n,clipped_gene_mae=rowmetrics(clip,y)["gene_mae"],**diagnostics(x,y-offset)))
                if view != "one_participant":
                    for t in np.linspace(0,1,101):
                        curves.append(dict(**background,family=family,slope=float(t),**rowmetrics(offset+t*x,y)))
    write("coefficients.tsv",coefficients)
    write("scores.tsv",scores)
    write("oracle.tsv",oracle)
    write("curves.tsv",curves)


if __name__ == "__main__":
    main()
