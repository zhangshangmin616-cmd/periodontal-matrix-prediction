"""Test a source-fitted amplitude while keeping every predicted panel mean fixed."""
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def read(path):
    with path.open(encoding="utf8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(name, rows):
    with (HERE/name).open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def inputs(rows, obs, forecast):
    ordered = sorted(rows, key=lambda r:r["gene"])
    assert len({r["gene"] for r in ordered}) == len(ordered) == 101
    p = np.array([float(r[forecast]) for r in ordered])
    y = np.array([float(r[obs]) for r in ordered])
    return p,y


def slope(deviations, residuals):
    x,y = deviations.ravel(), residuals.ravel()
    usable = x != 0
    if not np.any(usable):
        return 0.0
    ratio = y[usable]/x[usable]
    weights = abs(x[usable])
    order = np.argsort(ratio,kind="stable")
    half = weights.sum()/2
    root = ratio[order[np.searchsorted(np.cumsum(weights[order]), half)]]
    return max(0.0,root)


def main():
    training = [r for r in read(WORK/"stage13/external/hongkong_target_effects.tsv") if r["method"]=="tissue"]
    hk=[inputs([r for r in training if int(r["donor"])==j],"rank_observed","rank_predicted") for j in range(1,6)]
    p,y = np.stack([a for a,_ in hk]),np.stack([o for _,o in hk])
    deviations=p-p.mean(axis=1,keepdims=True)
    residuals=y-p.mean(axis=1,keepdims=True)
    fit=[]
    for omit in range(6):
        include=np.array([j!=omit for j in range(1,6)])
        a=slope(deviations[include],residuals[include])
        fit.append(dict(source_donor_omitted=omit, donor_pairs=int(sum(include)), amplitude=a,
                        source_gene_mae=abs(residuals[include]-a*deviations[include]).mean()))
    save("training.tsv",fit)
    observations=[r for r in read(HERE.parent/"transfer/genes.tsv") if r["method"]=="atlas_tissue"]
    scores=[]
    for target in dict.fromkeys(r["destination"] for r in observations):
        forecast,seen=inputs([r for r in observations if r["destination"]==target],"observed","forecast")
        m=forecast.mean()
        d=forecast-m
        candidate=[("common_shift",0.0),("original_atlas",1.0),("fit_in_Hong_Kong",fit[0]["amplitude"]),
                   ("destination_best_possible",slope(d,seen-m))]
        for label,a in candidate:
            pred=m+a*d
            assert abs(pred.mean()-m)<1e-12
            gain=abs(seen).mean()-abs(seen-pred).mean()
            scores.append(dict(destination=target, model=label, amplitude=a, targets=len(seen),panel_mean=m,
                               panel_error=abs(m-seen.mean()),gene_mae=abs(seen-pred).mean(),gain_over_zero=gain,
                               gain_over_common=abs(seen-m).mean()-abs(seen-pred).mean()))
        optimum=scores[-1]["gene_mae"]
        assert all(s["gene_mae"]>optimum-1e-12 for s in scores[-4:])
        # Piecewise linear loss attains the constrained minimum at 0 or a nonnegative knot.
        knots=(seen[d!=0]-m)/d[d!=0]
        choices=np.r_[0.,knots[knots>=0]]
        best=np.min(np.mean(abs((seen-m)[None,:]-choices[:,None]*d[None,:]),axis=1))
        assert abs(best-optimum)<1e-12
    all_d=deviations.ravel()
    all_r=residuals.ravel()
    cand=all_r[all_d!=0]/all_d[all_d!=0]
    grid=np.r_[0,cand[cand>=0]]
    assert abs(min(abs(all_r-q*all_d).mean() for q in grid)-fit[0]["source_gene_mae"])<1e-12
    save("scores.tsv",scores)
    print("Source fit",fit[0]["amplitude"],"; held-out loss recomputed for",len(scores),"rules and contrasts.")


if __name__=="__main__":
    main()
