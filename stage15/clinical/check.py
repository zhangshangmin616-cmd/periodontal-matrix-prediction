"""Verify errors from saved full-scale inputs and optima by convexity or enumeration."""
import csv
import gzip
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def records(name):
    with name.open(encoding="utf8") as stream:
        return list(csv.DictReader(stream,delimiter="\t"))


def same(left,right):
    np.testing.assert_allclose(left,right,atol=3e-11,rtol=1e-11)


def endpoints(x,y,b):
    predicted = b*x
    return np.abs(y-predicted).mean(), np.abs((y-predicted).mean(axis=1)).mean()


def verify_fit(x,y,b,expected):
    """Count the sign change in the derivative and compare every target breakpoint."""
    residual = b*x-y
    tolerance = 1e-11
    ties = np.abs(residual) < tolerance
    known = (np.sign(residual)*x)[~ties].sum()
    total = np.abs(x[ties]).sum()
    assert known+total >= -1e-7
    if b > 0:
        assert known-total <= 1e-7
    same(np.abs(residual).mean(),expected)
    if x.size<=1200:
        candidates = np.unique(np.r_[0., (y[x!=0]/x[x!=0]).ravel()])
        candidates = candidates[candidates>=0]
        losses = np.abs(candidates[:,None]*x.ravel()-y.ravel()).mean(axis=1)
        same(min(losses),expected)


def leuven_removals(genes,target):
    """Rank target counts against every mapped gene, then recompute arm means."""
    directory = ROOT/"stage5/zhao"
    mapping = {r["gene"]:r["entrez_id"] for r in records(directory/"feature_map.tsv")}
    ids = [mapping[g] for g in genes]
    arms = {}
    for row in records(directory/"joined_manifest.tsv"):
        run = row["run"]
        with gzip.open(directory/"counts"/(run+".tsv.gz"),"rt") as stream:
            header = next(stream)
            assert header.startswith("Geneid")
            counts = {k:float(v) for k,v in csv.reader(stream,delimiter="\t")}
        c = np.array([counts[k] for k in ids])
        sorted_counts = np.sort(c)
        lo = np.searchsorted(sorted_counts,c[target],side="left")
        hi = np.searchsorted(sorted_counts,c[target],side="right")-1
        rank = (lo+hi)*50./(len(ids)-1)
        arms.setdefault((row["batch"],row["group"]),{})[run] = rank
    outcomes = {}
    for day,batch in [("day7","4"),("day28","3")]:
        baseline = np.mean(list(arms[batch,"control"].values()),axis=0)
        for run in arms[batch,day]:
            retained = [v for key,v in arms[batch,day].items() if key!=run]
            outcomes["Leuven_"+day,run] = (np.mean(retained,axis=0)-baseline)[None,:]
    return outcomes


def main():
    with np.load(ROOT/"stage11/specificity/rank_changes.npz") as src:
        units,genes,kinds = src["units"],src["gene"],src["component"]
        full_y,full_p = src["observed"],src["fitted"]
    target = kinds=="other"
    outside = kinds=="outside_ecm"
    assert target.sum()==101 and outside.sum()==10589
    masks = {"Seoul_week3":np.char.startswith(units,"Kim_"),
             "Leuven_day7":units=="Zhao_day7", "Leuven_day28":units=="Zhao_day28"}
    slopes = {}
    for row in records(HERE/"coefficients.tsv"):
        who = masks[row["source"]]
        x,y = full_p[who][:,outside],full_y[who][:,outside]
        if row["estimand"]=="source_group_mean":
            x,y=x.mean(axis=0),y.mean(axis=0)
        signed = float(row["signed_slope"])
        verify_fit(x,y,signed,np.abs(y-signed*x).mean())
        constrained = float(row["nonnegative_slope"])
        same(constrained,max(0,signed))
        verify_fit(x,y,constrained,float(row["training_gene_mae"]))
        slopes[row["source"],row["estimand"]] = constrained

    def retrieve(row):
        who = masks[row["destination"]].copy()
        if row.get("unit","all")!="all":
            who &= units==row["unit"]
        if row.get("removed","").startswith("Kim_"):
            who &= units!=row["removed"]
        x,y = full_p[who][:,target],full_y[who][:,target]
        if row["estimand"]=="group_mean":
            x,y = np.mean(x,axis=0,keepdims=True),np.mean(y,axis=0,keepdims=True)
        return x,y

    scores = records(HERE/"scores.tsv")
    for row in scores:
        x,y=retrieve(row)
        rule=row["rule"]
        b=slopes[row["source"],rule] if rule.startswith("source_") else {"full":1,"quarter":.25,"zero":0}[rule]
        if row["assignment"]=="common_panel_shift":
            x=x.mean(axis=1,keepdims=True)
        error,panel=endpoints(x,y,b)
        same([error,panel,b],[float(row[k]) for k in ("gene_mae","panel_mae","slope")])
    removed=leuven_removals(genes,target)
    oracles=records(HERE/"oracle.tsv")+records(HERE/"deletions.tsv")
    for row in oracles:
        x,y=retrieve(row)
        if (row["destination"],row.get("removed","")) in removed:
            y=removed[row["destination"],row["removed"]]
        m=x.mean(axis=1,keepdims=True) if row["family"]=="preserve_mean" else np.zeros_like(x)
        b=float(row["nonnegative_optimum"])
        verify_fit(x-m,y-m,b,float(row["gene_mae"]))
        if "gain_vs_zero" in row:
            same(np.abs(y).mean()-float(row["gene_mae"]),float(row["gain_vs_zero"]))
    curves=records(HERE/"curves.tsv")
    for row in curves:
        x,y=retrieve(row)
        m=x.mean(axis=1,keepdims=True) if row["family"]=="preserve_mean" else 0.
        a,h=endpoints(x-m,y-m,float(row["slope"]))
        same([a,h],[float(row["gene_mae"]),float(row["panel_mae"])])
    print(f"Verified 6 source fits, {len(scores)} forecasts, {len(oracles)} exact outcome fits and {len(curves)} curve points.")


if __name__=="__main__":
    main()
