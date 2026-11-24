"""Independently validate target inputs, transferred coefficients, losses and identities."""
import csv
import gzip
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def rows(path):
    with path.open(encoding="utf8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def close(x, y, tol=3e-11):
    np.testing.assert_allclose(x, y, rtol=tol, atol=tol)


def check_subgradient(a, y, b):
    # At a weighted-median kink, the L1 derivative is an interval.
    residual = b*a-y
    at_kink = abs(residual) < 1e-12
    core = np.sum(a[~at_kink]*np.sign(residual[~at_kink]))
    radius = np.sum(abs(a[at_kink]))
    assert core+radius >= -1e-8
    assert b >= 0 and (b==0 or core-radius <= 1e-8)


def rank_columns(counts):
    cols=[]
    for x in counts.T:
        _,ix,n=np.unique(x,return_inverse=True,return_counts=True)
        q=100*(n.cumsum()-(n+1)/2)/(len(x)-1)
        cols.append(q[ix])
    return np.column_stack(cols)


def main():
    original = rows(WORK / "stage14/xian/target_changes.tsv")
    hkobs = {r["gene"]: float(r["hongkong"]) for r in original if r["scale"]=="log2"}
    xaobs = {(r["freq"],r["gene"]):float(r["xian"]) for r in original if r["scale"]=="log2"}
    old = rows(WORK / "stage12/culture/predictions.tsv")
    old_atlas = {
        "tissue":{r["gene"]:float(r["mixture"]) for r in old},
        "fibroblast":{r["gene"]:float(r["fibroblast"]) for r in old},
    }
    target_data = rows(HERE / "targets.tsv")
    assert len(target_data)==606
    for r in target_data:
        g,src,dst = r["gene"],r["source"],r["destination"]
        expected = hkobs[g] if dst=="HK" else xaobs[dst.split("_")[1],g]
        close(float(r["destination_observed"]),expected)
        close(float(r["source_observed"]),hkobs[g] if src=="HK" else (xaobs["0.1",g]+xaobs["0.5",g])/2)
        close(float(r["atlas"]),old_atlas[r["model"]][g])

    # Reconstruct all training effects and atlas inputs independently of score.py.
    with np.load(WORK / "stage13/external/paired_counts.npz") as hk:
        genes=hk["genes"]
        hcount=hk["counts"]
        htotal=hk["total"]
        hc=np.log2(1+hk["counts"]*1e6/hk["total"])
        hdy=(hc[:,5:]-hc[:,:5]).mean(axis=1)
    with np.load(WORK / "stage14/xian/xian_counts.npz") as x:
        assert np.array_equal(x["genes"],genes)
        xc=np.log2(1+x["counts"]*1e6/x["total"])
        xcount=x["counts"]
        xtotal=x["total"]
        arms=x["arm"]
        runs=x["run"].tolist()
        xdy=(xc[:,arms=="0.1"].mean(1)+xc[:,arms=="0.5"].mean(1))/2-xc[:,arms=="0"].mean(1)
    group={r["gene"]:r["ecm_component"] for r in rows(WORK / "stage5/reference/features.tsv")}
    eligible=np.array([group[g]=="outside_ecm" for g in genes])
    assert eligible.sum()==10588 and not any(g in hkobs for g in genes[eligible])
    with gzip.open(WORK / "stage8/axial_log2_multipliers.tsv.gz","rt",encoding="utf8") as fh:
        maps={r["gene"]:float(r["E3"]) for r in csv.DictReader(fh,delimiter="\t")}
    with np.load(WORK / "stage10/reference_counts.npz") as ref:
        ids={g:i for i,g in enumerate(ref["genes"])}
        idx=np.array([ids[g] for g in genes])
        k=list(ref["classes"]).index("fibroblast_like")
        after=1+ref["TP3D"][idx,k]/ref["TP3D_total"][k]*1e6
        before=1+ref["CP"][idx,k]/ref["CP_total"][k]*1e6
        a_fibro=np.log2(after)-np.log2(before)
    atlas={"tissue":np.array([maps[g] for g in genes]), "fibroblast":a_fibro}
    parameters = rows(HERE / "parameters.tsv")
    lookup={}
    for r in parameters:
        src,model,rule=r["source"],r["model"],r["rule"]
        c,b=float(r["intercept"]),float(r["slope"])
        lookup[src,model,rule]=c,b
        a=atlas[model][eligible]
        y={"HK":hdy,"XA":xdy}[src][eligible]
        assert int(r["training_genes"])==len(y)
        if rule=="LS_scale":
            close(b,max(0,float(np.linalg.lstsq(a[:,None],y,rcond=None)[0][0])))
        elif rule=="LAD_scale":
            check_subgradient(a,y,b)
        elif rule=="constant_LS":
            close([c,b],[y.mean(),0.])
        elif rule=="affine_LS":
            braw,craw=np.linalg.lstsq(np.column_stack((a,np.ones_like(a))),y,rcond=None)[0]
            close([c,b],[craw,braw] if braw>=0 else [y.mean(),0.])
        else:
            close([c,b],[0,{"zero":0,"quarter":.25,"full":1}[rule]])
    score_table=rows(HERE / "scores.tsv")
    for r in score_table:
        s,dst,model,rule=[r[k] for k in ("source","destination","model","rule")]
        mask={"target":np.array([g in hkobs for g in genes]),
              "all":np.ones(len(genes),bool),"outside_ecm":eligible}[r["scope"]]
        assert mask.sum()==int(r["genes"])
        y=(hdy if dst=="HK" else xc[:,arms==dst.split("_")[1]].mean(1)-xc[:,arms=="0"].mean(1))[mask]
        c,b=lookup[s,model,rule]
        p=c+b*atlas[model][mask]
        mae=np.mean(abs(p-y))
        flat=np.mean(abs(p.mean()-y))
        close([float(r[k]) for k in ("gene_mae","gain","panel_error","uniform_same_mean_mae","gain_over_uniform")],
              [mae,abs(y).mean()-mae,abs(p.mean()-y.mean()),flat,flat-mae])
    for r in rows(HERE / "oracle.tsv"):
        data=[q for q in target_data if q["destination"]==r["destination"] and q["model"]==r["model"]]
        a=np.array([float(q["atlas"]) for q in data])
        y=np.array([float(q["destination_observed"]) for q in data])
        b=float(r["oracle_LAD_scale"])
        check_subgradient(a,y,b)
        # Check the optimum against every admissible breakpoint, including zero.
        candidates=np.r_[0.,(y[a!=0]/a[a!=0])]
        feasible=candidates[candidates>=0]
        close(float(r["gene_mae"]),abs(a[None,:]*feasible[:,None]-y[None,:]).mean(1).min())
    for r in rows(HERE / "assignment.tsv"):
        src,dst,model,rule = [r[k] for k in ("source","destination","model","rule")]
        data=[q for q in target_data if q["source"]==src and q["destination"]==dst and q["model"]==model]
        c,b=lookup[src,model,rule]
        average=0.
        for q in data:
            matching=[s for s in data if s["target_decile"]==q["target_decile"]]
            mean_error=sum(abs(c+b*float(s["atlas"])-float(q["destination_observed"])) for s in matching)/len(matching)
            average+=mean_error/len(data)
        close(average,float(r["expected_shuffled_mae"]))
    donors=rows(HERE / "donors.tsv")
    deletions=rows(HERE / "deletions.tsv")
    assert len(parameters)==28 and len(score_table)==126 and len(donors)==70 and len(deletions)==340
    mask=np.array([g in hkobs for g in genes])
    for r in donors:
        donor=int(r["donor"])-1
        y=(hc[:,donor+5]-hc[:,donor])[mask]
        c,b=lookup["XA",r["model"],r["rule"]]
        close(float(r["gene_mae"]),abs(c+b*atlas[r["model"]][mask]-y).mean())
    source_cache={}
    for r in deletions:
        src,dst,model,rule,lab,omitted=[r[k] for k in ("source","destination","model","rule","omitted_lab","omitted")]
        if (lab,omitted) not in source_cache:
            hkkeep=np.ones(5,dtype=bool)
            xkeep=np.ones(len(arms),dtype=bool)
            if lab=="HK":
                hkkeep[int(omitted.split("_")[1])-1]=False
            else:
                xkeep[runs.index(omitted)]=False
            hkref=hc[:,5:][:,hkkeep].mean(1)-hc[:,:5][:,hkkeep].mean(1)
            xresponses=[xc[:,(arms==hz)&xkeep].mean(1)-xc[:,(arms=="0")&xkeep].mean(1) for hz in ("0.1","0.5")]
            source_cache[lab,omitted]={"HK":hkref, "XA":np.mean(xresponses,axis=0),
                                       "XA_0.1":xresponses[0], "XA_0.5":xresponses[1]}
        alt=source_cache[lab,omitted]
        a=atlas[model][eligible]
        ytrain=alt[src][eligible]
        c,b=float(r["intercept"]),float(r["slope"])
        if rule=="LS_scale":
            close(b,max(0,float(np.linalg.lstsq(a[:,None],ytrain,rcond=None)[0][0])))
        elif rule=="LAD_scale":
            check_subgradient(a,ytrain,b)
        elif rule=="affine_LS":
            br,cr=np.linalg.lstsq(np.column_stack([a,np.ones_like(a)]),ytrain,rcond=None)[0]
            close([c,b],[cr,br] if br>=0 else [ytrain.mean(),0.])
        else:
            close([c,b],[0.,{"quarter":.25,"full":1.}[rule]])
        y,p=alt[dst][mask],c+b*atlas[model][mask]
        mae=abs(p-y).mean()
        close([float(r[k]) for k in ("gene_mae","gain","gain_over_uniform")],
              [mae,abs(y).mean()-mae,abs(p.mean()-y).mean()-mae])
    for dst,lab,n in [("XA_0.1","XA",5),("XA_0.5","XA",6),("HK","XA",8),
                      ("XA_0.1","HK",5),("XA_0.5","HK",5),("HK","HK",5)]:
        subset=[r for r in deletions if r["destination"]==dst and r["omitted_lab"]==lab]
        assert len({r["omitted"] for r in subset})==n and len(subset)==10*n
    rank_rows=rows(HERE/"ranks.tsv")+rows(HERE/"rank_deletions.tsv")
    assert len(rank_rows)==62
    cached={}
    for r in rank_rows:
        src,dst,model,rule,deleted=[r[k] for k in ("source","destination","model","rule","omitted")]
        if (src,model,rule) not in cached:
            raw,total=(hcount,htotal) if src=="XA" else (xcount,xtotal)
            base=raw*1e6/total
            close(float(r["scale"]),lookup[src,model,rule][1])
            # Use the original published log ratios before exponentiation to match numerical ties.
            shifts=atlas[model] if model=="tissue" else np.log2(after/before)
            loaded=base*np.exp2(float(r["scale"])*shifts[:,None])
            cached[src,model,rule]=rank_columns(base)[mask],rank_columns(loaded)[mask]
        rbase,rnew=cached[src,model,rule]
        if dst=="HK":
            ctrl=np.arange(5)
            if deleted!="none":
                ctrl=np.delete(ctrl,int(deleted.split("_")[1])-1)
            load=ctrl+5
        else:
            keep=np.arange(len(arms))!=(-1 if deleted=="none" else runs.index(deleted))
            ctrl=np.where((arms=="0")&keep)[0]
            load=np.where((arms==dst.split("_")[1])&keep)[0]
        y=rbase[:,load].mean(1)-rbase[:,ctrl].mean(1)
        p=(rnew[:,ctrl]-rbase[:,ctrl]).mean(1)
        close([float(r[k]) for k in ("gene_mae","gain","panel_error")],
              [abs(p-y).mean(),abs(y).mean()-abs(p-y).mean(),abs(p.mean()-y.mean())])
    print("Verified 606 target records, 28 coefficient sets, 126 scores, 42 expected reassignments, "
          "six diagnostic optima, 70 donor scores, 340 single deletions and 62 rank-endpoint checks.")


if __name__=="__main__":
    main()
