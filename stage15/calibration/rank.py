"""Apply the calibrated log multipliers before ranking each untreated profile."""
import csv

import numpy as np

from score import HERE, evaluate, load_atlas, load_data, write


def percentile_profiles(cpm):
    result=np.empty(cpm.shape,dtype=float)
    for i,column in enumerate(cpm.T):
        ordered=np.sort(column)
        left,right=np.searchsorted(ordered,column,"left"),np.searchsorted(ordered,column,"right")
        result[:,i]=50*(left+right-1)/(len(ordered)-1)
    return result


def main():
    genes,target,_,h,x,arms,runs=load_data()
    atlas=load_atlas(genes)
    with (HERE/"parameters.tsv").open(encoding="utf8") as stream:
        params={(r["source"],r["model"],r["rule"]):float(r["slope"])
                for r in csv.DictReader(stream,delimiter="\t")}
    hong=percentile_profiles(h)[target]
    xian=percentile_profiles(x)[target]
    hk_change=hong[:,5:]-hong[:,:5]
    xc=x[:,arms=="0"]
    xctrl=xian[:,arms=="0"]
    main_rows,deletion_rows=[],[]
    for dst in ("XA_0.1","XA_0.5","HK"):
        source="XA" if dst=="HK" else "HK"
        if dst=="HK":
            control,ref=h[:,:5],hong[:,:5]
            loaded=hong[:,5:]
            identical=np.arange(5)
            deletions=[(f"donor_{j+1}",identical!=j,identical!=j) for j in range(5)]
        else:
            arm=dst.split("_")[1]
            control,ref=xc,xctrl
            loaded=xian[:,arms==arm]
            deletions=[(runs[j],np.flatnonzero(arms=="0")!=j,np.flatnonzero(arms==arm)!=j)
                       for j in np.flatnonzero((arms=="0")|(arms==arm))]
        omitted=[("none",np.ones(ref.shape[1],bool),np.ones(loaded.shape[1],bool))]
        for model,p in atlas.items():
            for rule in ("zero","full","quarter","LS_scale","LAD_scale"):
                b=params[source,model,rule]
                shift=percentile_profiles(control*np.exp2(b*p[:,None]))[target]-ref
                use=omitted+deletions if rule=="quarter" else omitted
                for deleted,keepctrl,keeploaded in use:
                    f=shift[:,keepctrl].mean(1)
                    if dst=="HK":
                        y=hk_change[:,keeploaded].mean(1)
                    else:
                        y=loaded[:,keeploaded].mean(1)-ref[:,keepctrl].mean(1)
                    data=dict(source=source,destination=dst,model=model,rule=rule,scale=b,
                              omitted=deleted,genes=int(target.sum()),**evaluate(y,f))
                    if deleted=="none":
                        main_rows.append(data)
                    else:
                        deletion_rows.append(data)
    write("ranks.tsv",main_rows)
    write("rank_deletions.tsv",deletion_rows)


if __name__=="__main__":
    main()
