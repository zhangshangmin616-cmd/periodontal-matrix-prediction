"""Reconstruct the original fixed-RNA-share rank predictor in an occlusal contrast."""
import csv
import gzip
import itertools
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
PRIOR = WORK/"stage5"
CLASSES = ("fibroblast_like","mural","immune","endothelial","glial","residual")


def read(path):
    opener = gzip.open if path.suffix==".gz" else open
    with opener(path,"rt",encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def write(path, rows):
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def metadata():
    with gzip.open(HERE/"series.gz","rt",encoding="latin1",newline="") as f:
        rows = list(csv.reader(f,delimiter="\t"))
    def get(key,prefix=""):
        return next(r[1:] for r in rows if len(r)>1 and r[0]=="!Sample_"+key and r[1].startswith(prefix))
    accessions=get("geo_accession")
    titles=get("title")
    # Matching blocks identified from extraction-date/individual prefixes in GEO
    # titles. The first matched block contributes 4 libraries, all others 2.
    def block(t):
        explicit=re.match(r'(?:G[12]S(\d)|S(\d+)G[12])',t)
        return "S"+next(x for x in explicit.groups() if x) if explicit else re.split(r'[_-]',t,1)[0]
    names = [block(t) for t in titles]
    numbered={n:f"pair{i+1:02d}" for i,n in enumerate(dict.fromkeys(names))}
    result=[]
    for s,n,arm,age,lab in zip(accessions,names,get("characteristics_ch1","treatment:"),
                               get("characteristics_ch1","age:"),get("description"),strict=True):
        result.append(dict(sample=s,unit=numbered[n],arm={"treatment: normal occlusal contact":"contact",
                           "treatment: no occlusal contact":"hypofunction"}[arm],age=age[5:],library=lab))
    assert len(result)==22 and len(numbered)==10
    for name in numbered.values():
        unit=[r for r in result if r["unit"]==name]
        assert len({r["age"] for r in unit})==1
        assert [r["arm"] for r in unit].count("contact")==len(unit)/2
    return result


def design():
    features=read(PRIOR/"zhao/feature_map.tsv")
    genes=[r["gene"] for r in features]
    targets=np.array([r["ecm_component"]=="other" for r in features])
    ids=[r["entrez_id"] for r in features]
    profile=read(PRIOR/"reference/profiles.tsv.gz")
    assert [r["gene"] for r in profile]==genes and len(genes)==10806 and sum(targets)==101
    def profiles(lib):
        return np.array([[float(r[lib+"/"+c]) for c in CLASSES] for r in profile])*1e6
    base=profiles("CP")
    # Average the two raw linear RNA profiles before forming the fold change.
    loaded=(profiles("PP7D")+profiles("TP7D"))/2
    loaded[:,-1]=base[:,-1]
    groups={r["cell_class"]:float(r["rna_weight"]) for r in read(PRIOR/"reference/groups.tsv") if r["library"]=="CP"}
    w=np.array([groups[c] for c in CLASSES])
    anchors=np.array([genes.index(r["gene"]) for r in read(PRIOR/"results/anchors.tsv")])
    a=(base[anchors]+1)/((base[anchors]@w+1)[:,None])
    return genes,targets,ids,base,loaded,w,anchors,a


def weights(x,calib,w,idx,a):
    y=np.clip((x[idx]+1)/(calib[idx]+1),.25,4)
    lhs=a.T@a/len(y)+.05*np.eye(len(w))
    rhs=a.T@y/len(y)+.05*w
    best=None
    best_loss=np.inf
    for k in range(1,len(w)+1):
        for subset in itertools.combinations(range(len(w)),k):
            ids=np.asarray(subset)
            gram=lhs[np.ix_(ids,ids)]
            system=np.block([[gram,np.ones((k,1))],[np.ones((1,k)),np.zeros((1,1))]])
            sol=np.linalg.solve(system,np.r_[rhs[ids],1])[:-1]
            if min(sol)<-1e-9:
                continue
            q=np.zeros(len(w))
            q[ids]=np.maximum(0,sol)
            q/=q.sum()
            loss=np.mean((a@q-y)**2)+.05*np.sum((q-w)**2)
            if loss<best_loss:
                best,best_loss=q,loss
    assert best is not None
    return best


def ranks(v,targets):
    ordered=np.sort(v)
    selected=v[targets]
    return (np.searchsorted(ordered,selected)+np.searchsorted(ordered,selected,side="right")-1)*50/(len(v)-1)


def read_counts(path,ids):
    with gzip.open(path,"rt",encoding="utf8",newline="") as stream:
        rows=csv.reader(stream,delimiter="\t")
        samples=next(rows)[1:]
        ordered={g:i for i,g in enumerate(ids)}
        counts=np.full((len(ids),len(samples)),np.nan)
        total=np.zeros(len(samples))
        for row in rows:
            v=np.asarray(row[1:],dtype=float)
            total+=v
            if row[0] in ordered:
                assert np.isnan(counts[ordered[row[0]],0])
                counts[ordered[row[0]]]=v
    assert np.isfinite(counts).all() and np.min(counts)>=0 and np.min(total)>0
    return dict(zip(samples,(counts*(1e6/total)).T,strict=True))


def measure(p,y):
    return np.array([abs(np.mean(p-y)),abs(np.mean(y)),np.abs(p-y).mean(),np.abs(y).mean()])


def rec(unit,p,y):
    em,ez,gm,gz=measure(p,y)
    return dict(unit=unit,n_targets=len(y),observed_mean=np.mean(y),model_mean=np.mean(p),
                model_panel_error=em,zero_panel_error=ez,panel_gain=ez-em,
                model_gene_error=gm,zero_gene_error=gz,gene_gain=gz-gm)


def main():
    info=metadata()
    write(HERE/"sample_units.tsv",info)
    genes,target,ids,base,loaded,w,idx,a=design()
    cpm=read_counts(HERE/"counts.gz",ids)
    assert set(cpm)=={r["sample"] for r in info}
    pools={}
    for r in info:
        pools.setdefault(r["unit"],dict(contact=[],hypofunction=[]))[r["arm"]].append(cpm[r["sample"]])
    weights_out=[]
    def analyse(units,save_weights=False):
        p,y,details=[],[],[]
        for u in units:
            others=[np.mean(pools[v]["hypofunction"],axis=0) for v in units if u!=v]
            calib=np.mean(others,axis=0)
            forecasts=[]
            baseline_samples=[r["sample"] for r in info if r["unit"]==u and r["arm"]=="hypofunction"]
            for sample,x in zip(baseline_samples,pools[u]["hypofunction"],strict=True):
                q=weights(x,calib,w,idx,a)
                if save_weights:
                    weights_out.append(dict(sample=sample,unit=u,**dict(zip(CLASSES,q,strict=True))))
                pred=x*(1+loaded@q)/(1+base@q)
                forecasts.append(ranks(pred,target)-ranks(x,target))
            p.append(np.mean(forecasts,axis=0))
            d={k:np.mean([ranks(x,target) for x in group],axis=0) for k,group in pools[u].items()}
            y.append(d["contact"]-d["hypofunction"])
            details.append(rec(u,p[-1],y[-1]))
        return np.asarray(p),np.asarray(y),details

    units=sorted(pools)
    p,y,details=analyse(units,save_weights=True)
    write(HERE/"weights.tsv",weights_out)
    group=rec("group_contrast",p.mean(axis=0),y.mean(axis=0))
    avg={k:(np.mean([r[k] for r in details]) if k!="unit" else "mean_of_paired_scores") for k in details[0]}
    write(HERE/"scores.tsv",details+[avg,group])
    rows=[]
    for u,pred,obs in zip(units,p,y,strict=True):
        for gene,px,yx in zip(np.asarray(genes)[target],pred,obs,strict=True):
            rows.append(dict(unit=u,gene=gene,predicted=px,observed=yx))
    write(HERE/"gene_changes.tsv",rows)
    deletions=[]
    for drop in units:
        pp,yy,dd=analyse([u for u in units if u!=drop])
        rr=rec(drop,pp.mean(axis=0),yy.mean(axis=0))
        deletions.append(dict(omitted=drop,participant_panel_gain=np.mean([d["panel_gain"] for d in dd]),
                              participant_gene_gain=np.mean([d["gene_gain"] for d in dd]),
                              group_panel_gain=rr["panel_gain"],group_gene_gain=rr["gene_gain"]))
    write(HERE/"deletions.tsv",deletions)
    print(avg,group)


if __name__=="__main__":
    main()
