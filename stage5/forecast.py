"""Freeze atlas-informed predictions using only untreated tissue measurements."""
import csv
import gzip
import hashlib
import itertools
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parent
REF, PAIR, OUT=(ROOT/x for x in ["reference","paired","results"])
COUNTS=PAIR/"GSE165988_raw_counts_GRCh38.p13_NCBI.tsv.gz"
SHEET=PAIR/"GSE165988_series_matrix.txt.gz"
CLASSES=["fibroblast_like","mural","immune","endothelial","glial","residual"]
SETTINGS=["primary","day3","pressure_only","tension_only","atlas_weights","all_profiles","state_only"]

def put(rows,name):
    pd.DataFrame(rows).to_csv(OUT/name,index=False,sep="\t",float_format="%.15g",lineterminator="\n")

def filehash(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()

def meta():
    with gzip.open(SHEET,"rt",encoding="utf-8") as stream:
        rows=list(csv.reader(stream,delimiter="\t"))
    def row(key): return next(r[1:] for r in rows if r and r[0]==key)
    gsm,title=row("!Sample_geo_accession"),row("!Sample_title")
    treatments=next(r[1:] for r in rows if r and r[0]=="!Sample_characteristics_ch1" and r[1].startswith("treatment:"))
    result=[]
    for s,t,condition in zip(gsm,title,treatments,strict=True):
        match=re.fullmatch(r"Sample (\d+)-(exp|con)\d*\s+\[.*\]",t)
        assert match and (match[2]=="con")==(condition=="treatment: untreated")
        result.append(dict(sample=s,patient=int(match[1]),arm="control" if match[2]=="con" else "force"))
    df=pd.DataFrame(result).set_index("sample")
    assert len(df)==29 and df.patient.nunique()==11 and df.index.is_unique
    return df

def load_bulk(samples,genes):
    # The forecast entry point selects only control columns on import.
    raw=pd.read_csv(COUNTS,usecols=["GeneID",*samples],sep="\t",dtype={"GeneID":str}).set_index("GeneID")
    assert raw.index.is_unique and raw.min().min()>=0
    totals=raw.sum()
    registry=pd.read_csv(PAIR/"hgnc_ids.tsv.gz",sep="\t",dtype=str)
    links=registry[["entrez_id","symbol"]].dropna().drop_duplicates()
    links=links[~links.entrez_id.duplicated(keep=False)&~links.symbol.duplicated(keep=False)]
    mapping=links.set_index("entrez_id").symbol
    raw.index=raw.index.map(mapping)
    raw=raw[raw.index.notna()&~raw.index.duplicated(keep=False)].loc[genes,samples]
    return raw.div(totals).mul(1.e6)

def score(values,idx):
    ranks=rankdata(values,method="average")-1
    return float(ranks[idx].mean()*100/(len(ranks)-1))

class Reference:
    def __init__(self):
        features=pd.read_csv(REF/"features.tsv",sep="\t").set_index("gene")
        prof=pd.read_csv(REF/"profiles.tsv.gz",sep="\t",index_col="gene")
        assert prof.index.equals(features.index)
        self.genes=features.index
        self.idx=np.flatnonzero(features.ecm_component.eq("other"))
        assert len(self.idx)==101 and len(self.genes)==10806
        groups=pd.read_csv(REF/"groups.tsv",sep="\t").set_index(["library","cell_class"])
        self.w={};self.expr={}
        for lib in ("CP","PP3D","PP7D","TP3D","TP7D"):
            self.w[lib]=groups.xs(lib).loc[CLASSES,"rna_weight"].to_numpy()
            self.expr[lib]=prof[[f"{lib}/{c}" for c in CLASSES]].to_numpy()*1e6
            assert np.isclose(self.w[lib].sum(),1)
        self.w0=self.w["CP"];self.base=self.expr["CP"]
        panel=json.loads((PAIR/"stage4_panel.json").read_text())["panels"]
        exclude={g for names in panel.values() for g in names}
        exclude.update(["CTHRC1","HBB","HBA1","HBA2","HBD","HBE1","HBG1","HBG2","HBQ1","HBZ","HBM"])
        eligible=features.ecm_component.eq("outside_ecm").to_numpy() & ~self.genes.isin(exclude) & ~self.genes.str.startswith(("MT-","RPS","RPL"))
        out=[]
        for i,c in enumerate(CLASSES):
            contrast=np.log2((self.base[:,i]+1)/(np.max(np.delete(self.base,i,axis=1),axis=1)+1))
            molecules=self.base[:,i]*groups.loc[("CP",c),"total_umi"]/1e6
            cand=np.flatnonzero(eligible&(molecules>=9.999999)&(contrast>0))
            best=cand[np.argsort(-contrast[cand],kind="stable")[:30]]
            for g in best:
                out.append(dict(gene=self.genes[g],cell_class=c,log2_specificity=contrast[g],atlas_umi=round(molecules[g])))
        self.anchors=pd.DataFrame(out)
        assert self.anchors.gene.is_unique and set(self.anchors.cell_class)==set(CLASSES)
        self.indices=self.genes.get_indexer(self.anchors.gene)
        b=self.base[self.indices]
        self.template=(b+1)/((b@self.w0)[:,None]+1)
        # Precompute constrained least-squares solutions for each of the 63 simplex faces.
        A=self.template.T@self.template/len(b)+.05*np.eye(len(CLASSES))
        self.faces=[]
        for k in range(1,len(CLASSES)+1):
            for subset in itertools.combinations(range(len(CLASSES)),k):
                ix=np.array(subset)
                mat=np.block([[A[np.ix_(ix,ix)],np.ones((k,1))],[np.ones((1,k)),np.zeros((1,1))]])
                self.faces.append((ix,np.linalg.inv(mat)))

    def coefficients(self,baseline,calibration):
        y=np.clip((baseline[self.indices]+1)/(calibration[self.indices]+1),.25,4)
        v=self.template.T@y/len(y)+.05*self.w0
        best=None;loss=np.inf
        for ix,inv in self.faces:
            sol=(inv@np.r_[v[ix],1])[:-1]
            if (sol< -1.e-9).any(): continue
            q=np.zeros(len(self.w0)); q[ix]=np.maximum(0,sol)
            q/=q.sum()
            ev=np.mean((self.template@q-y)**2)+.05*np.sum((q-self.w0)**2)
            if ev<loss:
                best,loss=q,ev
        assert best is not None
        rmse=np.mean((self.template@best-y)**2)**.5
        return best,rmse

    def predict(self,bulk,q,setting):
        if setting=="atlas_weights": q=self.w0
        sides={"day3":["PP3D","TP3D"],"pressure_only":["PP7D"],"tension_only":["TP7D"]}.get(setting,["PP7D","TP7D"])
        before=self.base@q
        comp=np.zeros(len(self.genes));state=np.zeros_like(comp)
        for lib in sides:
            postq=q*self.w[lib]/self.w0
            postq=postq/postq.sum()
            if setting=="state_only": postq=q
            comp+=self.base@postq
            within=self.expr[lib]@postq
            if setting!="all_profiles":
                within += (self.base[:,-1]-self.expr[lib][:,-1])*postq[-1]
            state+=within
        start=score(bulk,self.idx)
        return tuple(score(bulk*((v/len(sides)+1)/(before+1)), self.idx)-start for v in [comp,state])

def main():
    OUT.mkdir(exist_ok=True)
    ref=Reference(); ref.anchors.to_csv(OUT/"anchors.tsv",sep="\t",index=False)
    info=meta()
    cs=info.query("arm=='control'")
    bulk=load_bulk(cs.index.tolist(),ref.genes)
    patients=sorted(cs.patient.unique())
    means={i:bulk.loc[:,cs.index[cs.patient.eq(i)]].mean(axis=1).to_numpy() for i in patients}
    pred,inner,weights,tooth=[],[],[],[]
    for held in patients:
        train=[i for i in patients if i!=held]
        for target in [held,*train]:
            calibration=np.mean([means[j] for j in train if j!=target],axis=0)
            if target==held: settings=SETTINGS
            else: settings=["primary"]
            aggregate={s:[] for s in settings}
            for s in cs.index[cs.patient.eq(target)]:
                values=bulk[s].to_numpy()
                q,fit=ref.coefficients(values,calibration)
                if held==target:
                    weights.extend(dict(patient=held,sample=s,cell_class=c,weight=v,anchor_rmse=fit) for c,v in zip(CLASSES,q))
                for setting in settings:
                    dp=ref.predict(values,q,setting)
                    aggregate[setting].append(dp)
                    if target==held:
                        tooth.append(dict(patient=held,sample=s,setting=setting,baseline=score(values,ref.idx),composition=dp[0],state=dp[1]))
            for setting,values in aggregate.items():
                c,m=np.mean(values,axis=0)
                rec=dict(patient=target,setting=setting,composition=c,state=m)
                if held==target: pred.append(rec)
                else: inner.append(dict(held_out=held,**rec))
    put(pred,"predictions.tsv");put(inner,"inner_predictions.tsv");put(weights,"weights.tsv");put(tooth,"tooth_predictions.tsv")
    paths=[ROOT/"design.txt",ROOT/"forecast.py",REF/"profiles.tsv.gz",REF/"groups.tsv",REF/"features.tsv",*PAIR.iterdir(),
           *[OUT/f for f in ["anchors.tsv","predictions.tsv","inner_predictions.tsv","weights.tsv","tooth_predictions.tsv"]]]
    (OUT/"lock.json").write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
        sha256={p.relative_to(ROOT).as_posix():filehash(p) for p in paths}),indent=2))
    print("Frozen",len(pred),"patient/setting rows,",len(inner),"inner-fold predictions,",len(ref.anchors),"anchors.")

if __name__=="__main__":main()
