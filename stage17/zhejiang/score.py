"""Score the fixed gene panel in the single deposited Zhejiang stretch-control comparison."""
import csv
import gzip
import hashlib
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def rows(path):
    with path.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def save(file, records):
    with (HERE/file).open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=records[0],delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(records)


def counts(genes):
    mapping={r["gene"]:r["entrez_id"] for r in rows(WORK/"stage5/zhao/feature_map.tsv")}
    manifest={r["run"]:r for r in rows(HERE/"count_manifest.tsv")}
    state={r["sample_title"]:r["run_accession"] for r in rows(HERE/"samples.tsv")}
    x=[]
    for sample in ["normal PDLSCs","stretch-induced PDLSCs"]:
        run=state[sample]
        with gzip.open(HERE/"counts"/(run+".tsv.gz"),"rb") as f:
            raw=f.read()
        assert hashlib.md5(raw).hexdigest()==manifest[run]["md5"]
        it=csv.reader(raw.decode().splitlines(),delimiter="\t")
        assert next(it)==["Geneid",run+"_count"]
        values={k:float(v) for k,v in it}
        x.append(np.array([values[mapping[g]] for g in genes])*1e6/sum(values.values()))
    return np.array(x)


def fpkm(genes):
    lookup={}
    cnt=Counter()
    with gzip.open(HERE/"genes.fpkm.tsv.gz","rt",encoding="utf8") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            g=r["gene_short_name"]
            cnt[g]+=1
            lookup[g]=(float(r["PDLZC"]),float(r["J12h"]))
    retain=np.array([cnt[g]==1 for g in genes])
    return np.array([lookup[g] for g in genes[retain]]).T,retain


def main():
    spec=spec_from_file_location("frozen",WORK/"stage13/external/score.py")
    prior=module_from_spec(spec)
    spec.loader.exec_module(prior)
    with np.load(WORK/"stage11/specificity/rank_changes.npz") as z:
        genes,target=z["gene"],z["component"]=="other"
    assert len(genes)==10806 and sum(target)==101
    expr,keeprows=fpkm(genes)
    matrices=[("NCBI_counts",counts(genes),np.ones(len(genes),dtype=bool)),("author_FPKM",expr,keeprows)]
    scores,detail=[],[]
    for source,x,mask in matrices:
        ix=target[mask]
        syms=genes[mask]
        if source=="author_FPKM":
            missing=genes[target & ~mask].tolist()
            print("FPKM missing or multiple loci for:",missing)
        for model,shift in prior.atlas(syms).items():
            base=prior.ranks(x[0])
            y=(prior.ranks(x[1])-base)[ix]
            p=(prior.ranks(x[0]*np.exp2(shift))-base)[ix]
            d=p-y
            scores.append(dict(input=source,method=model,targets=int(ix.sum()),background=len(syms),
                               control_GSM="GSM2856840",loaded_GSM="GSM2856842",
                               observed_mean=y.mean(),predicted_mean=p.mean(),panel_abs_error=abs(d.mean()),
                               zero_panel_abs_error=abs(y.mean()),panel_gain=abs(y.mean())-abs(d.mean()),
                               gene_mae=abs(d).mean(),zero_gene_mae=abs(y).mean(),gene_gain=abs(y).mean()-abs(d).mean()))
            detail.extend(dict(input=source,method=model,gene=g,observed=yi,predicted=pi)
                          for g,yi,pi in zip(syms[ix],y,p,strict=True))
    save("scores.tsv",scores)
    save("targets.tsv",detail)
    save("fpkm_coverage.tsv",[dict(gene=g,unique_symbol_in_fpkm=int(k)) for g,k in zip(genes[target],keeprows[target],strict=True)])


if __name__=="__main__":
    main()
