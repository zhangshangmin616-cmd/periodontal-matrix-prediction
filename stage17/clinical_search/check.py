"""Verify the two endpoint scores from the raw tissue matrix and atlas molecule totals."""
import csv
import gzip
import hashlib
import re
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def tsv(path):
    with path.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def check_identity():
    a,b=[tsv(HERE/(p+".tsv")) for p in ("PRJNA698902","PRJNA782351")]
    assert len(a)==29 and len(b)==22
    for key in ["sample_accession","run_accession","experiment_accession","sample_alias","fastq_md5","read_count"]:
        assert not {r[key] for r in a}.intersection(r[key] for r in b)
    with gzip.open(WORK/"stage5/paired/GSE165988_series_matrix.txt.gz","rt",encoding="latin1") as f:
        previous=set(re.findall(r'TN[0-9]{4}[A-Z][0-9]{4}',next(line for line in f if line.startswith("!Sample_title"))))
    new={r["library"] for r in tsv(HERE/"sample_units.tsv")}
    assert len(new)==22 and len(previous)==29 and not new & previous
    path=HERE/"fastq_prefixes.tsv"
    if path.is_file():
        hashes=tsv(path)
        assert len({r["sha256_read_sequences"] for r in hashes})==len(hashes)==51
        assert {r["run"] for r in hashes}=={r["run_accession"] for r in a+b}


def tied_percentiles(x):
    _,inverse,counts=np.unique(x,return_inverse=True,return_counts=True)
    cum=np.cumsum(counts)
    return ((cum-1)-(counts-1)*.5)[inverse]*100/(len(x)-1)


def close(a,b):
    np.testing.assert_allclose(a,b,atol=5e-11,rtol=0)


def main():
    check_identity()
    metadata=tsv(HERE/"sample_units.tsv")
    feature=tsv(WORK/"stage5/zhao/feature_map.tsv")
    genes=[r["gene"] for r in feature]
    target=np.array([r["ecm_component"]=="other" for r in feature])
    ids=[r["entrez_id"] for r in feature]
    with gzip.open(HERE/"counts.gz","rt",encoding="utf8") as f:
        samples=f.readline().rstrip().split("\t")[1:]
        whole={parts[0]:np.fromiter(parts[1:],dtype=float) for line in f if (parts:=line.rstrip().split("\t"))}
    matrix=np.array([whole[g] for g in ids]).T
    cpm=matrix*1e6/np.sum(list(whole.values()),axis=0)[:,None]
    assert len(genes)==10806 and target.sum()==101
    sample_loc={s:i for i,s in enumerate(samples)}
    weights={r["sample"]:r for r in tsv(HERE/"weights.tsv")}
    with np.load(WORK/"stage10/reference_counts.npz",allow_pickle=False) as archive:
        positions={g:i for i,g in enumerate(archive["genes"])}
        indices=[positions[g] for g in genes]
        def cell(lib):
            return archive[lib][indices,:]*1e6/archive[lib+"_total"]
        cp,pp,tp=map(cell,("CP","PP7D","TP7D"))
        classes=list(archive["classes"])
        resid=classes.index("residual")
        w=archive["CP_total"]/archive["CP_total"].sum()
        pp[:,resid]=tp[:,resid]=cp[:,resid]
    anchors=np.array([genes.index(r["gene"]) for r in tsv(WORK/"stage5/results/anchors.tsv")])
    template=(cp[anchors]+1)/(1+cp[anchors]@w)[:,None]
    units={r["unit"] for r in metadata}
    data={r["sample"]:r for r in metadata}
    observed={s:tied_percentiles(matrix[i])[target] for i,s in enumerate(samples)}
    predicted={}
    for s,row in weights.items():
        unit=data[s]["unit"]
        other=[]
        for u in units-{unit}:
            other.append(np.mean([cpm[sample_loc[r["sample"]]] for r in metadata
                                  if r["unit"]==u and r["arm"]=="hypofunction"],axis=0))
        calib=np.mean(other,axis=0)
        x=cpm[sample_loc[s]]
        q=np.array([float(row[c]) for c in classes])
        close(q.sum(),1)
        assert (q>=0).all()
        response=np.clip((1+x[anchors])/(1+calib[anchors]),.25,4)
        gradient=(template.T@(template@q-response))/len(anchors)+.05*(q-w)
        active=q>1e-9
        close(gradient[active],np.repeat(gradient[active].mean(),active.sum()))
        assert np.all(gradient[~active]>=gradient[active].mean()-1e-10)
        pre=1+cp@q
        after=((1+pp@q)/pre*x+(1+tp@q)/pre*x)*.5
        predicted[s]=tied_percentiles(after)[target]-observed[s]
    saved={(r["unit"],r["gene"]):r for r in tsv(HERE/"gene_changes.tsv")}
    assert len(saved)==1010
    by_unit={r["unit"]:r for r in tsv(HERE/"scores.tsv")}
    for unit in sorted(units):
        arms={k:[r["sample"] for r in metadata if r["unit"]==unit and r["arm"]==k]
              for k in ("hypofunction","contact")}
        y=np.mean([observed[s] for s in arms["contact"]],axis=0)-np.mean([observed[s] for s in arms["hypofunction"]],axis=0)
        p=np.mean([predicted[s] for s in arms["hypofunction"]],axis=0)
        for gene,px,yx in zip(np.array(genes)[target],p,y,strict=True):
            close([float(saved[unit,gene][k]) for k in ("predicted","observed")],[px,yx])
        close([float(by_unit[unit][k]) for k in ("model_panel_error","zero_panel_error","model_gene_error","zero_gene_error")],
              [abs((p-y).mean()),abs(y.mean()),np.abs(p-y).mean(),np.abs(y).mean()])
    print("Verified distinct libraries, RNA-share optima, 1,010 predictions and observations, and paired loss endpoints.")


if __name__=="__main__":
    main()
