"""Verify target ranks directly by counting lower and tied background features."""
from collections import defaultdict
import csv
import gzip
import hashlib
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def read(path):
    with path.open(encoding="utf8") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def percentiles(x,ix):
    return np.array([(np.sum(x < x[j])+0.5*(np.sum(x==x[j])-1))*100/(len(x)-1) for j in ix])


def main():
    assert hashlib.sha256((HERE/"genes.fpkm.tsv.gz").read_bytes()).hexdigest()=="cdbfcc7fcfb7184e2ddc51a6a6ace236a1013e1bca695afc1d797d529a53c96c"
    with np.load(WORK/"stage8/mixture_inputs.npz") as raw:
        all_genes=list(raw["symbols"])
        panel=set(raw["symbols"][raw["idx"]])
    ids={r["gene"]:r["entrez_id"] for r in read(WORK/"stage5/zhao/feature_map.tsv")}
    runs=[("SRR6291210","GSM2856840","normal PDLSCs"),("SRR6291212","GSM2856842","stretch-induced PDLSCs")]
    samples={r["run_accession"]:r for r in read(HERE/"samples.tsv")}
    counts=[[],[]]
    for i,(run,gsm,title) in enumerate(runs):
        assert (samples[run]["sample_alias"],samples[run]["sample_title"])==(gsm,title)
        with gzip.open(HERE/"counts"/(run+".tsv.gz"),"rt") as f:
            tab={r["Geneid"]:float(r[run+"_count"]) for r in csv.DictReader(f,delimiter="\t")}
        counts[i]=[tab[ids[g]] for g in all_genes]
    groups=defaultdict(list)
    with gzip.open(HERE/"genes.fpkm.tsv.gz","rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            groups[r["gene_short_name"]].append((float(r["PDLZC"]),float(r["J12h"])))
    genef=[g for g in all_genes if len(groups[g])==1]
    fpkm=np.array([groups[g][0] for g in genef]).T
    with gzip.open(WORK/"stage8/axial_log2_multipliers.tsv.gz","rt") as f:
        tissue={r["gene"]:2**float(r["E3"]) for r in csv.DictReader(f,delimiter="\t")}
    with np.load(WORK/"stage10/reference_counts.npz") as ref:
        k=list(ref["classes"]).index("fibroblast_like")
        fibro=dict(zip(ref["genes"],(1+ref["TP3D"][:,k]*1e6/ref["TP3D_total"][k])/
                                    (1+ref["CP"][:,k]*1e6/ref["CP_total"][k]),strict=True))
    summ={(r["input"],r["method"]):r for r in read(HERE/"scores.tsv")}
    targets={(r["input"],r["method"],r["gene"]):r for r in read(HERE/"targets.tsv")}
    assert len(targets)==3*(101+80) and len(summ)==6
    for kind,genes,(before,after) in [("NCBI_counts",all_genes,np.asarray(counts)),("author_FPKM",genef,fpkm)]:
        ix=[i for i,g in enumerate(genes) if g in panel]
        start=percentiles(before,ix)
        obs=percentiles(after,ix)-start
        for name,ratio in [("zero",np.ones(len(genes))),("fibroblast",[fibro[g] for g in genes]),("tissue",[tissue[g] for g in genes])]:
            pred=percentiles(before*ratio,ix)-start
            effects=[targets[kind,name,genes[j]] for j in ix]
            for vec,field in [(obs,"observed"),(pred,"predicted")]:
                np.testing.assert_allclose(vec,[float(r[field]) for r in effects],atol=5e-12,rtol=0)
            d=pred-obs
            r=summ[kind,name]
            assert int(r["targets"])==len(ix) and int(r["background"])==len(genes)
            values=[abs(d.mean()),abs(obs.mean()),abs(d).mean(),abs(obs).mean()]
            np.testing.assert_allclose(values,[float(r[k]) for k in ("panel_abs_error","zero_panel_abs_error",
                                                                     "gene_mae","zero_gene_mae")],atol=5e-12,rtol=0)
    print("Independently verified 543 target rows and six pairs of error scores from archived GEO/NCBI inputs.")


if __name__=="__main__":
    main()
