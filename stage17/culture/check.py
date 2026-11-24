"""Independently reconstruct the FPKM rank errors from the primary worksheet and MINiML."""
from collections import Counter
import csv
import gzip
import hashlib
from pathlib import Path
import tarfile
from xml.etree import ElementTree
from openpyxl import load_workbook
import numpy as np

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def read(file):
    with (HERE/file).open(encoding="utf8") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def ranks(v):
    ordered=np.sort(v)
    mid=ordered.searchsorted(v,side="left")+ordered.searchsorted(v,side="right")-1
    return mid*(50/(len(v)-1))


def close(a,b):
    np.testing.assert_allclose(a,b,rtol=0,atol=5e-12)


def main():
    assert hashlib.sha256((HERE/"fpkm.xlsx").read_bytes()).hexdigest()=="96b167b8b00a228c193ed47b2e9a57090fcfd0106da1c132c22254ab9f4a6e89"
    with tarfile.open(HERE/"family.xml.tgz") as bundle:
        doc=ElementTree.parse(bundle.extractfile(bundle.getmembers()[0]))
    labels={s.findtext("{*}Accession"):(s.findtext("{*}Description").strip().splitlines()[-1],
                                        s.findtext("{*}Title")) for s in doc.getroot().findall("{*}Sample")}
    assert len(labels)==18
    for r in read("sample_map.tsv"):
        assert labels[r["gsm"]]==(r["worksheet_column"],r["title"])
    with np.load(WORK/"stage8/mixture_inputs.npz") as raw:
        symbols=raw["symbols"]
        panel=set(symbols[raw["idx"]])
    source=load_workbook(HERE/"fpkm.xlsx",read_only=True,data_only=True)
    cells=iter(source.active.values)
    headings=next(cells)
    samples=[(gsm,x) for gsm,(x,_) in labels.items() if x.startswith("T-") and "-SH-" not in x]
    cols=[headings.index(x) for gsm,x in samples]
    arm=np.array([int(x.split("-")[1]) for gsm,x in samples])
    matrix={}
    counts=Counter()
    for r in cells:
        counts[r[0]]+=1
        matrix[r[0]]=np.array([float(r[j]) for j in cols])
    source.close()
    genes=sorted(g for g in symbols if counts[g]==1)
    ix=np.array([i for i,g in enumerate(genes) if g in panel])
    assert len(genes)==10554 and len(ix)==101
    vals=np.array([matrix[g] for g in genes]).T
    with gzip.open(WORK/"stage8/axial_log2_multipliers.tsv.gz","rt") as f:
        weights={r["gene"]:float(r["E3"]) for r in csv.DictReader(f,delimiter="\t")}
    with np.load(WORK/"stage10/reference_counts.npz") as z:
        which={g:i for i,g in enumerate(z["genes"])}
        i=[which[g] for g in genes]
        k=list(z["classes"]).index("fibroblast_like")
        cp=z["CP"][i,k]/z["CP_total"][k]*1e6
        tp=z["TP3D"][i,k]/z["TP3D_total"][k]*1e6
    factors={"zero":np.ones(len(genes)),"tissue":np.exp2([weights[g] for g in genes]),
             "fibroblast":(1+tp)/(1+cp)}
    effects={(int(r["hours"]),r["method"],r["gene"]):r for r in read("targets.tsv")}
    assert len(effects)==606
    rows=read("scores.tsv")+read("deletions.tsv")
    assert len(rows)==42
    for r in rows:
        keep=np.array([gsm!=r.get("omitted_gsm") for gsm,_ in samples])
        control=vals[(arm==0)&keep]
        load=vals[(arm==int(r["hours"]))&keep]
        base=np.mean([ranks(x)[ix] for x in control],axis=0)
        y=np.mean([ranks(x)[ix] for x in load],axis=0)-base
        p=np.mean([ranks(x*factors[r["method"]])[ix] for x in control],axis=0)-base
        d=p-y
        actual=(y.mean(),p.mean(),abs(d.mean()),abs(y.mean()),abs(y.mean())-abs(d.mean()),
                abs(d).mean(),abs(y).mean(),abs(y).mean()-abs(d).mean())
        fields=("observed_mean","predicted_mean","panel_error","zero_panel_error","panel_gain",
                "gene_mae","zero_gene_mae","gene_gain")
        close(actual,[float(r[k]) for k in fields])
        if "omitted_gsm" not in r:
            pairs=[effects[int(r["hours"]),r["method"],genes[i]] for i in ix]
            close(y,[float(q["observed"]) for q in pairs])
            close(p,[float(q["prediction"]) for q in pairs])
    print("Verified 18 sample labels against GEO, 606 forecast rows, 6 full scores and 36 library-deletion scores.")


if __name__=="__main__":
    main()
