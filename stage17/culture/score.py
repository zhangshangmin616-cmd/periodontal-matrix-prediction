"""Score frozen atlas shifts on the GEO depositor's FPKM worksheet."""
from collections import Counter, defaultdict
import csv
import gzip
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import re
import tarfile
from xml.etree import ElementTree

import numpy as np
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def write(name, data):
    with (HERE/name).open("w",encoding="utf8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=data[0],delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def load():
    with gzip.open(HERE/"geo_matrix.txt.gz","rt",encoding="utf8") as handle:
        rows=list(csv.reader(handle,delimiter="\t"))
    geo=dict((r[0],r[1:]) for r in rows if r and r[0] in ("!Sample_title","!Sample_geo_accession"))
    with tarfile.open(HERE/"family.xml.tgz","r:gz") as bundle:
        member,=bundle.getmembers()
        root=ElementTree.parse(bundle.extractfile(member)).getroot()
    labels={s.findtext("{*}Accession"):s.findtext("{*}Description").strip().splitlines()[-1]
            for s in root.findall("{*}Sample")}
    assert len(labels)==18
    records=json.loads((HERE/"ena_runs.json").read_bytes())
    enaruns=defaultdict(list)
    for row in records:
        gsm=re.search(r"\bGSM\d+\b",row["experiment_title"])[0]
        enaruns[gsm].append(row["run_accession"])
    with (HERE/"fpkm.xlsx").open("rb") as data:
        book=load_workbook(data,read_only=True,data_only=True)
        assert len(book.worksheets)==1
        source=iter(book.active.values)
        header=next(source)
        assert header[0]=="GeneID" and len(header)==19
        mapping=[]
        for gsm,title in zip(geo["!Sample_geo_accession"],geo["!Sample_title"],strict=True):
            match=re.fullmatch(r"sh-(Control|PIEZO1), tension (0|24|48) h, rep ([123])",title)
            sh,time,rep=match.groups()
            col=("T-" if sh=="Control" else "T-SH-")+time+"-"+rep
            assert labels[gsm]==col
            mapping.append(dict(gsm=gsm,title=title,worksheet_column=col,order=header.index(col),
                                hours=int(time),short_hairpin=sh,replicate_label=rep,
                                runs="|".join(sorted(enaruns[gsm]))))
        assert len({r["order"] for r in mapping})==18
        selected=[r for r in mapping if r["short_hairpin"]=="Control"]
        ids=[r["order"] for r in selected]
        features=Counter()
        with np.load(WORK/"stage11/specificity/rank_changes.npz") as z:
            background=z["gene"]
            original_target=set(background[z["component"]=="other"])
        desired=set(background)
        values={}
        for row in source:
            symbol=row[0]
            features[symbol]+=1
            if symbol in desired:
                values[symbol]=[float(row[i]) for i in ids]
        book.close()
    assert len(desired)==10806 and len(original_target)==101
    unique=np.array([features[g]==1 for g in background])
    missing=set(background[~unique])
    if missing & original_target:
        raise ValueError("A frozen target is not uniquely measured: "+str(missing & original_target))
    genes=background[unique]
    target=np.array([g in original_target for g in genes])
    x=np.asarray([values[g] for g in genes])
    assert x.shape==(len(genes),9) and np.isfinite(x).all() and x.min()>=0
    write("sample_map.tsv",mapping)
    write("gene_coverage.tsv",[dict(gene=g,worksheet_rows=features[g],target=int(g in original_target),
                                    used=int(features[g]==1)) for g in background])
    return x,genes,target,selected,sum(features.values())


def main():
    spec=spec_from_file_location("frozen",WORK/"stage13/external/score.py")
    frozen=module_from_spec(spec)
    spec.loader.exec_module(frozen)
    x,genes,target,samples,source_features=load()
    times=np.array([s["hours"] for s in samples])
    ranks=np.stack([frozen.ranks(c)[target] for c in x.T],axis=1)
    metrics,detail,deleted=[],[],[]
    for name, shift in frozen.atlas(genes).items():
        pcols={j:frozen.ranks(x[:,j]*np.exp2(shift))[target]-ranks[:,j] for j in np.flatnonzero(times==0)}
        for time in (24,48):
            c0=np.flatnonzero(times==0)
            e0=np.flatnonzero(times==time)
            for omit in [None,*c0,*e0]:
                c=c0 if omit is None else c0[c0!=omit]
                e=e0 if omit is None else e0[e0!=omit]
                y=ranks[:,e].mean(axis=1)-ranks[:,c].mean(axis=1)
                p=np.mean([pcols[j] for j in c],axis=0)
                d=p-y
                record=dict(hours=time,method=name,gene_background=len(genes),targets=len(y),
                            control_libraries=len(c),loaded_libraries=len(e),observed_mean=y.mean(),
                            predicted_mean=p.mean(),panel_error=abs(d.mean()),zero_panel_error=abs(y.mean()),
                            panel_gain=abs(y.mean())-abs(d.mean()),gene_mae=np.abs(d).mean(),
                            zero_gene_mae=np.abs(y).mean(),gene_gain=np.abs(y).mean()-np.abs(d).mean())
                if omit is None:
                    metrics.append(record)
                    detail.extend(dict(hours=time,method=name,gene=g,observed=yi,prediction=pi)
                                  for g,yi,pi in zip(genes[target],y,p,strict=True))
                else:
                    deleted.append(dict(omitted_gsm=samples[omit]["gsm"],**record))
    write("scores.tsv",metrics)
    write("deletions.tsv",deleted)
    write("targets.tsv",detail)
    print(source_features,"worksheet rows;",len(genes),"uniquely mapped frozen genes;",target.sum(),"targets")


if __name__=="__main__":
    main()
