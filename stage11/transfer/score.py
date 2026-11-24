"""Score fixed Leuven carry-over forecasts with two prespecified effect definitions."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import numpy as np


HERE=Path(__file__).resolve().parent
COEFFICIENTS={'no_change':0.,'atlas_D3_D7':0.23238716799216566,'quarter':.25,'half':.5,'full':1.}
PROTOCOL='ce4624f2783c63d3bce00cb4dee1360aa7966616a4702903d829ee399cc5f273'


def main():
    args=argparse.ArgumentParser()
    args.add_argument('--outdir',type=Path,default=HERE)
    out=args.parse_args().outdir
    out.mkdir(parents=True,exist_ok=True)
    assert hashlib.sha256((HERE/'protocol.txt').read_bytes()).hexdigest()==PROTOCOL
    with (HERE/'specimen_metadata.tsv').open(encoding='utf8') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    with np.load(HERE/'specimen_cpm.npz',allow_pickle=False) as a:
        cpms,genes,selected=a['cpms'],a['genes'],a['target']
        assert list(a['runs'])==[r['run'] for r in rows]
    assert cpms.shape==(63,10806) and (cpms>=0).all() and np.isfinite(cpms).all()
    assert selected.sum()==101 and len(set(r['biosample'] for r in rows))==63
    group_defs=[('4','control',5),('4','day7',19),('3','control',22),('3','day28',17)]
    locs=[]
    for batch,arm,n in group_defs:
        ix=[i for i,r in enumerate(rows) if (r['batch'],r['group'])==(batch,arm)]
        assert len(ix)==n
        locs.append(ix)
    assert len({i for ix in locs for i in ix})==63
    contrasts={}
    summary=[]
    for name,aggregate in [('mean_log',lambda v:np.log2(v+1).mean(axis=0)),
                           ('log_mean',lambda v:np.log2(v.mean(axis=0)+1))]:
        pre7,post7,pre28,post28=[aggregate(cpms[ix]) for ix in locs]
        x,y=post7-pre7,post28-pre28
        contrasts[name+'_day7']=x
        contrasts[name+'_day28']=y
        for lab,keep in [('101',selected),('10806',np.ones(len(y),dtype=bool))]:
            ey,ex=y[keep],x[keep]
            zero=np.abs(ey).mean()
            quarter=np.abs(.25*ex-ey).mean()
            for method,coefficient in COEFFICIENTS.items():
                err=coefficient*ex-ey
                mae=np.abs(err).mean()
                summary.append(dict(estimand=name,panel=lab,n_genes=len(ey),method=method,coefficient=coefficient,
                                    mae=mae,gain_over_zero=zero-mae,gain_over_quarter=quarter-mae,bias=err.mean()))
    with (out/'scores.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,summary[0].keys(),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(summary)
    with gzip.GzipFile(filename=str(out/'gene_contrasts.tsv.gz'),mode='wb',mtime=0) as buffer:
        with io.TextIOWrapper(buffer,encoding='utf8',newline='') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n')
            w.writerow(['gene','in_101',*contrasts])
            for vals in zip(genes,selected,*contrasts.values(),strict=True):
                w.writerow(vals)
    for row in summary:
        print(row['estimand'],row['panel'],row['method'],round(float(row['mae']),6))


if __name__=='__main__':
    main()
