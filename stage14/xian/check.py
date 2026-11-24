"""Independently check every reported cross-study error and single-unit influence."""
import csv
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def tied_ranks(x):
    _,ids,n=np.unique(x,return_inverse=True,return_counts=True)
    return ((np.cumsum(n)-(n+1)/2)[ids])/(len(x)-1)*100


def read(path):
    with path.open(encoding='utf8') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def main():
    with np.load(HERE/'xian_counts.npz') as a,np.load(WORK/'stage13/external/paired_counts.npz') as b:
        x=a['counts']/a['total']*1e6
        h=b['counts']/b['total']*1e6
        target=a['target']
        arms=a['arm']
        runs=a['run']
    scores=read(HERE/'scores.tsv')
    checked=0
    deleted=[]
    for scale,fn in [('rank',tied_ranks),('log2',lambda col:np.log2(1+col))]:
        xv=np.stack([fn(c) for c in x.T],axis=1)
        hv=np.stack([fn(c) for c in h.T],axis=1)
        hd=hv[:,5:]-hv[:,:5]
        for freq in ['0.1','0.5']:
            dest='Xian_'+freq+'Hz'
            hmean=hd.mean(1)
            idx=np.flatnonzero(arms==freq)
            base=np.flatnonzero(arms=='0')
            xmean=xv[:,idx].mean(1)-xv[:,base].mean(1)
            for row in scores:
                if row['scale']!=scale or dest not in (row['source'],row['destination']):
                    continue
                mask=target if row['scope']=='target' else np.ones(len(target),bool)
                if row['source']=='HongKong':
                    pred,obs=hmean[mask],xmean[mask]
                else:
                    pred,obs=xmean[mask],hmean[mask]
                p={'quarter':.25*pred,'source_gene':pred,'source_mean':np.full_like(pred,pred.mean()),
                   'zero':np.zeros_like(pred)}[row['method']]
                assert abs(np.linalg.norm(p-obs,ord=1)/len(p)-float(row['gene_mae']))<1e-10
                if row['panel_error']:
                    assert abs(abs(p.mean()-obs.mean())-float(row['panel_error']))<1e-10
                checked+=1
            for channel,j in ([('source',k) for k in range(5)]+[('control',k) for k in base]+[('loaded',k) for k in idx]):
                p=np.delete(hd,j,1).mean(1) if channel=='source' else hmean
                y=(xv[:,idx[idx!=j] if channel=='loaded' else idx].mean(1)
                   - xv[:,base[base!=j] if channel=='control' else base].mean(1))
                err=(p-y)[target]
                deleted.append(dict(scale=scale,hz=freq,omitted=channel,
                                    unit=str(j+1) if channel=='source' else runs[j],
                                    gene_gain=abs(y[target]).mean()-abs(err).mean(),
                                    panel_gain=abs(y[target].mean())-abs(err.mean())))
    assert checked==len(scores)==64
    with (HERE/'deletion.tsv').open('w',encoding='utf8',newline='') as f:
        out=csv.DictWriter(f,list(deleted[0]),delimiter='\t',lineterminator='\n')
        out.writeheader();out.writerows(deleted)
    print('Independent reconstruction passed:',checked,'summary rows;',len(deleted),'deletion diagnostics.')


if __name__=='__main__':
    main()
