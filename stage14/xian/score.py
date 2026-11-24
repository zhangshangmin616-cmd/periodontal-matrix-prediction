"""Transport observed RNA-seq responses between the Hong Kong and Xi'an tensile cultures."""
import csv
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def ranks(x):
    ordered=np.sort(x)
    return (ordered.searchsorted(x)+ordered.searchsorted(x,side='right')-1)*50/(len(x)-1)


def save(name,rows):
    with (HERE/name).open('w',encoding='utf8',newline='') as f:
        out=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        out.writeheader()
        out.writerows(rows)


def main():
    with np.load(HERE/'xian_counts.npz') as dat:
        genes,tar,arm=dat['genes'],dat['target'],dat['arm']
        x=dat['counts']*1e6/dat['total']
    with np.load(WORK/'stage13/external/paired_counts.npz') as hk:
        assert np.array_equal(genes,hk['genes']) and np.array_equal(tar,hk['target'])
        h=hk['counts']*1e6/hk['total']
    assert arm.size==8 and h.shape==(10805,10) and tar.sum()==101
    results, panel_data=[],[]
    for scale,transform in [('log2',lambda v:np.log2(1+v)),('rank',lambda v:np.column_stack([ranks(a) for a in v.T]))]:
        hh,xx=transform(h),transform(x)
        mh=hh[:,:5].mean(1)
        mx=xx[:,arm=='0'].mean(1)
        dest_hk=hh[:,5:].mean(1)-mh
        for freq in ['0.1','0.5']:
            xa='Xian_'+freq+'Hz'
            dest_x=xx[:,arm==freq].mean(1)-mx
            for idx in np.where(tar)[0]:
                panel_data.append(dict(scale=scale,freq=freq,gene=genes[idx],hongkong=dest_hk[idx],xian=dest_x[idx]))
            for source,dest,p,y in [('HongKong',xa,dest_hk,dest_x),(xa,'HongKong',dest_x,dest_hk)]:
                for title,mask in [('target',tar),('background',np.ones(len(tar),bool))]:
                    g,observed = p[mask],y[mask]
                    for method,forecast in [('zero',np.zeros_like(g)),('source_mean',np.full_like(g,g.mean())),
                                            ('source_gene',g),('quarter',.25*g)]:
                        diff=forecast-observed
                        nominal=(scale=='rank' and title=='background')
                        results.append(dict(scale=scale,source=source,destination=dest,genes=len(g),scope=title,
                                            method=method,gene_mae=abs(diff).mean(),
                                            gene_gain=abs(observed).mean()-abs(diff).mean(),
                                            panel_error='' if nominal else abs(diff.mean()),
                                            panel_gain='' if nominal else abs(observed.mean())-abs(diff.mean())))
    save('scores.tsv',results)
    save('target_changes.tsv',panel_data)


if __name__=='__main__':
    main()
