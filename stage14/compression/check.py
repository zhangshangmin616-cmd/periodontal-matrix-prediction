"""Independent recomputation, plus a bounded rank-scale sensitivity."""
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rank_by_group(values):
    _,inv,n = np.unique(values,return_inverse=True,return_counts=True)
    mid = n.cumsum() - (n+1)/2
    return mid[inv]*100/(len(values)-1)


def main():
    with np.load(HERE/'counts.npz') as a:
        rank = np.stack([rank_by_group(c) for c in a['count']])
        log = np.log2(1+a['count']/a['library_total'][:,None]*1e6)
        mask = a['target']
        study,arm = a['accession'],a['arm']
    with (HERE/'scores.tsv').open(encoding='utf8') as stream:
        rows = list(csv.DictReader(stream,delimiter='\t'))
    assert len(rows)==120
    bounded=[]
    for r in rows:
        s,d,rule = r['source'],r['destination'],r['rule']
        x = rank if r['scale']=='rank' else log
        if r['scope']=='panel101':
            x=x[:,mask]
        def group_mean(site, condition):
            return np.sum(x[(study==site)&(arm==condition)],axis=0)/int(((study==site)&(arm==condition)).sum())
        dest_base=group_mean(d,'control')
        source=group_mean(s,'force')-group_mean(s,'control')
        if rule=='zero':
            p=np.zeros_like(source)
        elif 'uniform' in rule:
            p=np.repeat(source.sum()/len(source),len(source))
        else:
            p=source.copy()
        if rule.endswith('quarter'):
            p/=4
        y = group_mean(d,'force')-dest_base
        delta=y-p
        gene_err=sum(map(abs,delta))/len(delta)
        panel_err=abs(sum(delta)/len(delta))
        assert np.allclose([gene_err,panel_err],[float(r['gene_mae']),float(r['panel_abs_error'])],atol=1e-10,rtol=0)
        if r['scope']=='panel101' and r['scale']=='rank' and rule.startswith('source_gene_'):
            final = dest_base+p
            projected = np.clip(final,0,100)-dest_base
            gain = np.mean(np.abs(y))-np.mean(np.abs(y-projected))
            bounded.append(dict(source=s,destination=d,rule=rule,entries_outside=((final<0)|(final>100)).sum(),
                                gene_mae=np.mean(abs(y-projected)),gene_gain_vs_zero=gain,
                                panel_abs_error=abs(np.mean(y-projected))))
    with (HERE/'bounded.tsv').open('w',newline='',encoding='utf8') as st:
        writer=csv.DictWriter(st,fieldnames=list(bounded[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(bounded)
    print('Checked',len(rows),'cross-study scores and',len(bounded),'bounded comparators.')


if __name__=='__main__':
    main()
