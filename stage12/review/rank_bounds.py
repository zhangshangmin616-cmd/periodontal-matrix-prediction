"""Check the score-preserving uniform shift after enforcing marginal rank bounds."""
import csv
from pathlib import Path

import numpy as np
from identity import strata

HERE = Path(__file__).resolve().parent


def main():
    with np.load(HERE.parents[1]/'stage11/specificity/rank_changes.npz',allow_pickle=False) as d:
        idx=np.flatnonzero(d['component']=='other')
        genes=d['gene'][idx]
        labs=d['units'].astype(str)
        y=d['observed'][:,idx]
        p=d['fitted'][:,idx]
        bins=strata(d['cp_counts'],idx)
    with (HERE/'control_ranks.tsv').open(encoding='utf8') as f:
        source={(r['unit'],r['gene']):float(r['untreated_rank']) for r in csv.DictReader(f,delimiter='\t')}
    rank=np.array([[source[lab,g] for g in genes] for lab in labs])
    low,high=-rank,100-rank
    assert np.all((p >= low-1e-12)&(p <= high+1e-12))
    assert np.all((y >= low-1e-12)&(y <= high+1e-12))
    bounded=np.zeros_like(p)
    for i in range(len(labs)):
        target=p[i].mean()
        lo,hi=-100.,100.
        for _ in range(70):
            mid=(lo+hi)/2
            if np.clip(mid,low[i],high[i]).mean() < target:
                lo=mid
            else:
                hi=mid
        bounded[i]=np.clip((lo+hi)/2,low[i],high[i])
    np.testing.assert_allclose(bounded.mean(axis=1),p.mean(axis=1),atol=1e-12)
    assert np.all((bounded>=low)&(bounded<=high))
    # Sensitivity restricted to exchanges feasible in *every* clinical unit.
    # Build disjoint pairs in untreated-rank proximity; leave unmatched genes
    # fixed. An exact expectation follows from one independent fair swap/pair.
    candidates=[]
    for i,first in enumerate(genes):
        for j in range(i+1,len(genes)):
            if bins[i] != bins[j]:
                continue
            first_ok=np.all((low[:,i]<=p[:,j]+1e-12)&(high[:,i]>=p[:,j]-1e-12))
            second_ok=np.all((low[:,j]<=p[:,i]+1e-12)&(high[:,j]>=p[:,i]-1e-12))
            if first_ok and second_ok:
                distance=float(np.mean((rank[:,i]-rank[:,j])**2))
                candidates.append((distance,first,genes[j],i,j))
    used=set()
    partner=np.arange(len(genes))
    for _,_,_,i,j in sorted(candidates):
        if i not in used and j not in used:
            used.update([i,j]);partner[i]=j;partner[j]=i
    swapped=p[:,partner]
    np.testing.assert_allclose(swapped.mean(axis=1),p.mean(axis=1),atol=1e-12)
    assert np.all((swapped>=low-1e-12)&(swapped<=high+1e-12))
    rows=[]
    for cohort,mask in [('Seoul_week3',np.char.startswith(labs,'Kim_')),
                        ('Leuven_day7',labs=='Zhao_day7'),('Leuven_day28',labs=='Zhao_day28')]:
        scalar=p[mask].mean(axis=1)[:,None]
        rows.append([cohort,len(y[mask]),len(genes),np.mean(np.abs(p[mask]-y[mask])),
                    np.mean(np.abs(scalar-y[mask])),np.mean(np.abs(bounded[mask]-y[mask])),
                    int(np.sum((scalar<low[mask]-1e-12)|(scalar>high[mask]+1e-12))),
                    len(used)//2,0.5*np.mean(np.abs(swapped[mask]-y[mask])+np.abs(p[mask]-y[mask]))])
    with (HERE/'bounded_baselines.tsv').open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['cohort','n_units','genes','model_mae','uniform_mae','bounded_uniform_mae',
                    'uniform_unit_genes_outside_bounds','eligible_disjoint_pairs','expected_paired_swap_mae'])
        w.writerows(rows)
    print(rows)


if __name__=='__main__':
    main()
