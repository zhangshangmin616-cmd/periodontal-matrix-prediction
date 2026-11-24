"""Switch cell-group expression states in a fixed clinical RNA mixture model."""
from itertools import permutations
from math import comb
from pathlib import Path
import csv
import json

import numpy as np


HERE = Path(__file__).resolve().parent


def table(name, rows):
    with (HERE/name).open('w', newline='', encoding='utf8') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def ranks(x, target):
    """Percent midrank with exact ties among all genes."""
    ordered = np.sort(x)
    left = np.searchsorted(ordered, x[target], side='left')
    right = np.searchsorted(ordered, x[target], side='right')
    return (left + right - 1) * (50. / (len(x)-1))


def predict(d):
    base, pressure, tension = (d[k] for k in ('base','pressure','tension'))
    ngroups = len(d['classes'])-1
    masks = 2**ngroups
    sums = np.zeros((masks+1,len(d['units']),len(d['target'])))
    counts = np.zeros(len(d['units']),dtype=int)
    c, p, t = d['source_shares']
    np.testing.assert_allclose(d['q'].sum(axis=1),1.,atol=1e-10,rtol=0)
    named = np.arange(ngroups)
    assignments = {u:i for i,u in enumerate(d['units'])}
    max_mixture_error = 0.
    for unit,q,x in zip(d['unit_per_control'],d['q'],d['control_cpm'],strict=True):
        row = assignments[unit]
        counts[row] += 1
        before = base @ q
        old_rank = ranks(x, d['target'])
        delta = ((pressure[:,:ngroups]+tension[:,:ngroups])*.5 - base[:,:ngroups])*q[:ngroups]
        for mask in range(1,masks):
            active = (mask & (1 << named)) != 0
            after = before + delta[:,active].sum(axis=1)
            # Independent reconstruction replaces profiles on both sides and then mixes.
            kept = np.append(active,False)
            direct = ((np.where(kept,pressure,base) @ q)+(np.where(kept,tension,base) @ q))*.5
            max_mixture_error = max(max_mixture_error,float(np.max(abs(after-direct))))
            np.testing.assert_allclose(after,direct,atol=3e-10,rtol=2e-13)
            rank_change = ranks(x*((after+1)/(before+1)), d['target'])-old_rank
            sums[mask,row] += rank_change
        # Re-use the side-specific change in captured RNA share as comparator.
        wp=q*p/c
        wt=q*t/c
        after = ((base @ (wp/wp.sum()))+(base @ (wt/wt.sum())))*.5
        sums[masks,row] += ranks(x*((after+1)/(before+1)),d['target'])-old_rank
    return sums / counts[None,:,None], max_mixture_error


def names(groups, mask):
    k=len(groups)
    if mask==2**k:
        return 'composition'
    if mask==0:
        return 'zero'
    if mask==2**k-1:
        return 'full'
    ids=[j for j in range(k) if (mask & (1<<j))!=0]
    if len(ids)==k-1:
        return 'omit_'+groups[next(j for j in range(k) if j not in ids)]
    if len(ids)==1:
        return 'only_'+groups[ids[0]]
    return '+'.join(groups[j] for j in ids)


def main():
    with np.load(HERE/'inputs.npz', allow_pickle=False) as archive:
        d = {k:archive[k] for k in archive.files}
    p, mixture_error = predict(d)
    y,units,groups=d['observed'],d['units'],d['classes'][:-1].tolist()
    assert y.shape==(13,101) and len(d['genes'])==10806 and d['classes'][-1]=='residual'
    assert p.shape==(33,13,101) and not np.any(p[0])
    max_full=float(np.max(abs(p[-2]-d['historical_full'])))
    max_comp=float(np.max(abs(p[-1]-d['historical_composition'])))
    np.testing.assert_allclose(p[-2],d['historical_full'],atol=1e-9,rtol=0)
    np.testing.assert_allclose(p[-1],d['historical_composition'],atol=1e-9,rtol=0)
    metrics={'panel':abs((p-y[None,:,:]).mean(axis=2)), 'gene':abs(p-y[None,:,:]).mean(axis=2)}
    cohorts={'Kim_week3':np.flatnonzero(np.char.startswith(units,'Kim_')),
             'Zhao_day7':np.flatnonzero(units=='Zhao_day7'), 'Zhao_day28':np.flatnonzero(units=='Zhao_day28')}
    summary=[];unit_rows=[];phi_rows=[]
    for cohort, ix in cohorts.items():
        for mask in range(len(p)):
            model=names(groups,mask)
            for metric,absolute in metrics.items():
                e,z,full=absolute[mask,ix], absolute[0,ix], absolute[-2,ix]
                summary.append(dict(cohort=cohort,model=model,mask=mask,metric=metric,n_units=len(ix),
                    mean_observed=y[ix].mean(),mean_predicted=p[mask,ix].mean(),mae=e.mean(),zero_mae=z.mean(),
                    gain_over_zero=(z-e).mean(),mae_minus_full=(e-full).mean(),units_with_positive_gain=int(np.sum(z-e>1e-10))))
                if model.startswith(('only_','omit_')) or model in ('full','zero','composition'):
                    for unit,ee,zz in zip(ix,e,z,strict=True):
                        unit_rows.append(dict(unit=units[unit],model=model,metric=metric,
                            observed=y[unit].mean(),predicted=p[mask,unit].mean(),mae=ee,zero_mae=zz,gain_over_zero=zz-ee))
        for metric,absolute in metrics.items():
            gains=(absolute[0,ix][None,:]-absolute[:-1,ix]).mean(axis=1)
            phi=np.zeros(len(groups))
            for order in permutations(range(len(groups))):
                previous=0
                for g in order:
                    new=previous | (1<<g)
                    phi[g]+=(gains[new]-gains[previous])/120.
                    previous=new
            np.testing.assert_allclose(phi.sum(),gains[-1],atol=1e-10,rtol=0)
            for j,group in enumerate(groups):
                marg=[(mask,gains[mask|(1<<j)]-gains[mask]) for mask in range(len(gains)) if not mask & (1<<j)]
                check=sum(value/(len(groups)*comb(len(groups)-1,mask.bit_count())) for mask,value in marg)
                np.testing.assert_allclose(phi[j],check,atol=1e-10,rtol=0)
                phi_rows.append(dict(cohort=cohort,metric=metric,group=group,mean_ordered_gain=phi[j],
                    smallest_marginal_gain=min(v for _,v in marg),largest_marginal_gain=max(v for _,v in marg),
                    full_model_gain=gains[-1]))
    table('all_subsets.tsv',summary)
    table('ablation_summary.tsv',[r for r in summary if r['model'] in ('full','zero','composition') or r['model'].startswith(('omit_','only_'))])
    table('unit_sensitivity.tsv',unit_rows)
    table('order_averaged_contributions.tsv',phi_rows)
    (HERE/'audit.json').write_text(json.dumps(dict(controls=len(d['q']),units=len(units),genes=len(d['genes']),targets=y.shape[1],
        max_mixture_discrepancy=mixture_error,max_full_vs_stage11=max_full,max_composition_vs_stage5=max_comp,
        combinations=32,attribution_orders=120),indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    main()
