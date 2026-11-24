"""Evaluate forecasts against held-out, selectively recombined pseudo-tissues."""
import csv
import itertools
from pathlib import Path
from statistics import fmean

import numpy as np


HERE = Path(__file__).resolve().parent
SCENARIOS = ('unchanged', 'state_only', 'shares_only', 'both')
MODELS = (*SCENARIOS, 'other_side_state')
METRICS = ('gene_rank_mae', 'panel_rank_abs_error', 'gene_log2cpm1_mae')


def write(name, rows):
    with (HERE/name).open('w', newline='', encoding='utf8') as stream:
        w = csv.DictWriter(stream, rows[0].keys(), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def profile(counts, depths, folds):
    summed = counts[:, folds].sum(axis=1)
    depth = depths[:, folds].sum(axis=1)
    shares = depth / depth.sum(axis=1, keepdims=True)
    return summed / depth[:,:,None] * 1e6, shares


def aggregate(profiles, shares):
    return np.einsum('c,cg->g', shares, profiles)


def mixtures(profiles, q0, q1, side, wrong):
    baseline = profiles[0]
    response = profiles[side].copy()
    response[-1] = baseline[-1]
    opposite = profiles[wrong].copy()
    opposite[-1] = baseline[-1]
    return {'unchanged':aggregate(baseline,q0), 'state_only':aggregate(response,q0),
            'shares_only':aggregate(baseline,q1), 'both':aggregate(response,q1),
            'other_side_state':aggregate(opposite,q0)}


def ranks_for_panel(values, target):
    values = np.round(values,10)
    ordered = np.sort(values)
    subset = values[target]
    left = np.searchsorted(ordered,subset,side='left')
    right = np.searchsorted(ordered,subset,side='right')
    return 50 * (left+right-1) / (len(values)-1)


def evaluate():
    with np.load(HERE/'block_counts.npz', allow_pickle=False) as f:
        libs, target, counts, depths, cells = (f[k] for k in ('libraries','panel','counts','total_umi','cells'))
    assert list(libs)==['CP','PP3D','TP3D','PP7D','TP7D']
    assert counts.shape==(5,4,6,10806) and target.sum()==101 and cells.sum()==62490
    records=[]
    for folds in itertools.combinations(range(4),2):
        test=sorted(set(range(4))-set(folds))
        train_profiles, train_q=profile(counts,depths,list(folds))
        test_profiles, test_q=profile(counts,depths,test)
        for side,other in [(1,2),(2,1),(3,4),(4,3)]:
            starting_q=test_q[0]
            altered_q=starting_q * (train_q[side]/train_q[0])
            altered_q/=altered_q.sum()
            truth=mixtures(test_profiles,starting_q,test_q[side],side,other)
            m=mixtures(train_profiles,starting_q,altered_q,side,other)
            baseline=truth['unchanged']
            # Multiplying by one exactly for the no-change model preserves the zero-error null.
            pred={k:baseline*((1+v)/(1+m['unchanged'])) for k,v in m.items()}
            assert np.array_equal(pred['unchanged'],baseline)
            baseline_rank=ranks_for_panel(baseline,target)
            y_rank={k:ranks_for_panel(truth[k],target)-baseline_rank for k in SCENARIOS}
            p_rank={k:ranks_for_panel(v,target)-baseline_rank for k,v in pred.items()}
            for condition in SCENARIOS:
                y=y_rank[condition]
                for method in MODELS:
                    p=p_rank[method]
                    error=p-y
                    log_error=np.log2(1+pred[method][target])-np.log2(1+truth[condition][target])
                    records.append(dict(train_blocks=''.join(map(str,folds)),test_blocks=''.join(map(str,test)),
                                        library=libs[side],constructed=condition,forecast=method,
                                        observed_panel=y.mean(),predicted_panel=p.mean(),
                                        gene_rank_mae=np.mean(abs(error)),panel_rank_abs_error=abs(error.mean()),
                                        gene_log2cpm1_mae=np.mean(abs(log_error))))
    return records


def main():
    records=evaluate()
    write('scores_by_split.tsv',records)
    totals=[]
    for condition in SCENARIOS:
        for model in MODELS:
            selected=[r for r in records if r['forecast']==model and r['constructed']==condition]
            assert len(selected)==24
            totals.append(dict(constructed=condition,forecast=model,
                               **{key:fmean(r[key] for r in selected) for key in METRICS}))
    write('scores.tsv',totals)
    groups=[]
    for condition in SCENARIOS:
        for lib in ('PP3D','TP3D','PP7D','TP7D'):
            for method in MODELS:
                seq=[r for r in records if (r['constructed'],r['library'],r['forecast'])==(condition,lib,method)]
                groups.append(dict(constructed=condition,library=lib,forecast=method,
                                   **{key:fmean(r[key] for r in seq) for key in METRICS}))
    write('scores_by_library.tsv',groups)
    for row in totals:
        print(row['constructed'],row['forecast'],*[round(row[c],4) for c in METRICS])


if __name__=='__main__':
    main()
