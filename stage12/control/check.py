"""Independently replay all split-cell responses and check barcode membership."""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import itertools
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read_tsv(name):
    path=HERE/name
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf8',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def rank_all(x):
    values=np.round(x,10)
    _, inverse, multiplicity=np.unique(values,return_inverse=True,return_counts=True)
    midpoint=multiplicity.cumsum()-(multiplicity+1)/2
    return (100/(len(values)-1))*midpoint[inverse]


def mix(q, profile):
    result=np.zeros(profile.shape[1])
    for index,weight in enumerate(q):
        result+=weight*profile[index]
    return result


def profiles(data,den,selected):
    c=data[:,selected[0]].astype('float64')+data[:,selected[1]]
    d=den[:,selected[0]].astype('float64')+den[:,selected[1]]
    return 1e6*(c/d[:,:,None]),d/d.sum(axis=1)[:,None]


def verify_members(lib,groups,n):
    assigned=read_tsv('heldout_assignments.tsv.gz')
    assert len({(r['library'],r['barcode']) for r in assigned})==len(assigned)==62490
    cells=defaultdict(list)
    for row in assigned:
        cells[row['library'],row['cell_group']].append(row)
    for i,l in enumerate(lib):
        for j,g in enumerate(groups):
            group=cells[l,g]
            def key(row):
                return hashlib.sha256(('PDL-holdout-v1|'+l+'|'+row['barcode']).encode('utf8')).digest()
            for jrank,r in enumerate(sorted(group,key=key)):
                assert int(r['block'])==jrank%4
            counts=Counter(int(r['block']) for r in group)
            assert list(n[i,:,j])==[counts[b] for b in range(4)]


def main():
    with np.load(HERE/'block_counts.npz',allow_pickle=False) as data:
        counts,den,libs,groups,panel,n=(data[k] for k in ('counts','total_umi','libraries','groups','panel','cells'))
    verify_members(libs,groups,n)
    saved={}
    for r in read_tsv('scores_by_split.tsv'):
        key=tuple(r[k] for k in ('train_blocks','library','constructed','forecast'))
        assert key not in saved
        saved[key]=[float(r[k]) for k in ('gene_rank_mae','panel_rank_abs_error','gene_log2cpm1_mae')]
    error=0.
    for perm in itertools.combinations(range(4),2):
        train,qt=profiles(counts,den,perm)
        test,qo=profiles(counts,den,sorted(set(range(4))-set(perm)))
        for side in range(1,5):
            other=(side+1) if side%2 else (side-1)
            base=mix(qo[0],test[0])
            fitbase=mix(qo[0],train[0])
            truth={}
            pred={}
            for setting in ('unchanged','state_only','shares_only','both','other_side_state'):
                has_shares=setting in ('shares_only','both')
                has_state=setting in ('state_only','both','other_side_state')
                weights=qo[0].copy()
                if has_shares:
                    weights *= qt[side]/qt[0]
                    weights /= weights.sum()
                trainmix=train[0].copy()
                if has_state:
                    sel=other if setting=='other_side_state' else side
                    trainmix[:-1]=train[sel,:-1]
                pred[setting]=base*((1+mix(weights,trainmix))/(1+fitbase))
                if setting!='other_side_state':
                    testmix=test[0].copy()
                    if has_state:
                        testmix[:-1]=test[side,:-1]
                    truth[setting]=mix(qo[side] if has_shares else qo[0],testmix)
            real_ranks={k:rank_all(v)[panel] for k,v in truth.items()}
            model_ranks={k:rank_all(v)[panel] for k,v in pred.items()}
            for scenario,ys in truth.items():
                for forecast,ps in pred.items():
                    delta=model_ranks[forecast]-real_ranks[scenario]
                    output=[abs(delta).mean(),abs(delta.mean()),abs(np.log2((ps[panel]+1)/(ys[panel]+1))).mean()]
                    key=(''.join(map(str,perm)),libs[side],scenario,forecast)
                    difference=max(abs(np.asarray(saved.pop(key))-output))
                    error=max(error,difference)
    assert not saved and error<1e-8,error
    totals=read_tsv('scores.tsv')
    splits=read_tsv('scores_by_split.tsv')
    for t in totals:
        part=[r for r in splits if (r['constructed'],r['forecast'])==(t['constructed'],t['forecast'])]
        assert len(part)==24
        for metric in ['gene_rank_mae','panel_rank_abs_error','gene_log2cpm1_mae']:
            np.testing.assert_allclose(float(t[metric]),np.mean([float(r[metric]) for r in part]),rtol=1e-12)
    print('Cross-checked 62,490 assignments, 480 model/scenario comparisons and all summaries. Max numeric difference:',error)


if __name__=='__main__':
    main()
