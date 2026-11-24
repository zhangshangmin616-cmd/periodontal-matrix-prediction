"""Independent scalar checks of target scores and four random draws in every condition."""
import csv
from pathlib import Path
from statistics import fmean

import numpy as np


HERE = Path(__file__).resolve().parent
DRAW_IDS = [0,111,4999,9999]


def scalar_errors(y, p, sel):
    panel = [[], []]
    genes = [[], []]
    for yi, pi in zip(y, p, strict=True):
        outcomes, prediction = yi[sel], pi[sel]
        panel[0].append(abs(fmean(outcomes)))
        panel[1].append(abs(fmean(outcomes)-fmean(prediction)))
        genes[0].append(fmean(map(abs,outcomes)))
        genes[1].append(fmean(abs(a-b) for a,b in zip(outcomes, prediction, strict=True)))
    return {name: tuple(map(fmean,lists)) for name,lists in [('panel',panel),('gene',genes)]}


def main():
    with np.load(HERE / 'rank_changes.npz', allow_pickle=False) as archive:
        y,p,cp,components,units = (archive[k] for k in ('observed','fitted','cp_counts','component','units'))
    expected = np.load(HERE/'matched_gains.npz',allow_pickle=False)
    with (HERE/'comparison.tsv').open(encoding='utf8') as stream:
        outputs = {(r['cohort'],r['metric']): r for r in csv.DictReader(stream, delimiter='\t')}
    assert np.allclose(p.mean(axis=1),0,atol=1e-11) and np.allclose(y.mean(axis=1),0,atol=1e-11)
    ordered = np.sort(cp)
    rank = (np.searchsorted(ordered,cp,side='left') + np.searchsorted(ordered,cp,side='right')-1) / 2
    strata = np.clip(np.floor(rank*10/(len(cp)-1)),0,9)
    target = np.flatnonzero(components=='other')
    assert (len(target),sum(components=='outside_ecm'))==(101,10589)
    sampled = {i: [] for i in DRAW_IDS}
    rng = np.random.Generator(np.random.PCG64(20261124))
    for b in range(10):
        pool = np.flatnonzero((components=='outside_ecm') & (strata==b))
        n = int(sum(strata[target]==b))
        if n:
            for i in range(10000):
                draw = rng.choice(pool, n, replace=False)
                if i in sampled:
                    sampled[i].extend(draw)
    for key, rr in [('Kim_week3',np.char.startswith(units,'Kim_')),
                    ('Zhao_day7',units=='Zhao_day7'),('Zhao_day28',units=='Zhao_day28')]:
        for metric,(z,m) in scalar_errors(y[rr],p[rr],target).items():
            record=outputs[key,metric]
            np.testing.assert_allclose([z,m],[float(record['zero_mae']),float(record['model_mae'])], rtol=0,atol=2e-12)
        for i,sel in sampled.items():
            for metric,(z,m) in scalar_errors(y[rr],p[rr],sel).items():
                np.testing.assert_allclose(z-m,expected[key+'_'+metric][i],rtol=0,atol=2e-12)
    print('Verified all six fixed-panel losses and 24 matched-set gains by independent scalar arithmetic.')


if __name__=='__main__':
    main()
