"""Compare measured culture effects with atlas predictions; bound unreported errors."""
import csv
import gzip
import itertools
import math
from pathlib import Path
from statistics import fmean

import numpy as np


HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def save(name,rows):
    with (HERE/name).open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,rows[0].keys(),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    with (HERE/'culture_readouts.tsv').open(encoding='utf8') as f:
        measured=list(csv.DictReader(f,delimiter='\t'))
    with gzip.open(WORK/'stage8/axial_log2_multipliers.tsv.gz','rt') as f:
        mix={r['gene']:float(r['E3']) for r in csv.DictReader(f,delimiter='\t') if r['in_panel']=='True'}
    with np.load(WORK/'stage10/reference_counts.npz',allow_pickle=False) as f:
        ix={g:i for i,g in enumerate(f['genes'])}
        k=list(f['classes']).index('fibroblast_like')
        tp=f['TP3D'][:,k]*(1.e6/f['TP3D_total'][k])
        cp=f['CP'][:,k]*(1.e6/f['CP_total'][k])
        fb=np.log2((1+tp)/(1+cp))
    assert set(mix)=={r['gene'] for r in measured} and len(mix)==len(measured)==101
    predictions=[]
    for row in measured:
        g=row['gene']
        y=float(row['reported_log2']) if row['reported_log2'] else None
        raw_delta=None
        concordance=None
        if row['control1']:
            controls=[int(row['control'+str(j)]) for j in range(1,4)]
            tensile=[int(row['tension'+str(j)]) for j in range(1,4)]
            raw_delta=sum(tensile)-sum(controls)
            concordance=int(raw_delta*y>0)
        predictions.append(dict(gene=g,observed=y,zero=0.,fibroblast=float(fb[ix[g]]),
                                mixture=mix[g],count_sum_difference=raw_delta,raw_direction_concordant=concordance))
    save('predictions.tsv',predictions)
    scores=[]
    seen=[r for r in predictions if r['observed'] is not None]
    unseen=[r for r in predictions if r['observed'] is None]
    for method,baseline in [('fibroblast','zero'),('mixture','zero'),('fibroblast','mixture')]:
        gains=[abs(r[baseline]-r['observed'])-abs(r[method]-r['observed']) for r in seen]
        width=sum(abs(r[method]-r[baseline]) for r in unseen)
        scores.append(dict(forecast=method,comparator=baseline,reported=len(seen),not_reported=len(unseen),
                           selected_forecast_mae=fmean(abs(r[method]-r['observed']) for r in seen),
                           selected_comparator_mae=fmean(abs(r[baseline]-r['observed']) for r in seen),
                           selected_mean_gain=fmean(gains),selected_forecast_closer=sum(x>0 for x in gains),
                           selected_ties=sum(x==0 for x in gains),
                           selected_forecast_sign_correct=sum(r[method]*r['observed']>0 for r in seen),
                           min_full_panel_gain=(sum(gains)-width)/len(predictions),
                           max_full_panel_gain=(sum(gains)+width)/len(predictions)))
        # Both sharp limits can be reached by assigning each unknown outcome exactly
        # its forecast or comparator, independently of the other missing genes.
        for target,limit in [(method,scores[-1]['max_full_panel_gain']),(baseline,scores[-1]['min_full_panel_gain'])]:
            trial=gains+[abs(r[baseline]-r[target])-abs(r[method]-r[target]) for r in unseen]
            assert math.isclose(fmean(trial),limit,abs_tol=1e-12)
    save('scores.tsv',scores)
    print(scores)
    print('counts agree',sum(r['raw_direction_concordant'] or 0 for r in seen),
          'of',sum(r['raw_direction_concordant'] is not None for r in seen))


if __name__=='__main__':
    main()
