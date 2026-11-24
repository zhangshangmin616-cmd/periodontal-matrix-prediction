"""Independently check the exported effects and all 20 forecast scores."""
import csv
import gzip
import math
from pathlib import Path
from statistics import fmean
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    with (HERE/'specimen_metadata.tsv').open(encoding='utf8') as f:
        meta=list(csv.DictReader(f,delimiter='\t'))
    with np.load(HERE/'specimen_cpm.npz',allow_pickle=False) as a:
        values=a['cpms']
        genes=a['genes']
        target=a['target']
    with gzip.open(HERE/'gene_contrasts.tsv.gz','rt',encoding='utf8') as f:
        contrasts=list(csv.DictReader(f,delimiter='\t'))
    assert [r['gene'] for r in contrasts]==list(genes)
    sets=[[i for i,r in enumerate(meta) if r['batch']==b and r['group']==g]
          for b,g in [('4','day7'),('4','control'),('3','day28'),('3','control')]]
    for j,r in enumerate(contrasts):
        group_values=[[values[i,j] for i in sample] for sample in sets]
        for key,estimate in [('mean_log',lambda v:fmean(math.log2(1+e) for e in v)),
                             ('log_mean',lambda v:math.log2(1+fmean(v)))]:
            y7,c7,y28,c28=map(estimate,group_values)
            assert abs(float(r[key+'_day7'])-(y7-c7))<1e-11 and abs(float(r[key+'_day28'])-(y28-c28))<1e-11
    with (HERE/'scores.tsv').open(encoding='utf8') as f:
        scores=list(csv.DictReader(f,delimiter='\t'))
    assert len(scores)==20
    for row in scores:
        group=contrasts if row['panel']=='10806' else [r for j,r in enumerate(contrasts) if target[j]]
        terms=[abs(float(row['coefficient'])*float(r[row['estimand']+'_day7'])-float(r[row['estimand']+'_day28'])) for r in group]
        assert abs(fmean(terms)-float(row['mae']))<1e-14
    print('All 43,224 gene-effect values and 20 summary errors verified independently against 63 specimens.')


if __name__=='__main__':
    main()
