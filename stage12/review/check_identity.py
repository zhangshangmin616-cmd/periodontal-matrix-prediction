"""Independently reconstruct the score-preserving controls with scalar sums."""
import csv
from math import fsum
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rows(file):
    with (HERE/file).open(encoding='utf8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    with np.load(HERE.parents[1]/'stage11/specificity/rank_changes.npz', allow_pickle=False) as d:
        chosen = np.flatnonzero(d['component'] == 'other')
        observed = d['observed'][:, chosen]
        fitted = d['fitted'][:, chosen]
        counts = d['cp_counts']
        labels = list(d['units'])
    n = len(chosen)
    # Exact midrank decile assignment with integer counts and no fractional rounding.
    def decile(c):
        below = int(np.sum(counts < c))
        equal = int(np.sum(counts == c))
        return min(9, (5*(2*below+equal-1))//(len(counts)-1))
    strata = [decile(c) for c in counts[chosen]]
    ids = {'Seoul_week3':[i for i,key in enumerate(labels) if key.startswith('Kim_')],
           'Leuven_day7':[labels.index('Zhao_day7')], 'Leuven_day28':[labels.index('Zhao_day28')]}
    for r in rows('gene_identity.tsv'):
        idx=ids[r['study']]
        groups = strata if r['control']=='matched_atlas_decile' else [0]*n
        choices = [[k for k in range(n) if groups[j]==groups[k]] for j in range(n)]
        expected = fsum(fsum(abs(float(fitted[i,k]-observed[i,j])) for k in choices[j])/len(choices[j])
                        for j in range(n) for i in idx)/(n*len(idx))
        assert abs(expected-float(r['assignment_expected_mae'])) < 1e-11
    for r in rows('baselines.tsv'):
        idx=ids[r['study']]
        residuals=[]
        for i in idx:
            mean=fsum(map(float,fitted[i]))/n
            def estimate(j):
                return {'unchanged':0.,'state_only':fitted[i,j],'same_panel_mean':mean}[r['per_unit_model']]
            residuals.append([float(observed[i,j]-estimate(j)) for j in range(n)])
        panel = fsum(abs(fsum(v))/n for v in residuals)/len(idx)
        gene = fsum(abs(x) for v in residuals for x in v)/(n*len(idx))
        assert abs(panel-float(r['panel_mae'])) < 1e-11
        assert abs(gene-float(r['gene_mae'])) < 1e-11
    print('Checked six assignment expectations and nine matched panel/gene calculations.')


if __name__=='__main__':
    main()
