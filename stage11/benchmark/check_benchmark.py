"""Independently evaluate each direct reconstruction by scalar summation."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def table(path):
    with path.open(encoding='utf8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, default=WORK/'stage10/reference_counts.npz')
    parser.add_argument('--features', type=Path, default=WORK/'stage5/reference/features.tsv')
    parser.add_argument('--results', type=Path, default=HERE)
    args = parser.parse_args()
    material = table(args.features)
    expected = {(r['model'],r['comparison'],r['subset']):r for r in table(args.results/'gene_metrics.tsv')
                if r['layer']=='direct' and r['scale']=='log2_cpm1'}
    with np.load(args.counts, allow_pickle=False) as data:
        names = data['genes'].tolist()
        lookup = {g:i for i,g in enumerate(names)}
        bg = [lookup[r['gene']] for r in material]
        set_ids = {'panel': [idx for idx, row in zip(bg,material) if row['ecm_component']=='other'],
                   'background':bg, 'all_unique':list(range(len(names)))}
        cp_rows, cp_depth = data['CP'].tolist(), data['CP_total'].tolist()
        true_states = {}
        work = {}
        # Independently walk the raw integer vectors and accumulate the mixtures.
        for side, prefix in [('pressure','PP'),('tension','TP')]:
            early, late = prefix+'3D', prefix+'7D'
            early_depth, late_depth = data[early+'_total'].tolist(), data[late+'_total'].tolist()
            early_rows, late_rows = data[early].tolist(), data[late].tolist()
            total_cp, total_early, total_late = sum(cp_depth), sum(early_depth), sum(late_depth)
            for model in ('unchanged','state_only','composition_only','state_and_composition','measured_day3'):
                work[model,side] = []
            true_states[side] = []
            for c, e, l in zip(cp_rows,early_rows,late_rows):
                b = 1e6*sum(c)/total_cp
                true_states[side].append(math.log2((1+1e6*sum(l)/total_late)/(1+b)))
                state = c[-1]/total_cp + math.fsum(cp_depth[j]*e[j]/(total_cp*early_depth[j]) for j in range(5))
                composition = math.fsum(c[j]*early_depth[j]/(cp_depth[j]*total_early) for j in range(6))
                combined = sum(e[:-1])/total_early + c[-1]*early_depth[-1]/(total_early*cp_depth[-1])
                abundances = [b, state*1e6, composition*1e6, combined*1e6, 1e6*sum(e)/total_early]
                for model, a in zip(('unchanged','state_only','composition_only','state_and_composition','measured_day3'),abundances):
                    work[model,side].append(math.log2((1+a)/(1+b)))
        errors = []
        for (model,comparison,subset), row in expected.items():
            selected = set_ids[subset]
            if comparison in true_states:
                prediction, real = work[model,comparison], true_states[comparison]
                differences = [abs(prediction[i]-real[i]) for i in selected]
            else:
                direction = -1 if comparison.startswith('swapped') else 1
                differences = [abs(direction*(work[model,'pressure'][i]-work[model,'tension'][i])
                                   - (true_states['pressure'][i]-true_states['tension'][i])) for i in selected]
            errors.append(abs(math.fsum(differences)/len(selected)-float(row['mae'])))
        assert len(errors)==60 and max(errors)<1e-11
    print(json.dumps(dict(direct_MAE_values=len(errors),max_scalar_vs_vector_difference=max(errors))))


if __name__=='__main__':
    main()
