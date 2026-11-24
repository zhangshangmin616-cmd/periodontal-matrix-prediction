"""Evaluate atlas transfer and day-3 persistence on the same unfiltered gene sets."""
import csv
import gzip
import json
import math
from pathlib import Path
from statistics import fmean

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def load_contrast(condition):
    folder = WORK / 'stage8' / 'hotm'
    path = folder / f'deg_C0_{condition}.json'
    meta = json.loads(path.with_name(path.stem + '_meta.json').read_text(encoding='utf8'))
    assert f'Select Group 1:\nC0\nSelect Group 2:\n{condition}\n' in meta['label']
    data = json.loads(path.read_text(encoding='utf8'))
    assert data['header'][1:4] == ['gene_name', 'baseMean', 'log2FoldChange']
    result = {r[1]: r[3] for r in data['rows'] if r[3] is not None and math.isfinite(r[3])}
    assert len(result) == len(data['rows']) == meta['page']['recordsDisplay']
    return result


def write_table(name, rows):
    with (HERE / name).open('w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, rows[0].keys(), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    with gzip.open(WORK / 'stage8' / 'axial_log2_multipliers.tsv.gz', 'rt') as f:
        universe = list(csv.DictReader(f, delimiter='\t'))
    groups = {label: load_contrast(label) for label in ['E3', 'I3', 'E10', 'I10']}
    alleffects = []
    missing = []
    for item in universe:
        gene = item['gene']
        absent = [key for key, rows in groups.items() if gene not in rows]
        if absent:
            missing.append(dict(gene=gene, in_panel=item['in_panel'], unavailable=','.join(absent)))
            continue
        alleffects.append(dict(gene=gene, in_panel=item['in_panel'], atlas_E=float(item['E10']),
                               atlas_I=float(item['I10']), **{c: v[gene] for c,v in groups.items()}))
    assert len(alleffects) + len(missing) == len(universe) == 10806
    write_table('paired_effects.tsv', alleffects)
    write_table('unavailable.tsv', missing)
    scores = []
    for label in ['panel101', 'background10806']:
        subset = [r for r in alleffects if label == 'background10806' or r['in_panel'] == 'True']
        def get(name):
            return np.array([r[name] for r in subset], dtype=float)
        for arm, opposite in [('E', 'I'), ('I', 'E'), ('I_minus_E', None)]:
            if opposite:
                observed = get(arm+'10')
                atlas = get('atlas_'+arm)
                atlas_flipped = get('atlas_'+opposite)
                prior = get(arm+'3')
                prior_flipped = get(opposite+'3')
            else:
                observed = get('I10') - get('E10')
                atlas = get('atlas_I') - get('atlas_E')
                prior = get('I3') - get('E3')
                atlas_flipped = -atlas
                prior_flipped = -prior
            methods = {'zero': np.zeros_like(observed), 'matching_atlas': atlas, 'swapped_atlas': atlas_flipped,
                       'matching_day3': prior, 'swapped_day3': prior_flipped,
                       'average_day3': (prior+prior_flipped)/2}
            for name, pred in methods.items():
                error = np.abs(pred-observed)
                gain = np.abs(observed)-error
                close = np.abs(gain) < 1e-12
                scores.append(dict(geneset=label, outcome=arm, model=name, genes=len(subset),
                                   mae=fmean(error), bias=fmean(pred-observed), gain_vs_zero=fmean(gain),
                                   wins=int(np.sum((gain>0) & ~close)), ties=int(np.sum(close)),
                                   losses=int(np.sum((gain<0) & ~close))))
    write_table('scores.tsv', scores)
    print(len(alleffects), 'common background genes;', sum(r['in_panel']=='True' for r in alleffects),
          'common panel genes')
    for row in scores:
        if row['outcome'] in ['I','E'] and row['model'] in ['zero','matching_day3','matching_atlas']:
            print(row['geneset'],row['outcome'],row['model'],round(row['mae'],4))


if __name__ == '__main__':
    main()
