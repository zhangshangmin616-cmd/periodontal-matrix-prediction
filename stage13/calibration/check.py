"""Independent scalar verification of the primary calibration and the observed contrast."""
import csv
from itertools import combinations
from pathlib import Path
from statistics import mean, median


here = Path(__file__).parent
def read(name):
    with (here / name).open(encoding='utf8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    rows = read('specimens.tsv')
    scores = {r['run']: float(r['score']) for r in rows}
    arms = {arm: [scores[r['run']] for r in rows if r['batch'] == '4' and r['group'] == arm]
            for arm in ['control','day7']}
    real = mean(arms['day7']) - mean(arms['control'])
    groups = [[r['run'] for r in rows if r['batch'] == '3' and r['group'] == 'control' and r['arch'] == arch]
              for arch in ['mandibular','maxillary']]
    values = []
    for mand in combinations(groups[0], 3):
        for maxil in combinations(groups[1], 2):
            a, b = set(mand), set(maxil)
            comparison = (8 * mean(scores[r] for r in groups[0] if r not in a) +
                          11 * mean(scores[r] for r in groups[1] if r not in b)) / 19
            values.append(comparison - mean(scores[r] for r in a | b))
    expected = read('calibration.tsv')[0]
    assert len(values) == int(expected['allocations']) == 7920
    assert abs(real-float(expected['true_day7_ref'])) < 1e-10
    assert abs(median(abs(v) for v in values) - float(expected['median_absolute'])) < 1e-10
    assert sum(abs(v) >= abs(real)-1e-10 for v in values) == int(expected['n_at_least_observed']) == 1411
    print('Primary untreated calibration and all 7,920 assignments independently checked.')


if __name__ == '__main__':
    main()
