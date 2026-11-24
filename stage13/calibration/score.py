"""Recompute periodontal tissue ranks and enumerate fixed jaw-matched contrasts."""
import csv
import gzip
import itertools
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
STAGE5 = HERE.parents[1] / 'stage5'
ZHAO = STAGE5 / 'zhao'
ARCHES = ('mandibular', 'maxillary')


def records(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_table(name, rows):
    with (HERE / name).open('w', encoding='utf8', newline='') as stream:
        out = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        out.writeheader()
        out.writerows(rows)


def score_run(run, ids, targets):
    with gzip.open(ZHAO / 'counts' / (run + '.tsv.gz'), 'rt') as stream:
        header = next(stream).rstrip().split('\t')
        assert header[0] == 'Geneid'
        raw = dict(line.rstrip().split('\t') for line in stream)
    counts = np.array([float(raw[identifier]) for identifier in ids])
    assert len(counts) == 10806 and np.isfinite(counts).all()
    ordered = np.sort(counts)
    selected = counts[targets]
    left = ordered.searchsorted(selected, side='left')
    right = ordered.searchsorted(selected, side='right')
    rank0 = (left + right - 1) / 2.0
    return float(np.mean(rank0) * 100 / (len(counts) - 1))


def enumerate_differences(rows, control_counts, test_weights, baseline_weights):
    values = [np.array([r['score'] for r in rows if r['arch'] == arch]) for arch in ARCHES]
    contrasts = []
    for arr, n, w1, w0 in zip(values, control_counts, test_weights, baseline_weights):
        config = np.asarray(list(itertools.combinations(range(len(arr)), n)))
        control_sums = np.sum(arr[config], axis=1)
        held_mean = (arr.sum() - control_sums) / (len(arr) - n)
        contrasts.append(w1 * held_mean - w0 * control_sums / n)
    return np.add.outer(*contrasts).ravel()


def main():
    mapping = records(ZHAO / 'feature_map.tsv')
    ids = [r['entrez_id'] for r in mapping]
    targets = [i for i, r in enumerate(mapping) if r['ecm_component'] == 'other']
    assert len(targets) == 101 and len(set(ids)) == len(ids)
    history = {r['run']: float(r['score']) for r in records(HERE / 'archived_scores.tsv')}
    samples = []
    for row in records(ZHAO / 'joined_manifest.tsv'):
        score = score_run(row['run'], ids, targets)
        assert abs(score - history[row['run']]) < 1e-10
        samples.append(dict(run=row['run'], batch=row['batch'], group=row['group'],
                            arch=row['arch'], score=score))
    write_table('specimens.tsv', samples)
    day7 = [r for r in samples if r['batch'] == '4']
    controls7 = [r for r in day7 if r['group'] == 'control']
    controls28 = [r for r in samples if r['batch'] == '3' and r['group'] == 'control']
    loaded7 = [r for r in day7 if r['group'] == 'day7']
    measured = np.mean([r['score'] for r in loaded7]) - np.mean([r['score'] for r in controls7])
    assert np.isclose(measured, 0.6286071311532143)
    equal = (.5, .5)
    treated = (8/19, 11/19)
    baseline = (.6, .4)
    tasks = [
        ('other_batch_controls_observed_jaw_mix', controls28, (3,2), treated, baseline, 7920),
        ('other_batch_controls_unweighted', controls28, (3,2), (7/17,10/17), baseline, 7920),
        ('other_batch_controls_equal_jaw', controls28, (3,2), equal, equal, 7920),
        ('same_batch_five_controls_equal_jaw', controls7, (1,1), equal, equal, 6),
        ('day7_within_jaw_relabelling', day7, (3,2), treated, baseline, 12870),
    ]
    summaries = []
    distributions = {}
    for name, rows, n, test_weights, base_weights, count in tasks:
        values = enumerate_differences(rows, n, test_weights, base_weights)
        assert len(values) == count
        # Use the same weighting on the real contrast as on the calibration.
        def arm_mean(arm, w):
            return sum(a * np.mean([r['score'] for r in arm if r['arch'] == arch])
                       for a, arch in zip(w, ARCHES))
        reference_shift = arm_mean(loaded7, test_weights) - arm_mean(controls7, base_weights)
        hits = int(np.sum(np.abs(values) >= abs(reference_shift) - 1.e-10))
        signed = np.quantile(values, [.025,.5,.975])
        summaries.append(dict(comparison=name, allocations=count, true_day7_ref=reference_shift,
                              median_absolute=float(np.median(abs(values))), signed_q025=signed[0],
                              signed_q50=signed[1], signed_q975=signed[2], minimum=values.min(),
                              maximum=values.max(), n_at_least_observed=hits, fraction=hits/count))
        distributions[name] = values
    assert min(abs(distributions['day7_within_jaw_relabelling']-measured)) < 1.e-10
    write_table('calibration.tsv', summaries)
    np.savez_compressed(HERE/'distributions.npz', **distributions)


if __name__ == '__main__':
    main()
