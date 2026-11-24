"""Compare panel-level and gene-level loss against expression-matched non-ECM sets."""
import argparse
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SEED = 20261124
DRAW_COUNT = 10_000


def table(path, rows):
    with path.open('w', newline='', encoding='utf8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def matched_sets(cp, components):
    """Assign ties to the same decile based on the midrank of integer atlas counts."""
    _, inverse, count = np.unique(cp, return_counts=True, return_inverse=True)
    end = np.cumsum(count)
    midrank = (end - (count + 1) / 2) / (len(cp) - 1)
    strata = np.minimum(9, (10 * midrank[inverse]).astype(int))
    target = np.flatnonzero(components == 'other')
    rng = np.random.Generator(np.random.PCG64(SEED))
    match = np.empty((DRAW_COUNT, len(target)), dtype=np.int32)
    start = 0
    bins = []
    for b in range(10):
        n = int(np.sum(strata[target] == b))
        pool = np.flatnonzero((strata == b) & (components == 'outside_ecm'))
        if n:
            assert len(pool) >= n
            for row in match:
                row[start:start+n] = rng.choice(pool, size=n, replace=False)
        bins.append(dict(decile=b+1, target_genes=n, candidate_genes=len(pool),
                         atlas_cp_min=int(cp[strata == b].min()), atlas_cp_max=int(cp[strata == b].max())))
        start += n
    assert start == 101 and np.all(components[match] == 'outside_ecm')
    assert np.all(np.diff(np.sort(match, axis=1),axis=1) > 0)
    return target, match, bins


def losses(y, p, sets):
    """Return separate panel and gene MAEs; units always receive equal weight."""
    dy = y[:, sets].mean(axis=2)
    dp = p[:, sets].mean(axis=2)
    panel_zero = np.abs(dy).mean(axis=0)
    panel_model = np.abs(dy-dp).mean(axis=0)
    by_gene_zero = np.abs(y).mean(axis=0)[sets].mean(axis=1)
    by_gene_model = np.abs(y-p).mean(axis=0)[sets].mean(axis=1)
    return dict(panel=(panel_zero, panel_model), gene=(by_gene_zero, by_gene_model))


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument('--input', type=Path, default=HERE/'rank_changes.npz')
    cli.add_argument('--outdir', type=Path, default=HERE)
    args = cli.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    with np.load(args.input, allow_pickle=False) as a:
        units, observed, fitted, cp, components = (a[k] for k in ('units','observed','fitted','cp_counts','component'))
        assert observed.shape == fitted.shape == (13,10806)
        genes=a['gene']
    target, match, bins = matched_sets(cp, components)
    table(args.outdir/'matching_deciles.tsv', bins)
    studies = {'Kim_week3': np.flatnonzero(np.char.startswith(units, 'Kim_')),
               'Zhao_day7': np.flatnonzero(units == 'Zhao_day7'), 'Zhao_day28': np.flatnonzero(units == 'Zhao_day28')}
    assert [len(v) for v in studies.values()] == [11,1,1]
    results, unit_rows, gene_rows, all_draws = [], [], [], {}
    for study, idx in studies.items():
        y, p = observed[idx], fitted[idx]
        actual = losses(y, p, target[None,:])
        for local, unit in enumerate(idx):
            per_unit = losses(y[local:local+1], p[local:local+1], target[None,:])
            for metric,(z,m) in per_unit.items():
                unit_rows.append(dict(unit=units[unit], metric=metric, zero_mae=z[0], model_mae=m[0], gain=z[0]-m[0]))
        for g in target:
            ez, em = np.abs(y[:,g]).mean(), np.abs(y[:,g]-p[:,g]).mean()
            gene_rows.append(dict(cohort=study, gene=genes[g], zero_mae=ez, model_mae=em, gain=ez-em))
        random_gains = {metric: [] for metric in actual}
        for offset in range(0,len(match),128):
            for metric,(z,m) in losses(y,p,match[offset:offset+128]).items():
                random_gains[metric].append(z - m)
        for metric,(z,m) in actual.items():
            g = (z-m)[0]
            null = np.concatenate(random_gains[metric])
            all_draws[study+'_'+metric] = null
            lo, med, hi = np.quantile(null, [.025,.5,.975])
            results.append(dict(cohort=study, metric=metric, units=len(idx), panel_genes=len(target),
                                model_mae=float(m[0]), zero_mae=float(z[0]), gain=float(g),
                                reference_draws=len(match), reference_q025=lo, reference_median=med, reference_q975=hi,
                                reference_positive=int(np.sum(null > 1e-12)),
                                reference_equal_or_higher=int(np.sum(null >= g-1e-12))))
    table(args.outdir/'comparison.tsv', results)
    table(args.outdir/'unit_errors.tsv', unit_rows)
    table(args.outdir/'target_gene_errors.tsv', gene_rows)
    np.savez_compressed(args.outdir/'matched_gains.npz', **all_draws)
    # Conditional on the published profiles and completed cross-fitted forecasts.
    # Keep the two outcomes together when resampling the same 11 individuals.
    kim = studies['Kim_week3']
    rng = np.random.Generator(np.random.PCG64(SEED + 1))
    resamples = rng.integers(len(kim), size=(50_000, len(kim)))
    patients = []
    for metric in ('panel', 'gene'):
        gains = np.array([r['gain'] for r in unit_rows if r['unit'].startswith('Kim_') and r['metric'] == metric])
        lo, hi = np.quantile(gains[resamples].mean(axis=1), [.025, .975])
        patients.append(dict(metric=metric, participants=len(gains), improve=int(sum(gains>1e-12)),
                             worsen=int(sum(gains < -1e-12)), ties=int(sum(np.abs(gains)<=1e-12)),
                             mean_gain=gains.mean(), bootstrap_resamples=len(resamples), seed=SEED+1,
                             percentile_025=lo, percentile_975=hi))
    table(args.outdir/'kim_patient_sensitivity.tsv', patients)
    for r in results:
        print(r['cohort'], r['metric'], *(round(r[k],4) for k in ('model_mae','zero_mae','gain','reference_median')),
              r['reference_equal_or_higher'])


if __name__ == '__main__':
    main()
