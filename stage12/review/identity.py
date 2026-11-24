"""Evaluate the gene assignment while holding every panel-average forecast fixed."""
import argparse
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE.parents[1] / 'stage11/specificity/rank_changes.npz'
SEED = 20261126
DRAWS = 10000


def save_table(path, rows):
    with path.open('w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, rows[0], delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def strata(counts, chosen):
    _, inverse, group_size = np.unique(counts, return_inverse=True, return_counts=True)
    rank0 = (np.cumsum(group_size) - (group_size + 1) / 2)[inverse]
    decile = np.minimum(9, (10 * rank0 / (len(counts) - 1)).astype(int))
    return decile[chosen]


def generate_permutations(bins, rng):
    base = np.arange(len(bins))
    permutations = np.tile(base, (DRAWS, 1))
    for label in np.unique(bins):
        members = base[bins == label]
        for row in permutations:
            row[members] = rng.permutation(members)
    assert np.all(np.sort(permutations, axis=1) == base)
    assert np.all(bins[permutations] == bins)
    return permutations


def expected_mae(y, p, bins):
    # Every source-to-destination pair in the same stratum is equally probable.
    total = 0.
    for label in np.unique(bins):
        mask = bins == label
        distances = np.abs(y[:, mask, None] - p[:, None, mask])
        total += distances.sum() / mask.sum()
    return total / y.size


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    ap.add_argument('--out', type=Path, default=HERE)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with np.load(args.input, allow_pickle=False) as d:
        idx = np.flatnonzero(d['component'] == 'other')
        labels = d['units'].astype(str)
        y = d['observed'][:, idx]
        p = d['fitted'][:, idx]
        bins = strata(d['cp_counts'], idx)
    assert y.shape == p.shape == (13, 101)
    units = {'Seoul_week3': np.char.startswith(labels, 'Kim_'),
             'Leuven_day7': labels == 'Zhao_day7', 'Leuven_day28': labels == 'Zhao_day28'}
    rng = np.random.Generator(np.random.PCG64(SEED))
    baselines, results = [], []
    for study, mask in units.items():
        fitted = p[mask]
        observed = y[mask]
        scalar = np.broadcast_to(fitted.mean(axis=1)[:, None], fitted.shape)
        for name, forecast in [('unchanged', np.zeros_like(fitted)), ('state_only', fitted),
                               ('same_panel_mean', scalar)]:
            error = forecast - observed
            baselines.append(dict(study=study, per_unit_model=name, n_units=mask.sum(),
                                  genes=len(idx), panel_mae=np.abs(error.mean(axis=1)).mean(),
                                  gene_mae=np.abs(error).mean()))
        assert abs(baselines[-1]['panel_mae'] - baselines[-2]['panel_mae']) < 1e-12

    for control, groups in [('matched_atlas_decile', bins), ('unrestricted', np.zeros_like(bins))]:
        perms = generate_permutations(groups, rng)
        # Shared label permutations for all patients and times retain cross-unit
        # dependence between a gene's forecasts. Construct loss draws in batches.
        all_draws = np.empty((len(labels), DRAWS))
        for i in range(0, DRAWS, 200):
            draw = perms[i:i + 200]
            pp = p[:, draw]
            np.testing.assert_allclose(pp.mean(axis=2) - p.mean(axis=1)[:, None], 0, atol=1e-12)
            all_draws[:, i:i + len(draw)] = np.abs(pp - y[:, None, :]).mean(axis=2)
        for study, mask in units.items():
            observed, forecast = y[mask], p[mask]
            exact = expected_mae(observed, forecast, groups)
            draws = all_draws[mask].mean(axis=0)
            original = np.abs(observed - forecast).mean()
            lo, med, hi = np.quantile(draws, [.025, .5, .975])
            results.append(dict(study=study, control=control, genes=len(idx), n_units=mask.sum(),
                                cell_atlas_decile_sizes=','.join(map(str, np.unique(groups, return_counts=True)[1])),
                                actual_gene_mae=original, assignment_expected_mae=exact,
                                actual_gain_vs_expected=exact-original,
                                sampled_mae_q025=lo,sampled_mae_q50=med,sampled_mae_q975=hi,
                                sampled_as_good_or_better=int(np.sum(draws<=original+1e-12)),
                                sampled_mean_mae=draws.mean(), draws=len(draws)))
            assert abs(draws.mean()-exact) < max(0.002, draws.std() * 5/np.sqrt(DRAWS))
    save_table(args.out / 'baselines.tsv', baselines)
    save_table(args.out / 'gene_identity.tsv', results)


if __name__ == '__main__':
    main()
