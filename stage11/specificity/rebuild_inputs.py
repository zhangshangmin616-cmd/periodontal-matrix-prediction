"""Rebuild gene-level rank shifts from the stage 5 count files and forecasts."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata


HERE = Path(__file__).resolve().parent
WORK = HERE.parent.parent
STAGE5 = WORK / 'stage5'
sys.path.insert(0, str(STAGE5))
from forecast import CLASSES, OUT, Reference, load_bulk, meta
from external import layout, one_run


def ranks(x):
    return (rankdata(x, method='average') - 1) * (100.0 / (len(x) - 1))


def within_group_shift(ref, values, weights):
    before = ref.base @ weights
    after = np.zeros_like(before)
    for side in ['PP7D', 'TP7D']:
        side_mean = ref.expr[side] @ weights
        side_mean += (ref.base[:, -1] - ref.expr[side][:, -1]) * weights[-1]
        after += side_mean
    transferred = values * ((after * 0.5 + 1) / (before + 1))
    return ranks(transferred) - ranks(values)


def get_weights(name, identifier):
    saved = pd.read_csv(OUT / name, sep='\t')
    return saved.pivot(index=identifier, columns='cell_class', values='weight')[CLASSES]


def main():
    ref = Reference()
    info = meta()
    bulk = load_bulk(info.index.tolist(), ref.genes)
    weights = get_weights('weights.tsv', 'sample')
    kim_pred = pd.read_csv(OUT / 'predictions.tsv', sep='\t').query("setting=='state_only'").set_index('patient')
    kim_obs = pd.read_csv(OUT / 'patient_errors.tsv', sep='\t').query("model=='zero'").set_index('patient')
    labels, fitted, observed = [], [], []

    def append(label, controls, treated, get_values, weight_table, expect_p, expect_y):
        c = {key: get_values(key) for key in controls}
        control_ranks = np.mean([ranks(x) for x in c.values()], axis=0)
        y = np.mean([ranks(get_values(k)) for k in treated], axis=0) - control_ranks
        p = np.mean([within_group_shift(ref, x, weight_table.loc[k].to_numpy()) for k, x in c.items()], axis=0)
        np.testing.assert_allclose([p[ref.idx].mean(), y[ref.idx].mean()], [expect_p, expect_y], rtol=0, atol=1e-9)
        labels.append(label)
        fitted.append(p)
        observed.append(y)

    for donor, group in info.groupby('patient'):
        append(f'Kim_{donor}', group.index[group.arm == 'control'], group.index[group.arm == 'force'],
               lambda key: bulk[key].to_numpy(), weights, kim_pred.loc[donor, 'state'], kim_obs.loc[donor, 'observed'])
    runs, ids = layout(ref)
    z_weights = get_weights('external_weights.tsv', 'run')
    saved = pd.read_csv(OUT / 'external_errors.tsv', sep='\t').query("model=='state_only' and analysis=='batch'").set_index('time')
    for batch, time in [(4, 'day7'), (3, 'day28')]:
        group = runs[runs.batch == batch]
        append(f'Zhao_{time}', group.index[group.group == 'control'], group.index[group.group == time],
               lambda key: one_run(key, ids), z_weights, saved.loc[time, 'predicted'], saved.loc[time, 'observed'])
    features = pd.read_csv(STAGE5 / 'reference/features.tsv', sep='\t').set_index('gene').loc[ref.genes]
    with np.load(WORK / 'stage10/reference_counts.npz', allow_pickle=False) as source:
        positions = pd.Index(source['genes']).get_indexer(ref.genes)
        assert (positions >= 0).all()
        cp_counts = source['CP'][positions].sum(axis=1)
    np.savez_compressed(HERE / 'rank_changes.npz', gene=np.asarray(ref.genes, dtype=str),
                        component=features.ecm_component.to_numpy(dtype=str), cp_counts=cp_counts,
                        units=np.asarray(labels, dtype=str), fitted=np.stack(fitted), observed=np.stack(observed))
    print('Verified 13 complete paired or group contrasts against the frozen aggregate endpoints.')


if __name__ == '__main__':
    main()
