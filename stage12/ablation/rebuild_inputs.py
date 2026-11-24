"""Gather public-data inputs and reference results for an offline cell-group ablation."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
STAGE5 = WORK / 'stage5'
sys.path.insert(0, str(STAGE5))
from forecast import CLASSES, OUT, Reference, load_bulk, meta
from external import layout, one_run


def weights(filename, identifier):
    frame = pd.read_csv(OUT/filename, sep='\t')
    return frame.pivot(index=identifier, columns='cell_class', values='weight')[CLASSES]


def main():
    ref = Reference()
    old = pd.read_csv(OUT/'gene_deltas_exploratory.tsv', sep='\t', dtype={'unit':str})
    old_ids = np.where(old.cohort == 'Kim', 'Kim_', 'Zhao_')
    old['unit_name'] = [a+b for a,b in zip(old_ids, old.unit, strict=True)]
    with np.load(WORK/'stage11/specificity/rank_changes.npz', allow_pickle=False) as previous:
        np.testing.assert_array_equal(previous['gene'], ref.genes)
        np.testing.assert_array_equal(ref.idx, np.flatnonzero(previous['component']=='other'))
        units = previous['units']
        previous_full = previous['fitted'][:,ref.idx]
        observed = previous['observed'][:,ref.idx]
    comp = old.pivot(index='unit_name', columns='gene', values='composition').loc[units,ref.genes[ref.idx]]
    old_observed = old.pivot(index='unit_name', columns='gene', values='observed').loc[units,ref.genes[ref.idx]]
    np.testing.assert_allclose(old_observed, observed, atol=1e-11, rtol=0)
    controls, mapping, q_rows, values = [], [], [], []
    def append(k, unit, qt, x):
        controls.append(k)
        mapping.append(str(unit))
        q_rows.append(qt.loc[k].to_numpy())
        values.append(x)
    info = meta().query("arm == 'control'")
    cpm = load_bulk(info.index.tolist(), ref.genes)
    w = weights('weights.tsv', 'sample')
    for k, r in info.iterrows():
        append(k, f'Kim_{r.patient}', w, cpm[k].to_numpy())
    runs, ids = layout(ref)
    z = weights('external_weights.tsv', 'run')
    for k,r in runs.query("group == 'control'").iterrows():
        append(k, 'Zhao_day7' if r.batch==4 else 'Zhao_day28', z, one_run(k,ids))
    assert set(mapping) == set(units)
    np.savez_compressed(HERE/'inputs.npz', classes=np.array(CLASSES), units=units,
                        sources=np.array(controls), unit_per_control=np.array(mapping),
                        genes=ref.genes.to_numpy(dtype=str), target=ref.idx,
                        base=ref.base, pressure=ref.expr['PP7D'], tension=ref.expr['TP7D'],
                        source_shares=np.stack([ref.w0, ref.w['PP7D'], ref.w['TP7D']]),
                        q=np.array(q_rows), control_cpm=np.array(values), observed=observed,
                        historical_full=previous_full, historical_composition=comp.to_numpy())
    print('Prepared',len(controls),'controls in',len(units),'clinical units')


if __name__ == '__main__':
    main()
