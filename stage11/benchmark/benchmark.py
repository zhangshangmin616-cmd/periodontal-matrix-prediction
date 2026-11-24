"""Retrospective time-forward ablations of cell-group RNA changes; NumPy only."""
import argparse
import csv
from fractions import Fraction
import hashlib
from itertools import groupby
import json
from math import lcm
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
SIDES = {'pressure': 'PP', 'tension': 'TP'}
MODELS = ('unchanged', 'state_only', 'composition_only', 'state_and_composition', 'measured_day3')
LAYERS = ('direct', 'baseline_transfer')


def midpercentile(values):
    ranks = np.empty(len(values))
    offset = 0
    for _, entries in groupby(sorted(enumerate(values), key=lambda x: x[1]), key=lambda x: x[1]):
        indices = [index for index, _ in entries]
        ranks[indices] = 100 * (offset + (len(indices)-1)/2) / (len(values)-1)
        offset += len(indices)
    return ranks


def rows(path):
    with path.open(encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, records):
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, records[0], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def observed(counts, depths, background):
    totals = counts.sum(axis=1)
    abundance = totals * (1.e6 / depths.sum())
    return np.log2(1 + abundance), midpercentile(totals[background].tolist())


def predict(model, cp_counts, cp_depths, d3_counts, d3_depths, background):
    # Use only CP and day 3. Heterogeneous groups contribute RNA mass fractions.
    weights = d3_depths if model in ('composition_only', 'state_and_composition', 'measured_day3') else cp_depths
    use_day3 = np.zeros(len(weights), dtype=bool)
    if model in ('state_only', 'state_and_composition'):
        use_day3[:-1] = True
    if model == 'measured_day3':
        use_day3[:] = True
    counts = np.where(use_day3, d3_counts, cp_counts)
    depths = np.where(use_day3, d3_depths, cp_depths)
    abundance = np.sum(counts * (1.e6 / depths) * (weights / weights.sum()), axis=1)
    pooled_counts = cp_counts.sum(axis=1)
    pooled_depth = int(cp_depths.sum())
    untreated = pooled_counts * (1.e6 / pooled_depth)
    # In the alternative transfer rule, a smoothed relative response multiplies
    # *unadjusted* CP abundance (the stage-10 rank convention).
    transferred = untreated * (1 + abundance) / (1 + untreated)
    # An exact rational construction ensures predicted ties have common midranks.
    common = lcm(*(int(depth) for depth in depths))
    factors = [1_000_000 * int(w) * (common // int(d)) for w, d in zip(weights, depths)]
    denominator = int(weights.sum()) * common
    numerators = [sum(int(x)*a for x, a in zip(row, factors)) for row in counts[background]]
    exact_transfer = [Fraction(int(n)*(denominator+m), pooled_depth + 1_000_000*int(n))
                      for n, m in zip(pooled_counts[background], numerators)]
    assert np.allclose(abundance[background], np.asarray([k / denominator for k in numerators]),
                       rtol=1e-13, atol=1e-13)
    return {'direct': (np.log2(1 + abundance), midpercentile(numerators)),
            'baseline_transfer': (np.log2(1 + transferred), midpercentile(exact_transfer))}


def metrics(y, p):
    # Changes smaller than roundoff are treated as zero for direction counting.
    def signs(x):
        return np.where(np.abs(x) <= 1.e-12, 0, np.sign(x))
    sy, sp = signs(y), signs(p)
    changed = sy != 0
    right = (sp == sy) & changed
    opposite = (sp == -sy) & changed
    abstain = (sp == 0) & changed
    assert np.all(right | opposite | abstain | ~changed)
    return dict(n=len(y), observed_change=int(changed.sum()), correct=int(right.sum()),
                opposite=int(opposite.sum()), abstain_on_changed=int(abstain.sum()),
                both_zero=int((~changed & (sp == 0)).sum()),
                predicted_change_on_observed_zero=int((~changed & (sp != 0)).sum()),
                mae=float(np.abs(p - y).mean()), zero_mae=float(np.abs(y).mean()),
                gain=float((np.abs(y) - np.abs(p-y)).mean()),
                direction_fraction=float(right.sum()/changed.sum()) if changed.any() else '')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, default=WORK/'stage10/reference_counts.npz')
    parser.add_argument('--features', type=Path, default=WORK/'stage5/reference/features.tsv')
    parser.add_argument('--earlier', type=Path, default=WORK/'stage10/time_transfer.tsv')
    parser.add_argument('--external-calibration', type=Path, default=WORK/'stage5/results/calibration.tsv')
    parser.add_argument('--out', type=Path, default=HERE)
    a = parser.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    features = rows(a.features)
    with np.load(a.counts, allow_pickle=False) as z:
        names = z['genes'].tolist()
        idx = {g: i for i, g in enumerate(names)}
        assert len(names) == len(idx) == 37433
        background = np.array([idx[r['gene']] for r in features])
        panel = np.array([r['ecm_component'] == 'other' for r in features])
        assert len(background) == 10806 and panel.sum() == 101
        assert z['classes'].tolist() == ['fibroblast_like','mural','immune','endothelial','glial','residual']
        cp = z['CP']
        depths = z['CP_total']
        base = observed(cp, depths, background)
        # Complete all day-3 forward predictions before accessing the day-7 arrays.
        computed = {side: {model: predict(model, cp, depths, z[key+'3D'], z[key+'3D_total'], background)
                           for model in MODELS} for side, key in SIDES.items()}
        truth = {side: observed(z[key+'7D'], z[key+'7D_total'], background) for side, key in SIDES.items()}
        for side, key in SIDES.items():
            remeasured = observed(z[key+'3D'], z[key+'3D_total'], background)
            for actual, tested in [(base, computed[side]['unchanged'][layer]) for layer in LAYERS] + [
                (remeasured, computed[side]['measured_day3']['direct'])]:
                assert np.allclose(actual[0], tested[0], rtol=1e-13, atol=1e-13)
                assert np.array_equal(actual[1], tested[1])
    slices = [(0, 'log2_cpm1', {'all_unique': slice(None), 'background': background, 'panel': background[panel]}),
              (1, 'percentile', {'background': slice(None), 'panel': panel})]
    # The external per-person nested fits used the target panel's tissue score.
    # Both older model variants have the same median multiplier. This sensitivity
    # applies it uniformly to every *predicted endpoint change* without refitting.
    previous = rows(a.external_calibration)
    medians = {float(np.median([float(r['multiplier']) for r in previous if r['model']==m]))
               for m in ('composition', 'state')}
    assert len(medians)==1
    damping = medians.pop()
    rows_gene, rows_panel, rows_damped = [], [], []
    for layer in LAYERS:
        for model in MODELS:
            forecasts = {s: computed[s][model][layer] for s in SIDES}
            for i, scale, sets in slices:
                y_pressure, y_tension = [truth[s][i] - base[i] for s in SIDES]
                p_pressure, p_tension = [forecasts[s][i] - base[i] for s in SIDES]
                comparisons = [('pressure', y_pressure, p_pressure), ('tension', y_tension, p_tension),
                               ('pressure_minus_tension', y_pressure-y_tension, p_pressure-p_tension),
                               ('swapped_pressure_minus_tension', y_pressure-y_tension, p_tension-p_pressure)]
                for side, y, p in comparisons:
                    details = dict(layer=layer, model=model, comparison=side, scale=scale)
                    for subset, mask in sets.items():
                        rows_gene.append(dict(**details, subset=subset, **metrics(y[mask], p[mask])))
                        reduced = metrics(y[mask], damping*p[mask])
                        rows_damped.append(dict(**details, subset=subset, endpoint='gene_wise',
                            factor=damping, mae=reduced['mae'], zero_mae=reduced['zero_mae'], gain=reduced['gain']))
                    if scale == 'percentile':
                        real, forecast = float(y[panel].mean()), float(p[panel].mean())
                        rows_panel.append(dict(**details, observed=real, predicted=forecast,
                                               absolute_error=abs(real-forecast), zero_error=abs(real),
                                               gain=abs(real)-abs(real-forecast)))
                        rows_damped.append(dict(**details, subset='panel', endpoint='mean_percentile',
                            factor=damping, mae=abs(real-damping*forecast), zero_mae=abs(real),
                            gain=abs(real)-abs(real-damping*forecast)))
    summaries = []
    for layer in LAYERS:
        for model in MODELS:
            for _, scale, subsets in slices:
                for subset in subsets:
                    selected = [r for r in rows_gene if r['comparison'] in SIDES and r['scale']==scale
                                and r['subset']==subset and r['model']==model and r['layer']==layer]
                    assert len(selected)==2
                    summaries.append(dict(layer=layer, model=model, scale=scale, subset=subset,
                        **{k:float(np.mean([r[k] for r in selected])) for k in ['mae','zero_mae','gain','direction_fraction']}))
    check = {r['contrast']: r for r in rows(a.earlier)}
    deviations = [abs(r['predicted']-float(check[r['comparison']]['state'])) for r in rows_panel
                  if r['comparison'] in SIDES and r['model']=='state_only' and r['layer']=='baseline_transfer']
    assert len(deviations)==2 and max(deviations)<1.e-10
    write(a.out/'gene_metrics.tsv', rows_gene)
    write(a.out/'two_side_mean.tsv', summaries)
    write(a.out/'panel_rank.tsv', rows_panel)
    write(a.out/'external_damping.tsv', rows_damped)
    audit = dict(input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [a.counts,a.features,a.earlier,a.external_calibration]},
                 unique=len(names), background=len(background), panel=int(panel.sum()),
                 largest_stage10_panel_state_deviation=max(deviations),
                 measured_days={s: k+'7D' for s,k in SIDES.items()},
                 prediction_days={s: ['CP',k+'3D'] for s,k in SIDES.items()}, numpy=np.__version__)
    (a.out/'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')


if __name__=='__main__':
    main()
