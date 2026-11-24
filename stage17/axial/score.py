"""Compare gene and panel-mean errors for archived, unfiltered HOTM contrasts."""
import csv
import gzip
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGES = HERE.parents[1]


def scores(p, y):
    pred, obs = np.mean(p), np.mean(y)
    panel_zero, panel_err = abs(obs), abs(obs-pred)
    gene_zero, gene_err = np.mean(abs(y)), np.mean(abs(y-p))
    return dict(predicted_mean=pred, observed_mean=obs, zero_panel_error=panel_zero,
                model_panel_error=panel_err, panel_gain=panel_zero-panel_err,
                zero_gene_mae=gene_zero, model_gene_mae=gene_err, gene_gain=gene_zero-gene_err)


def main():
    with gzip.open(STAGES/'stage8/axial_log2_multipliers.tsv.gz', 'rt', encoding='utf8') as handle:
        panel = [r for r in csv.DictReader(handle, delimiter='\t') if r['in_panel']=='True']
    assert len({r['gene'] for r in panel}) == len(panel) == 101
    summary = []
    for cond in ['E3', 'I3', 'E10', 'I10']:
        data = json.loads((STAGES/f'stage8/hotm/deg_C0_{cond}.json').read_bytes())
        columns = {name:i for i,name in enumerate(data['header'])}
        effects = {r[columns['gene_name']]:r[columns['log2FoldChange']] for r in data['rows']}
        assert len(effects) == len(data['rows'])
        p = np.array([float(r[cond]) for r in panel])
        y = np.array([effects.get(r['gene']) for r in panel], dtype=float)
        keep = np.isfinite(y)
        m = scores(p[keep], y[keep])
        gain_sum = (abs(y[keep]) - abs(p[keep]-y[keep])).sum()
        gain_radius = abs(p[~keep]).sum()
        if np.all(keep):
            panel_lo = panel_hi = m['panel_gain']
        else:
            panel_lo, panel_hi = -abs(np.mean(p)), abs(np.mean(p))
        observed = np.flatnonzero(keep)
        n_robust = 0
        for i in observed:
            mask = keep.copy()
            mask[i] = False
            d = scores(p[mask], y[mask])
            n_robust += (d['panel_gain']>0) and (d['gene_gain']<0)
        record = dict(condition=cond, atlas_day=3 if cond.endswith('3') else 7,
                      n_total=len(p),n_measured=int(keep.sum()),
                      missing='|'.join(r['gene'] for i,r in enumerate(panel) if not keep[i]))
        record.update({'measured_'+key:value for key,value in m.items()})
        record.update(full_panel_gain_lo=panel_lo,full_panel_gain_hi=panel_hi,
                      full_gene_gain_lo=(gain_sum-gain_radius)/len(p),
                      full_gene_gain_hi=(gain_sum+gain_radius)/len(p),
                      measured_gene_deletions_with_reversal=n_robust)
        summary.append(record)
    with (HERE/'scores.tsv').open('w',newline='',encoding='utf8') as fp:
        w = csv.DictWriter(fp, fieldnames=list(summary[0]), delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(summary)


if __name__ == '__main__':
    main()
