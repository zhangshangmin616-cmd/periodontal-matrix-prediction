"""Score published human axial-force group contrasts against frozen atlas predictors."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DAYS = (1, 3, 10)
NAMES = ('zero', 'atlas', 'reversed_atlas')


def read_direct(source, first, second):
    path = source / f'deg_{first}_{second}.json'
    data = json.loads(path.read_bytes())
    assert data['header'][1:4] == ['gene_name', 'baseMean', 'log2FoldChange']
    values = {row[1]: float(row[3]) for row in data['rows']}
    assert len(values) == len(data['rows']) and np.isfinite(list(values.values())).all()
    return values, (len(values), hashlib.sha256(path.read_bytes()).hexdigest())


def rank(v):
    unique, inverse, counts = np.unique(v, return_inverse=True, return_counts=True)
    center = counts.cumsum() - (counts + 1) / 2
    return center[inverse]


def metrics(genes, outcomes, prediction, method, group, day):
    y = np.asarray([outcomes[g] for g in genes])
    p = np.asarray([prediction[g] for g in genes], dtype=float)
    err = p - y
    def corr(a, b):
        return float(np.corrcoef(a, b)[0, 1]) if a.std() and b.std() else None
    return dict(day=day, group=group, method=method, n=len(genes),
                y_mean=float(y.mean()), p_mean=float(p.mean()), mean_bias=float(err.mean()),
                mae=float(abs(err).mean()), zero_mae=float(abs(y).mean()),
                rmse=float(np.mean(err ** 2) ** .5),
                gain_vs_zero=float((abs(y) - abs(err)).mean()),
                direction_correct=int((p*y>0).sum()), direction_opposite=int((p*y<0).sum()),
                direction_tie=int((p*y==0).sum()), pearson=corr(y, p), spearman=corr(rank(y), rank(p)),
                per_gene_model_closer=int((abs(err)<abs(y)).sum()),
                centered_mae=float(np.abs((p-np.median(p))-(y-np.median(y))).mean()),
                centered_zero_mae=float(np.abs(y-np.median(y)).mean()))


def write(path, rows):
    with path.open('w', newline='', encoding='utf8') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=HERE/'source')
    parser.add_argument('--atlas', type=Path, default=HERE.parents[1]/'stage8/axial_log2_multipliers.tsv.gz')
    parser.add_argument('--replay', action='store_true', help='Use bundled per-gene table instead of network-source tables.')
    args = parser.parse_args()
    with gzip.open(args.atlas, 'rt', encoding='utf8') as stream:
        frozen = {r['gene']: r for r in csv.DictReader(stream, delimiter='\t')}
    all_genes = list(frozen)
    panel = [g for g in all_genes if frozen[g]['in_panel']=='True']
    assert len(all_genes)==10806 and len(panel)==101
    cache = HERE/'direct_log2fc.tsv.gz'
    if args.replay:
        with gzip.open(cache, 'rt', encoding='utf8') as stream:
            compact = list(csv.DictReader(stream, delimiter='\t'))
        direct = {day: {r['gene']: float(r[f'day{day}']) for r in compact if r[f'day{day}']} for day in DAYS}
    else:
        direct = {}
        manifest = []
        for day in DAYS:
            fwd, fm = read_direct(args.source, f'E{day}', f'I{day}')
            rev, rm = read_direct(args.source, f'I{day}', f'E{day}')
            assert fwd.keys()==rev.keys() and all(fwd[g]==-rev[g] for g in fwd)
            direct[day] = fwd
            for label, (n, h) in zip(['E_to_I','I_to_E'], [fm, rm]):
                manifest.append(dict(day=day, order=label, n=n, sha256=h))
        write(HERE/'source_manifest.tsv', manifest)
        with cache.open('wb') as dest, gzip.GzipFile(fileobj=dest, mode='wb', mtime=0) as gz:
            gz.write(b'gene\tday1\tday3\tday10\n')
            for g in all_genes:
                gz.write(('\t'.join([g]+[str(direct[d].get(g,'')) for d in DAYS])+'\n').encode())
    predictions = {}
    for day in DAYS:
        atlas_day = 3 if day==1 else day
        scale = 1/3 if day==1 else 1
        predictions[day] = {g: (float(r['I'+str(atlas_day)])-float(r['E'+str(atlas_day)]))*scale for g,r in frozen.items()}
    rows = []
    bounds = []
    common = set.intersection(*[set(direct[d]) for d in DAYS])
    for day in DAYS:
        same = predictions[day]
        observed = direct[day]
        for group, group_list in [('fixed_background',all_genes), ('fixed_panel',panel)]:
            missing = sorted(set(group_list)-observed.keys())
            eligible = [g for g in group_list if g in observed]
            gain_sum = sum(abs(observed[g])-abs(observed[g]-same[g]) for g in eligible)
            width = sum(abs(same[g]) for g in missing)
            bounds.append(dict(day=day,group=group,n_all=len(group_list),n_missing=len(missing),
                               missing_genes='|'.join(missing),lower_gain=(gain_sum-width)/len(group_list),
                               upper_gain=(gain_sum+width)/len(group_list)))
            selections = [('available',eligible), ('all_days',[g for g in group_list if g in common])]
            if day==10:
                selections.append(('d3_d10', [g for g in eligible if g in direct[3]]))
            for selection, genes in selections:
                candidates = {'zero':dict.fromkeys(genes,0.),'atlas':same,'reversed_atlas':{g:-same[g] for g in genes}}
                if day==10 and selection in ['all_days', 'd3_d10']:
                    candidates['empirical_day3']=direct[3]
                for method, pred in candidates.items():
                    rows.append(metrics(genes, observed, pred, method, group+'_'+selection, day))
    write(HERE/'direct_metrics.tsv', rows)
    write(HERE/'missing_bounds.tsv', bounds)
    for r in rows:
        if r['group'].endswith('_available') and r['method']=='atlas' or r['method']=='empirical_day3':
            print(r['day'],r['group'],r['method'],r['n'],r['mae'],r['zero_mae'],r['gain_vs_zero'],r['spearman'])


if __name__ == '__main__':
    main()
