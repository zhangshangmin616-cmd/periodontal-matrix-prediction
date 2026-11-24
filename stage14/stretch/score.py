"""Cross-study prediction of tensile effects from public culture measurements."""
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
EXTERNAL = WORK/'stage13'/'external'


def read(path):
    with path.open(encoding='utf8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def save(name, rows):
    with (HERE/name).open('w', encoding='utf8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def tied_ranks(values):
    _, index, counts = np.unique(values, return_inverse=True, return_counts=True)
    return (np.cumsum(counts) - 0.5 * (counts + 1))[index]


def predictions(source):
    return {'zero':np.zeros_like(source), 'uniform':np.full_like(source,np.mean(source)),
            'full':source, 'quarter':0.25*source}


def loss(y, p):
    return dict(gene_mae=np.mean(abs(y-p)), panel_abs_error=abs(np.mean(y-p)),
                gene_gain_vs_zero=np.mean(abs(y))-np.mean(abs(y-p)),
                panel_gain_vs_zero=abs(np.mean(y))-abs(np.mean(y-p)))


def load_hong_kong():
    with np.load(EXTERNAL/'paired_counts.npz',allow_pickle=False) as h:
        genes = list(h['genes'])
        target = set(h['genes'][h['target']])
        cpm = h['counts'] * 1.e6 / h['total']
    assert len(genes)==len(set(genes))==10805 and len(target)==101
    effects = np.log2(1+cpm[:,5:])-np.log2(1+cpm[:,:5])
    return {g:e for g,e in zip(genes,effects)}, target


def run_pair(hk, osaka, genes, gate, scope, outputs):
    rows, detail, relationships = outputs
    subjects = np.stack([hk[g] for g in genes])
    x = subjects.mean(axis=1)
    y = np.array([osaka[g] for g in genes])
    same = (x*y)>0
    nonzero = (x*y)!=0
    labels = dict(gate=gate, scope=scope, genes=len(genes))
    relationships.append({**labels,'spearman':np.corrcoef(tied_ranks(x),tied_ranks(y))[0,1],
                          'concordant':int(same.sum()),'nonzero_in_both':int(nonzero.sum()),
                          'hk_mean':x.mean(),'osaka_mean':y.mean()})
    for direction, source, dest in [('HK12_to_OS48', x, y), ('OS48_to_HK12', y, x)]:
        for method,p in predictions(source).items():
            rows.append({**labels,'direction':direction,'method':method,**loss(dest,p)})
            if direction == 'OS48_to_HK12':
                for i in range(subjects.shape[1]):
                    detail.append({**labels,'method':method,'hk_donor':i+1,**loss(subjects[:,i],p)})


def run_vienna(hk, target):
    info = read(WORK/'stage12'/'culture'/'culture_readouts.tsv')
    published = {r['gene']:float(r['reported_log2']) for r in info if r['reported_log2']}
    assert {r['gene'] for r in info} == target and len(published)==40
    genes = sorted(target)
    x = np.array([hk[g].mean() for g in genes])
    present = np.array([g in published for g in genes])
    y = np.array([published[g] for g in genes if g in published])
    rows = []
    for method,p in predictions(x).items():
        available = p[present]
        known_gain = np.sum(np.abs(y) - np.abs(y-available))
        interval_width = np.sum(np.abs(p[~present]))
        rows.append(dict(direction='HK12_to_Vienna72', method=method, measured=present.sum(), total=len(genes),
                         source_mean=x.mean(), **loss(y,available),
                         full_panel_gene_gain_min=(known_gain-interval_width)/len(genes),
                         full_panel_gene_gain_max=(known_gain+interval_width)/len(genes)))
    save('vienna_bounds.tsv', rows)


def main():
    hk,target = load_hong_kong()
    chip = read(EXTERNAL/'chip_effects.tsv')
    rows, detail, relationships, common_targets = [], [], [], []
    for gate in ('positive_signal','all_published_values'):
        eligible = [r for r in chip if r['rule']==gate]
        osaka = {r['gene']:float(r['log2_change']) for r in eligible}
        assert len(osaka)==len(eligible)
        common = sorted(hk.keys() & osaka.keys())
        targets = sorted(target & set(common))
        assert len(targets)==(84 if gate=='positive_signal' else 95)
        for gene in targets:
            row = dict(gate=gate,gene=gene,osaka=osaka[gene],hong_kong_mean=hk[gene].mean())
            row.update({f'hong_kong_{i+1}':v for i,v in enumerate(hk[gene])})
            common_targets.append(row)
        for scope,genes in [('target',targets), ('background',common)]:
            run_pair(hk, osaka, genes, gate, scope, (rows,detail,relationships))
    save('common_target_effects.tsv',common_targets)
    save('group_scores.tsv',rows)
    save('donor_scores.tsv',detail)
    save('associations.tsv',relationships)
    run_vienna(hk,target)


if __name__ == '__main__':
    main()
