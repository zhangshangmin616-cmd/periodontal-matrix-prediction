"""Sensitivity and baseline-expression-matched panel analysis of tensile cultures."""
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def dataset():
    with np.load(WORK/'stage13/external/paired_counts.npz',allow_pickle=False) as dat:
        counts, totals = dat['counts'], dat['total']
        genes, target = dat['genes'], dat['target']
    with (WORK/'stage5/reference/features.tsv').open(encoding='utf8') as inp:
        group = {r['gene']:r['ecm_component'] for r in csv.DictReader(inp,delimiter='\t')}
    assert len(genes)==10805 and target.sum()==101 and len(set(genes))==len(genes)
    assert np.array_equal(target,[group[g]=='other' for g in genes])
    outside = np.array([group[g]=='outside_ecm' for g in genes])
    return genes, target, outside, counts, counts * 1e6/totals


def ranks(column):
    ordered = np.sort(column)
    return 50*(ordered.searchsorted(column) + ordered.searchsorted(column,side='right') - 1)/(len(column)-1)


def save(name, rows):
    with (HERE/name).open('w',encoding='utf8',newline='') as f:
        out = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        out.writeheader()
        out.writerows(rows)


def reference_panels(control, target, outside):
    mean_control = control.mean(axis=1)
    edges = np.quantile(mean_control,np.arange(1,10)/10)
    decile = np.searchsorted(edges,mean_control,side='right')
    rng = np.random.default_rng(20261124)
    draws = []
    strata = []
    for i in range(10):
        source = np.flatnonzero(outside & (decile==i))
        take = int(np.sum(target & (decile==i)))
        strata.append(dict(stratum=i,min_mean_cpm=float(mean_control[decile==i].min()),
                           max_mean_cpm=float(mean_control[decile==i].max()),eligible=len(source),target=take))
        if take:
            draws.append(np.stack([rng.choice(source,take,replace=False) for _ in range(10000)],axis=0))
    return np.hstack(draws),strata


def main():
    genes,target,outside,counts,cpm = dataset()
    c = cpm[:,:5]
    draws,strata = reference_panels(c,target,outside)
    save('strata.tsv',strata)
    detected = c >= 1
    rows, match, draws_out = [], [], {}
    eligibility = [dict(gene=g, at_least_1cpm_controls=int(np.sum(detected[i])),
                        positive_libraries=int(np.sum(counts[i]>0)))
                   for i,g in enumerate(genes) if target[i]]
    save('target_detectability.tsv',eligibility)
    for scale, matrix in [('log2',np.log2(1+cpm)),('rank',np.stack([ranks(a) for a in cpm.T],axis=1))]:
        y = matrix[:,5:]-matrix[:,:5]
        p = (y.sum(axis=1,keepdims=True)-y)/4
        for donor in range(5):
            train_detect = np.delete(detected,donor,axis=1)
            order = np.lexsort((genes[target],-abs(p[target,donor])))
            high = np.zeros(len(genes),bool)
            high[np.flatnonzero(target)[order[:21]]] = True
            groups = {'panel':target, 'panel_control_3of5':target & (detected.sum(1)>=3),
                      'panel_training_control_3of4':target & (train_detect.sum(1)>=3),
                      'panel_top21_training_effect':high,'panel_other80_training_effect':target & ~high,
                      'background':np.ones(len(genes),bool),'outside_panel':~target,'outside_ecm':outside}
            for name,mask in groups.items():
                study = y[mask,donor]
                pred = p[mask,donor]
                for method,vector in [('zero',np.zeros_like(pred)),('uniform_training_mean',pred.mean()),
                                      ('training_gene_mean',pred),('training_gene_quarter',pred*0.25)]:
                    diff = vector-study
                    panel = '' if scale=='rank' and name=='background' else float(abs(diff.mean()))
                    rows.append(dict(scale=scale, scope=name,donor=donor+1,method=method,genes=int(mask.sum()),
                                     gene_mae=float(abs(diff).mean()),panel_mae=panel))
        gdiff = (abs(y)-abs(p-y)).mean(1)
        gene = gdiff[draws].mean(1)
        panel = np.empty(len(draws))
        for start in range(0,len(draws),200):
            block = draws[start:start+200]
            panel[start:start+len(block)] = (abs(y[block].mean(1))-abs((p-y)[block].mean(1))).mean(1)
        t_gene = gdiff[target].mean()
        t_panel = (abs(y[target].mean(0))-abs((p-y)[target].mean(0))).mean()
        match.append(dict(scale=scale, draws=len(draws), target_gene_gain=t_gene,target_panel_gain=t_panel,
                          gene_positive=int((gene>0).sum()),panel_negative=int((panel<0).sum()),
                          positive_gene_negative_panel=int(((gene>0)&(panel<0)).sum()),
                          gene_at_least_target=int((gene>=t_gene).sum()),panel_at_most_target=int((panel<=t_panel).sum()),
                          gene_q025=np.quantile(gene,.025),gene_median=np.median(gene),gene_q975=np.quantile(gene,.975),
                          panel_q025=np.quantile(panel,.025),panel_median=np.median(panel),panel_q975=np.quantile(panel,.975)))
        draws_out[scale+'_gene_gain'] = gene
        draws_out[scale+'_panel_gain'] = panel
    summary = []
    for scale,scope,method in dict.fromkeys((r['scale'],r['scope'],r['method']) for r in rows):
        selected = [r for r in rows if (r['scale'],r['scope'],r['method'])==(scale,scope,method)]
        zero = [r for r in rows if (r['scale'],r['scope'],r['method'])==(scale,scope,'zero')]
        pg = np.array([r['panel_mae'] for r in zero]) if zero[0]['panel_mae']!='' else None
        pm = np.array([r['panel_mae'] for r in selected]) if pg is not None else None
        g = np.array([r['gene_mae'] for r in selected])
        z = np.array([r['gene_mae'] for r in zero])
        summary.append(dict(scale=scale,scope=scope,method=method,min_genes=min(r['genes'] for r in selected),
                            max_genes=max(r['genes'] for r in selected),gene_mae=np.mean(g),gene_gain=np.mean(z-g),
                            donors_gene_gain=int(np.sum(z>g)),panel_mae='' if pg is None else np.mean(pm),
                            panel_gain='' if pg is None else np.mean(pg-pm)))
    save('donors.tsv',rows)
    save('summary.tsv',summary)
    save('matched.tsv',match)
    np.savez_compressed(HERE/'reference_gains.npz',**draws_out)


if __name__=='__main__':
    main()
