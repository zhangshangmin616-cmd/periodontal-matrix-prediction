"""Transfer each study's observed mean vector into each of the other studies."""
import csv
from itertools import permutations
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def save(name, rows):
    with (HERE/name).open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def ranks(values):
    ordered = np.sort(values)
    rank0 = (np.searchsorted(ordered,values,side='left') + np.searchsorted(ordered,values,side='right') - 1)/2
    return 100*rank0/(len(values)-1)


def main():
    with np.load(HERE/'counts.npz') as src:
        gene, target, study, arm = [src[key] for key in ('gene','target','accession','arm')]
        counts, totals = src['count'], src['library_total']
    cpm = counts*1e6/totals[:,None]
    summaries, effects = [], []
    for scale,samples in [('rank',np.stack([ranks(c) for c in counts])),('log2',np.log2(1+cpm))]:
        delta = {}
        for s in np.unique(study):
            x = samples[(study==s)&(arm=='force')].mean(axis=0) - samples[(study==s)&(arm=='control')].mean(axis=0)
            delta[s] = x
            effects.extend(dict(study=s,scale=scale,gene=g,observed=y) for g,y in zip(gene[target],x[target]))
        for s,d in permutations(delta,2):
            for scope,mask in [('panel101',target),('all10806',np.ones(len(gene),dtype=bool))]:
                x,y = delta[s][mask],delta[d][mask]
                scalar = np.full_like(x,x.mean())
                for rule,p in [('zero',np.zeros_like(x)),('source_gene_full',x),('source_gene_quarter',.25*x),
                               ('source_uniform_full',scalar),('source_uniform_quarter',.25*scalar)]:
                    err = np.mean(abs(y-p))
                    mean_err=abs(np.mean(y-p))
                    summaries.append(dict(source=s,destination=d,scale=scale,scope=scope,genes=len(x),rule=rule,
                                          source_mean=x.mean(),destination_mean=y.mean(),gene_mae=err,
                                          panel_abs_error=mean_err,gene_gain_vs_zero=np.mean(abs(y))-err,
                                          panel_gain_vs_zero=abs(y.mean())-mean_err))
    save('scores.tsv',summaries)
    save('target_effects.tsv',effects)


if __name__=='__main__':
    main()
