"""Independent check of model errors and the sign pattern in matched reference panels."""
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def read(path):
    with path.open(encoding='utf8', newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def main():
    with np.load(WORK/'stage13/external/paired_counts.npz') as source:
        cpm = source['counts'].astype(float)/source['total']*1e6
        genes, target = source['genes'], source['target']
    reference = {r['gene']:r['ecm_component'] for r in read(WORK/'stage5/reference/features.tsv')}
    outside = np.array([reference[g]=='outside_ecm' for g in genes])
    matrix = np.empty(cpm.shape)
    for k, column in enumerate(cpm.T):
        _,inv,n = np.unique(column,return_inverse=True,return_counts=True)
        matrix[:,k] = 100 * (np.cumsum(n) - (n+1)/2)[inv]/(len(n[inv])-1)
    table = {(r['scale'],r['scope'],int(r['donor']),r['method']): r for r in read(HERE/'donors.tsv')}
    checkrows = 0
    for scale,m in [('rank',matrix),('log2',np.log2(1+cpm))]:
        errors = []
        for i in range(5):
            train = [j for j in range(5) if j!=i]
            ys = m[:,5+i]-m[:,i]
            prediction = np.mean([m[:,5+j]-m[:,j] for j in train],axis=0)
            masks = {'panel':target,'panel_control_3of5':target & (np.sum(cpm[:,:5]>=1,axis=1)>=3),
                     'panel_training_control_3of4':target & (np.sum(cpm[:,train]>=1,axis=1)>=3),
                     'background':np.ones(len(genes),bool), 'outside_ecm':outside}
            for scope,ix in masks.items():
                for method,p in [('zero', np.zeros(np.sum(ix))),('training_gene_mean',prediction[ix]),
                                 ('training_gene_quarter',.25*prediction[ix]),('uniform_training_mean',prediction[ix].mean())]:
                    vals = table[scale,scope,i+1,method]
                    obs = ys[ix]
                    assert int(vals['genes']) == len(obs)
                    assert abs(np.mean(abs(obs-p))-float(vals['gene_mae'])) < 1e-10
                    if vals['panel_mae']:
                        assert abs(abs(np.mean(obs)-np.mean(p))-float(vals['panel_mae'])) < 1e-10
                    checkrows += 1
        if scale=='rank':
            delta = m[:,5:]-m[:,:5]
            fit = np.column_stack([np.mean(np.delete(delta,i,axis=1),axis=1) for i in range(5)])
            # Rebuild reference memberships without the main script's stratum assignment shortcut.
            means = np.array([np.mean(x) for x in cpm[:,:5]])
            edges = np.quantile(means,np.arange(1,10)/10)
            strata = np.sum(means[:,None]>=edges[None,:],axis=1)
            sampler = np.random.default_rng(20261124)
            panels = np.empty((10000,101),dtype=np.int32)
            start = 0
            for s in range(10):
                width = np.count_nonzero(target & (strata==s))
                avail = np.flatnonzero(outside & (strata==s))
                if width:
                    for dst in panels:
                        dst[start:start+width] = sampler.choice(avail,size=width,replace=False)
                start+=width
            assert start==101 and outside[panels].all()
            with np.load(HERE/'reference_gains.npz') as original:
                # All 10,000 reference sets checked, one donor at a time.
                benefit = np.zeros(10000)
                mean_benefit = np.zeros(10000)
                for donor in range(5):
                    target_values = delta[panels,donor]
                    fitted = fit[panels,donor]
                    benefit += (abs(target_values).mean(1) - abs(target_values-fitted).mean(1))/5
                    mean_benefit += (abs(target_values.mean(1)) - abs(target_values.mean(1)-fitted.mean(1)))/5
                np.testing.assert_allclose(benefit,original['rank_gene_gain'], atol=1e-10)
                np.testing.assert_allclose(mean_benefit,original['rank_panel_gain'], atol=1e-10)
    assert checkrows==200
    print('Checked 200 donor-specific model results and all 10,000 reference panels on the rank scale.')


if __name__ == '__main__':
    main()
