"""Independent group replacement and tie-ranking check of the headline ablations."""
import csv
from pathlib import Path
from statistics import fmean

import numpy as np


HERE = Path(__file__).resolve().parent


def target_percentiles(x, target):
    _, group, n = np.unique(x,return_inverse=True,return_counts=True)
    midrank = np.cumsum(n)-(1+n)/2
    return midrank[group[target]]*100/(len(x)-1)


def main():
    with (HERE/'unit_sensitivity.tsv').open(encoding='utf8') as f:
        units = {(r['unit'],r['model'],r['metric']):r for r in csv.DictReader(f,delimiter='\t')}
    with np.load(HERE/'inputs.npz',allow_pickle=False) as z:
        data={k:z[k] for k in z.files}
    base,pressure,tension = [data[k] for k in ('base','pressure','tension')]
    classes=data['classes'][:-1].tolist()
    count=0
    for i, unit in enumerate(data['units']):
        ys=data['observed'][i]
        selectors={'full':classes,'zero':[]}
        for group in classes:
            selectors['only_'+group]=[group]
            selectors['omit_'+group]=[k for k in classes if k!=group]
        where=np.flatnonzero(data['unit_per_control']==unit)
        for name,chosen in selectors.items():
            shifts=[]
            for row in where:
                q=data['q'][row]
                x=data['control_cpm'][row]
                subsets=[k in chosen for k in data['classes']]
                current=np.zeros_like(x)
                start=np.zeros_like(x)
                for j,use in enumerate(subsets):
                    start += q[j]*base[:,j]
                    current += q[j]*((pressure[:,j]+tension[:,j])*.5 if use else base[:,j])
                prediction=x*((current+1)/(start+1))
                shifts.append(target_percentiles(prediction,data['target'])-target_percentiles(x,data['target']))
            estimated=np.mean(shifts,axis=0)
            for metric in ['panel','gene']:
                prior=units[unit,name,metric]
                error=abs(fmean(estimated)-fmean(ys)) if metric=='panel' else fmean(abs(estimated-ys))
                np.testing.assert_allclose([fmean(estimated),error],[float(prior['predicted']),float(prior['mae'])],
                                           rtol=0,atol=1e-8)
                count+=1
    print('Verified',count,'per-unit results using independent class replacement and midranks.')


if __name__=='__main__':
    main()
