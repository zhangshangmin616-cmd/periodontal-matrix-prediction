"""Post hoc checks of two small gains when individual specimens or jaw weights change."""
import csv
from pathlib import Path
import numpy as np


HERE=Path(__file__).resolve().parent
from score import COEFFICIENTS


def main():
    with (HERE/'specimen_metadata.tsv').open(encoding='utf8') as f:
        meta=list(csv.DictReader(f,delimiter='\t'))
    with np.load(HERE/'specimen_cpm.npz',allow_pickle=False) as a:
        expression,genes=a['cpms'],a['target']
    groups=[[i for i,r in enumerate(meta) if (r['batch'],r['group']) == key]
            for key in [('4','control'),('4','day7'),('3','control'),('3','day28')]]
    alpha=COEFFICIENTS['atlas_D3_D7']

    def evaluate(label, selected, equal_arch):
        rows=[]
        for name,initial,final in [('mean_log',np.log2(expression[:,genes]+1),lambda v:v),
                                   ('log_mean',expression[:,genes],lambda v:np.log2(v+1))]:
            means=[]
            for ix in selected:
                if equal_arch:
                    jaws={meta[i]['arch'] for i in ix}
                    assert jaws=={'maxillary','mandibular'}
                    value=np.mean([initial[[i for i in ix if meta[i]['arch']==jaw]].mean(axis=0)
                                   for jaw in sorted(jaws)],axis=0)
                else:
                    value=initial[ix].mean(axis=0)
                means.append(final(value))
            x=means[1]-means[0]
            y=means[3]-means[2]
            m=np.abs(y-alpha*x).mean()
            rows.append(dict(case=label,estimand=name,model_mae=m,
                             gain_over_zero=np.abs(y).mean()-m,gain_over_quarter=np.abs(y-.25*x).mean()-m))
        return rows

    rows=evaluate('all',groups,False)+evaluate('equal_jaw',groups,True)
    for i,r in enumerate(meta):
        keep=[[j for j in ix if j!=i] for ix in groups]
        rows.extend(evaluate('omit_'+r['run'],keep,False))
    with (HERE/'sensitivity.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys(),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    for scale in ['mean_log','log_mean']:
        removed=[r for r in rows if r['case'].startswith('omit_') and r['estimand']==scale]
        for key in ['gain_over_zero','gain_over_quarter']:
            g=[r[key] for r in removed]
            jaw=next(r[key] for r in rows if r['case']=='equal_jaw' and r['estimand']==scale)
            print(scale,key,'deletions with gain',sum(x>0 for x in g),'range',min(g),max(g),'equal jaw gain',jaw)


if __name__=='__main__':
    main()
