"""RNA checks with an explicit definition of the measurement for each comparison."""
from pathlib import Path
import csv
import gzip
import itertools
import json
import math
from statistics import fmean

root=Path(__file__).resolve().parent

def save(rows,path):
    with (root/path).open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,rows[0],delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)
def main():
    with gzip.open(root/'state_multipliers.tsv.gz','rt') as f:
        pred={r['gene']:{k:float(r[k]) for k in ['PP7D','TP7D']} for r in csv.DictReader(f,delimiter='\t')}
    side=[('pressure','PP7D'),('tension','TP7D')]
    # Between-side transcriptional contrast is on the article's ACTB-normalized scale.
    # Reported non-significant differences have unknown numerical direction.
    direction={'TGFB1':None,'TNF':1,'IL10':-1,'COL1A1':-1,'MMP1':1,'TIMP1':-1,'BGLAP':-1,'TNFSF11':1,'TNFRSF11B':-1}
    rows=[]
    def entry(g,comparison,value,sign):
        return dict(gene=g,comparison=comparison,prediction=value,known_direction=sign,
            sign_check='unresolved' if sign is None else 'correct' if sign*value>0 else 'opposite' if sign*value<0 else 'zero',
            # Excess MAE bounds from knowing only the sign of a nonzero observation.
            smallest_excess=-sign*value if sign is not None else -abs(value),largest_excess=abs(value))
    for g,s in direction.items():
        normalized=[pred[g][lib]-pred['ACTB'][lib] for _,lib in side]
        rows.append(entry(g,'pressure_minus_tension',normalized[0]-normalized[1],s))
        for (name,_),value in zip(side,normalized,strict=True):
            observed_sign=None if name=='pressure' and g in ('COL1A1','TIMP1') else 1
            rows.append(entry(g,name+'_minus_control',value,observed_sign))
    save(rows,'garlet07.tsv')
    checks={}
    for name,subset in [('side_difference',[r for r in rows if r['comparison']=='pressure_minus_tension']),
                        ('sides_vs_control',[r for r in rows if r['comparison']!='pressure_minus_tension'])]:
        checks[name]={'total':len(subset),'correct':sum(r['sign_check']=='correct' for r in subset),
            'opposite':sum(r['sign_check']=='opposite' for r in subset),'unresolved':sum(r['sign_check']=='unresolved' for r in subset),
            'excess_range_from_marginal_signs':[fmean(r[key] for r in subset) for key in ('smallest_excess','largest_excess')]}
    # Figure 4B, Otero. Audit thick horizontal lines against the downloaded PDF.
    axis=(479,485);c4=(424,430);t4=(450,456)
    def ratio(z,c,t):return math.log2((z-c)/(z-t))
    obs=ratio(fmean(axis),fmean(c4),fmean(t4))
    spread=[ratio(*v) for v in itertools.product(*[[a-.5,b+.5] for a,b in (axis,c4,t4)])]
    p=(pred['TNFSF11']['PP7D']-pred['TNFRSF11B']['PP7D'])-(pred['TNFSF11']['TP7D']-pred['TNFRSF11B']['TP7D'])
    assert min(spread)>p>0
    checks['Otero_Fig4B']={'dose_oz':4,'observed_log2_C_over_T_ratio_marks':obs,'pixel_interval':[min(spread),max(spread)],
        'prediction':p,'model_error':abs(obs-p),'zero_error':abs(obs),'gain':obs-abs(obs-p),
        'dose7':'tension mark coincides with baseline; log ratio unresolved'}
    (root/'clinical_summary.json').write_text(json.dumps(checks,indent=2)+'\n')
    # Supplementary non-site-specific samples. Mix expression before GAPDH calibration.
    zhu=json.loads((root/'raw/measurement_audit.json').read_text())['measurements']
    other=[]
    for g,measurement in zhu.items():
        y=measurement['approximate_log2_ratio_of_plotted_group_bars']
        def forecast(a):
            num=(1-a)*2**pred[g]['TP7D']+a*2**pred[g]['PP7D']
            den=(1-a)*2**pred['GAPDH']['TP7D']+a*2**pred['GAPDH']['PP7D']
            return math.log2(num/den)
        vals=[forecast(a) for a in [0,.5,1]]
        other.append(dict(gene=g,observed=y,pressure_0=vals[0],pressure_half=vals[1],pressure_1=vals[2],
            model_error_half=abs(vals[1]-y),zero_error=abs(y),closer_half='state' if abs(vals[1]-y)<abs(y) else 'zero'))
    save(other,'zhu_all_targets.tsv')
    print(checks,other)
if __name__=='__main__':main()
