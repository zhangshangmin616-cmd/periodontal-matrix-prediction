"""Check external fold anchors and score paired, probe-keyed public array exports."""
import csv
import gzip
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def read(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf8',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def write(name, rows):
    with (HERE/name).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def average_ranks(x):
    sorted_x=np.sort(x)
    return .5*(sorted_x.searchsorted(x,'left')+sorted_x.searchsorted(x,'right'))


def rho(x,y):
    return float(np.corrcoef(average_ranks(x),average_ranks(y))[0,1])


def main():
    records=read(HERE/'probes.tsv.gz')
    assert len({r['probe'] for r in records})==len(records)==30336
    x=np.array([[float(r[k]) for k in ['logA','logB','ch1A','ch2A','ch1B','ch2B','depositedA','depositedB']] for r in records])
    assert np.isfinite(x).all()
    keys=[r['symbol_by_entrez'] for r in records]
    symbol=np.array([r['platform_symbol'] for r in records])
    valid_signal=(x[:,2:6]>0).all(axis=1)
    anchors=[]
    for r in read(HERE/'figure1a.tsv'):
        p=np.where(symbol==r['symbol'])[0]
        if len(p):
            ratio=float(2**(x[p,:2].mean()))
            anchors.append(dict(symbol=r['symbol'],platform_probes=len(p),reported=r['printed_fold'],
                                scanner_fold=ratio,within_one_percent=int(abs(ratio/float(r['printed_fold'])-1)<.01)))
        else:
            anchors.append(dict(symbol=r['symbol'],platform_probes=0,reported=r['printed_fold'],
                                scanner_fold='',within_one_percent=''))
    write('orientation.tsv',anchors)
    assert np.count_nonzero(x[:,0]==x[:,6])==963 and np.array_equal(x[:,1],x[:,7])
    # Anchors are independent of the fixed panel and rule. Record multiple/no-probe symbols too.
    assert sum(r['within_one_percent']==1 for r in anchors if r['platform_probes']==1)>=3
    audit=[]
    genes=[]
    forecasts=[]
    with np.load(WORK/'stage13/external/paired_counts.npz') as data:
        hong=data['genes'].astype(str)
        is_target=dict(zip(hong,data['target'].tolist()))
        abund=np.log2(1+data['counts']*(1.e6/data['total']))
        paired=(abund[:,5:]-abund[:,:5])
        hk=dict(zip(hong,paired.mean(axis=1)))
    for rule,mask in [('positive_both_scans',valid_signal),('all_mapped_probes',np.ones(len(x),dtype=bool))]:
        audit.append(dict(rule=rule,platform_probes=int(mask.sum()),rho_between_scans=rho(x[mask,0],x[mask,1])))
        grouped={}
        for i in np.where(mask)[0]:
            if keys[i] in hk:
                grouped.setdefault(keys[i],[]).append(i)
        group_rows=[]
        for key,ix in sorted(grouped.items()):
            a,b=x[ix,:2].mean(axis=0)
            group_rows.append(dict(rule=rule,gene=key,target=int(is_target[key]),probes=len(ix),
                                   hk=hk[key],array_A=float(a),array_B=float(b),both=float((a+b)/2)))
        genes.extend(group_rows)
        for scope,selected in [('targets',[r for r in group_rows if r['target']]),('background',group_rows)]:
            h=np.array([r['hk'] for r in selected])
            for arr in ('array_A','array_B','both'):
                o=np.array([r[arr] for r in selected])
                for direction,source,dest in [('HongKong_to_Osaka',h,o),('Osaka_to_HongKong',o,h)]:
                    for label,p in [('zero',np.zeros_like(source)),('uniform',np.full_like(source,source.mean())),
                                    ('full',source),('quarter',source*.25)]:
                        err=p-dest
                        forecasts.append(dict(rule=rule,scope=scope,array=arr,direction=direction,
                                              method=label,genes=len(dest),rho=rho(source,dest),
                                              gene_mae=float(np.abs(err).mean()),panel_mae=abs(float(err.mean())),
                                              gain_vs_zero=float(np.mean(abs(dest))-np.mean(abs(err)))))
    write('scanner_qc.tsv',audit)
    write('gene_effects.tsv',genes)
    write('scores.tsv',forecasts)


if __name__=='__main__':
    main()
