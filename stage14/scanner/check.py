"""Independently recompute the joint-signal target results without vectorized scoring."""
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path
from statistics import fmean

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def table(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf8',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def main():
    targets={r['gene'] for r in table(WORK/'stage5/reference/features.tsv') if r['ecm_component']=='other'}
    hk=defaultdict(list)
    for r in table(WORK/'stage13/external/hongkong_target_effects.tsv'):
        if r['method']=='zero':
            hk[r['gene']].append(float(r['log2_observed']))
    assert set(hk)==targets and all(len(v)==5 for v in hk.values())
    o=defaultdict(list)
    for r in table(HERE/'probes.tsv.gz'):
        if r['symbol_by_entrez'] in targets and min(float(r[k]) for k in ['ch1A','ch2A','ch1B','ch2B'])>0:
            o[r['symbol_by_entrez']].extend([float(r['logA']),float(r['logB'])])
    genes=sorted(o)
    assert len(genes)==78
    means=[list(map(fmean,[data[g] for g in genes])) for data in [hk,o]]
    look={ (r['direction'],r['method']):r for r in table(HERE/'scores.tsv') if
          r['array']=='both' and r['scope']=='targets' and r['rule']=='positive_both_scans'}
    for direction, src, dst in [('HongKong_to_Osaka',*means),('Osaka_to_HongKong',*reversed(means))]:
        for method, p in [('full',src),('quarter',[.25*q for q in src]),('uniform',[fmean(src)]*len(src)),('zero',[0.]*len(src))]:
            residual=[q-y for q,y in zip(p,dst)]
            row=look[direction,method]
            assert math.isclose(fmean(abs(a) for a in residual),float(row['gene_mae']),abs_tol=1e-12)
            assert math.isclose(abs(fmean(residual)),float(row['panel_mae']),abs_tol=1e-12)
    print('Verified the two raw-array co-detection forecasts in both directions.')


if __name__=='__main__':
    main()
