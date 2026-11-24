"""Independently check the observed 84/95-gene comparisons and missing-effect bounds."""
import csv
from collections import defaultdict
from math import isclose
from pathlib import Path
from statistics import mean


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def read(path):
    with path.open(encoding='utf8') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def near(x, y):
    assert isclose(float(x),float(y),rel_tol=1e-9,abs_tol=1e-11), (x,y)


def prediction(source, method):
    if method == 'zero':
        return [0.] * len(source)
    if method == 'uniform':
        return [mean(source)] * len(source)
    return [v*(0.25 if method == 'quarter' else 1.) for v in source]


def check(y, p, r):
    mae = mean(abs(obs-pred) for obs,pred in zip(y,p))
    panel = abs(mean(y)-mean(p))
    near(mae,r['gene_mae'])
    near(panel,r['panel_abs_error'])
    near(mean(map(abs,y))-mae,r['gene_gain_vs_zero'])
    near(abs(mean(y))-panel,r['panel_gain_vs_zero'])


def main():
    hk_rows=read(WORK/'stage13/external/hongkong_target_effects.tsv')
    hk=defaultdict(dict)
    for r in hk_rows:
        if r['method']=='zero':
            hk[r['gene']][int(r['donor'])]=float(r['log2_observed'])
    chip=read(WORK/'stage13/external/chip_effects.tsv')
    group=[r for r in read(HERE/'group_scores.tsv') if r['scope']=='target']
    units=[r for r in read(HERE/'donor_scores.tsv') if r['scope']=='target']
    for r in group + units:
        y_o={x['gene']:float(x['log2_change']) for x in chip if x['rule']==r['gate'] and x['target']=='1'}
        common=sorted(y_o.keys() & hk.keys())
        assert len(common)==int(r['genes'])
        dest_osaka=r.get('direction')=='HK12_to_OS48'
        if 'hk_donor' in r:
            h=[hk[g][int(r['hk_donor'])] for g in common]
        else:
            h=[mean(hk[g].values()) for g in common]
        o=[y_o[g] for g in common]
        src,dst=(h,o) if dest_osaka else (o,h)
        check(dst, prediction(src,r['method']),r)
    bounds=read(HERE/'vienna_bounds.tsv')
    v={r['gene']:float(r['reported_log2']) for r in read(WORK/'stage12/culture/culture_readouts.tsv') if r['reported_log2']}
    genes=sorted(hk)
    for r in bounds:
        p=prediction([mean(hk[g].values()) for g in genes],r['method'])
        for side in ('min','max'):
            # Complete the missing readings so the triangle-inequality limit is attained.
            obs=[v[g] if g in v else (forecast if side=='max' else 0.) for g,forecast in zip(genes,p)]
            gain=mean(abs(y)-abs(y-f) for y,f in zip(obs,p))
            near(gain,r['full_panel_gene_gain_'+side])
    print(f'Independently verified {len(group)} group and {len(units)} donor endpoint rows and eight sharp bounds.')


if __name__ == '__main__':
    main()
