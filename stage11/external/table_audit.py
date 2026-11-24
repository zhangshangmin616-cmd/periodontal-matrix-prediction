"""Compare same-day I-E coefficients with the arithmetic difference of I-C and E-C."""
from pathlib import Path
import argparse
import json
import numpy as np
from score_direct import write


HERE = Path(__file__).resolve().parent


def table(file):
    blob = json.loads(file.read_bytes())
    assert blob['header'][1:4] == ['gene_name','baseMean','log2FoldChange']
    return {r[1]:r for r in blob['rows']}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--work', type=Path, default=HERE.parents[1])
    args = parser.parse_args()
    rows = []
    direct = {}
    for day in (1,3,10):
        c0dir = args.work / ('stage9/day1' if day==1 else 'stage8/hotm')
        direction = table(HERE/'source'/f'deg_E{day}_I{day}.json')
        direct[day] = direction
        c0i = table(c0dir/f'deg_C0_I{day}.json')
        c0e = table(c0dir/f'deg_C0_E{day}.json')
        shared = sorted(set(direction) & set(c0i) & set(c0e))
        y,a,b = [np.asarray([d[g][3] for g in shared]) for d in [direction,c0i,c0e]]
        subtract = a-b
        same = sum(direction[g][2] == c0i[g][2] == c0e[g][2] for g in shared)
        assert same == len(shared)
        rows.append(dict(day=day,genes=len(shared),exactly_equal_basemean=same,
                         same_direction=int(np.sum(np.sign(y)==np.sign(subtract))),
                         direct_median_abs=float(np.median(abs(y))),
                         indirect_median_abs=float(np.median(abs(subtract))),
                         median_abs_discrepancy=float(np.median(abs(y-subtract))),
                         max_abs_discrepancy=float(max(abs(y-subtract))),
                         correlation=float(np.corrcoef(y,subtract)[0,1])))
    across = set.intersection(*(set(direct[d]) for d in direct))
    assert len(across)==14981 and all(direct[1][g][2]==direct[3][g][2]==direct[10][g][2] for g in across)
    write(HERE/'table_audit.tsv',rows)


if __name__=='__main__':
    main()
