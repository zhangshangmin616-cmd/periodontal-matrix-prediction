"""Recover baseline rank bounds for the 13 clinical comparisons."""
import csv
from pathlib import Path
import sys

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
STAGE5 = HERE.parents[1]/'stage5'
sys.path.insert(0, str(STAGE5))
from forecast import Reference, load_bulk, meta
from external import layout, one_run


def main():
    ref=Reference()
    info=meta()
    info=info[info.arm=='control']
    bulk=load_bulk(info.index.tolist(),ref.genes)
    inputs={}
    def panel_ranks(x):
        return (rankdata(x,method='average')[ref.idx]-1)*100/(len(x)-1)
    for donor, entries in info.groupby('patient'):
        inputs[f'Kim_{donor}']=np.mean([panel_ranks(bulk[s].to_numpy()) for s in entries.index],axis=0)
    runs, genes=layout(ref)
    for batch, day in [(4,7),(3,28)]:
        entries=runs[(runs['group']=='control') & (runs.batch==batch)]
        inputs[f'Zhao_day{day}']=np.mean([panel_ranks(one_run(s,genes)) for s in entries.index],axis=0)
    with (HERE/'control_ranks.tsv').open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['unit','gene','untreated_rank'])
        for unit, ranks in inputs.items():
            for g,x in zip(ref.genes[ref.idx],ranks,strict=True):
                w.writerow([unit,g,x])


if __name__=='__main__':
    main()
