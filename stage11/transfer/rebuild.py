"""Extract four study groups of whole-run CPM values in the common gene universe."""
import csv
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE5=HERE.parents[1]/'stage5'
sys.path.insert(0,str(STAGE5))
from forecast import Reference
from external import layout,one_run


def main():
    ref=Reference()
    runs, ids=layout(ref)
    cpms=[]
    rows=[]
    for run,r in runs.iterrows():
        cpms.append(one_run(run,ids))
        rows.append({'run':run,'biosample':r['sample'],'group':r['group'],'batch':r['batch'],'arch':r['arch']})
    selected=np.zeros(len(ref.genes),dtype=bool)
    selected[ref.idx]=True
    np.savez_compressed(HERE/'specimen_cpm.npz',genes=np.array(ref.genes,dtype=str),target=selected,
                        runs=np.array(runs.index,dtype=str),cpms=np.stack(cpms))
    with (HERE/'specimen_metadata.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys(),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print('Extracted',len(rows),'samples and',len(ref.genes),'genes.')


if __name__=='__main__':
    main()
