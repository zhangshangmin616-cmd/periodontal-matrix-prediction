"""Calculate conditional paired bootstrap intervals without simulation noise."""
import csv
import itertools
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    with (HERE/"scores.tsv").open(encoding="utf8",newline="") as f:
        rows=[r for r in csv.DictReader(f,delimiter="\t") if r["unit"].startswith("pair")]
    n=len(rows)
    assert n==10
    gains=np.array([[float(r[k]) for k in ("panel_gain","gene_gain")] for r in rows])
    outcomes=[]
    multiplicities=[]
    for bars in itertools.combinations(range(2*n-1),n-1):
        counts=np.diff([-1,*bars,2*n-1])-1
        outcomes.append(counts@gains/n)
        multiplicities.append(math.factorial(n)//math.prod(math.factorial(c) for c in counts))
    assert sum(multiplicities)==n**n and len(outcomes)==math.comb(2*n-1,n-1)
    outcomes=np.array(outcomes)
    mass=np.asarray(multiplicities,dtype="int64")
    with (HERE/"uncertainty.tsv").open("w",encoding="utf8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["endpoint","n_matched_units","mean_gain","lower_pct_2.5","upper_pct_97.5","resampling"])
        for i,endpoint in enumerate(["panel_mean","mean_individual_gene"]):
            order=np.argsort(outcomes[:,i])
            cum=np.cumsum(mass[order])
            lo,hi=np.searchsorted(cum,[n**n//40,n**n*39//40])
            w.writerow([endpoint,n,gains[:,i].mean(),outcomes[order[lo],i],outcomes[order[hi],i],"ten paired units, predictions held fixed"])


if __name__=="__main__":
    main()
