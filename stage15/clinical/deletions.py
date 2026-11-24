"""Influence of loaded Leuven specimens and Seoul participants on attainable gains."""
import importlib.util
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
STAGES = HERE.parents[1]
spec = importlib.util.spec_from_file_location("magnitude_analysis", HERE/"score.py")
magnitude = importlib.util.module_from_spec(spec)
spec.loader.exec_module(magnitude)
sys.path.insert(0,str(STAGES/"stage13/transfer"))
from deletions import clinical_ranks, observed_leuven


def evaluate(destination, view, removed, p, y):
    rows = []
    for family,m in [("amplitude",np.zeros_like(p)),
                     ("preserve_mean",p.mean(axis=1,keepdims=True))]:
        slope = max(0., magnitude.identify(p-m, y-m))
        forecast = m + slope*(p-m)
        mae = np.abs(y-forecast).mean()
        rows.append(dict(destination=destination, estimand=view, removed=removed,
            family=family, nonnegative_optimum=slope, gene_mae=mae,
            gain_vs_zero=np.abs(y).mean()-mae, gain_vs_anchor=np.abs(y-m).mean()-mae))
    return rows


def main():
    with np.load(STAGES/"stage11/specificity/rank_changes.npz") as z:
        genes, units, target = z["gene"],z["units"],z["component"]=="other"
        p,y = z["fitted"][:,target],z["observed"][:,target]
    ranks = clinical_ranks(genes,target)
    rows = []
    for day in ("day7","day28"):
        at = units == "Zhao_"+day
        full,removed = observed_leuven(ranks,day)
        np.testing.assert_allclose(full,y[at][0],rtol=0,atol=1e-10)
        for role,run,contrast in removed:
            if role == "Leuven_exposed_specimen":
                rows.extend(evaluate("Leuven_"+day,"group_mean",run,p[at],contrast[None,:]))
    seoul = np.flatnonzero(np.char.startswith(units,"Kim_"))
    for i in seoul:
        kept = seoul[seoul!=i]
        rows.extend(evaluate("Seoul_week3","participant_first",units[i],p[kept],y[kept]))
        rows.extend(evaluate("Seoul_week3","group_mean",units[i],p[kept].mean(axis=0,keepdims=True),
                             y[kept].mean(axis=0,keepdims=True)))
    magnitude.write("deletions.tsv",rows)


if __name__ == "__main__":
    main()
