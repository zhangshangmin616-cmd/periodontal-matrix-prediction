"""Partial identification using only the published direction of TGFB1 on each side."""
import csv
import gzip
import numpy as np
from score import HERE, WORK, write


def main():
    with np.load(WORK/"stage11/specificity/rank_changes.npz",allow_pickle=False) as q:
        panel=q["gene"][q["component"]=="other"]
    with gzip.open(WORK/"stage10/state_multipliers.tsv.gz","rt",encoding="utf8") as f:
        pred={r["gene"]:r for r in csv.DictReader(f,delimiter="\t") if r["gene"] in panel}
    assert len(pred)==len(panel)==101
    coeff={}
    for side in ("PP7D","TP7D"):
        coeff[side]=np.array([float(pred[g][side])-float(pred["ACTB"][side]) for g in panel])
    outputs=[]
    for contrast,x,known_up in [("pressure_vs_untreated",coeff["PP7D"],True),
                                ("tension_vs_untreated",coeff["TP7D"],True),
                                ("pressure_vs_tension",coeff["PP7D"]-coeff["TP7D"],False)]:
        lower=-np.abs(x)
        upper=np.abs(x)
        target=list(panel).index("TGFB1")
        # y>0 permits the full gain interval for a positive forecast; when p<0
        # this single gene's gain equals p for every value in that half-line.
        if known_up and x[target]<0:
            upper[target]=x[target]
        assert x[list(panel).index("ACTB")]==0.
        outputs.append(dict(contrast=contrast,targets=len(x),biological_targets_with_reported_assay=1,
                            target_predictions_relative_to="ACTB",
                            TGFB1_prediction=x[target],TGFB1_min_gain=lower[target],TGFB1_max_gain=upper[target],
                            full101_min_gain=np.mean(lower),full101_max_gain=np.mean(upper)))
    write(HERE/"bounds.tsv",outputs)


if __name__=="__main__":
    main()
