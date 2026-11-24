"""Score the probe-reconciled Osaka 48-hour tensile culture array."""
import csv
import numpy as np
from score import HERE, atlas, save


def main():
    with (HERE/"chip_effects.tsv").open(encoding="utf8") as stream:
        rows = list(csv.DictReader(stream,delimiter="\t"))
    out = []
    for rule in ("positive_signal","all_published_values"):
        usable = [r for r in rows if r["rule"]==rule]
        target = np.array([r["target"]=="1" for r in usable])
        y = np.array([float(r["log2_change"]) for r in usable])
        for name, prediction in atlas([r["gene"] for r in usable]).items():
            for scope,sel in (("target",target),("background",np.ones(len(y),dtype=bool))):
                err = prediction[sel]-y[sel]
                out.append(dict(sample="GSM474198",rule=rule,scope=scope,method=name,
                             genes=len(err),gene_mae=np.mean(abs(err)),panel_abs_error=abs(np.mean(err))))
    save("chip_scores.tsv",out)


if __name__ == "__main__":
    main()
