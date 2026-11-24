"""Recalculate cross-study clinical transfer after individual one-unit deletions."""
import csv
import gzip

import numpy as np

from score import WORK, average_ranks, write


def read_rows(path):
    with path.open(encoding="utf8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def clinical_ranks(genes, targets):
    root = WORK / "stage5/zhao"
    by_symbol = {row["gene"]:row["entrez_id"] for row in read_rows(root / "feature_map.tsv")}
    entrez = [by_symbol[gene] for gene in genes]
    assert len(entrez) == len(set(entrez)) == 10806
    groups = {}
    for row in read_rows(root / "joined_manifest.tsv"):
        run = row["run"]
        with gzip.open(root / "counts" / (run + ".tsv.gz"), "rt", encoding="utf8") as stream:
            records = csv.reader(stream, delimiter="\t")
            header = next(records)
            assert header[0] == "Geneid" and len(header) == 2
            counts = {key:float(value) for key,value in records}
        values = np.array([counts[key] for key in entrez])
        ranks = (average_ranks(values) - 1) * (100 / (len(genes) - 1))
        # The CPM factor is constant within a specimen and preserves all ranks.
        groups.setdefault((row["batch"], row["group"]), []).append((run, ranks[targets]))
    return groups


def observed_leuven(groups, day):
    batch, treated = {"day7":("4","day7"), "day28":("3","day28")}[day]
    arms = {"control":groups[batch,"control"], "exposed":groups[batch,treated]}
    means = {arm:np.mean([x for _,x in records],axis=0) for arm,records in arms.items()}
    original = means["exposed"] - means["control"]
    dropped = []
    for arm, records in arms.items():
        sign = 1 if arm == "exposed" else -1
        for run, x in records:
            effect = original + sign * (means[arm] - x) / (len(records) - 1)
            dropped.append(("Leuven_" + arm + "_specimen", run, effect))
    return original, dropped


def main():
    with np.load(WORK / "stage11/specificity/rank_changes.npz", allow_pickle=False) as d:
        genes, units = d["gene"], d["units"]
        target = d["component"] == "other"
        y = d["observed"][:,target]
    groups = clinical_ranks(genes, target)
    seoul_rows = np.flatnonzero(np.char.startswith(units, "Kim_"))
    seoul = y[seoul_rows]
    seoul_drop = [(units[i], np.delete(seoul, j, axis=0)) for j,i in enumerate(seoul_rows)]
    detail = []
    def evaluate(source,destination,dataset,unit,prediction,observed):
        for model, multiplier in [("clinical_full",1.), ("clinical_quarter",.25)]:
            errors = np.abs(observed - prediction * multiplier)
            detail.append(dict(source=source,destination=destination,model=model,exclusion_type=dataset,
                               removed=unit,gene_mae=errors.mean(),gene_gain_vs_zero=np.abs(observed).mean()-errors.mean()))
    for day in ("day7", "day28"):
        label = "Leuven_" + day
        leuven, leuven_drop = observed_leuven(groups,day)
        np.testing.assert_allclose(leuven, y[np.flatnonzero(units=="Zhao_"+day)[0]], atol=1e-10,rtol=0)
        for who, rest in seoul_drop:
            evaluate("Seoul_week3",label,"Seoul_source_participant",who,rest.mean(axis=0),leuven)
        for arm, who, rest in leuven_drop:
            evaluate("Seoul_week3",label,arm,who,seoul.mean(axis=0),rest)
        if day == "day7":
            for who, rest in seoul_drop:
                evaluate(label,"Seoul_week3","Seoul_destination_participant",who,leuven,rest)
            for arm, who, rest in leuven_drop:
                evaluate(label,"Seoul_week3",arm,who,rest,seoul)
    factors = ["source","destination","model","exclusion_type"]
    keys = sorted({tuple(row[k] for k in factors) for row in detail})
    summary = []
    for key in keys:
        rows = [r for r in detail if tuple(r[k] for k in factors)==key]
        gains = np.array([r["gene_gain_vs_zero"] for r in rows])
        summary.append(dict(**dict(zip(factors,key)),deletions=len(gains),improved=int((gains>0).sum()),
                            worst_gain=gains.min(), best_gain=gains.max()))
    write("deletions_by_unit.tsv",detail)
    write("deletions_summary.tsv",summary)


if __name__ == "__main__":
    main()
