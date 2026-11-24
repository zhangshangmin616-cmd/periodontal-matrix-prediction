"""Rebuild 24 h selected targets from the public XLS (optional) and score all sources."""
import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

HERE = Path(__file__).resolve().parent
SHA256 = "1f13a1c05d167e44735a9286474f4ea96274bbb3287084005b97f4210d599d84"


def write_table(name, rows):
    with (HERE / name).open("w", newline="", encoding="utf8") as handle:
        out = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def read_table(path):
    with path.open(encoding="utf8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract(workbook, genes):
    import xlrd
    if hashlib.sha256(workbook.read_bytes()).hexdigest() != SHA256:
        raise ValueError("Public supplement checksum differs")
    sheet = xlrd.open_workbook(str(workbook)).sheet_by_name("B.after 24h compression")
    rows = []
    for i in range(2, sheet.nrows):
        for column, direction in [(0, "up"), (3, "down")]:
            gene = str(sheet.cell_value(i, column)).strip()
            if gene in genes:
                signed_fold = float(sheet.cell_value(i, column + 1))
                assert (signed_fold > 0) == (direction == "up") and abs(signed_fold) >= 1
                rows.append(dict(gene=gene, xls_row=i+1, block=direction, signed_fold=signed_fold))
    write_table("reported_targets.tsv", rows)


def evaluate(source, rule, predicted, observed):
    covered = list(observed)
    n = len(predicted)
    assert n == 101
    zero = mean(abs(y) for y in observed.values())
    mae = mean(abs(observed[g] - predicted[g]) for g in covered)
    reported_gain = (zero - mae) * len(covered)
    unknown_reach = sum(abs(p) for g,p in predicted.items() if g not in observed)
    return dict(source=source, rule=rule, reported_n=len(covered), total_n=n,
                reported_mae=mae, reported_zero_mae=zero, reported_gain=zero-mae,
                full_gain_lower=(reported_gain-unknown_reach)/n,
                full_gain_upper=(reported_gain+unknown_reach)/n)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, help="Re-extract matched rows using xlrd 2.0.2")
    args=parser.parse_args()
    sources=defaultdict(dict)
    for row in read_table(HERE.parent / "target_effects.tsv"):
        if row["scale"] == "log2":
            sources[row["study"]][row["gene"]] = float(row["observed"])
    genes = set.intersection(*(set(data) for data in sources.values()))
    assert len(sources) == 3 and len(genes) == 101
    if args.workbook:
        extract(args.workbook, genes)
    logs=defaultdict(list)
    for row in read_table(HERE / "reported_targets.tsv"):
        fold = float(row["signed_fold"])
        logs[row["gene"]].append(math.copysign(math.log2(abs(fold)),fold))
    observed={g:mean(values) for g,values in logs.items()}
    sources["equal_mean_three"] = {g:mean(src[g] for src in sources.values()) for g in sorted(genes)}
    records = [evaluate("zero", "zero", dict.fromkeys(genes, 0), observed)]
    for name, source in sources.items():
        level = mean(source.values())
        for weight in [1, 0.25]:
            gene_specific = {g:y*weight for g,y in source.items()}
            scalar = dict.fromkeys(genes, level*weight)
            records.append(evaluate(name, f"gene_x{weight}", gene_specific, observed))
            records.append(evaluate(name, f"uniform_x{weight}", scalar, observed))
    write_table("scores.tsv", records)


if __name__ == "__main__":
    main()
