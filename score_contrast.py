"""Score a prediction for one contrast of the article's fixed 101-gene panel."""
import argparse
import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXAMPLE = ROOT / "stage19/figures/Figure_1_source_data.csv"


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as file:
        return list(csv.DictReader(file))


def metrics(predicted, observed):
    residuals = [p - y for p, y in zip(predicted, observed, strict=True)]
    return abs(math.fsum(residuals) / len(residuals)), math.fsum(map(abs, residuals)) / len(residuals)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, help="CSV with gene, observed and a forecast column, one contrast per file")
    parser.add_argument("--forecast", default="predicted", help="Column containing the predicted rank change")
    args = parser.parse_args()
    try:
        rows = read_rows(args.file)
        genes = [r["gene"].strip() for r in rows]
        panel = {r["gene"] for r in read_rows(EXAMPLE)}
        if len(genes) != 101 or len(panel) != 101 or set(genes) != panel:
            raise ValueError("Use exactly one row per gene from the fixed 101-gene panel.")
        y = [float(r["observed"]) for r in rows]
        p = [float(r[args.forecast]) for r in rows]
        if not all(map(math.isfinite, y + p)):
            raise ValueError("All observed and predicted values must be finite.")
    except (KeyError, ValueError, OSError) as exc:
        parser.error(str(exc))
    print("forecast\tpanel_error_pp\tmean_gene_error_pp")
    for rule, values in (("zero", [0.0] * len(y)), (args.forecast, p)):
        panel_error, gene_error = metrics(values, y)
        print(f"{rule}\t{panel_error:.6f}\t{gene_error:.6f}")
    closer = sum(abs(a - b) < abs(b) for a, b in zip(p, y, strict=True))
    farther = sum(abs(a - b) > abs(b) for a, b in zip(p, y, strict=True))
    print(f"Genes closer than zero: {closer}; farther: {farther}; ties: {len(y)-closer-farther}.")


if __name__ == "__main__":
    main()
