"""Verify paired outcomes and donor-held-out forecasts by direct independent calculation."""
import csv
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def table(filename):
    with (HERE / filename).open(encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def percentile(values):
    sorted_values = np.sort(values)
    lower = np.searchsorted(sorted_values, values, side="left")
    upper = np.searchsorted(sorted_values, values, side="right")
    return (lower + upper - 1) * 50 / (len(values)-1)


def main():
    with np.load(HERE/"paired_counts.npz", allow_pickle=False) as matrix:
        genes = matrix["genes"]
        count = matrix["counts"]
        library_depth = matrix["total"]
        selected = matrix["target"]
    assert len(genes)==10805 and selected.sum()==101 and (count.sum(axis=0) <= library_depth).all()
    assert [int(r["all_feature_total"]) for r in table("sample_map.tsv")] == list(library_depth)
    p = {r["gene"]:r for r in table("../../stage12/culture/predictions.tsv")}
    selected_genes = list(genes[selected])
    assert set(selected_genes) == set(p)
    target_cpm = count * (1.e6/library_depth)
    scales = {"rank":np.stack([percentile(col)[selected] for col in target_cpm.T]),
              "log2":np.log2(1+target_cpm[selected].T)}
    effects = {s:vals[5:]-vals[:5] for s,vals in scales.items()}
    by_gene = {(int(r["donor"]),r["gene"],r["method"]):r for r in table("hongkong_target_effects.tsv")}
    lodo = {(r["scale"],int(r["donor"]),r["method"]):r for r in table("donor_transfer_units.tsv")}
    ranges = {(int(r["donor"]),r["method"]):r for r in table("donor_transfer_bounds.tsv")}
    for i in range(5):
        others=[j for j in range(5) if j!=i]
        for s, delta in effects.items():
            transfer = sum(delta[j] for j in others) / len(others)
            for mode,value in [("zero",np.zeros(101)),("gene_mean",transfer),
                               ("gene_mean_quarter",transfer/4),("panel_mean",np.full(101,transfer.mean()))]:
                err = value - delta[i]
                row = lodo[s,i+1,mode]
                np.testing.assert_allclose([float(row["gene_mae"]),float(row["panel_abs_error"])],
                                           [np.mean(np.abs(err)),np.abs(err.mean())],atol=1e-9)
                if s=="rank":
                    baseline = scales["rank"][i]
                    clipped = np.maximum(-baseline, np.minimum(100-baseline, value))
                    error = clipped-delta[i]
                    bounded = ranges[i+1,mode]
                    np.testing.assert_allclose([float(bounded["gene_mae"]),float(bounded["panel_abs_error"])],
                                               [np.mean(np.abs(error)),np.abs(error.mean())],atol=1e-9)
            for j,gene in enumerate(selected_genes):
                for label,prev in [("fibroblast","fibroblast"),("tissue","mixture"),("zero","zero")]:
                    row=by_gene[i+1,gene,label]
                    np.testing.assert_allclose(float(row[s+"_observed"]),delta[i,j],atol=1e-10)
                    if s=="log2":
                        np.testing.assert_allclose(float(row["log2_predicted"]),float(p[gene][prev]),atol=1e-10)
    print("Verified paired measurements, historical target multipliers, all 40 leave-one-donor-out scores and bounded rank checks.")


if __name__ == "__main__":
    main()
