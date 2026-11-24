"""Leave one donor out in the Hong Kong tensile culture experiment."""
import numpy as np

from score import HERE, ranks, save


def main():
    with np.load(HERE / "paired_counts.npz", allow_pickle=False) as data:
        target = data["target"]
        cpm = data["counts"] * 1.e6 / data["total"]
    rows, bounds = [], []
    for scale, value in [("log2", np.log2(1. + cpm)),
                         ("rank", np.column_stack([ranks(col) for col in cpm.T]))]:
        delta = (value[:, 5:] - value[:, :5])[target]
        for i in range(5):
            others = np.delete(delta, i, axis=1)
            assert others.shape == (101, 4)
            learned = others.mean(axis=1)
            for label, pred in [("zero", np.zeros_like(learned)), ("panel_mean", np.full_like(learned,learned.mean())),
                                ("gene_mean_quarter", learned*0.25), ("gene_mean", learned)]:
                y = delta[:, i]
                err = pred-y
                rows.append(dict(scale=scale, donor=i+1, method=label, genes=len(y),
                                 gene_mae=np.mean(abs(err)), panel_abs_error=abs(err.mean()),
                                 observed_mean=y.mean(), predicted_mean=pred.mean()))
                if scale == "rank":
                    untreated = value[target, i]
                    before = untreated + pred
                    clipped = np.clip(before, 0., 100.) - untreated
                    clipped_error = clipped-y
                    bounds.append(dict(donor=i+1,method=label,targets=len(y),
                                       below_zero=int((before<0).sum()),above_100=int((before>100).sum()),
                                       gene_mae=np.mean(abs(clipped_error)),panel_abs_error=abs(clipped_error.mean())))
    summary = []
    for (scale, method) in dict.fromkeys((r["scale"],r["method"]) for r in rows):
        subset = [r for r in rows if r["scale"]==scale and r["method"]==method]
        base = [r for r in rows if r["scale"]==scale and r["method"]=="zero"]
        summary.append(dict(scale=scale, method=method, donors=len(subset), genes=101,
                            gene_mae=np.mean([r["gene_mae"] for r in subset]),
                            panel_abs_error=np.mean([r["panel_abs_error"] for r in subset]),
                            better_than_zero=sum(r["gene_mae"] < z["gene_mae"] for r,z in zip(subset,base))))
    save("donor_transfer_scores.tsv", summary)
    save("donor_transfer_units.tsv", rows)
    save("donor_transfer_bounds.tsv", bounds)


if __name__ == "__main__":
    main()
