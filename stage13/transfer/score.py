"""Transfer clinical gene responses or an atlas scaling factor between studies."""
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def write(name, rows):
    with (HERE / name).open("w", newline="", encoding="utf8") as stream:
        table = csv.DictWriter(stream, rows[0], delimiter="\t", lineterminator="\n")
        table.writeheader()
        table.writerows(rows)


def average_ranks(values):
    _, indices, counts = np.unique(values, return_inverse=True, return_counts=True)
    return (np.cumsum(counts) - (counts - 1) / 2)[indices]


def main():
    with np.load(WORK / "stage11/specificity/rank_changes.npz", allow_pickle=False) as data:
        units = data["units"]
        component = data["component"]
        targets = component == "other"
        calibration = component == "outside_ecm"
        genes = data["gene"][targets]
        observed = data["observed"]
        atlas = data["fitted"]
    assert len(units) == 13 and len(genes) == 101 and calibration.sum() == 10589
    groups = {"Seoul_week3": np.flatnonzero(np.char.startswith(units, "Kim_")),
              "Leuven_day7": np.flatnonzero(units == "Zhao_day7"),
              "Leuven_day28": np.flatnonzero(units == "Zhao_day28")}
    with (WORK / "stage12/review/control_ranks.tsv").open(encoding="utf8") as stream:
        baselines = {(r["unit"], r["gene"]):float(r["untreated_rank"])
                     for r in csv.DictReader(stream, delimiter="\t")}
    untreated = np.array([[baselines[unit, gene] for gene in genes] for unit in units])
    analysis = [("Seoul_week3", "Leuven_day7"), ("Seoul_week3", "Leuven_day28"),
                ("Leuven_day7", "Seoul_week3")]
    summary, unit_values, coefficients, forecasts, associations = [], [], [], [], []
    for source, destination in analysis:
        train = groups[source]
        target = groups[destination]
        assert not set(train).intersection(target)
        x = atlas[train][:, calibration]
        y = observed[train][:, calibration]
        slope = float(np.sum(x * y) / np.sum(x * x))
        calibrated = np.clip(slope, 0., 1.)
        coefficients.append(dict(source=source, destination=destination, source_units=len(train),
                                 calibration_genes=x.shape[1], unconstrained_slope=slope,
                                 bounded_slope=calibrated))
        source_means = observed[train][:, targets].mean(axis=0)
        endpoint = observed[target][:, targets]
        destination_atlas = atlas[target][:, targets]
        predictions = {
            "zero": np.zeros_like(endpoint),
            "clinical_uniform": np.full_like(endpoint, source_means.mean()),
            "clinical_full": np.broadcast_to(source_means, endpoint.shape),
            "clinical_quarter": np.broadcast_to(.25 * source_means, endpoint.shape),
            "atlas_full": destination_atlas,
            "atlas_quarter": .25 * destination_atlas,
            "atlas_source_calibrated": calibrated * destination_atlas,
        }
        target_average = endpoint.mean(axis=0)
        rho = np.corrcoef(average_ranks(source_means), average_ranks(target_average))[0, 1]
        associations.append(dict(source=source, destination=destination, genes=len(genes),
                                 source_panel_change=source_means.mean(),
                                 destination_panel_change=target_average.mean(),
                                 spearman_rho=rho,
                                 gene_direction_agreements=int((source_means * target_average > 0).sum()),
                                 gene_nonzero_in_both=int(((source_means!=0)&(target_average!=0)).sum())))
        for i, idx in enumerate(target):
            for k, gene in enumerate(genes):
                forecasts.append(dict(source=source, destination=destination, unit=units[idx],
                    gene=gene, untreated_rank=untreated[idx,k], observed_change=endpoint[i,k],
                    **{name:p[i,k] for name,p in predictions.items()}))
        for name, p in predictions.items():
            feasible = np.clip(p, -untreated[target], 100 - untreated[target])
            changed = np.abs(feasible - p) > 1e-10
            for bound, values in [("original", p), ("bounded", feasible)]:
                residual = values - endpoint
                gene_mae = np.abs(residual).mean(axis=1)
                panel = np.abs(residual.mean(axis=1))
                zero_gene = np.abs(endpoint).mean(axis=1)
                zero_panel = np.abs(endpoint.mean(axis=1))
                for i, idx in enumerate(target):
                    unit_values.append(dict(source=source, destination=destination, unit=units[idx],
                            model=name, bounds=bound, gene_mae=gene_mae[i], panel_abs_error=panel[i]))
                summary.append(dict(source=source, destination=destination, model=name, bounds=bound,
                        n_destination_units=len(target), genes=len(genes),
                        gene_mae=gene_mae.mean(), panel_mae=panel.mean(),
                        gene_gain_vs_zero=(zero_gene-gene_mae).mean(),
                        panel_gain_vs_zero=(zero_panel-panel).mean(),
                        units_gene_better_than_zero=int((gene_mae<zero_gene-1e-12).sum()),
                        entries_requiring_clipping=int(changed.sum())))
    write("coefficients.tsv", coefficients)
    write("associations.tsv", associations)
    write("forecasts.tsv", forecasts)
    write("unit_metrics.tsv", unit_values)
    write("scores.tsv", summary)


if __name__ == "__main__":
    main()
