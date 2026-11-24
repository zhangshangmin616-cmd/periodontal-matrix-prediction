"""Apply the unchanged stage 13 atlas models to another tensile culture."""
import csv
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def write(name, rows):
    with (HERE / name).open("w", encoding="utf8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    spec = spec_from_file_location("atlas13", WORK / "stage13/external/score.py")
    frozen = module_from_spec(spec)
    spec.loader.exec_module(frozen)
    with np.load(HERE.parent / "xian/xian_counts.npz", allow_pickle=False) as data:
        genes, target, arm, runs = [data[k] for k in ("genes", "target", "arm", "run")]
        cpm = data["counts"] * 1e6 / data["total"]
    assert cpm.shape == (10805, 8) and target.sum() == 101
    assert sorted(arm) == ["0"] * 3 + ["0.1"] * 2 + ["0.5"] * 3
    models = frozen.atlas(genes)
    ranks = np.stack([frozen.ranks(column)[target] for column in cpm.T], axis=1)
    observed = {"log2": np.log2(1+cpm[target]), "rank": ranks}
    full, deletions, effects = [], [], []
    for method, shift in models.items():
        forecast_rank = np.zeros_like(ranks)
        for j in np.flatnonzero(arm == "0"):
            forecast_rank[:, j] = frozen.ranks(cpm[:, j] * np.exp2(shift))[target] - ranks[:, j]
        for freq in ["0.1", "0.5"]:
            relevant = np.flatnonzero((arm == "0") | (arm == freq))
            for deleted in [None, *relevant]:
                keep = relevant[relevant != deleted] if deleted is not None else relevant
                control = keep[arm[keep] == "0"]
                loaded = keep[arm[keep] == freq]
                forecasts = {"log2": shift[target], "rank": forecast_rank[:, control].mean(axis=1)}
                for scale, x in observed.items():
                    y = x[:, loaded].mean(axis=1) - x[:, control].mean(axis=1)
                    p = forecasts[scale]
                    delta = p-y
                    result = dict(frequency_Hz=freq, scale=scale, method=method, genes=len(y),
                                  control_runs=len(control), loaded_runs=len(loaded), gene_mae=np.mean(abs(delta)),
                                  gain_vs_zero=np.mean(abs(y))-np.mean(abs(delta)),
                                  panel_abs_error=abs(np.mean(delta)), panel_gain_vs_zero=abs(y.mean())-abs(delta.mean()))
                    if deleted is None:
                        full.append(result)
                        effects.extend(dict(frequency_Hz=freq, scale=scale, method=method, gene=g,
                                            observed=o, prediction=i) for g,o,i in zip(genes[target], y,p))
                    else:
                        deletions.append(dict(omitted_run=runs[deleted], **result))
    write("scores.tsv", full)
    write("deletions.tsv", deletions)
    write("targets.tsv", effects)


if __name__ == "__main__":
    main()
