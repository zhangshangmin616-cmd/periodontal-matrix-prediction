"""Check physical percentile bounds for the additive cross-laboratory rank forecasts."""
import csv
import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def load_module(name, file):
    spec = importlib.util.spec_from_file_location(name, file)
    code = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(code)
    return code


def table(path):
    with path.open(encoding="utf8", newline="") as source:
        return list(csv.DictReader(source, delimiter="\t"))


def main():
    transfer = load_module("transfer", HERE.parent / "transfer/score.py")
    original = load_module("frozen", WORK / "stage13/external/score.py")
    with np.load(WORK / "stage13/external/paired_counts.npz", allow_pickle=False) as inputs:
        genes, mask = inputs["genes"], inputs["target"]
    target = genes[mask]
    multipliers = np.exp2(original.atlas(genes)["tissue"])
    # Independent midpoint calculation. Zero-based midranks use denominator N-1.
    def target_ranks(values):
        ordered = np.sort(values)
        q = values[mask]
        return (ordered.searchsorted(q) + ordered.searchsorted(q, side="right")-1)*50/(len(values)-1)

    effects = table(HERE.parent / "transfer/genes.tsv")
    rows = []
    for cohort, ids, loaded, cpm in transfer.destination_inputs(genes, mask):
        x = np.stack([target_ranks(v) for v in cpm])
        control = x[~loaded]
        baseline = x[~loaded].mean(axis=0)
        y = x[loaded].mean(axis=0) - baseline
        atlas = np.stack([target_ranks(v*multipliers) for v in cpm[~loaded]]) - control
        for method, source in [("atlas_common", "atlas_tissue"), ("atlas_tissue", "atlas_tissue"),
                               ("hong_kong_common", "hong_kong_common"), ("hong_kong_gene", "hong_kong_gene")]:
            selected = {r["gene"]:r for r in effects if r["destination"] == cohort and r["method"] == source}
            p = np.array([float(selected[g]["forecast"]) for g in target])
            np.testing.assert_allclose(y, [float(selected[g]["observed"]) for g in target], atol=2e-12, rtol=0)
            if source == "atlas_tissue":
                np.testing.assert_allclose(p, atlas.mean(axis=0), atol=2e-12, rtol=0)
                by_library = atlas.copy()
                if method == "atlas_common":
                    # Take the atlas panel mean within each control library.
                    by_library = np.broadcast_to(by_library.mean(axis=1)[:,None], control.shape)
                    p[:] = p.mean()
            else:
                by_library = np.broadcast_to(p[None,:], control.shape)
            np.testing.assert_allclose(p, by_library.mean(axis=0), atol=2e-12, rtol=0)
            panel_gain = lambda z: abs(y.mean())-abs(np.mean(z-y))
            gene_gain = lambda z: np.mean(abs(y))-np.mean(abs(z-y))
            for level, b, preds in [("each_control_library",control,by_library),
                                     ("mean_control_profile", baseline[None,:],p[None,:])]:
                end = b + preds
                clipped = (np.clip(end,0,100)-b).mean(axis=0)
                violating = (end < -1e-10) | (end > 100+1e-10)
                outside = violating.any(axis=0)
                rows.append(dict(destination=cohort, method=method, clipping=level, measured_targets=len(y),
                                 background=len(genes), genes_outside=int(outside.sum()),
                                 library_gene_values_outside=int(violating.sum()),
                                 markers_outside="|".join(target[outside]), min_final_rank=end.min(),
                                 max_final_rank=end.max(), raw_panel_gain=panel_gain(p),
                                 clipped_panel_gain=panel_gain(clipped), raw_gene_gain=gene_gain(p),
                                 clipped_gene_gain=gene_gain(clipped)))
    assert len(rows) == 24
    with (HERE/"rank_limits.tsv").open("w", encoding="utf8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
