"""Extract all gene-level observations and predictions behind the three figures."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def read(path):
    with (WORK / path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def save(number, rows):
    with (HERE / f"Figure_{number}_source_data.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def unit(accession, condition, scored_unit, baseline_libraries, loaded_libraries, background):
    return dict(accession=accession, condition=condition, scored_unit=scored_unit,
                baseline_libraries=baseline_libraries, loaded_libraries=loaded_libraries,
                rank_background=background, panel_genes=101, unit="expression percentile points")


def fig1():
    rows = [r for r in read("stage13/transfer/forecasts.tsv") if r["destination"] == "Leuven_day7"]
    shared = mean(float(r["atlas_full"]) for r in rows)
    assert len(rows) == len({r["gene"] for r in rows}) == 101
    meta = unit("PRJNA1136245", "Leuven day 7", "within-batch specimen-group contrast", 5, 19, 10806)
    save(1, [dict(**meta, gene=r["gene"], observed=float(r["observed_change"]), zero=0.0,
                  atlas=float(r["atlas_full"]), atlas_shared_mean=shared) for r in rows])


def fig2():
    rows = []
    for r in read("stage17/wuhan_counts/targets.tsv"):
        if r["model"] == "tissue" and r["scale"] == "rank":
            meta = unit("GSE312900", f"Wuhan {r['hours']} h", "group contrast", 3, 3, 10806)
            rows.append(dict(**meta, prediction_source="day-3 tension atlas", donor="", gene=r["gene"],
                             observed=float(r["observed"]), zero=0.0, predicted=float(r["predicted"])))
    for r in read("stage17/zhejiang/targets.tsv"):
        if r["method"] == "tissue" and r["input"] == "NCBI_counts":
            meta = unit("GSE106887", "GSE106887 12 h", "single library per arm", 1, 1, 10806)
            rows.append(dict(**meta, prediction_source="day-3 tension atlas", donor="", gene=r["gene"],
                             observed=float(r["observed"]), zero=0.0, predicted=float(r["predicted"])))
    donor_changes = defaultdict(dict)
    for r in read("stage18/transfer/hong_kong_donor_changes.tsv"):
        donor_changes[r["gene"]][int(r["donor"])] = float(r["observed_rank_change"])
    assert len(donor_changes) == 101 and all(len(values) == 5 for values in donor_changes.values())
    for donor in range(1, 6):
        meta = unit("GSE173891", "Hong Kong 12 h", "one held-out donor pair", 1, 1, 10805)
        for g, vals in sorted(donor_changes.items()):
            rows.append(dict(**meta, prediction_source="mean of the other four donor changes", donor=donor,
                             gene=g, observed=vals[donor], zero=0.0,
                             predicted=mean(v for d, v in vals.items() if donor != d)))
    assert len(rows) == 808
    save(2, rows)


def fig3():
    lookup = defaultdict(dict)
    for r in read("stage18/transfer/genes.tsv"):
        key = r["destination"], r["gene"]
        lookup[key][r["method"]] = float(r["forecast"])
        lookup[key]["observed"] = float(r["observed"])
    amplitude = float(read("stage18/deviations/training.tsv")[0]["amplitude"])
    rows = []
    for condition in ("Wuhan_24h", "Wuhan_48h", "GSE106887_12h"):
        wuhan = condition.startswith("Wuhan")
        accession, nlib = ("GSE312900", 3) if wuhan else ("GSE106887", 1)
        meta = unit(accession, condition, "group contrast" if wuhan else "single library per arm", nlib, nlib, 10805)
        values = sorted((g, v) for (c, g), v in lookup.items() if c == condition)
        assert len(values) == 101
        m = mean(v["atlas_tissue"] for _, v in values)
        for gene, data in values:
            rows.append(dict(**meta, paired_Hong_Kong_source_donors=5, gene=gene, observed=data["observed"],
                             zero=data["zero"], atlas=data["atlas_tissue"], atlas_shared_mean=m,
                             centred_amplitude=amplitude,
                             atlas_calibrated=m + amplitude * (data["atlas_tissue"] - m),
                             measured=data["hong_kong_gene"], measured_shared_mean=data["hong_kong_common"]))
    save(3, rows)


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
