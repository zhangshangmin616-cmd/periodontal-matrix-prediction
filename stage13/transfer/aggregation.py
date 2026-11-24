"""Compare three transports on the same group-average clinical rank contrast."""
from collections import defaultdict
from statistics import mean

from check import rows
from score import write


def main():
    entries = defaultdict(lambda:defaultdict(list))
    for row in rows("forecasts.tsv"):
        entries[row["source"],row["destination"]][row["gene"]].append(row)
    output = []
    for (source,destination), groups in entries.items():
        for model in ["zero","clinical_uniform","clinical_full","clinical_quarter",
                      "atlas_full","atlas_quarter","atlas_source_calibrated"]:
            errors = defaultdict(list)
            for gene, group in groups.items():
                obs = mean(float(r["observed_change"]) for r in group)
                baseline = mean(float(r["untreated_rank"]) for r in group)
                unbounded = mean(float(r[model]) for r in group)
                clipped_group = min(100.,max(0.,baseline+unbounded))-baseline
                clipped_units = mean(min(100.,max(0.,float(r["untreated_rank"])+float(r[model])))
                                     -float(r["untreated_rank"]) for r in group)
                for label,p in [("original",unbounded),("bound_unit_then_mean",clipped_units),
                                ("mean_then_bound_group",clipped_group)]:
                    errors[label].append(p-obs)
            for label,values in errors.items():
                output.append(dict(source=source,destination=destination,genes=len(values),
                                   model=model,bounds=label,group_gene_mae=mean(map(abs,values)),
                                   group_panel_abs_error=abs(mean(values))))
    write("group_metrics.tsv",output)


if __name__ == "__main__":
    main()
