"""Build one machine-readable source file per supplementary table from the public analysis archive."""
import argparse
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
args = argparse.ArgumentParser()
args.add_argument("--archive-root", type=Path, default=HERE.parents[1], help="Directory containing stage12 through stage18")
args.add_argument("--output", type=Path, default=HERE / "data")
opt = args.parse_args()
opt.output.mkdir(parents=True, exist_ok=True)

def read(rel):
    with (opt.archive_root / rel).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def save(n, rows):
    with (opt.output / f"Supplementary_Table_{n}.csv").open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def only(rows, **match):
    found = [r for r in rows if all(r[k] == str(v) for k, v in match.items())]
    assert len(found) == 1, match
    return found[0]

side_labels = dict(original="Original assignment", side_labels_exchanged="Pressure and tension exchanged", same_on_both_sides="Mean profile assigned to both sides")
rows = []
for r in read("stage16/sides/metrics.tsv"):
    rows.append(dict(assignment=side_labels[r["assignment"]], **{k:v for k,v in r.items() if k!="assignment"}))
save(1, rows)

rows = []
def score_row(test, forecast, scale, panel0, panel, gene0, gene, units):
    return dict(test=test,forecast=forecast,scale=scale,n_targets=101,scoring_units=units,
                zero_panel_error=panel0, forecast_panel_error=panel,zero_gene_mae=gene0,forecast_gene_mae=gene)
mix = read("stage12/control/scores.tsv")
for cond, label in [("state_only","Expression changed"),("shares_only","RNA shares changed"),("both","Both changed")]:
    base=only(mix,constructed=cond,forecast="unchanged")
    pred=only(mix,constructed=cond,forecast=cond)
    rows.append(score_row("Constructed mixture",label,"rank points",base["panel_rank_abs_error"],pred["panel_rank_abs_error"],base["gene_rank_mae"],pred["gene_rank_mae"],"4 side-time conditions x 6 overlapping cell splits, 5 original libraries"))
don=read("stage13/external/donor_transfer_scores.tsv")
for scale, label in [("rank","rank points"),("log2","log2(CPM + 1)")]:
    base=only(don,scale=scale,method="zero")
    pred=only(don,scale=scale,method="gene_mean")
    rows.append(score_row("Hong Kong, held-out donor","Four-donor gene mean",label,base["panel_abs_error"],pred["panel_abs_error"],base["gene_mae"],pred["gene_mae"],"5 donor pairs, overlapping four-donor training folds"))
clinical=read("stage13/transfer/scores.tsv")
for dest,label in [("Leuven_day7","Leuven day 7"),("Leuven_day28","Leuven day 28")]:
    for method,predictor in [("clinical_full","Seoul gene mean"),("clinical_quarter","Quarter Seoul gene mean")]:
        base=only(clinical,source="Seoul_week3",destination=dest,model="zero",bounds="original")
        pred=only(clinical,source="Seoul_week3",destination=dest,model=method,bounds="original")
        rows.append(score_row(label,predictor,"rank points",base["panel_mae"],pred["panel_mae"],base["gene_mae"],pred["gene_mae"],"11 Seoul participants, 1 destination group contrast"))
save(2,rows)

rows=[]
for t in [24,48]:
    for key,label in [("rank","Count rank"),("fpkm","FPKM rank"),("log2","log2(CPM + 1)")]:
        if key=="fpkm":
            r=only(read("stage17/culture/scores.tsv"),method="tissue",hours=t)
            row=score_row(f"Wuhan {t} h","Tension atlas",label,r["zero_panel_error"],r["panel_error"],r["zero_gene_mae"],r["gene_mae"],"3 loaded and 3 baseline libraries")
            background=r["gene_background"]
        else:
            r=only(read("stage17/wuhan_counts/scores.tsv"),model="tissue",scale=key,hours=t)
            row=score_row(f"Wuhan {t} h","Tension atlas",label,r["zero_panel_abs_error"],r["panel_abs_error"],r["zero_gene_mae"],r["gene_mae"],"3 loaded and 3 baseline libraries")
            background=r["background"]
        rows.append(dict(**row,background=background))
save(3,rows)

destnames={"Wuhan_24h":"Wuhan 24 h","Wuhan_48h":"Wuhan 48 h","GSE106887_12h":"GSE106887 12 h"}
ends=["panel_gain","gene_gain"]
deletions=read("stage18/transfer/deletions.tsv")
bounds=read("stage18/bounds/rank_limits.tsv")
rows=[]
for d,dest in destnames.items():
    b=only(bounds,destination=d,method="hong_kong_gene",clipping="each_control_library")
    for kind,label in [("source_donor","Source donor"),("destination_library","Destination library")]:
        chosen=[r for r in deletions if r["destination"]==d and r["method"]=="hong_kong_gene" and r["deletion"]==kind]
        if not chosen: continue
        interval={key:[float(r[key]) for r in chosen] for key in ends}
        row=dict(destination=dest,unit_omitted=label,n_deletions=len(chosen),targets=101,background=10805)
        for k,a in interval.items():
            row.update({k+"_min":min(a),k+"_max":max(a),k+"_positive":sum(v>0 for v in a)})
        rows.append(dict(**row,all_units_raw_panel_gain=b["raw_panel_gain"],all_units_clipped_panel_gain=b["clipped_panel_gain"],all_units_raw_gene_gain=b["raw_gene_gain"],all_units_clipped_gene_gain=b["clipped_gene_gain"],out_of_range_genes=b["markers_outside"],n_library_gene_values_outside=b["library_gene_values_outside"]))
save(4,rows)

scores=read("stage18/deviations/scores.tsv")
cal=read("stage18/deviations/training.tsv")
rows=[]
for d,dest in destnames.items():
    zero=only(read("stage18/transfer/scores.tsv"),method="zero",destination=d)
    rows.append(dict(destination=dest,targets=101,background=10805,zero_gene_mae=zero["gene_mae"],atlas_gene_mae=only(scores,destination=d,model="original_atlas")["gene_mae"],source_calibrated_gene_mae=only(scores,destination=d,model="fit_in_Hong_Kong")["gene_mae"],constant_atlas_gene_mae=only(scores,destination=d,model="common_shift")["gene_mae"],fixed_atlas_mean=only(scores,destination=d,model="common_shift")["panel_mean"],source_amplitude=only(cal,source_donor_omitted=0)["amplitude"],source_leave_one_out_min=min(float(r["amplitude"]) for r in cal if r["source_donor_omitted"]!="0"),source_leave_one_out_max=max(float(r["amplitude"]) for r in cal if r["source_donor_omitted"]!="0"),destination_best_amplitude=only(scores,destination=d,model="destination_best_possible")["amplitude"]))
save(5,rows)

paired=read("stage17/clinical_search/scores.tsv")
uncertainty=read("stage17/clinical_search/uncertainty.tsv")
rows=[]
for r in paired:
    row={**r}
    for ep,prefix in [("panel_mean","panel"),("mean_individual_gene","gene")]:
        limits=only(uncertainty,endpoint=ep)
        row[prefix+"_gain_conditional_95pct_lo"]=limits["lower_pct_2.5"] if r["unit"]=="mean_of_paired_scores" else ""
        row[prefix+"_gain_conditional_95pct_hi"]=limits["upper_pct_97.5"] if r["unit"]=="mean_of_paired_scores" else ""
    rows.append(row)
save(6,rows)

rows=[]
xian=read("stage14/xian/scores.tsv")
for freq in ["0.1","0.5"]:
    for s,label in [("rank","rank points"),("log2","log2(CPM + 1)")]:
        pred=only(xian,source="HongKong",destination=f"Xian_{freq}Hz",scope="target",method="source_gene",scale=s)
        zero=only(xian,source="HongKong",destination=f"Xian_{freq}Hz",scope="target",method="zero",scale=s)
        row=score_row(f"Xi'an {freq} Hz","Hong Kong gene mean",label,zero["panel_error"],pred["panel_error"],zero["gene_mae"],pred["gene_mae"],"3 shared controls, "+("2" if freq=="0.1" else "3")+" loaded libraries")
        rows.append(dict(**row,full_panel_gain_lo=pred["panel_gain"],full_panel_gain_hi=pred["panel_gain"],full_gene_gain_lo=pred["gene_gain"],full_gene_gain_hi=pred["gene_gain"]))
for r in read("stage17/axial/scores.tsv"):
    cond=r["condition"]
    row=score_row("Chengdu "+("extrusion " if cond[0]=="E" else "intrusion ")+cond[1:]+" d","Atlas side proxy","log2 fold change",r["measured_zero_panel_error"],r["measured_model_panel_error"],r["measured_zero_gene_mae"],r["measured_model_gene_mae"],"1 published cohort contrast")
    row["n_targets"]=r["n_measured"]
    rows.append(dict(**row,**{k:r[k] for k in ("full_panel_gain_lo","full_panel_gain_hi","full_gene_gain_lo","full_gene_gain_hi")}))
save(7,rows)
