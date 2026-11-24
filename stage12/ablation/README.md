# Which cell-group expression changes affect the clinical forecast?

This retrospective analysis keeps all cells and RNA shares in the model. It lets the day-7 profile replace the unloaded profile in any subset of the five named cell groups, or restores individual groups to their unloaded profiles. The residual group always retains its unloaded profile. The independent clinical readouts come from the paired Seoul series and the Leuven 7- and 28-day comparisons. Each clinical baseline and every predicted expression vector retain the original 10,806-gene background and 101-gene endpoint. [Force atlas](https://doi.org/10.1177/00220345261491310) · [Seoul clinical study](https://doi.org/10.3390/jcm9020358) · [Leuven clinical study](https://doi.org/10.1111/jre.13352)

## Clinical results

Gain is the absolute error for zero change minus that for the indicated model, in within-sample percentile points. In Seoul, errors were calculated within each participant and averaged equally over 11 people. Each Leuven analysis has one tissue-group contrast. Leuven at day 7 uses five untreated specimens and 19 loaded specimens. Day 28 uses 22 untreated specimens and 17 loaded specimens.

| State changes allowed | Seoul 3 weeks, panel gain | Leuven 7 days, panel gain | Leuven 28 days, panel gain | Leuven 7 days, gene-level gain |
|---|---:|---:|---:|---:|
| All five groups | +0.072 | +0.224 | −0.363 | −0.363 |
| All except fibroblast-like | +0.063 | +0.346 | −0.437 | −0.223 |
| All except mural | +0.038 | −0.066 | +0.014 | −0.401 |
| All except immune | +0.079 | +0.252 | −0.381 | −0.137 |
| All except endothelial | +0.054 | +0.235 | −0.371 | −0.308 |
| All except glial | +0.060 | +0.148 | −0.296 | −0.243 |
| Mural alone | +0.045 | +0.310 | −0.397 | +0.099 |
| Change RNA shares alone | −0.053 | −0.399 | −0.072 | −0.403 |

The Leuven day-7 observed shift was +0.629. Restoring only the mural profile changed the prediction from +0.224 to −0.066, taking absolute error from 0.405 to 0.695. At day 28, restoring the same group changed the forecast from +0.363 to −0.014, against an observed −0.141. Error fell from 0.504 to 0.127. All five single omissions still left gene-level error above the zero baseline in both Leuven comparisons and on average in Seoul. Mural alone improved day-7 gene-level accuracy by 0.099 points and worsened it at day 28 by 0.331 and in Seoul by 0.342. All configurations were reported, and no forecast was refitted to prefer a particular subset.

Because the expression mixture is ranked after transfer and scored by absolute error, removing a group from the full model and adding it to the all-off model can have different effects. We scored all 32 subsets and averaged each group's marginal gain across all 120 possible activation orders. Across the 16 possible backgrounds of the other groups, activating the mural profile improved day-7 panel accuracy by 0.276 to 0.319 points and worsened day-28 panel accuracy by 0.368 to 0.407 points. The averaged contributions from the five groups sum to the gain of the full model over zero. These analyses establish dependence on a component of this prediction rule. They do not resolve cell-intrinsic regulation or ascertain which populations drove the actual tissue change.

## What the mural group represents

The mural gate required detection of at least two of RGS5, NOTCH3, MYH11 and KCNJ8, with no second qualifying marker panel. It accounted for 46.9% of captured RNA in the source untreated library. The five historical controls in the Leuven day-7 batch had fitted mural RNA shares from 40.5% to 53.0%. Both RNA share and the inferred expression change influence a group's leverage. The direct sensitivity analysis does not distinguish their relative biological contributions.

For an annotation check prompted by this ablation result, we used the previously tabulated, donor-resolved cross-classification of 1,590 barcodes in an independent healthy periodontal dataset. All 91 cells called mural by the strict gate had Perivascular or SMC labels in the published integrated atlas. Among 385 cells assigned those labels, the gate recovered 91 (23.6%). Eighty-six of the 91 recovered cells came from one donor and five from a second donor. No strict mural calls were recovered in the remaining three donors, who collectively had 188 reference-labelled mural cells. This measures agreement with author labels in this dataset and gives no estimate of accuracy in the loaded source atlas. [Original tissue study](https://doi.org/10.1016/j.isci.2021.102405) · [Integrated cell atlas](https://doi.org/10.3389/fcell.2022.1021459) · [Label release, Mendeley Data v1, CC BY 4.0](https://doi.org/10.17632/2kskdknngb.1)

## Reproduction

Run `python score.py`, `python check.py` and `python check_markers.py` in this directory. The first two require NumPy. The marker check uses the Python standard library. The three commands work offline using the delivered `inputs.npz` and `reference_label_cross_table.tsv`.

- `inputs.npz` contains 41 control-specimen CPM vectors and their frozen RNA weights, source profiles and RNA shares, the 101 target indices, gene names, unit and control identifiers, observed target shifts, and the historical full-model and composition predictions. Active states alter only within-group expression. The original composition model is reproduced separately using source changes in each side's RNA share.
- `all_subsets.tsv` reports the 32 subsets and the separate composition benchmark for each clinical unit of analysis and both metrics. `ablation_summary.tsv` retains the all-on, all-off, individual omissions, single-group and composition cases. `unit_sensitivity.tsv` gives the patient-level and batch-level errors underlying these summaries. `order_averaged_contributions.tsv` reports the exact averaged allocations and minimum and maximum gains across backgrounds. `audit.json` records numerical cross-checks.
- `check.py` independently recomputes the all-on, all-off, single-group and individual-omission models with explicit class summation and direct tied-rank calculation. All 312 forecast/error rows agree within 10⁻⁸ percentile points. `check_markers.py` independently sums the healthy-data donor cross-table, copied from stage4/validation/results/author_gate_crosstab.tsv. The reference project matched source labels to healthy original counts by donor ID and barcode before applying the identical marker rule.
- Optional rebuilding from source counts uses `rebuild_inputs.py`, with NumPy, pandas and SciPy. It expects the accompanying stage5 source files and results, including stage5/results/gene_deltas_exploratory.tsv, and stage11/specificity/rank_changes.npz at their original relative paths. The previously published stage5 and stage11 local research bundles also contain those files. Rebuilding does not refit clinical weights or use loaded expression for group selection. All rank shifts for the full model exactly reproduce the preceding stage11 calculation, and the composition comparator differs from the stage5 gene effects by less than 5 × 10⁻¹⁴ percentile points.
