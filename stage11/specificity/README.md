# External expression-matched gene-set analysis

The aggregate score and gene-level scores use percentile points among the same 10,806 genes within each sample. The aggregate loss takes an absolute value after averaging 101 predicted changes and 101 observed changes separately. The gene-level loss first computes 101 absolute prediction errors and then averages them. The paired Kim analysis computes each loss within a participant before weighting all 11 participants equally. The Zhao data supply one group contrast per time and archive batch. Patient pairing is unavailable there. The observed markers have no day-7 pressure/tension labels in these cohorts. Forecasts use equal pressure and tension contributions and the atlas day-7 cell-group expression changes at every target time.

## Results

| Contrast | Aggregate change MAE, model / zero | Individual-gene MAE, model / zero | Matched sets with equal or greater aggregate gain |
|---|---:|---:|---:|
| Kim, 3 weeks, 11 paired participants | 1.812 / 1.884 | 6.249 / 5.554 | 271 / 10,000 |
| Zhao, day 7, 19 exposed and 5 control specimens | 0.405 / 0.629 | 6.536 / 6.173 | 434 / 10,000 |
| Zhao, day 28, 17 exposed and 22 control specimens | 0.504 / 0.141 | 3.299 / 2.511 | 7,720 / 10,000 |

The untreated atlas contains 10,589 eligible non-ECM genes in the fixed background. The 101 targets, 32 formation-only genes, 51 degradation-only genes and 33 genes belonging to both groups are excluded from that comparison pool. Sets preserve the 101 targets' counts per expression decile in the untreated atlas. This comparison does not match genes for co-expression, variance or function. Its proportions describe relative performance among these sets and cannot measure patient-level significance or establish that the result is specific to ECM biology.

In Kim, 5 participants improved and 6 worsened on aggregate error, with mean gain 0.072 points. A descriptive 95% bootstrap interval was −0.011 to +0.163. Individual-gene error increased in all 11 participants. Mean gain was −0.696 (interval −0.893 to −0.496). Each of 50,000 draws resampled 11 people with replacement and kept both metrics paired by person. Predictions and source weights remained fixed. The interval describes the variation across those completed predictions. It does not propagate uncertainty from estimating the source atlas, estimating each participant's weights using other controls, estimating individual treatment effects or choosing the scoring rule.

## Offline reproduction

From this directory run `python score.py` and `python check.py` using Python 3.10 or above with NumPy. rank_changes.npz contains all 13 × 10,806 predicted and observed rank-change pairs plus gene labels, pathway categories and untreated atlas UMI totals. score.py regenerates comparison.tsv, unit_errors.tsv, target_gene_errors.tsv, kim_patient_sensitivity.tsv, matching_deciles.tsv and matched_gains.npz. check.py independently reproduces six observed score comparisons and 24 randomly matched score comparisons by scalar arithmetic.

For source-count reconstruction, put the stage 5 research package contents in work/stage5, the stage 10 archive's stage10/reference_counts.npz in work/stage10, and this directory in work/stage11/specificity. Run `python rebuild_inputs.py` with NumPy, pandas and SciPy. It reads the frozen stage 5 per-control RNA weights, the published count files, metadata and profiles, then verifies the aggregate prediction and observed change of each of the 13 units against the earlier stage 5 outputs (tolerance 10⁻⁹ percentile points). No patient pairing is inferred in Zhao. The 63 samples there contribute only to their own time and batch group contrasts.

The source atlas is GSE327337. The three-week paired cohort is GSE165988. The Zhao count files derive from PRJNA1136245. [Atlas](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327337) · [Kim](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165988) · [Zhao](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1136245)
