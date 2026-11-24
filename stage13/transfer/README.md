# Can empirical clinical responses forecast the matrix genes in another study?

The 101 targets are the fixed measurable subset of the Reactome extracellular matrix organisation panel outside collagen formation and matrix degradation. All changes and errors refer to within-specimen expression rank, measured in percentile points with the same 10,806-gene reference. Seoul provides 11 paired participant changes after three weeks. Leuven provides separate batch-matched, between-group contrasts at 7 and 28 days. The genes, source, times, baselines and historical cell-atlas predictions are the same as in the earlier article. This cross-study benchmark is retrospective. [Seoul primary study](https://doi.org/10.3390/jcm9020358) · [Leuven primary study](https://doi.org/10.1111/jre.13352) · [Force atlas](https://doi.org/10.1177/00220345261491310)

### Cohort-level forecast performance

Each row compares a forecast with the gene-wise **mean observed rank change in the destination cohort**. The Seoul vector averages its 11 participants equally. Gene error takes the absolute difference for each gene after that averaging. The direct clinical predictor carries the source cohort's mean for that gene across studies. The two quarter-strength predictors multiply their own full-scale forecasts by 0.25. The uniform clinical predictor gives all targets the same value, equal to the source panel mean.

| Source → destination | Zero | Uniform clinical | Direct clinical | Quarter clinical | Atlas | Quarter atlas |
|---|---:|---:|---:|---:|---:|---:|
| Seoul at 3 weeks → Leuven at 7 days | 6.173 | 6.166 | 5.978 | 5.998 | 6.536 | 6.181 |
| Seoul at 3 weeks → Leuven at 28 days | 2.511 | 2.531 | 2.767 | 2.389 | 3.299 | 2.499 |
| Leuven at 7 days → Seoul at 3 weeks | 2.098 | 2.239 | 5.978 | 2.343 | 3.573 | 2.293 |

The same quarter-strength clinical rule improved gene-level accuracy against zero at both Leuven times. It also outperformed the quarter-strength atlas rule in both comparisons. In the reciprocal test, zero was more accurate. The unscaled absolute distance between the two cohort vectors is 5.978 **in either direction by definition**. The two zero-change errors are different because the measured magnitudes differ across the cohorts. A change in gain across directions cannot establish biological directionality. The rank correlations across the 101 observed changes were 0.272 (Seoul versus Leuven day 7) and 0.374 (Seoul versus Leuven day 28).

The original patient-first analysis gives every paired Seoul participant a separate error before averaging. When Leuven supplied the predictor, its unscaled and quarter-strength errors by this measure were 8.701 and 5.926, versus 5.554 for zero. After forcing each gene's predicted final rank to stay within that participant's 0–100 range the errors were 7.255 and 5.750. In the common group-level comparison, bounding each participant's prediction before averaging gave 4.546 and 2.168 versus 2.098 for zero. Clipping after pooling gave 4.823 and 2.192. The forward quarter-strength clinical predictions required no clipping. The gene-wise clinical average specifies a scoring comparator and does not reconstruct a complete library-wide rank order.

The direct clinical forecast and the uniform source-panel forecast return the same panel mean in each destination when unbounded. In Leuven day 7, aggregate absolute error was 0.558 for both, 0.405 for the atlas and 0.629 for zero. Thus the atlas better approximated the panel average than the unscaled historical human data in that comparison, while the empirical human gene forecast had smaller gene-level errors. When fitting a single nonnegative attenuation of atlas effects in the source cohort, the unconstrained slopes across 10,589 non-ECM genes were −0.143 (Seoul) and −0.143 (Leuven day 7). Both constrained slopes were zero. These source fits address one parameterisation of the atlas forecast.

### One-unit deletion influence

Deleting one Seoul source participant or one Leuven untreated or treated specimen at a time retained a positive forward quarter-strength gain at all 35 deletions at day 7 and all 50 at day 28. The minimum gains across deletions were 0.107 and 0.033 points. The advantage of the unscaled clinical predictor at day 7 reversed slightly after omitting one Seoul source participant (−0.003 points). In the reciprocal comparison the patient-first negative quarter-strength gain persisted at all 35 deletions. Frequencies refer to overlapping recalculations. They are not intervals across prospective patients, estimates of treatment effects or additional validation cohorts. Leuven specimen identities cannot be aggregated by donor.

### Reproduction

The module requires Python 3.10 or newer and NumPy. From the extracted project root run the following commands in order.

```text
python stage13/transfer/score.py
python stage13/transfer/aggregation.py
python stage13/transfer/deletions.py
python stage13/transfer/check.py
```

`score.py` reads `stage11/specificity/rank_changes.npz` and `stage12/review/control_ranks.tsv`. It calculates the source slope on non-matrix genes, builds all 1,313 held-study target-gene forecast records, and saves original and bounded scores. `aggregation.py` scores cohort averages, including the reciprocal group-matched sensitivity prompted by the earlier unit-first analysis. `deletions.py` reads `stage5/zhao/joined_manifest.tsv`, `stage5/zhao/feature_map.tsv`, all 63 `stage5/zhao/counts/*.tsv.gz` and the frozen rank changes. It first verifies the two Leuven target contrasts reconstructed from the count files. All score scripts run offline. `check.py` independently re-evaluates target assignments, source slopes, unit-first and group-first errors, and the deletion summaries.
