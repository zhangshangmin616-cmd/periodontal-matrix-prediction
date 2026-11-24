# Time-forward reconstruction in GSE327337

The archive contains one distinguishable library per condition. The unloaded library (CP) and each day-3 library supply the predictors. The two day-7 libraries supply the outcomes. A previous stage examined these outcomes, which makes this comparison a retrospective internal analysis. Aggregate counts came from 62,490 retained barcodes across five libraries. No patient identifiers link the retained barcodes to individual donors. [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327337) · [Atlas publication](https://doi.org/10.1177/00220345261491310)

## Inputs, predictions and estimands

For gene g, a group profile equals the gene's UMI count in that group divided by the sum of UMIs for **all** assayed features in the group, multiplied by one million. Group RNA shares sum to one within a library. The weighted sum of group profiles is the modeled abundance M. The counts file retains 37,433 uniquely named genes and the full library totals, including UMIs assigned to symbols with ambiguous names. Capture efficiency and per-cell RNA content can change the RNA shares alongside the number of captured cells.

| Prediction | RNA shares | Group profiles |
|---|---|---|
| Unchanged | CP | CP |
| State only | CP | day 3 in five named cell groups, CP in the residual group |
| Composition only | day 3 | CP |
| State and composition | day 3 | day 3 in five named groups, CP in the residual group |
| Measured day 3 | day 3 | day 3 in all groups, equivalent to ordinary whole-tissue persistence in the direct reconstruction |

The five named groups correspond to fibroblast-like, mural, immune, endothelial and glial cells. For untreated abundance B, the direct reconstruction uses M. The earlier transfer rule uses B × (1 + M)/(1 + B). For prediction and observation we then compute log₂(1 + abundance) minus log₂(1 + B). The direct reconstruction is the natural within-atlas comparison because it estimates tissue CPM from the six profiles. The alternative rule was originally intended to transfer a regularized fold change onto a separate tissue baseline. It gives zero abundance whenever untreated abundance is zero. Both layers use the same untreated and observed measurements. The transfer multiplier by itself equals the direct change when expressed on the log₂ scale. Its subsequent application to B changes the abundance that enters the rank calculation.

For the percentile score, each gene receives its midrank on a 0–100 scale in a common background of 10,806 genes. Integer or rational values enforce exact ties. The 101 targets form the fixed “other ECM” set from the earlier project. Their mean percentile shift is one scalar per side. Gene-wise absolute error gives equal weight to the absolute error of every gene; the error of the average percentile allows positive and negative errors to cancel. The mean percentile across the entire ranked background is constant by construction, so the background rank endpoint here measures gene-wise rank displacement. Two-side means weight pressure and tension equally. The additional pressure-minus-tension comparison removes the common CP term. Reversing the predictions' side labels tests sign sensitivity descriptively.

## Results

| Day-3 prediction | 10,806-gene mean absolute error | 101-gene mean absolute error | 101 mean-percentile score error, direct | 101 mean-percentile score error, transfer |
|---|---:|---:|---:|---:|
| Unchanged | 0.243 | 0.384 | 1.184 | 1.184 |
| State only | 0.369 | 0.554 | 1.219 | 1.130 |
| Composition only | 0.311 | 0.651 | 2.037 | 1.989 |
| State and composition | 0.422 | 0.774 | 2.631 | 2.472 |
| Measured day 3 | 0.545 | 0.923 | 2.428 | 2.236 |

The first two columns use changes in directly reconstructed log₂(CPM+1). The last two use percentage points of mean rank change. Each column averages the two side-specific absolute errors. Across 37,433 unique genes, the corresponding direct expression errors were 0.188, 0.266, 0.227, 0.280 and 0.384. This set includes undetected genes.

Within-group state updates improved the panel-mean rank only under the earlier transfer rule. Under direct reconstruction both sides had negative observed changes, −1.053 and −1.316 percentage points. The state-only model predicted +0.247 and −0.177. The illustrative 0.2-fold amplitude sensitivity uses the median of the pre-existing patient-wise calibration coefficients in Kim, separately 0.2 for the two model variants in that earlier dataset. We applied the same factor to every predicted endpoint change in the present analysis. The older fit used a panel-score endpoint and was not a calibration for these per-gene outcomes. This exploratory sensitivity was added after the undamped summary. Under direct reconstruction, state-only gene-wise errors at 0.2 were 0.390 for the 101 targets and 0.246 for the background, compared with the unchanged values of 0.384 and 0.243.

For the direct pressure-minus-tension expression contrast, the composition-only model achieved a 101-target error of 0.379 versus 0.416 under no difference and 0.503 after exchanging the day-3 side labels. Among 100 targets with a measured side difference, it placed 70 in the correct direction, 29 in the opposite direction, and predicted zero for one. Across the 10,806-gene background the corresponding errors were 0.298, 0.312 and 0.339. For the same 101-target contrast the composition-only change in *mean percentile* had the opposite sign to the observed change (+0.263 observed versus −0.086 predicted). These estimates describe this atlas and show how the endpoint and contrast affect a model comparison. They do not establish a causal composition effect or between-person accuracy.

## Reproduce

An exploratory depth check kept the day-3 composition-only prediction fixed and thinned the two observed day-7 libraries without replacement. With both sides measured at 7,816,091 UMIs, the per-gene 101-target side-contrast MAE gain over zero was positive in all 200 draws, with median 0.0371 and 2.5th–97.5th percentiles 0.0316–0.0439. The gain also remained positive in all 200 draws at half this depth. Each draw reuses the same original material. These percentiles describe variation under technical subsampling and cannot estimate variation across donors.

Keep stage5, stage10 and stage11 beside one another. Python 3.10+ and NumPy are sufficient. Run from a terminal:

```
python stage11/benchmark/benchmark.py
python stage11/benchmark/check_benchmark.py
python stage11/benchmark/depth_sensitivity.py
```

- `gene_metrics.tsv` reports the number of genes with nonzero observed change, correct direction, opposite direction, abstentions on a change, ties and MAE for each side or cross-side contrast. Direction accuracy is correct divided by all nonzero observations. Numerical differences smaller than or equal to 10⁻¹² count as zero in direction tallies.
- `two_side_mean.tsv` reports equally weighted means of the two side-specific gene-level metrics. `panel_rank.tsv` gives the 101-gene scalar endpoint by side and its side difference. `external_damping.tsv` reports the same endpoint comparisons after multiplying predictions of change by 0.2.
- `audit.json` holds the SHA-256 hashes for the integer counts, calibration, rank-check table and features file. The checks reproduce the earlier frozen state panel result within 10⁻¹⁰ points and reproduce the raw day-3 profile under direct reconstruction. A second implementation independently walks the group counts and checks 60 direct gene-wise MAE values, with a maximum observed discrepancy of 2.22 × 10⁻¹⁶.
