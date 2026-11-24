# Human periodontal ligament cells under tensile loading

## Hong Kong paired RNA sequencing

Lin et al. exposed primary human periodontal ligament cells from each of five healthy donors to 12% tensile strain for 12 hours. Cells were obtained from the middle third of the root surface. GEO supplies the authors' unfiltered matrix with 60,661 rows, five untreated columns (A1–A5) and five matched loaded columns (B1–B5). Each sample's GEO title and description were checked against the count column. Before restricting features, we used all rows to calculate library depth. Version-free Ensembl identifiers with a unique HGNC mapping linked 10,805 of the 10,806 atlas background genes and all 101 panel genes. The missing background identifier in this join was RUNX2 (ENSG00000124813). Some present target genes have zero counts in some libraries. They were retained. [Original study](https://doi.org/10.1093/ejo/cjab043) · [Series and matrix](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173891)

The frozen day-3 tension-side atlas supplied the fixed-share tissue and fibroblast-like log₂ multipliers. The log-scale analysis compared these multipliers directly with paired differences of log₂(CPM + 1). For the rank-scale analysis we multiplied each donor's untreated culture CPM by the corresponding linear multiplier and re-ranked it across the same mapped genes. We fitted no coefficient to the donor outcomes. All comparisons are gene-by-gene within donors and equal-weight averages across donors. Twelve hours and three days represent different windows of response, and the in vitro cells were not selected with the atlas marker panel.

| Forecast for the 101 targets | Gene MAE, log₂ | Gene MAE, rank points | Mean-rank panel error, points |
|---|---:|---:|---:|
| Zero change | 0.169 | 1.928 | 0.288 |
| Day-3 tissue multiplier | 0.358 | 3.780 | 0.476 |
| Day-3 fibroblast-like multiplier | 0.627 | 6.573 | 0.799 |
| Mean of that gene's changes in the four other donors | 0.138 | 1.626 | 0.341 |
| One quarter of that four-donor mean | 0.150 | 1.740 | 0.285 |
| Mean across the panel and four other donors, applied to all genes | 0.169 | 1.956 | 0.341 |

The last three rules are retrospective, within-study leave-one-donor-out comparisons. The first three forecasts were specified before assessing the Hong Kong expression measurements. The four-donor gene mean improved individual-gene accuracy in all five held-out donors, whereas its panel-mean rank error was higher on average than that of zero. Under the four-donor mean, four of 505 predicted final rank percentiles fell below zero. Clipping the final prediction of each gene to its marginal 0–100 range yielded 1.618 individual-gene MAE and 0.333 panel MAE. Bounds on individual genes do not impose a joint transcriptome rank distribution. This follow-up establishes within-study predictability of an early culture response. It is not validation at an anatomically specified day-7 clinical root surface.

## Osaka 48-hour two-channel array

The GSE19118 culture series used cyclic stretching at 0.5 Hz for 48 hours. Its two two-channel arrays had opposite sample dye labels. In GSM474198, all 30,336 probe-keyed log₂ measurements matched the normalized log-ratio column from the associated scanner file. In GSM474118, only 963 matched when joined on probe ID. A few individual normalized values match in a different row order, but we could not establish a reliable measurement map for this array. Consequently the calculation is descriptive for one GEO array. It does not provide donor-resolved replication. [Original study](https://doi.org/10.1074/jbc.M109.097303) · [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19118)

For GSM474198 the submitted numerator was the loaded sample. We kept probes with a unique, current Entrez-to-HGNC match and positive normalized background-subtracted median intensity in both channels. Averaging eligible probes per gene retained 84/101 targets and 7,929 background genes. The target gene MAEs were 0.884 (zero), 0.938 (tissue), and 1.096 (fibroblast-like) log₂ units. Including finite deposited values from weak/missing channel signal expanded coverage to 95 targets and preserved this ordering. The remaining genes were not imputed. This microarray measures fluorescence ratios. The scanner's normalized ratios and the RNA-based smoothed fold estimates have different measurement properties.

## Reproduction

From this directory with NumPy installed, run `python score.py`, `python donor_transfer.py`, `python score_chip.py`, and `python check.py`. The scripts use the two compact local inputs `paired_counts.npz` and `chip_effects.tsv` plus the archive's `stage10/reference_counts.npz` and `stage8/axial_log2_multipliers.tsv.gz`. The independent numerical check uses `stage12/culture/predictions.tsv` to verify the frozen atlas multipliers for each target. The Hong Kong tables report results by donor and scale, including full target-level forecasts. The complete-background mean rank change is zero by the rank definition, so its signed mean and panel error carry no accuracy information. All-background individual-gene errors remain informative.

For optional reconstruction from the public records, run `python fetch.py`, then `python extract.py` and `python extract_chip.py`. The fetcher uses NCBI's GEO text service and FTP HTTPS host and validates SHA-256 for every input. Source mapping uses `stage5/paired/hgnc_ids.tsv.gz` and `stage5/reference/features.tsv`. Do not interpret barcode, probe, array, gene, or repeat cross-validation predictions as additional donors.
