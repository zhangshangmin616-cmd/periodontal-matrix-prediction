# Transfer of a fixed orthodontic-force predictor to natural occlusal contact

GSE189288 provides a separate human PDL RNA dataset with a natural difference in loading. Chang et al. enrolled ten patients with contact between opposing premolars on one side and absent contact on the other. RNA was obtained directly from the root middle third. The cell-culture and knockdown experiments in the paper were separate. The paper reports 34 teeth, 17 of each occlusal status. The archive contains 22 sequencing libraries, 11 per status, with some sample titles containing two tooth positions. The 34 teeth and the 22 libraries do not constitute 34 or 22 independent patients. The duration and magnitude of natural occlusal loading are unknown. [Author methods](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2021.796274/full) · [GEO study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189288)

## Clinical comparison and independence

The experimental unit was a pair of occlusal states within a person. The primary article describes the pairing, but the GEO characteristics fields do not include a separate donor-number field. We reconstructed ten matched units from the sample-title prefixes and verified that age and conditions were internally consistent in each. Four early units use the shared extraction-date prefix, one uses the shared name prefix and the remaining five use an explicit S code. One unit has two libraries on each side. Nine have one per side. In the released table units are labelled pair01 to pair10. Laboratory identifiers and GEO accessions permit checking the grouping. The two similarly aged 14-year-old units have distinct prefixes. The score first averages libraries within each unit and state, then gives each unit equal weight. [GEO sample metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189288/matrix/GSE189288_series_matrix.txt.gz)

The previously used three-week orthodontic trial at Yonsei and this observational study had distinct exposures and ethics approvals, respectively 2-2017-0028 and 2-2018-0002. Comparing the respective 29 and 22 libraries showed no shared accession for runs, experiments or BioSamples, no shared provider library identifier, read count, ENA file-MD5 pair, or digest of the first 20 forward-read sequences. Enrollment windows overlapped in March–May 2018. No cross-study participant key is publicly reported, so these checks establish distinct library sets and do not prove that every person was unique across the two studies. [Orthodontic clinical trial](https://doi.org/10.3390/jcm9020358) · [Orthodontic project](https://www.ebi.ac.uk/ena/browser/view/PRJNA698902) · [Occlusal project](https://www.ebi.ac.uk/ena/browser/view/PRJNA782351)

## Scores

The direction, forecast, gene panel, background and both loss functions were recorded before scoring. The contrast is contact minus hypofunction. The model uses the original day-7 changes in expression within cell groups, keeps baseline RNA shares fixed, and averages the two root-side forecasts in linear abundance. For each hypofunction baseline, all other participants' hypofunction profiles supply the calibration. The observed endpoint is a within-person difference between anatomical specimens. It does not measure longitudinal change, causal benefit of restoring contact, or a seven-day orthodontic intervention. The rank background contains the full 10,806 locked genes, including all 101 targets.

| Error in expression percentile points | Fixed model | Zero change | Gain over zero |
|---|---:|---:|---:|
| Panel-mean absolute error, averaged across 10 matched units | 1.0664 | 1.0633 | −0.0031 |
| Per-gene absolute error, averaged across genes and units | 4.7964 | 4.2650 | −0.5313 |
| Panel-mean error of the group contrast | 0.3314 | 0.2375 | −0.0939 |
| Per-gene error of the group contrast | 2.7791 | 1.9928 | −0.7863 |

Six matched units improved on the panel average, and one improved on individual genes. Leaving out a unit and recalibrating the remaining baselines gave panel gains between −0.034 and +0.027, and per-gene gains between −0.593 and −0.458. A supplementary percentile bootstrap over matched units gave conditional 95% intervals of −0.100 to +0.093 for panel gain and −0.730 to −0.318 for per-gene gain. The bootstrap holds fitted predictions fixed. We enumerated the 92,378 possible resampling multiplicity vectors and weighted them by their multinomial counts. The pooled scores first average participants' gene changes and then calculate each endpoint, so they are a different estimand. The paired averages show no detectable panel-level advantage; they do not establish equivalence within a prespecified margin.

## Other human in vivo force studies

Broad English, Chinese and Portuguese searches of original experimental force and general mechanical loading studies recovered two close clinical candidates whose released outcomes have previously been audited. Ratanasereeprasert et al. measured longitudinal human PDL in an extraction series; the published quantitative RNA table reports selected genes and baseline controls are third molars. In Spitz et al., the 331 reported RNA effects in the conventional orthodontic arm are selected by effect size and adjusted significance. The paper's supplement and the investigator's thesis describe the same observations. The Chengdu extrusion/intrusion pilot, Leuven orthodontic cohort and Seoul three-week trial are already part of the project. No new complete intervention dataset was established by this search. [Ratanasereeprasert et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10829649/) · [Spitz et al.](https://doi.org/10.1007/s00784-021-04178-y) · [Chengdu pilot](https://doi.org/10.1111/ocr.70103)

## Reproduction

The working raw-count matrix is the [NCBI-generated GRCh38.p13 file](https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE189288&format=file&file=GSE189288_raw_counts_GRCh38.p13_NCBI.tsv.gz) (SHA-256 b3883a9f5e88775aadc473afd9905f6c6ff90d583bd442d31b24bd51c450c492). Source links for all input files are in `sources.py`. Run the following from the extracted cumulative archive with Python 3.10+ and NumPy.

```text
python stage17/clinical_search/score.py
python stage17/clinical_search/uncertainty.py
python stage17/clinical_search/check.py
```

`check.py` independently constructs midranks, recomputes response predictions from the integer atlas count matrices, checks optimality conditions for each baseline RNA-share estimate, and verifies the per-unit observations and scores. The verification is offline. Optional live recovery is available through `sources.py --outdir PATH` for the GEO/ENA sources and `sources.py --reads --outdir PATH` for a hash of 20 forward reads from each ENA run. These fingerprints check for duplicate files and are not a test of personal identity. Human read sequences are not stored in the fingerprint table.
