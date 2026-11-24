# Original atlas quality control

The analysis used five human periodontal ligament libraries from GEO series GSE327337. The 15 input names, public download URLs and original SHA-256 checksums are recorded in `inputs_manifest.tsv`. The expression matrices themselves are available from the public GEO deposit. The per-library thresholding script and doublet classification script are included here. Precomputed downstream inputs for the manuscript's clinical and culture analyses are supplied with the three accompanying archives.

To rerun this optional preprocessing, place the 15 compressed files listed in the manifest under `stage4/input/`. Preserve all filenames. From the extracted archive root run:

```text
python stage4/prepare_counts.py
Rscript --vanilla stage4/doublets/run_doublets.R CP stage4/cache stage4/doublets
```

Repeat the R command for PP3D, PP7D, TP3D and TP7D. The Python command uses NumPy, pandas and SciPy. The original R session used R 4.5.3, Matrix 1.7-4, SingleCellExperiment 1.32.0 and scDblFinder 1.24.10. The per-library R session records and package version list accompany the scripts. Read and write access to the extracted directory is needed. A full matrix rerun uses substantially more memory and disk space than the worked scoring example. Across the five original libraries, 77,638 barcodes were available, 71,070 passed count-quality filters and 62,490 remained after doublet removal. This optional source-data preprocessing was not repeated as part of the submission packaging run.
