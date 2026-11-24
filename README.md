# Averaging matrix gene responses masks prediction errors in mechanically loaded human periodontal ligament

Shangmin Zhang | ORCID [0009-0003-1846-9333](https://orcid.org/0009-0003-1846-9333)

Analysis software and worked example accompanying the manuscript. This repository contains the code, example data and intermediate inputs included in Supplementary Software 1. The complete Supplementary Data 1 and 2 archives are supplied with the manuscript and are not mirrored here.

This software evaluates a fixed panel of 101 matrix genes. For each contrast it calculates the absolute error of the mean signed change and the mean absolute error across individual genes. The values are changes in within-sample expression percentiles unless another scale is explicitly identified. Positive gain over a forecast of zero change means lower error.

## Run the clinical example

Download or clone this repository to a writable directory and open a terminal there. With Python 3.12 installed, the example runs immediately using the Python standard library. No package installation, network connection, GPU or specialised hardware is needed. Use this command from the repository root.

```text
python score_contrast.py stage19/figures/Figure_1_source_data.csv --forecast atlas
```

The input contains one row per target gene for the Leuven day-7 group contrast. The model and zero-change forecasts should give these values (percentile points).

```text
forecast     panel_error_pp   mean_gene_error_pp
zero         0.628607         6.172725
atlas        0.404912         6.536164
Genes closer than zero: 41; farther: 60; ties: 0.
```

## Score your own observations

Prepare a UTF-8 CSV with the columns `gene,observed,predicted`, one row for each of the 101 gene symbols in the example. Each file should contain one contrast. Supply signed loaded-minus-reference changes in percentile points. Use the same gene background for observations and predictions when ranking the two profiles. Keep full precision in the CSV. The scorer checks panel membership, duplicate or missing rows and nonfinite values.

```text
python score_contrast.py your_contrast.csv --forecast predicted
```

The script scores already computed rank changes. To change the model inputs or regenerate predictions, see the study-specific scripts and their input files below. Means across individuals in paired designs are calculated after each individual's errors have been scored.

## Rebuild the reported results

For full reconstruction, obtain Supplementary_Data_1.zip and Supplementary_Data_2.zip supplied with the manuscript, and extract both into the repository root, allowing the stage folders to merge. The three original ZIPs have disjoint file paths. Supplementary Data 1 holds clinical inputs and reference profiles. Supplementary Data 2 holds additional tissue and culture benchmarks. This repository retains the worked example data, analysis intermediates and file manifest originally supplied in Supplementary Software 1. The two complete data archives are manuscript attachments; the underlying source datasets can also be obtained from the public accessions listed in the manuscript and study-specific documentation. Access to manuscript attachments follows the manuscript's review or publication arrangements. There is no single fetch command for rebuilding all archived inputs from their public sources.

The SHA-256 manifest describes the merged software and two data archives. The standalone example above runs with this repository alone; the full reconstruction commands below require both data archives.

From that directory, `python reproduce.py` checks SHA-256 digests and reconstructs the data for the three main figures and seven supplementary tables using the archived numerical inputs. For the principal cross-laboratory culture comparison, use the following commands to rebuild from archived counts and to redraw all three figures.

```text
python -m pip install -r requirements-figures.txt
python reproduce.py --from-counts --figures
```

Outputs include `stage19/figures/Figure_1.pdf` to `Figure_3.pdf`, their PNG and CSV counterparts, and seven CSVs in `stage19/supp/data`. Other input matrices supplied with the software are processed gene-level or tissue-level data. The extraction, scoring and plotting steps above can all run offline after the packages have been installed. For a repeated integrity check after running the numerical pipeline in a different environment, use a fresh extraction because text representations of floating-point values can differ in their final digits.

## Environment and timing

The original software archive records the following environment and timing information. The independent check performed for this repository is documented separately in `REPRODUCTION_CHECK.txt`; that check did not rerun plotting.

The example and reconstruction commands were verified on 64-bit Windows 10, build 19045, using Python 3.12.14 on an Intel Core Ultra 9 275HX workstation with 16 GB nominal RAM. The numerical and plotting workflow used NumPy 2.3.5 and Matplotlib 3.11.2. With Python preinstalled, extracting the software ZIP took about 0.22 seconds and the clinical example took about 0.05 seconds including Python startup. The full count-based reproduction and plotting run took about 4.5 seconds with its dependencies installed. Detailed timing measurements are provided in `VERIFICATION.txt`. Runtime depends on storage, hardware and whether the plotting library has already built its font cache. Dependency installation time on another computer depends on package caches and the network. The small demonstration requires only archive extraction with an existing Python installation.

Optional atlas preprocessing and its exact software versions are documented in `stage4/ATLAS_QC.md`. These steps are separate from the offline demo and are unnecessary for reproducing the figures from the archived inputs.

## Access and reuse

Code is publicly accessible. No public reuse licence has been assigned to the custom analysis scripts. Public access does not grant a general licence to reuse or redistribute the code.

Third-party data, annotations, publications and software retain the terms of their original providers. Public source accessions are listed in the manuscript and in the study-specific documentation.

## Analysis files
- `stage5/` contains the atlas reference, target panel, clinical inputs and tissue forecast implementation. The atlas quality-control and reference data retain the source library labels. The article reports a secondary analysis, and this archive uses processed source material at the stage indicated by each module.
- `stage11/specificity/rank_changes.npz` contains the frozen clinical predictions and observed changes. `stage12/review/` scores zero change, the atlas and a shift shared across genes. `stage12/control/` evaluates the constructed mixtures of held-out cells.
- `stage13/transfer/` contains the transfer of measured changes between clinical cohorts. `stage13/external/` contains the five donor-paired tensile cultures, their public count matrix and donor-level analyses.
- `stage14/` contains the six-hour tensile comparisons and the compression studies. `stage15/` contains the clinical ranking and calibration analyses. `stage16/sides/` evaluates which side-specific predictions are distinguishable in pooled tissue.
- `stage17/wuhan_counts/` and `stage17/zhejiang/` contain the two external stretch series. The 24- and 48-hour Wuhan contrasts share the same baseline samples. `stage17/clinical_search/` contains the occlusal-contact comparison, and `stage17/axial/` the axial-loading comparison.
- `stage18/transfer/`, `stage18/deviations/` and `stage18/bounds/` score the measured Hong Kong transfers, the centred calibration and the bounded predictions. `stage18/check/` independently reconstructs the shared-background forecasts and checks them against the public counts.
- `stage19/figures/` contains data extraction and plotting code. `stage19/supp/` contains the source tables and their reconstruction script. The Source_Data.zip file contains copies organised by figure and table for direct inspection.
