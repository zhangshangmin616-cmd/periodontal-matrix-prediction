# Frozen day-3 tensile atlas forecasts in the Xi'an six-hour cultures

This analysis asks whether retaining the measured cell-group expression responses but fixing cell-group RNA shares, or retaining only the fibroblast-like group's responses, improves on predicting no change. It evaluates the identical two day-3 human tensile-side atlas predictors previously used for the Hong Kong culture. Xi'an primary human periodontal ligament stem cells were cultured under 12% equibiaxial cyclic stretch for 6 hours. The three unloaded libraries were shared across the 0.1 and 0.5 Hz comparisons. The public record identifies two and three loaded libraries, respectively. There is no recoverable donor-to-run pairing. Bulk culture measures a group-level response and does not resolve an individual fibroblast's response. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9259061/) · [Public project](https://www.ebi.ac.uk/ena/browser/view/PRJNA665587)

| Endpoint, across 101 genes | Frequency | Zero change | Fixed-share tissue | Fibroblast-like |
|---|---:|---:|---:|---:|
| Mean absolute error, log₂(CPM+1) | 0.1 Hz | 0.270 | 0.440 | 0.706 |
| Mean absolute error, log₂(CPM+1) | 0.5 Hz | 0.627 | 0.680 | 0.868 |
| Mean absolute error, percentile points | 0.1 Hz | 2.718 | 4.389 | 6.832 |
| Mean absolute error, percentile points | 0.5 Hz | 7.279 | 7.936 | 9.763 |

Both models raised gene-level error relative to zero at both frequencies and on both scales. Every relevant single-run deletion retained that result (five overlapping deletions for 0.1 Hz and six for 0.5 Hz). The absolute error in mean panel response sometimes fell, for both models at 0.1 Hz on the rank scale and at 0.5 Hz on the log scale. The comparison evaluates a cross-time, tissue-to-culture forecast. It cannot establish performance for human day-7 root-face tissue or identify whether the transfer error is due to tissue context, strain, cell population, duration or assay.

Run `python stage14/atlas_xian/score.py` then `python stage14/atlas_xian/check.py`. The analysis reuses the stage 13 model definition and the processed Xi'an input. Independent checks use the frozen stage 12 per-gene predictions, the original reference counts and tissue shifts, and the separately generated Xi'an target measurements. Outputs include 12 full comparisons, 66 descriptive single-deletion comparisons, and 1,212 per-target observed and predicted values.
