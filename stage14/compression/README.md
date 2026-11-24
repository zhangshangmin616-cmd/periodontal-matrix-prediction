# Transfer across 24-hour static compression experiments

These three centre-specific human periodontal ligament culture experiments previously contributed to the measurable universe or assessment of extracellular matrix subpathways. All 101 prespecified targets and all 10,806 background genes were present in all three matrices. The six directed comparisons are retrospective and reuse three group-level contrasts. There is one reference vector per study.

| Centre and accession | Data used |
|---|---|
| Bangkok, GSE112122 | Three independently identified donors, two samples per donor (static arm and its matched control), 24 h of continuous 1.5 g/cm² compression, under serum-free conditions. Samples GSM3058133–GSM3058138. The intermittent arm and its separate control are excluded. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112122) · [Primary data article](https://doi.org/10.1016/j.dib.2019.104553) |
| Beijing, GSE109167 | Human periodontal ligament stem cells, two unloaded libraries and two loaded libraries after 24 h of 1.5 g/cm² compression. GSM2934477–78 are the controls and GSM2934479–80 the loaded samples. The published study specifies two RNA replicates per arm. We calculated a group contrast only. NCBI's series-level human gene counts provided all four profiles. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109167) · [Primary study](https://doi.org/10.1038/s41419-020-02963-3) |
| Xi'an, GSE347200 | Primary human periodontal ligament stem cells, three biological replicates per arm, 24 h under a coverslip plus 2.0 g/cm², with a coverslip-only control. GSM10047288–90 correspond to C_1–C_3 and GSM10047291–93 to A_1–A_3. The GEO description identifies the sample arms and no donor-resolved pairing was assumed. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347200) · [Primary study](https://doi.org/10.3390/antiox15101292) |

| Source → destination | Target gene error, zero | Target gene error, full | Target gene error, quarter | Panel error, zero | Panel error, full |
|---|---:|---:|---:|---:|---:|
| Bangkok → Beijing | 5.070 | 8.421 | 5.582 | 2.599 | 2.263 |
| Bangkok → Xi'an | 2.375 | 6.200 | 2.868 | 0.857 | 1.193 |
| Beijing → Bangkok | 5.502 | 8.421 | 5.875 | 0.336 | 2.263 |
| Beijing → Xi'an | 2.375 | 5.768 | 2.784 | 0.857 | 3.456 |
| Xi'an → Bangkok | 5.502 | 6.200 | 5.524 | 0.336 | 1.193 |
| Xi'an → Beijing | 5.070 | 5.768 | 5.182 | 2.599 | 3.456 |

The table reports percentile points. The complete source gene response gave larger gene-level errors than zero in all six transfers; the fixed quarter response did so as well. Both results persisted on the log₂(CPM+1) scale and with the rank predictions clipped to their marginal feasible range in the destination. Panel-mean error nevertheless decreased in one full transfer and in two quarter transfers on the original rank scale. On the whole 10,806-gene background, quarter transfer gave lower gene-level error than zero in two of six directions on each scale. These outcomes concern the specific genes, endpoint and protocols evaluated.

Run `python stage14/compression/score.py` and `python stage14/compression/check.py` with NumPy to reproduce all tables from the compact count matrices. The independent checker recomputes all 120 comparisons using a second implementation of tied ranks and errors. It also tests bounds for 12 gene-specific rank predictions in a sensitivity analysis added after the main scores were known.

To reconstruct `counts.npz`, download the three public sources listed in `sources.json`, verify SHA-256 and run `python stage14/compression/extract.py --raw-dir DIRECTORY`. All three downloads were rechecked against the recorded hashes. The extractor uses the pre-existing `stage5/paired/hgnc_ids.tsv.gz` registry and `stage5/reference/features.tsv` panel and background. It retains total mapped gene counts for normalisation, including contributions outside the common background. Archive files and output values do not require sending or requesting any information externally.

## Audit of the other Beijing accession

GSE116787 contains one compressed and one control library. The matching original paper specifies 12 h at 2 g/cm². The public GEO page omits the duration, its title and contributors match the paper, and a later paper explicitly links the RNA-seq data to this accession. ENA shows SRR7499436–37 and SAMN09631441–42. In GSE109167 the identifiers are SRR6470156–59 and SAMN08361218–21. All FASTQ MD5 values differ between these studies. Donor identity across the accessions remains unresolved. GSE116787 already contributed to the stage 3 common-gene intersection and had published-in-package pathway summaries in `culture_module_summary.tsv`. It therefore cannot supply a new 24 h evaluation. [Original experiment](https://doi.org/10.1093/ejo/cjy057) · [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116787) · [Later use of the accession](https://pmc.ncbi.nlm.nih.gov/articles/PMC13104472/) · [ENA 12 h project](https://www.ebi.ac.uk/ena/browser/view/PRJNA480245) · [ENA 24 h project](https://www.ebi.ac.uk/ena/browser/view/PRJNA429831)

An external 24 h culture destination, limited to two of the frozen genes in the thresholded public table, is documented in `external3d/README.md`.
