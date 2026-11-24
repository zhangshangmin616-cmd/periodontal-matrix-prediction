# Regional search for a second human day-7 periodontal ligament RNA cohort

Search completed on 25 November 2026 (Asia/Shanghai). This extension examined DDBJ Search and Chinese-, Japanese- and Korean-language literature, thesis, and grant sources. The question was whether a different human cohort had day-7 RNA measurements from separately sampled compression and tension regions of the periodontal ligament after clinical loading. Testing the full fixed panel requires 101 gene measurements. Testing cell-group-specific change further requires cell identities or a corresponding physically isolated population.

## Repository search

Three terms were queried independently against each of four DDBJ Search record types, without an organism filter. All pages were downloaded and unique accessions were reconciled against the reported totals. Counts below represent each query before combining overlapping results. Cross-links between BioProject and SRA often refer to the same samples. [DDBJ Search API](https://ddbj.nig.ac.jp/search/api/openapi.json)

| DDBJ record type | orthodontic | “periodontal ligament” | “tooth movement” |
| --- | ---: | ---: | ---: |
| GEA | 0 | 0 | 0 |
| JGA Study | 0 | 0 | 0 |
| SRA Study | 68 | 89 | 21 |
| BioProject | 79 | 153 | 25 |

There were 217 distinct BioProject accessions and 149 distinct SRA Study accessions. Positive controls for the two empty search routes returned 75 GEA records and 192 JGA study records for “cancer”. The queries recovered the known Nanjing human side-specific atlas, the Leuven human clinical bulk benchmark, and the previously audited Japanese rat PDL time course. All 366 metadata records and their query membership are retained for review.

## Resolution of close matches

* Tsurumi University's publicly catalogued microarray project describes microdissection of molar PDL during tooth movement. The final report specifies that rat molars were exposed to force for 6 hours or 5 days. [Grant report](https://kaken.nii.ac.jp/en/report/KAKENHI-PROJECT-17592153/175921532006kenkyu_seika_hokoku_gaiyo/) · [CiNii catalog](https://ci.nii.ac.jp/ncid/BA84764690)
* Kikuta and colleagues measured a panel gene, TGFB1. Day-7 pressure-side staining came from rats. Human expression was examined in compressed cells in culture. [Original study](https://www.jstage.jst.go.jp/article/ijoms/19/3/19_193/_article/-char/en)
* A study of BMP signaling in human PDL measured another panel gene, BMP4, by qPCR in cells subjected to 24 hours of tensile strain. Its clinical arm sampled teeth after 3–8 weeks for immunohistology. [Original study, Methods and Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC11640950/)
* Taiwan's national thesis catalogue provides the original account of 24- and 48-hour cyclic stretching of human PDL. A Chosun University thesis labels pressure and tension surfaces but measured proteins in rats at days 1, 3, and 5. [Taiwan thesis](https://ndltd.ncl.edu.tw/cgi-bin/gs32/gsweb.cgi?o=dnclcdr&s=id%3D%22095NTU05594015%22.&searchmode=basic) · [Korean thesis](https://oak.chosun.ac.kr/handle/2020.oak/6532)
* A Guangzhou study sampled human gingival crevicular fluid after one week. Its experiments on isolated human PDL cells tested a cytokine treatment. A separate, previously checked Nanjing human qPCR study compared loaded and contralateral control premolars. Its reported human target outcomes PIEZO1, IL1B, TNF and ZBP1 have zero overlap with the fixed panel. [Guangzhou study](https://www.ejpd.eu/wp-content/pdf/EJPD_2021_22_02_09.pdf) · [Nanjing study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13239690/)

No second day-7 human pressure–tension PDL RNA cohort suitable for assessing the fixed 101 genes was identified by this regional search. No eligible new donor-level contrasts were available to calculate gene coverage or prediction error on an independent cohort. Ineligible studies may measure genes in the panel in other conditions. The search depends on indexed terminology and public disclosures. The Korean, Taiwanese and Japanese literature searches were targeted primary-source checks and do not specify complete national holdings.

## Audit files

* `ddbj_queries.tsv` preserves the query URLs, counts and page totals. `ddbj_records.tsv` preserves the deduplicated titles, organisms, accessions and query membership. `ddbj_controls.tsv` records the GEA and JGA positive controls.
* `sources.tsv` records the clinical timing, sampling unit, assay and decision with a primary source. `discovery.tsv` records representative regional search terms and primary records checked.
* Run `python check.py` for offline verification of this snapshot. Run `python ddbj_search.py` to retrieve the current DDBJ index. Live totals may change.
