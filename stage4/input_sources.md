# 数据来源与文库结构

本轮下载 GSE327337 的 15 个已处理文件，每个文库各有一个特征表、条形码表和基因 × 条形码整数计数矩阵。下载的字节数均与 GEO 文件清单一致，SHA-256 记录于 `inputs_manifest.tsv`。全部文库使用相同的 37,487 个唯一 Ensembl 标识，对应 37,460 个不同的提交者基因名，重复命名涉及 27 个符号。条形码表在每个文库内只含一列。论文报告其自身质控后分析了 65,032 个细胞；本轮独立质控的起点是提交矩阵中的 77,638 个条形码。[GEO 文件清单](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327337/suppl/filelist.txt) · [原文](https://doi.org/10.1177/00220345261491310)

| 文库 | 组织条件 | GEO 样本 | 条形码数 |
|---|---|---|---:|
| CP | 对照 | GSM9653908 | 16,629 |
| PP3D | 压力侧第 3 天 | GSM9653909 | 23,289 |
| PP7D | 压力侧第 7 天 | GSM9653910 | 17,775 |
| TP3D | 张力侧第 3 天 | GSM9653911 | 11,583 |
| TP7D | 张力侧第 7 天 | GSM9653912 | 8,362 |

GEO 的样本信息说明使用 BD Rhapsody 捕获细胞、BGI T7 测序，合并唯一分子标识符计数后得到 GRCh38 表达矩阵。ENA 在同一 BioProject 下记录了五个测序运行，各自对应上述一个 GEO 样本。每个条件均只有一个可区分的文库，未从公开记录确定各文库的供者人数和是否混池。[GEO 样本信息](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327337/matrix/GSE327337_series_matrix.txt.gz) · [ENA 项目](https://www.ebi.ac.uk/ena/browser/view/PRJNA1450544)

GEO 文库级元数据将五个文库的细胞类型都填为“Goblet cell”。原论文报告 15 个细胞群，文库级标签无法用于逐细胞分类。GEO 文件清单和论文数据可用性声明未指向带标签的 RDS、h5ad 或条形码到受试者的映射。SAGE 列出一份补充 DOCX，本次访问返回 HTTP 403。对 Figshare 的资源 DOI 与题名、DataCite 的登录号及 GitHub 公开仓库的核查，未找到可核验的作者注释对象。图中的细胞分类采用本轮独立标志基因定义。每库的 barcode 只用于归属本库，任何跨库比较均保留文库维度。[作者数据可用性声明与补充材料入口](https://journals.sagepub.com/doi/10.1177/00220345261491310)
