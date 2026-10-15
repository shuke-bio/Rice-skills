# 水稻数据库资源导航（v3.0）

> 本包各分析模块的**外部数据库配套清单**：何时用哪个库、怎么接入。
> 网址可能随版本更新变化——若失效请按库名检索；重要分析建议同时核对多个来源。

## 1. 参考基因组与基因注释

| 数据库 | 网址 | 内容 | 用途 |
|---|---|---|---|
| **RAP-DB** | rapdb.dna.affrc.go.jp | IRGSP-1.0（日本晴）基因组 + 注释；**本包内置其 MSU↔RAP-DB ID 映射快照** | 基因 ID 体系基准；`assets/annotations/` 已内置 |
| **RGAP / MSU** | rice.uga.edu | MSU7 注释（水稻注释计划，已迁至 UGA）；**IRGSP-1.0 R7 数据下载：rice.uga.edu/download_osa1r7.shtml** | 老文献基因号（LOC_Os）核对；基因组/注释/重复注释批量下载 |
| **Ensembl Plants FTP** | jun2026-plants.ensembl.org/info/data/ftp | Ensembl Plants 全量数据 FTP（FASTA/GTF/VCF/比对） | 脚本化批量下载参考数据（wget/curl 友好） |
| **Gramene** | gramene.org | 比较基因组学 + 通路 + 多样性门户 | 跨物种比较、通路注释 |
| **Ensembl Plants** | plants.ensembl.org/Oryza_sativa | 基因组浏览器 + BioMart 批量导出 | 批量序列/注释下载（BioMart） |
| **Phytozome** | phytozome-next.jgi.doe.gov | 多物种基因组门户 | 同源基因序列获取 |

> 本包 `assets/annotations/rapdb_to_msu_map.csv`（34,133 对 ID 映射）用于模块 1/7 的注释转换；
> 更新用 `python scripts/utils/fetch_rice_annotations.py --force`。

## 2. 群体遗传与 GWAS 数据库

| 数据库 | 网址 | 内容 | 用途 |
|---|---|---|---|
| **RiceVarMap2** | ricevarmap.ncpgr.cn | 4,726 份材料的变异图谱与单倍型-表型关联 | **候选 SNP→单倍型分析**；查自己位点的已知关联 |
| **RFGB v2** | rmbreeding.cn | 3K-RGP 基因型-表型育种平台（含 **Rice3KGS** 基因组选择模块） | 单倍型/育种值查询；GS 模型参考 |
| **Rice SNP-Seek** | snp-seek.irri.org | 3,000 份 Rice Diversity Panel 变异与表型 | 3K 材料位点核对、多样性查询 |
| **3K-RGP 原始数据** | 各镜像站 | 3,010 份重测序（Nature 2018） | 群体遗传参考面板下载 |
| **RPAN（3K Rice DB）** | cgm.sjtu.edu.cn | 3K 泛基因组浏览器（PAV/SV 展示） | 模块 11 泛基因组结果的在线复核 |

**典型用法（与本包衔接）**：

- 模块 7/11 GWAS 命中的 lead SNP → RiceVarMap2 查单倍型与已知表型关联；
- 模块 11 PAV/稀有基因 → RPAN 看该基因在 3K 材料中的存在分布；
- 模块 8 群体结构结论 → SNP-Seek 的群体分类（indica/japonica/aus…）交叉验证。

## 3. 功能基因与表达数据库

| 数据库 | 网址 | 内容 | 用途 |
|---|---|---|---|
| **funRiceGenes** | funricegenes.github.io | 已克隆功能基因（>6,000）与克隆状态 | 候选基因**是否已克隆/已知功能**速查 |
| **RiceXPro** | ricexpro.dna.affrc.go.jp | 全生育期/组织表达谱（芯片） | 候选基因表达模式；模块 1 结果的对照 |
| **RiceFREND** | ricefrend.dna.affrc.go.jp | 共表达网络 | 候选基因的共表达模块（对照模块 1 WGCNA） |
| **IC4R / RED** | ic4r.org | 表达、变异、通路整合库 | 综合查询 |
| **Expression Atlas** | ebi.ac.uk/gxa/plants | 多物种 RNA-seq 表达图谱 | 跨条件表达查询 |
| **Rice eFP Browser** | bar.utoronto.ca | 电子荧光图像表达浏览 | 可视化表达模式（演示/制图） |
| **PlantTFDB** | planttfdb.gao-lab.org | 转录因子分类与结合基序；**本包内置水稻 TF 列表快照（1,862 TF / 56 家族）** | 模块 1 的 TF 分析、motif 注释 |

## 3b. 蛋白注释与 ID 映射（v3.1 补充）

| 数据库 | 网址 | 内容 | 用途 |
|---|---|---|---|
| **UniProt 水稻参考 proteome** | uniprot.org/uniprotkb/Q67UK9/entry（单条目）；**UP000059680** = 日本晴参考蛋白组 | UniProt ID ↔ 基因名 ↔ RAP-DB/MSU 交叉引用 | 蛋白组模块（模块 14）的 ID 映射源头；`protein_gene_map.csv` 的权威来源 |
| **InterPro**（EBI） | ebi.ac.uk/interpro/protein/UniProt/Q2QQV9/entry/InterPro | 蛋白结构域/家族/GO（InterProScan 结果） | 蛋白功能注释；差异蛋白的结构域富集；与 STRING/PPI 互补 |

## 4. 表观基因组与单细胞/空间资源（v3.0）

| 资源 | 网址 | 内容 | 用途 |
|---|---|---|---|
| **RiceENCODE** | ricencode.cn | 水稻表观基因组（ATAC、组蛋白修饰、甲基化） | 模块 15 表观组结果的公共数据对照 |
| **RSTA（Rice Spatio-Temporal Atlas）** | github.com/RiceSpatioTemporalAtlas/RSTA | 全生命周期时空图谱（78 万单核 / 202 细胞类型）+ **RICE scGPT 模型** | 模块 17/18 单细胞与空间分析的参考资源与基础模型（详见 `docs/spatial_usage.md`） |
| **单细胞图谱（多团队）** | 各论文附件 | 根尖/叶片/花药等组织 scRNA-seq | 细胞类型 marker 参考 |

## 5. 通路与代谢资源

| 数据库 | 网址 | 内容 | 用途 |
|---|---|---|---|
| **KEGG** | genome.jp/kegg | 通路图谱（水稻代码 `osa`） | 模块 1/2/3 富集分析（本包已集成） |
| **ROAD（Rice Oligonucleotide Array Database）** | roadv2.khu.ac.kr/rice/kegg/ | 水稻 KEGG 通路面板 + 芯片/表达整合 | 通路可视化与基因-通路映射（KEGG osa 的水稻本地化视图） |
| **PlantCyc / RiceCyc** | plantcyc.org | 植物代谢通路 | 代谢组通路补充注释 |
| **HMDB / MoNA** | hmdb.ca | 代谢物光谱与分类 | 代谢物鉴定核对（模块 2） |
| **MapMan** | mapman.gabipd.org | 植物功能分类与通路绘图 | 代谢/转录数据的 MapMan 注释 |

## 5b. 富集与 motif 在线分析工具（v3.1 补充）

| 工具 | 网址 | 内容 | 用途 |
|---|---|---|---|
| **agriGO v2** | systemsbiology.cau.edu.cn/agriGOv2（水稻粳稻 SpeciseID=18） | 植物专用 GO/SEA/PageGap 富集（水稻内置背景） | 本包富集结果的**在线交叉验证**；无 R 环境时的替代 |
| **PlantPAN 4.0** | plantpan.itps.ncku.edu.tw/plantpan4 | 植物启动子分析：TF/motif 数据库 + 启动子扫描 + 共表达调控对 | `01_promoter_kmer.py` 检出的 top k-mer 的**已知 TF motif 匹配**（替代 TOMTOM 的植物专用方案） |

## 6. 常用检索策略（省时提示）

| 场景 | 步骤 |
|---|---|
| 快速了解一个未知基因 | funRiceGenes 查克隆状态 → RAP-DB 看结构与表达 → RiceXPro/eFP 看模式 → RiceFREND 看共表达伙伴 |
| 候选 SNP 找基因 | RAP-DB 按坐标检索 → 本包 `scripts/gwas/04_candidate_loci.R` 批量求交 → funRiceGenes 过滤已克隆 |
| 位点已知关联 | RiceVarMap2 / SNP-Seek 输入基因或区域 → 看全部已知关联与单倍型 |
| 泛基因组核对 | 基因 ID → RPAN / 3K 泛基因组 → 存在/缺失与结构变异 |

> 提示：所有在线数据库的结果以**复核为目的**；本包的本地产出（GWAS/eQTL/PAV 等）应作为主结果，
> 在线库用于交叉验证与已知信息补充。
