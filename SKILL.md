---
name: rice-skills
description: 水稻（Oryza sativa）多组学分析 Skill（全组学平台）。覆盖转录组、代谢组及联合分析（数据质控、DESeq2/edgeR、GO/KEGG/Reactome 富集、GSEA、WGCNA、转录因子与 PPI；OPLS-DA/VIP、代谢通路、ROC/随机森林、代谢网络；MOFA/DIABLO/O2PLS）；基因组学（变异检测 BWA→GATK→SnpEff、GWAS 三引擎、Meta-GWAS、群体遗传、泛基因组 SV/PAV/Pan-GWAS、eQTL/mQTL QTL 映射）；蛋白质组（limma 差异蛋白、火山/热图）；表观基因组（WGBS/DMR、ATAC/ChIP 差异区域）；表型组（BLUP/遗传力、光谱指数、GS 指引）；单细胞与空间组学（Scanpy/Squidpy 流程、RSTA/RICE scGPT 集成指引）；引物与 CRISPR 靶点设计；配套分子实验协议参考。全部图表遵循 Nature 级绘图规范（低饱和色盲友好配色、Arial/Helvetica 7-9 pt、PDF/SVG/TIFF ≥300 dpi 矢量输出）。当用户提到水稻转录组、水稻代谢组、水稻GWAS、群体结构/群体遗传、泛基因组、eQTL、蛋白质组、表观组/甲基化/ATAC、表型/BLUP、单细胞/scRNA、空间转录组/Stereo-seq、引物设计、CRISPR 靶点、多组学联合分析、WGCNA、GSEA、DESeq2、OPLS-DA、MOFA、DIABLO、差异表达、富集分析等关键词时触发。
version: 3.2.0
license: MIT
species: Oryza sativa (水稻, Japonica)
trigger:
  - 水稻转录组
  - 水稻代谢组
  - 多组学联合分析
  - 转录组 代谢组 整合
  - WGCNA
  - GSEA
  - DESeq2
  - 差异表达分析
  - OPLS-DA
  - MOFA
  - DIABLO
  - 富集分析
  - KEGG 通路
  - 火山图
  - Nature 绘图
  - 水稻GWAS
  - 群体结构
  - 群体遗传
  - Fst
  - LD 衰减
  - 曼哈顿图
  - 系统发育树
  - 亲缘矩阵
  - VCF
  - PLINK
  - 引物设计
  - CRISPR 靶点设计
  - sgRNA
  - 基因敲除
  - 载体构建
  - 蛋白纯化
  - 亚细胞定位
  - GUS 染色
  - TAIL-PCR
  - 转录因子功能研究
  - 变异检测
  - 泛基因组
  - PAV
  - 结构变异
  - Pan-GWAS
  - eQTL
  - Meta-GWAS
  - 蛋白质组
  - 差异蛋白
  - 表观基因组
  - 甲基化
  - ATAC
  - ChIP-seq
  - DMR
  - 表型组
  - BLUP
  - 光谱指数
  - NDVI
  - 单细胞
  - scRNA-seq
  - 空间转录组
  - Stereo-seq
  - 单倍型
  - 共定位
  - GBLUP
  - 基因组选择
  - AMMI
  - GGE biplot
  - XP-EHH
  - 选择信号
  - KA/KS
  - 基因家族进化
  - 启动子
  - motif
  - k-mer 富集
  - miRNA
  - 可变剪接
---

# Rice-skills — 水稻全组学与基因组学分析 Skill

> 面向水稻（*Oryza sativa*）的**全组学分析全流程工具包（v3.0）**：转录组、代谢组及联合分析、
> 基因组学（变异检测/GWAS/群体遗传/泛基因组/QTL/Meta-GWAS）、蛋白质组、表观基因组、表型组、
> 单细胞与空间组学、引物与 CRISPR 靶点设计，输出符合 Nature 系列期刊规范的出版级图表。
> 物种可编辑：将 `org.Osativa.eg.db` / KEGG 代码 `osa` 等替换为目标物种即可复用整套流程。

## 1. 适用边界

本 Skill 是**面向全体水稻组学研究者的通用工具链**：数据格式通用、参数可配置，不绑定任何实验室、
平台或测序公司；仓库内示例数据均为合成数据（性状名为演示占位，可直接替换为你自己的数据）。

| 维度 | 支持 | 不支持 |
|---|---|---|
| 物种 | 水稻（默认 Japonica；可改） | —（改 3 处物种参数即可迁移） |
| 组学 | 转录组（bulk RNA-seq count 矩阵）、代谢组（广靶/靶向代谢物丰度矩阵）、两者联合、**基因组变异（变异检测/GWAS/群体遗传/泛基因组/eQTL，v2.0–v3.0）**、**蛋白质组（v3.0）**、**表观基因组（WGBS/ATAC/ChIP，v3.0）**、**表型组（BLUP/光谱/GS 指引，v3.0）**、**单细胞与空间组学（v3.0）**、**引物与 CRISPR 靶点设计（v2.1）** | 泛基因组图构建为工具链模板（重计算依赖外部服务器）；Hi-C 三维基因组 |
| 比较设计 | 两组比较、多组两两比较、多因素设计（~ batch + group）、多点表型（BLUP 后） | 专用时间序列模型（可用 group=时间点近似；单细胞有 PAGA 轨迹） |
| 平台 | R ≥ 4.3（主分析栈）+ Python ≥ 3.10（质控/预处理/可视化） | 单细胞/空间需独立 Scanpy/Squidpy 环境（文档指引） |

## 2. 输入数据格式

三个标准输入（CSV，UTF-8，逗号分隔，首行表头）：

**① 转录组 count 矩阵** `counts.csv` — 行=基因，列=样本，整数计数：

```
gene_id,CTRL_1,CTRL_2,CTRL_3,SALT_1,SALT_2,SALT_3
LOC_Os01g01010,1520,1380,1602,4210,3980,4150
LOC_Os01g01019,88,102,76,231,205,198
```

**② 代谢组丰度矩阵** `metabolites.csv` — 行=代谢物，列=样本；推荐附带注释列
`metabolite_name` / `class`（化学分类）/ `kegg_id`：

```
metabolite_id,metabolite_name,class,kegg_id,CTRL_1,CTRL_2,CTRL_3,SALT_1,SALT_2,SALT_3
M0001,Flavonoid-3-glucoside,Flavonoids,C12271,15200,14100,15900,38900,41200,37600
```

**③ 样本元数据** `metadata.csv` — 首列 `sample_id`（必须与矩阵列名完全一致），
**必含 `group` 列**（分组），可选 `batch`、`variety` 及 WGCNA 性状列：

```
sample_id,group,batch,variety,tiller_number,Na_content
CTRL_1,control,1,Nipponbare,12.4,0.31
SALT_1,salt,1,Nipponbare,8.1,1.87
```

> 样本 ID 大小写敏感、不含空格；矩阵样本列顺序无需与元数据一致（脚本自动对齐）。
> 缺失值：count 矩阵不允许；代谢组矩阵允许空值（由 `01_preprocess.py` 填充）。

## 3. 工作流总览

```
counts.csv ──┐
             ├─► [模块1 转录组] QC过滤 ─► DESeq2差异表达 ─► GO/KEGG富集 ─► GSEA ─► WGCNA ─► TF/PPI ─► 出版图
metadata.csv ┘
metabolites.csv ──► [模块2 代谢组] 预处理 ─► OPLS-DA/VIP ─► 通路富集 ─► 分类统计 ─► 生物标志物 ─► 代谢网络 ─► 出版图
                     │
                     ▼
        [模块3 联合分析] MOFA / DIABLO / O2PLS ─► 基因-代谢物相关性 ─► 联合通路 ─► 多组学网络 ─► 整合图
```

三模块可独立运行，也可由 `run_pipeline.sh` 一条命令串行执行。

**v3.0 全组学扩展**（各模块独立运行，均带合成数据 demo 与独立文档）：

```
FASTQ ─► [10 变异检测] BAM ─► VCF ─► 过滤 ─► SnpEff 注释 ─┬─► [7 GWAS] ─► [13 Meta-GWAS]
                                                          ├─► [8 群体遗传]
                                                          ├─► [11 泛基因组] SV/PAV ─► Pan-GWAS
                                                          └─► [12 QTL 映射] eQTL/mQTL/pQTL
蛋白矩阵 ─► [14 蛋白质组] 过滤/归一化 ─► 差异蛋白(limma) ─► 火山/热图
甲基化矩阵 / peak 计数 ─► [15 表观组] DMR / DAR ─► 火山/热图
多点表型 / 光谱 ─► [16 表型组] BLUP+遗传力 / 植被指数 ─► GWAS / GS
单细胞 / 空间数据 ─► [17/18] Scanpy / Squidpy 流程 ─► RSTA / RICE scGPT 集成
```

### 模块 7：GWAS 全链（scripts/gwas/，v2.0 新增）

| # | 脚本 | 步骤 | 工具 | 关键输出 |
|---|---|---|---|---|
| 01 | `01_plink_qc.sh` | VCF→PLINK 转换 + MAF/缺失率/HWE 质控 + LD 剪枝 + 剂量导出（模板，需 PLINK） | PLINK | PLINK fileset、`.prune.in`、`.raw` |
| 02 | `02_plink_to_matrix.R` | PLINK `.raw/.map` → 标准剂量矩阵 | R | `genotype.csv`、`snp_info.csv` |
| 03 | `03_association.R` | GWAS 关联（三引擎：simple 内置 OLS-GLM / rMVP FarmCPU·MLM / GAPIT），PCA 协变量 | 可选 rMVP(推荐)/GAPIT | 全量结果表、显著 SNP 表、PCA 得分 |
| 04 | `04_candidate_loci.R` | 显著 SNP→位点聚类→候选区间→与 GFF3 注释求交 | R | `gwas_loci.csv`、`candidate_genes.csv` |
| 05 | `05_plot_gwas.R` | 出版级 Manhattan / QQ / 组合总览图 | ggplot2+patchwork | 三类图 × 3 格式 |

### 模块 8：群体遗传（scripts/popgen/，v2.0 新增）

| # | 脚本 | 步骤 | 关键输出 |
|---|---|---|---|
| 01 | `01_pca_structure.R` | PCA 群体结构 + VanRaden 亲缘矩阵（可作 GWAS 协变量） | 得分/矩阵/碎石图/散点/总览图 |
| 02 | `02_ld_decay.R` | LD 衰减曲线 + 半衰减距离（r²=0.2 交点，指导候选窗口） | `ld_decay.csv`、衰减曲线图 |
| 03 | `03_selection_stats.R` | Hudson Fst（两组）+ 核苷酸多样性 π（窗口化） | Fst/π 表、全基因组扫描图、总览图 |
| 04 | `04_admixture.sh` | ADMIXTURE 祖先成分估计（K 扫描 + CV 误差，模板，Linux/WSL） | `.Q/.P` 矩阵、`cv_errors.txt` |
| 05 | `05_plot_admixture.R` | 祖先成分堆积图 + CV 曲线（任意 `.Q` 矩阵均可绘图） | K 图 × N、CV 曲线图 |
| 06 | `06_phylogeny.R` | IBS 距离 + NJ/UPGMA 树（Newick 导出，类别着色） | `tree.nwk`、系统树图 |

### 模块 9：设计与验证（scripts/design/，v2.1 新增）

| # | 脚本 | 步骤 | 关键输出 |
|---|---|---|---|
| 01 | `01_sgrna_design.py` | SpCas9 靶点设计：标准方案七原则全部代码化（CDS 优先 / GC 40-80% / 禁连续≥4T / 骨架配对≥7bp 剔除 / 中上游偏好 / BsaI 位点检查 / 序列内脱靶统计）+ **双靶点自动配对**（间隔 ≥100 bp）+ **pYLCRISPR/Cas9 构建引物自动生成**（gRT#+/U6#T#-、SP-L1/SP-R、通用引物表） | 候选表（含逐条通过/剔除原因）、推荐双靶点、构建引物、靶点分布图 |
| 02 | `02_primer_design.py` | 引物设计（primer3 引擎）：qPCR 跨内含子（逐结合点定向）/ CRISPR 靶点验证（±flank）/ 通用克隆；Tm、GC、发夹、同源/异源二聚体完整质量指标 | 引物表、质量报告、覆盖图 |

配套实验参考（三份文档互补）：
- `docs/lab_protocols.md`：靶点原则、双靶点构建流程、CTAB/TPS 提取、编辑体鉴定与解码（整理自公开发表的 pYLCRISPR/Cas9 协议并注明出处）；
- `docs/protocols_molecular.md`：引物设计评估、同源重组构建、RNAi 载体、蛋白纯化与 SDS-PAGE、烟草亚细胞定位、LUC 互补、mhiTAIL-PCR、GUS 染色、BioEdit 序列分析；
- `docs/research_workflow_tf.md`：转录因子功能研究全流程路线图（背景调查 → 功能验证 → 表达模式 → 分子机理，含代际操作要点与周期估算）。

### 模块 10–18（v3.0 全组学扩展）

| 模块 | 目录 | 核心脚本 | 说明 |
|---|---|---|---|
| **10 · 变异检测与注释** | `scripts/variants/` | 01_bwa_align.sh → 02_gatk_call.sh → 03_vcf_filter.sh → 04_snpeff_annotate.sh | FASTQ→BAM→VCF（GATK/bcftools）→硬过滤→SnpEff 注释；重工具模板式调用（用户自装，脚本自动探测+安装指引） |
| **11 · 泛基因组** | `scripts/pangenome/` | 01_pangenome_graph.sh（minigraph/cactus/PGGB）· 02_sv_detect.sh（minimap2+SYRI）· 03_pav_analysis.R · 04_pan_gwas.R | 图构建→SV→PAV（核心/可变/私有 + 性状关联 + Manhattan）→Pan-GWAS（复用 GWAS 引擎） |
| **12 · QTL 映射（eQTL/mQTL/pQTL）** | `scripts/qtl/` | 01_qtl_mapping.R · 02_colocalization.R · 03_coloc_bayesian.R · 04_qtl_region_snps.py | 自研 cis/trans 引擎（向量化 OLS + 保守 BH）、cis-距离图、Manhattan、每特征最佳命中；共定位筛查/PP4；**QTL/BSA 区间 SNP 提取 + 基因区域注释（区间 VCF + 三表 Excel）** |
| **13 · Meta-GWAS** | `scripts/gwas/06_meta_gwas.R` | 多研究逆方差合并 + Cochran Q/I² + 随机效应 + 森林图；配 `docs/rice_databases.md` 水稻数据库导航 |
| **14 · 蛋白质组** | `scripts/proteomics/` | 01_preprocess.py · 02_diff_protein.R（limma/Welch）· 03_plot_proteome.R | 过滤（≥70%）/归一化/批次/填充→差异蛋白→火山+热图 |
| **15 · 表观基因组** | `scripts/epigenomics/` | 01_methylation_analysis.R（DMR）· 02_atac_chip_analysis.R（DAR/DHR） | 区域甲基化/peak 计数→差异检验→火山+热图；上游 Bismark/MACS2 模板见文档 |
| **16 · 表型组** | `scripts/phenomics/` | 01_phenotype_blup.R · 02_spectral_indices.R | MAD 异常标记 + 环境调整 BLUP + 遗传力；7 种植被指数；PlantCV/GS 指引 |
| **17 · 单细胞** | `scripts/single_cell/` | 01_scanpy_pipeline.py | Scanpy 全流程模板（QC→聚类→marker→打分→PAGA）；RSTA/RICE scGPT 集成指引 |
| **18 · 空间组学** | `scripts/spatial/` | 01_spatial_pipeline.py | Stereo-seq/Visium→表达域→Moran's I SVG→邻域富集；RSTA 时空图谱对接 |

> v3.0 各模块均配**合成数据生成器**（`scripts/utils/make_demo_*.py`）与独立文档
> （`docs/variants_usage.md`、`pangenome_usage.md`、`qtl_usage.md`、`proteomics_usage.md`、
> `epigenomics_usage.md`、`phenomics_usage.md`、`single_cell_usage.md`、`spatial_usage.md`、
> `rice_databases.md`），可独立验证各模块。

### v3.1 进阶分析扩展（v3.1 新增）

| # | 脚本 | 能力 | 验证 |
|---|---|---|---|
| 1 | `scripts/gwas/07_haplotype.R` | 候选位点**单倍型分析**：EM 相位推断 + 主单倍型分组 + ANOVA/Kruskal 表型检验 + 箱线图；`--mode snp` 等位剂量分组 | demo 因果位点 p=2.45e-12（EM）/ 3.49e-28（剂量） |
| 2 | `scripts/qtl/02_colocalization.R` | **GWAS×QTL 共定位筛查**：双轨迹叠加图、lead 距离、|Z| 相关三级判定（colocalized/partial/distinct） | 构造数据 r=1 完美回收 |
| 3 | `scripts/phenomics/03_gs_gblup.R` | **GBLUP 基因组选择**：VanRaden 亲缘 + 岭回归解 + 自动 λ + 重复 k 折 CV 准确率 | demo r=0.31±0.11 |
| 4 | `scripts/phenomics/04_ammi_gge.R` | **AMMI/GGE biplot**：SVD 互作分解、which-won-where 多边形、ASV 稳定性表 | 3 环境 IPC 秩 2 数学自洽 |
| 5 | `scripts/popgen/07_xpehh.R` | **XP-EHH 选择扫描**：精确单倍型 EHH + iEHH 积分 + 窗口 z 标准化 + 候选窗口 | demo 亚群分化 |z|=33 |
| 6 | `scripts/evolution/01_ka_ks.py` | **KA/KS 基因家族进化**：BLOSUM62 NW 比对（或外部比对）+ pal2nal + Nei-Gojobori/JC | 模拟 dN/dS≈0.11 → 中位 0.195 |
| 7 | `scripts/tf/01_promoter_kmer.py` | **启动子 k-mer 富集**：链感知启动子提取 + Fisher/BH + top k-mer 条形图（top k-mer 可与 PlantTFDB motif 比对） | 注入 motif 排名 #1（FDR 3e-9） |
| 8 | `docs/small_rna_splicing_usage.md` | **miRNA 与可变剪接流程模板**（miRDeep2/rMATS/SUPPA2 命令 + 水稻资源 + 与模块衔接） | 文档级 |

### v3.2 进阶分析扩展（v3.2 新增）

| # | 脚本 | 能力 | 验证 |
|---|---|---|---|
| 1 | `scripts/qtl/03_coloc_bayesian.R` | **Bayesian 共定位（PP4）**：coloc.abf 自实现（近似 Bayes 因子 + H0-H4 后验，零依赖），SNP 级共享因果后验，对接 02 筛查结果 | demo 共享因果 PP4=0.994 lead=SNP017；不同因果 PP3>PP4 |
| 2 | `scripts/gwas/08_2d_association.R` | **2D-GWAS 成对性状关联**：性状 A 极端分组的 B 关联对比（PATOWAS 方法整理），条件特异性位点判定 + 双向散点 | demo 植入上位性 SNP030 以 +A 方向回收，假阳性≈期望值 |
| 3 | `scripts/integration/07_grn_boost.py` | **TF 调控网络推断**：GRNBoost2 式随机梯度提升（sklearn 实现，arboreto 免安装），按 target 排序边表 | demo 三模块重建、诱饵 TF 排除 |
| 4 | `scripts/phenomics/05_power_analysis.R` | **功效与样本量**：两群 t 检验功效/n、GWAS 位点检测功效（非中心 χ²）与 80% 功效所需 PVE、GS 训练群大小（Daetwyler） | 与 R 内置 power.t.test 对拍一致（0.4657） |
| 5 | `scripts/proteomics/04_uniprot_mapping.py` | **UniProt 蛋白→基因映射自动化**：下载 UP000059680 TSV→解析 MSU/RAP-DB 双风格 locus（MSU 优先、版本号剥离、去重）→ 输出 `protein_gene_map.csv`（模块 1 富集衔接）；缓存 + `--tsv-file`/`--insecure` 离线兜底 | demo 自检（LOC_ 归一化/MSU 优先/unmatched 检测）通过 |
| 6 | `scripts/proteomics/05_interpro_domains.py` | **InterProScan TSV 解析 + 差异蛋白结构域富集**：超几何过度表示 + BH FDR + 条形图；`--background` 限定检测蛋白全集 | demo 植入结构域 top 排名（fold 4.1，FDR 1.2e-30） |
| 7 | `scripts/evolution/02_mirna_targets.py` | **miRNA 靶基因预测**：psRNATarget 式种子/3' 配对打分本地实现（G:U 半配对、cleavage/inhibition 分级），零依赖 | demo 完美位点=cleavage(0 分)、种子错配=inhibition；真实 FASTA 冒烟通过 |
| 8 | `scripts/tf/02_motif_tf_match.py` | **motif→TF 匹配**：MEME PWM 库解析 + log-odds 双链滑窗匹配（衔接 01 k-mer 输出）；TF 排名；后续 TOMTOM/PlantPAN 定量 | demo 3 consensus 命中对应 PWM、2 噪声 k-mer 未匹配 |

### 模块 1：转录组单独分析（scripts/transcriptomics/）

| # | 脚本 | 步骤 | 工具 | 关键输出 |
|---|---|---|---|---|
| 01 | `01_qc_filter.py` | 数据质控与过滤：低表达基因过滤、样本 PCA、样本相关性、文库大小检查 | pandas/numpy/matplotlib | `counts_filtered.csv`、QC 报告与图（PDF/SVG/TIFF） |
| 02 | `02_deseq2.R` | 差异表达：两组或多组两两比较，log2FC / pvalue / padj；支持 `~ batch + group` | DESeq2（edgeR 可选） | `de_all_*.csv`、`de_sig_*.csv`、VST 归一化矩阵、离散度诊断图 |
| 03 | `03_enrichment.R` | 功能富集：GO BP/CC/MF、KEGG（Reactome 可选） | clusterProfiler + `org.Osativa.eg.db` | 富集结果表、点图/气泡图 |
| 04 | `04_gsea.R` | GSEA：排序基因列表 + 自定义 GMT（水稻无 MSigDB，提供 GMT 构建说明） | fgsea / clusterProfiler::GSEA | NES/p/padj 表、富集曲线 |
| 05 | `05_wgcna.R` | 共表达网络：软阈值筛选、模块识别、模块-性状关联、Hub 基因、模块富集 | WGCNA | 模块-性状热图、hub 基因列表、模块特征基因 |
| 06 | `06_tf_analysis.R` | 转录因子：PlantTFDB/iTAK 注释整合，差异表达 TF 筛选与家族统计 | R + 注释文件 | `de_TF.csv`、TF 家族统计图 |
| 07 | `07_plot_transcriptome.R` | 转录组出版图汇编：火山图、MA 图、热图、富集总览 | ggplot2 + plot_theme.R | 多面板组合图 |

### 模块 2：代谢组单独分析（scripts/metabolomics/）

| # | 脚本 | 步骤 | 工具 | 关键输出 |
|---|---|---|---|---|
| 01 | `01_preprocess.py` | 预处理：缺失值汇总与填充（半最小值/KNN）、总和/log/Pareto 归一化、异常样本检查 | pandas/scikit-learn | 预处理矩阵、QC 图（缺失值、归一化前后 PCA/箱线） |
| 02 | `02_oplsda.R` | 差异代谢物：PCA、OPLS-DA、VIP>1 且 p<0.05（可调）、置换检验 | ropls | 得分图、S-plot、VIP 表 |
| 03 | `03_pathway.R` | 代谢通路富集：KEGG osa（MetaboAnalystR / KEGGREST） | MetaboAnalystR / clusterProfiler | 通路富集表、气泡图 |
| 04 | `04_metabolite_classification.R` | 化学分类统计：黄酮/生物碱/脂质等差异代谢物分布 | R | 分类统计表与图 |
| 05 | `05_biomarker.R` | 生物标志物：ROC/AUC、随机森林重要性、组合筛选 | pROC / randomForest | ROC 曲线、重要性排序 |
| 06 | `06_metabolite_network.R` | 代谢网络：代谢物-代谢物相关性（\|r\|>0.8）、枢纽识别 | igraph | 网络边表/点表、网络图 |
| 07 | `07_plot_metabolome.R` | 代谢组出版图汇编：OPLS-DA 得分图、VIP、分类组成、ROC | ggplot2 + plot_theme.R | 多面板组合图 |

### 模块 3：多组学联合分析（scripts/integration/，仅转录组+代谢组）

| # | 脚本 | 步骤 | 工具 | 关键输出 |
|---|---|---|---|---|
| 01 | `01_mofa.R` | 数据整合：MOFA 潜在因子识别、方差分解 | MOFA2 | 因子得分/载荷、方差解释图 |
| 02 | `02_diablo.R` | 判别整合：DIABLO 多组学判别 + O2PLS（`--method o2pls`） | mixOmics | 联合得分图、载荷 circos、性能评估 |
| 03 | `03_correlation.R` | 相关性：差异基因 × 差异代谢物 Pearson/Spearman，\|r\|>0.8 且 p<0.05 | R | 相关性对表、网络边表、相关性热图 |
| 04 | `04_joint_pathway.R` | 联合通路：基因与代谢物共同富集的 KEGG 通路 | clusterProfiler + KEGGREST | 联合通路表、通路整合图 |
| 05 | `05_network_integration.R` | 多组学网络：WGCNA 模块 × 代谢物模块/枢纽整合 | WGCNA + igraph | 模块关联表、枢纽基因/代谢物、多层网络图 |
| 06 | `06_plot_integration.R` | 整合出版图：联合得分图、桑基图、Circos、UpSet | ggplot2/ggalluvial/circlize | 多面板组合图 |

## 4. Nature 级绘图规范（混合后端，全包统一）

- **后端策略**：R 分析脚本内使用 ggplot2/ComplexHeatmap；Python 质控/预处理脚本使用
  matplotlib。两端共用 `config/nature_palette.json` 调色板与
  `scripts/utils/plot_theme.R` / `plot_theme.py` 主题，**视觉完全一致**。
- **配色**：低饱和、色盲友好（NMI pastel 风格）；语义色固定：上调红 `#CE6A6B`、
  下调蓝 `#5B84AE`、不显著灰 `#C9C4BE`、高亮琥珀 `#E8A33D`；各面板配色统一。
- **字体**：Arial（或度量兼容的 Helvetica）；正文 8 pt、标题 9 pt 加粗、刻度/图例 7 pt；
  线宽 0.5 pt；PDF 文本可编辑可嵌入。
- **输出**：每张图同时输出 **PDF + SVG + TIFF（300 dpi）**；单栏 89 mm、1.5 栏 120 mm、
  双栏 180 mm。
- **绘图逻辑**：每张图对应明确科学结论；无意义装饰一律省略；图注科学含义写在生成脚本注释中。
- **主题实现**：自行实现，未复制任何第三方主题源码（含 `natureskills/nature-figure`）。

## 5. 内置 Nature 图表与写作规范（第三方集成：nature-skills）

本仓库在 `third_party/nature-skills/` 内置了 [nature-skills](https://github.com/xiaolai/nature-skills)
（MIT License, Copyright (c) 2026 Yuan Yizhe；原始声明见 `third_party/nature-skills/LICENSE`，
溯源见 `NOTICE.md`）的**九个 skill** 源码，按需加载：

| 目录 | 能力 | 加载时机 |
|---|---|---|
| `skills/nature-figure/` | Nature 级图表规范与 matplotlib 模板（配色/字体/导出策略/多面板信息层级/图件契约） | "Nature 级作图/出版图规范/多面板组合" |
| `skills/nature-polishing/` | 英文论文润色（25 条规则：句长、时态、hedging、引用规范、英式拼写） | "润色/翻译为 Nature 风格英文" |
| `skills/nature-writing/` | 论文各章节写作（摘要/引言/方法/结果/讨论/结论/审稿自查） | "写论文/写摘要/写引言" |
| `skills/nature-academic-search/` | 学术检索（PubMed/MEDLINE 查询、DOI 提取、RIS/BibTeX 导出） | "查文献/找参考文献/导出引用" |
| `skills/nature-citation/` | 引用格式规范与转换 | "引用格式/BibTeX" |
| `skills/nature-data/` | 数据呈现规范（图表中的数据组织） | "数据如何呈现/数据规范" |
| `skills/nature-paper2ppt/` | 论文一键转 Nature 风格中文汇报 PPTX | "论文转 PPT/组会汇报" |
| `skills/nature-reader/` | 全文中英对照精读辅助（图表感知） | "精读论文/对照阅读" |
| `skills/nature-response/` | 审稿意见逐条回复信工作流 | "审稿回复/返修/response letter" |

说明：
- 仅复制各 skill 的源码与参考资料；上游 `assets/`（约 30 MB 演示资产，含 figures4papers
  第三方素材）**未包含**，需要时请从上游获取。
- 与本包自研绘图主题（`config/nature_palette.json` + `scripts/utils/plot_theme.*`）相互独立：
  本包分析脚本默认使用自研主题（R + Python 双后端，输出 PDF/SVG/TIFF）；
  nature-figure 作为"规范参考与模板库"按需加载，两者不冲突。
- 二次分发本项目时请保留 `third_party/nature-skills/LICENSE` 与 `NOTICE.md`。

此外，`third_party/` 还内置 5 个第三方技能库（均为精选子集，含 LICENSE 与 NOTICE 溯源；
数量与内容详见 `docs/skill_routing.md` 库速查）：

| 目录 | 内容 | 加载时机 |
|---|---|---|
| `third_party/bioskills/` | 354 个技能（GPTomics/bioSkills，MIT，41 类）：首轮 10 类——差异表达、代谢组、多组学整合、通路分析、数据可视化、表达矩阵、变异检测、群体遗传、**引物设计、CRISPR 筛选**；二轮 31 类——**单细胞(17)、空间转录组(12)、ChIP/ATAC/甲基化/Hi-C(43)、蛋白组(9)、剪接/small RNA/表转录组(20)、GRN(6)、因果基因组学(11)、比较基因组(13)、长读长(9)、基因组组装/注释/区间/编辑(29)、CNV/相位(15)、序列/比对/QC/定量(65)、数据库(15)、ML/报告/实验设计(17)** | 需要具体环节的专家经验时：批次校正、时间序列 DE、XCMS/MS-DIAL 上游、脂质组、整合设计与数据协调、Circos/UpSet/多面板绘图、引物/敲除设计、单细胞 QC-聚类-轨迹-通信、ChIP/ATAC 峰 calling、WGCNA→GRN、孟德尔随机化/共定位、基因组注释与 CRISPR 编辑等 |
| `third_party/kdense-scientific-skills/` | 28 个技能（K-Dense-AI，MIT，活跃维护）：bulk-rnaseq、pydeseq2、pathway-enrichment、scientific-visualization、statistical-analysis + 单细胞/空间（anndata、scanpy、scvi-tools、scvelo）、GRN 推断（arboreto、networkx）、数据库 API（bioservices、gget）、质谱（matchms、pyopenms）、表观组（deeptools）、序列与进化（biopython、etetoolkit、phylogenetics）、变异处理（pysam）、统计与 ML（statistical-power、statsmodels、scikit-learn、umap-learn）、实验支撑（mageck、primer-design、experimental-design、cellprofiler） | 需要上游流程规范（FASTQ→counts、nf-core）、统计方法选择、出版级可视化与"可辩护分析"规范时；单细胞分析（模块 16/17）、GRN 推断、UniProt 映射、质谱鉴定、CRISPR 筛选/引物设计时。注：GPL 库（bioservices/etetoolkit 等）仅收 MIT 文档，库本体运行时从 PyPI/bioconda 安装 |
| `third_party/omicverse-skills/` | 61 个技能（omicverse 官方，MIT 文档；⚠️ 驱动的主包 omicverse 为 GPL-3.0，仅收文档，主包运行时 pip/conda 安装）：单细胞全流程（预处理/注释/批次整合/聚类/轨迹/velocity/SCENIC/CytoTRACE/CNMF）、**单细胞基础模型（scGPT/Geneformer/scFoundation/UCE/CellPLM）**、参考映射与跨模态注释迁移、细胞通信（CellPhoneDB/LIANA/Scissor）、单细胞→空间映射、bulk WGCNA、**代谢组多因子通路**、脂质组、PPI（StringDB）、bulk 反卷积、微生物组（16S DADA2/比较）、数据 IO/统计/绘图/导出 | 模块 16 单细胞/空间分析的 Python 路线参考（尤其 scGPT 类基础模型与 RSTA 衔接）；代谢组多因子通路分析；根际微生物组候选扩展；注意区分"方法思路"与 omicverse 具体 API |
| `third_party/cell-skills/` | 12 个技能（yrui-cmd/Cell + Cell_CNS_Figure，MIT）：科研工作流全套——选题/研究计划/综述、英文论文润色、审稿回复、投稿整理、科研可视化代码、数据图、CNS 风格绘图、PPT 矢量转换 | 论文全流程写作与投稿辅助；注意：CNS 绘图走外部付费服务（额度制），PPT 相关技能需 Windows + Office |
| `third_party/sciagent-skills/` | 64 个技能（jaechang-hits/SciAgent-Skills，CC BY 4.0）：多组学整合（MOFA+/通路/PPI）、数据可视化、生物统计、RNA-seq 全链（DESeq2/PyDESeq2/Salmon/gseapy）、变异检测（bcftools/GATK/PLINK2/SnpEff）、**分子生物学（sgRNA 设计指南/质粒注释/RNA 二级结构）**、QC、科研写作与七刊图规范（Nature/Cell/Science/Lancet/NEJM/PNAS/eLife） | 需要跨领域方法参考、分刊图规范或写作辅助时 |

方法参考笔记：`docs/references/patowas-2d-association.md`——群体水平 2D 关联作图（GWAS×组学）
的方法导航，说明何时从"组内比较"切换到"群体关联"路线（上游 PATOWAS 为 GPL-3.0 独立程序，未打捆）。

配套实验协议（供实验设计参考）：
- `docs/protocols_molecular.md`——引物设计与评估 9 条原则、同源重组构建、RNAi（pYLRNAi.5 类）、
  蛋白表达纯化（GST/MBP/His）与 SDS-PAGE 配方、烟草亚细胞定位、LUC 互补（split-LUC + LUC/Ren）、
  mhiTAIL-PCR、GUS 染色、BioEdit 序列分析；
- `docs/research_workflow_tf.md`——转录因子功能研究全流程路线图（背景调查/功能/表达/机理 + 周期估算）。

配套数据与工具：仓库**内置官方注释快照** `assets/annotations/`（RAP-DB MSU↔RAP-DB ID 转换表
34,133 对、PlantTFDB 水稻 TF 列表 1,862 个 TF / 56 家族、53 个 TF 家族 GSEA 基因集），
开箱即用、无需联网；需要更新时运行 `scripts/utils/fetch_rice_annotations.py --force`
（自动从官方源刷新；网络受限时支持 `--insecure` 或"浏览器下载 + 传本地路径"）。

### 5.1 分析模块 → 第三方技能路由（528 个技能的入口）

执行某模块分析前，AI 助手应**先加载对应技能**获取专家经验（版本兼容性、常见陷阱、
参数选择、可辩护分析规范）。完整路由表见 **`docs/skill_routing.md`**（104 条映射，
全部路径已核验）。常用路由速查：

| 模块（自研脚本） | 首选技能路径（third_party/ 下） |
|---|---|
| 转录组上游（FASTQ→counts） | `kdense-scientific-skills/skills/bulk-rnaseq/` + `bioskills/read-qc/`、`read-alignment/`、`rna-quantification/` |
| 差异表达 | `bioskills/differential-expression/` + `kdense-scientific-skills/skills/pydeseq2/` |
| 代谢组（预处理/多元/通路） | `bioskills/metabolomics/` + `omicverse-skills/skills/bulk-metabol-*/` + `kdense-scientific-skills/skills/matchms/` |
| 多组学联合（MOFA/DIABLO） | `bioskills/multi-omics-integration/` + `sciagent-skills/skills/systems-biology-multiomics/mofaplus-multi-omics/` |
| GWAS 与共定位/精细定位 | `bioskills/population-genetics/` + `bioskills/causal-genomics/`（colocalization/fine-mapping/MR，v3.2 路径）+ `kdense-scientific-skills/skills/statistical-power/` |
| 引物/CRISPR 设计与筛选 | `bioskills/genome-engineering/`、`crispr-screens/`、`primer-design/` + `kdense-scientific-skills/skills/mageck/`、`primer-design/` |
| 泛基因组/组装/注释 | `bioskills/comparative-genomics/`、`genome-assembly/`、`genome-annotation/`、`long-read-sequencing/` |
| 蛋白质组 + UniProt 映射 | `bioskills/proteomics/` + `kdense-scientific-skills/skills/pyopenms/`、`bioservices/` |
| 表观组（ChIP/ATAC/甲基化/m6A/Hi-C） | `bioskills/chip-seq/`、`atac-seq/`、`methylation-analysis/`、`epitranscriptomics/`、`hi-c-analysis/` + `kdense-scientific-skills/skills/deeptools/` |
| 单细胞（含 scGPT 类基础模型） | `bioskills/single-cell/` + `kdense-scientific-skills/skills/anndata/`、`scanpy/`、`scvi-tools/` + `omicverse-skills/skills/single-cell-foundation-model/` |
| 细胞通信 | `sciagent-skills/skills/systems-biology-multiomics/cellchat-cell-communication/` + `omicverse-skills/skills/single-cell-liana-communication/` |
| 空间组学 | `bioskills/spatial-transcriptomics/` + `omicverse-skills/skills/single-to-spatial-mapping/` |
| GRN/共表达进阶（v3.2） | `bioskills/gene-regulatory-networks/` + `kdense-scientific-skills/skills/arboreto/`（GRNBoost2/GENIE3）+ `omicverse-skills/skills/bulk-wgcna-analysis/` |
| 进化/KA-KS/TF | `bioskills/comparative-genomics/` + `kdense-scientific-skills/skills/biopython/`、`etetoolkit/`、`phylogenetics/` |
| 出版图规范与论文 | `bioskills/data-visualization/` + `nature-skills/skills/` + `sciagent-skills/skills/scientific-writing/`（七刊图规范） |

## 6. 参数配置

所有脚本通过命令行参数配置（`--help` 查看），常用可调项：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--padj-cutoff` | 0.05 | 差异显著性阈值 |
| `--log2fc-cutoff` | 1 | \|log2FC\| 阈值 |
| `--vip-cutoff` | 1 | OPLS-DA VIP 阈值 |
| `--cor-cutoff` | 0.8 | 基因-代谢物相关性 \|r\| 阈值 |
| `--min-counts` / `--min-samples` | 10 / 3 | 低表达基因过滤（counts ≥ 10 且至少在 3 个样本表达） |
| `--seed` | 42 | 随机种子（全流程固定，保证可重复） |

完整参数表见 `docs/usage.md`。

## 7. 依赖

**R ≥ 4.3**：DESeq2、edgeR、clusterProfiler、org.Osativa.eg.db、ReactomePA（可选）、
fgsea、WGCNA、MOFA2、mixOmics、ropls、MetaboAnalystR、imputeLCMD、ComplexHeatmap、
ggplot2、ggrepel、igraph、ggraph、pROC、randomForest、UpSetR、ggalluvial、circlize、
jsonlite、dplyr、readr、tidyr；ragg（可选，高质量 TIFF）。

**Python ≥ 3.10**：numpy、pandas、scipy、scikit-learn、matplotlib、seaborn、pillow；
gseapy、networkx（可选）。

**基因组学模块（v2.0）**：R 侧核心依赖 `ape`（系统发育）、`data.table`、`matrixStats`、
`patchwork`（组合图）；可选增强：`rMVP`（推荐，Apache-2.0，`install.packages("rMVP")`）、
`GAPIT`；可选二进制：PLINK（VCF 转换，conda 或官方静态包）、ADMIXTURE（Linux/WSL）。
核心计算（PCA/LD/Fst/π/NJ/GWAS-simple）全部为自实现，缺任一二进制亦可运行。

**设计模块（v2.1）**：`primer3-py`（`pip install primer3-py`，GPL-2.0，仅作外部引擎调用）、
matplotlib（出图，缺失时优雅跳过）。sgRNA 设计引擎为纯自研，零外部依赖。

**v3.0 模块（全部可选增强，缺失时脚本给出指引并优雅退出）**：
- 变异检测：`bwa`/`samtools`/`gatk4`/`bcftools`/`snpeff`（conda bioconda，模板式调用不捆绑）；
- 泛基因组：`minigraph`/`cactus`/`pggb`/`minimap2`/`syri`（重计算，服务器运行）；PAV/Pan-GWAS 为纯自研零依赖；
- QTL 映射（`scripts/qtl/`）：纯自研零依赖；
- Meta-GWAS：纯自研零依赖（ggplot2 出图）；
- 蛋白质组：核心流程零额外依赖；`limma` 装上后自动启用（推荐）；
- 表观组：核心流程零额外依赖（t/wilcoxon + BH）；
- 表型组：纯自研（BLUP/遗传力/光谱指数）；GS 建议 `rrBLUP`（外部安装）；
- 单细胞/空间：独立环境 `scanpy` + `squidpy`（文档指引；未装时脚本给出安装命令）。

一键安装：`conda env create -f environment.yml`（详见 `docs/installation.md`）。

## 8. 快速开始

```bash
# 1) 安装环境
conda env create -f environment.yml && conda activate rice-skills

# 2) 一条命令跑通完整流程（自定义数据）
bash run_pipeline.sh --counts counts.csv --metadata metadata.csv \
    --metabolites metabolites.csv --outdir results --stages all

# 3) 或分模块运行（示例）
Rscript scripts/transcriptomics/02_deseq2.R \
    --counts results/transcriptome/01_qc/counts_filtered.csv \
    --metadata metadata.csv --contrast salt_vs_control \
    --outdir results/transcriptome/02_de

# 4) 测试流程自检（含工具层冒烟测试）
bash tests/test_pipeline.sh
```

基因组学模块快速体验（零外部依赖，合成数据全流程）：

```bash
python scripts/utils/make_demo_genomics.py --outdir demo_data/genomics
Rscript scripts/gwas/03_association.R --genotype demo_data/genomics/genotype.csv \
    --snp-info demo_data/genomics/snp_info.csv \
    --phenotype demo_data/genomics/phenotype.csv --trait salt_tolerance_index
# 完整用法（GWAS 全链 + 群体遗传 6 步）见 docs/genomics_usage.md
```

设计模块快速体验（v2.1，引物 + CRISPR 靶点设计）：

```bash
python scripts/design/01_sgrna_design.py --fasta demo_data/genomics/demo_gene.fa \
    --gff demo_data/genomics/demo_gene.gff3 --outdir results/design/01_sgrna
python scripts/design/02_primer_design.py --mode qpcr \
    --fasta demo_data/genomics/demo_gene.fa --gff demo_data/genomics/demo_gene.gff3 \
    --outdir results/design/02_primers
# 完整用法见 docs/design_usage.md；实验流程参考见 docs/lab_protocols.md
```

## 9. 输出目录结构

```
results/
├── transcriptome/01_qc/            # 过滤后矩阵 + QC 图 + run_info.txt
├── transcriptome/02_de/            # 差异表达表 + 火山图/MA 图
├── transcriptome/03_enrichment/    # GO/KEGG 富集表 + 点图
├── transcriptome/04_gsea/          # GSEA 结果 + 富集曲线
├── transcriptome/05_wgcna/         # 模块-性状热图 + hub 基因
├── transcriptome/06_tf/            # 差异表达 TF
├── transcriptome/07_figures/       # 出版级组合图
├── metabolome/01_preprocess/ ... 07_figures/
├── integration/01_mofa/ ... 06_plot_integration/
└── pipeline_log.txt                # 全流程运行日志
```

每个输出目录均含 `run_info.txt`（运行时间、参数、关键包版本、随机种子），保证可追溯。

基因组学模块输出（详见 docs/genomics_usage.md）：
`results/gwas/{01_qc,03_association,04_candidates,05_figures}/`、
`results/popgen/{01_pca,02_ld,03_selection,04_admixture,05_admix,06_phylogeny}/`。
设计模块输出：`results/design/{01_sgrna,02_primers_*}/`（详见 docs/design_usage.md）。

v3.0 模块输出：
`results/variants/{01_align,02_call,03_filter,04_annot}/`、
`results/pangenome/{01_graph,02_sv,03_pav,04_pan_gwas}/`、
`results/qtl/01_*/`、`results/gwas/06_meta/`、
`results/proteomics/{01_preprocess,02_diff,03_plots}/`、
`results/epigenome/{01_dmr,02_dar}/`、
`results/phenomics/{01_blup,02_indices}/`、
`results/single_cell/01_*/`、`results/spatial/01_*/`（各配独立文档）。

## 10. 注意事项（notes）

1. **数据格式**：见第 2 节；样本 ID 必须与元数据完全一致，否则脚本会给出明确报错与修复建议。
2. **内存需求**：WGCNA 与 MOFA 是大头。经验值：2 万基因 × 100 样本建议 ≥ 16 GB RAM；
   WGCNA 已用 `blockwiseModules` 分块计算；基因-代谢物相关性按块计算。
3. **可重复性**：所有脚本固定 `seed=42`，并在 `run_info.txt` 记录；DESeq2/WGCNA 等
   结果受版本影响，升级包后建议重跑并对比。
4. **物种迁移**：换物种只需修改 ①注释库（如 `org.Osativa.eg.db` → 目标物种）②KEGG 代码
   （`osa`）③TF 注释来源三处，流程不变。
5. **网络依赖**：clusterProfiler/KEGGREST 的富集分析需要联网访问 KEGG/GO 数据库；
   离线环境请预先下载基因集 GMT 并在 `04_gsea.R` 中本地运行。
6. **水稻无 MSigDB**：GSEA 请使用自定义 GMT（脚本内附从 KEGG osa 构建 GMT 的说明）。
7. **示例数据**：本仓库默认不附示例数据；可用
   `python scripts/utils/make_demo_data.py --outdir demo_data` 生成小型合成数据自检流程
   （合成数据仅供测试，勿用于生物学结论）。
