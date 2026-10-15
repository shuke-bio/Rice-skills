# 第三方技能路由表（分析模块 → 技能映射）

> **目的**：`third_party/` 内置 6 个库共 **528 个技能**。本表把每个分析模块（自研脚本）
> 映射到对应的第三方技能路径，让 AI 助手在执行该模块分析前**先加载对应技能**获取
> 专家经验（版本兼容性、常见陷阱、参数选择、可辩护分析规范）。
>
> **用法**：执行某模块分析 → 在下表找到该模块行 → 加载列出的 `SKILL.md`（及其 `references/`）
> → 按其规范指导分析。技能是方法论参考，不替代自研脚本；两者冲突时以自研脚本 +
> `docs/usage.md` 为准。技能清单与许可溯源见各库 `NOTICE.md`。
>
> 路径均相对 `third_party/`；所有条目已核对实际存在（2026-10-15）。

---

## 库速查（6 个库一句话定位）

| 库 | 定位 | 规模 |
|---|---|---|
| `bioskills/` | 生信全领域专家经验（两层：`类别/技能/`） | 41 类 354 个 |
| `kdense-scientific-skills/` | "可辩护分析"规范 + 单细胞/GRN/质谱/统计功效 | 28 个 |
| `omicverse-skills/` | OmicVerse 工作流指引（单细胞基础模型/代谢多因子/微生物组；主包 GPL-3.0 运行时装） | 61 个 |
| `sciagent-skills/` | 多组学/CellChat/分子生物学/分刊图规范/科研写作 | 64 个 |
| `nature-skills/` | Nature 级图表/写作/润色/检索/审稿回复 | 9 个 |
| `cell-skills/` | Cell 风格科研全流程（选题/绘图/PPT/投稿） | 12 个 |

---

## 模块 1 转录组（`scripts/transcriptomics/`）

| 场景 | 技能路径 |
|---|---|
| FASTQ→counts 上游全流程规范（nf-core/STAR/Salmon/featureCounts、链特异性） | `kdense-scientific-skills/skills/bulk-rnaseq/` |
| FASTQ 质控/接头/UMI | `bioskills/read-qc/`（fastp-workflow、adapter-trimming、quality-filtering、quality-reports、rnaseq-qc、umi-processing、contamination-screening） |
| 比对 | `bioskills/read-alignment/`（hisat2-alignment、star-alignment、bwa-alignment、bowtie2-alignment） |
| SAM/BAM 操作 | `bioskills/alignment-files/`（sam-bam-basics、alignment-sorting、alignment-indexing、bam-statistics、duplicate-handling 等 11 个） |
| 定量 | `bioskills/rna-quantification/`（featurecounts-counting、alignment-free-quant、tximport-workflow、count-matrix-qc） |
| 差异表达 | `bioskills/differential-expression/`（deseq2-basics、edger-basics、batch-correction、timeseries-de、de-results、de-visualization）+ `kdense-scientific-skills/skills/pydeseq2/` |
| 序列读取/操作基础 | `bioskills/sequence-io/`（9）、`bioskills/sequence-manipulation/`（7）+ `kdense-scientific-skills/skills/biopython/` |

## 模块 2 代谢组（`scripts/metabolomics/`）

| 场景 | 技能路径 |
|---|---|
| 上游预处理（XCMS/MS-DIAL） | `bioskills/metabolomics/`（xcms-preprocessing、msdial-preprocessing、normalization-qc） |
| 质谱谱图匹配/注释 | `kdense-scientific-skills/skills/matchms/` + `bioskills/metabolomics/metabolite-annotation/` |
| 质谱数据处理 | `kdense-scientific-skills/skills/pyopenms/` |
| 多因子/通路 | `omicverse-skills/skills/bulk-metabol-pathway-multifactor/` + `bioskills/metabolomics/pathway-mapping/` |
| 多元统计（OPLS-DA 类） | `omicverse-skills/skills/bulk-metabol-multivariate/` + `bioskills/metabolomics/statistical-analysis/` |
| 靶向/脂质组/同位素 | `bioskills/metabolomics/`（targeted-analysis、lipidomics、isotope-tracing）+ `omicverse-skills/skills/bulk-metabol-untargeted-lipidomics/` |

## 模块 3 多组学联合（`scripts/integration/`）

| 场景 | 技能路径 |
|---|---|
| MOFA/DIABLO 整合 | `bioskills/multi-omics-integration/`（mofa-integration、mixomics-analysis、data-harmonization、integration-design、similarity-network）+ `sciagent-skills/skills/systems-biology-multiomics/mofaplus-multi-omics/` |
| 多组学方法选择 | `sciagent-skills/skills/systems-biology-multiomics/omics-analysis-guide/` + `omicverse-skills/skills/single-multiomics/` |
| 关联网络 | `kdense-scientific-skills/skills/networkx/` + `sciagent-skills/skills/systems-biology-multiomics/string-database-ppi/` |

## 模块 7 GWAS（`scripts/gwas/`）

| 场景 | 技能路径 |
|---|---|
| 关联检验/PLINK/群体结构 | `bioskills/population-genetics/`（association-testing、plink-basics、population-structure、rare-variant-association） |
| 线性/混合模型 | `kdense-scientific-skills/skills/statsmodels/` |
| GWAS 图（Manhattan/QQ/LocusZoom） | `bioskills/data-visualization/manhattan-qq-locuszoom/` |
| 基因型填充/单倍型相位 | `bioskills/phasing-imputation/`（genotype-imputation、reference-panels、imputation-qc、haplotype-phasing） |
| 样本量与功效 | `kdense-scientific-skills/skills/statistical-power/` + `bioskills/experimental-design/`（power-analysis、sample-size） |
| **共定位/精细定位/MR（v3.2 升级路径）** | `bioskills/causal-genomics/`（colocalization-analysis、fine-mapping、mendelian-randomization、effector-gene-prioritization、transcriptome-wide-association、genetic-correlation、heritability-partitioning、mediation-analysis 等 11 个） |
| 单倍型分析 | `bioskills/population-genetics/linkage-disequilibrium/` + `bioskills/phasing-imputation/haplotype-phasing/` |

## 模块 8 群体遗传（`scripts/popgen/`）

| 场景 | 技能路径 |
|---|---|
| PCA/Admixture/结构 | `bioskills/population-genetics/population-structure/` + `bioskills/data-visualization/dimensionality-reduction-plots/` |
| LD 衰减 | `bioskills/population-genetics/linkage-disequilibrium/` |
| 选择分析（含 XP-EHH） | `bioskills/population-genetics/selection-statistics/` |
| allele 层面计算 | `bioskills/population-genetics/scikit-allel-analysis/` |
| 进化树可视化 | `kdense-scientific-skills/skills/etetoolkit/`、`skills/phylogenetics/` |

## 模块 9 引物与 CRISPR（`scripts/design/`）

| 场景 | 技能路径 |
|---|---|
| 引物设计/qPCR/特异性 | `bioskills/primer-design/`（primer-basics、qpcr-primers、primer-specificity、primer-validation）+ `kdense-scientific-skills/skills/primer-design/`（含热力学/特异性脚本与报告模板） |
| sgRNA 设计 | `bioskills/genome-engineering/`（grna-design、base-editing-design、off-target-prediction、hdr-template-design、prime-editing-design）+ `sciagent-skills/skills/molecular-biology/sgrna-design-guide/` |
| CRISPR 筛选分析 | `bioskills/crispr-screens/`（15 个：mageck-analysis、library-design、crispresso-editing、hit-calling、screen-qc、perturb-seq-analysis 等）+ `kdense-scientific-skills/skills/mageck/` |
| 质粒注释 | `sciagent-skills/skills/molecular-biology/plannotate-plasmid-annotation/` |
| 实验设计与功效 | `kdense-scientific-skills/skills/experimental-design/` + `bioskills/experimental-design/`（5 个） |

## 模块 10 泛基因组（`scripts/pangenome/`）

| 场景 | 技能路径 |
|---|---|
| 泛基因组/PAV/正交基因 | `bioskills/comparative-genomics/`（pangenome-analysis、ortholog-inference、synteny-analysis、whole-genome-alignment、gene-family-evolution、positive-selection 等 13 个） |
| 组装 | `bioskills/genome-assembly/`（hifi-assembly、long-read-assembly、short-read-assembly、scaffolding、assembly-qc、assembly-polishing、contamination-detection） |
| 长读长 | `bioskills/long-read-sequencing/`（basecalling、clair3-variants、structural-variants、isoseq-analysis、medaka-polishing 等 9 个） |
| 结构变异 | `bioskills/variant-calling/structural-variant-calling/` |
| 基因组注释 | `bioskills/genome-annotation/`（eukaryotic-gene-prediction、functional-annotation、repeat-annotation、annotation-transfer、ncrna-annotation 等 7 个） |

## 模块 11 QTL（`scripts/qtl/`）

| 场景 | 技能路径 |
|---|---|
| QTL/共定位 | `bioskills/causal-genomics/colocalization-analysis/`、`fine-mapping/` |
| 区间运算/候选基因筛选 | `bioskills/genome-intervals/`（bed-file-basics、gtf-gff-handling、interval-arithmetic、overlap-significance、proximity-operations 等 8 个） |
| 基因注释查询 | `bioskills/database-access/`（biomart-queries、ensembl-rest、entrez-search、entrez-link 等 15 个） |

## 模块 13 蛋白质组（`scripts/proteomics/`）

| 场景 | 技能路径 |
|---|---|
| 蛋白组全流程（DIA/肽段/定量/QC/PTM） | `bioskills/proteomics/`（9 个：data-import、dia-analysis、peptide-identification、protein-inference、quantification、differential-abundance、ptm-analysis、proteomics-qc、spectral-libraries） |
| 质谱底层处理 | `kdense-scientific-skills/skills/pyopenms/` |
| UniProt 蛋白映射（v3.2 自动化路径） | `kdense-scientific-skills/skills/bioservices/`（含 `scripts/batch_id_converter.py`）+ `bioskills/database-access/uniprot-access/` |
| 蛋白结构延伸（未收 structural-biology，见路线图） | 暂用 `bioskills/database-access/` 查询类 |

## 模块 14 表观组（`scripts/epigenomics/`）

| 场景 | 技能路径 |
|---|---|
| ChIP-seq/CUT&RUN | `bioskills/chip-seq/`（12 个：peak-calling、differential-binding、motif-analysis、chromatin-state-segmentation、cut-and-run-tag、super-enhancers、spike-in-normalization 等） |
| ATAC-seq | `bioskills/atac-seq/`（12 个：atac-peak-calling、differential-accessibility、footprinting、enhancer-gene-linking、nucleosome-positioning、single-cell-atac 等） |
| 甲基化/WGBS/DMR | `bioskills/methylation-analysis/`（10 个：bismark-alignment、methylation-calling、dmr-detection、methylkit-analysis、differential-cpg-testing 等） |
| 信号矩阵/富集 | `kdense-scientific-skills/skills/deeptools/` |
| m6A/RNA 修饰 | `bioskills/epitranscriptomics/`（5 个）+ `bioskills/long-read-sequencing/nanopore-methylation/` |
| 3D 基因组（Hi-C） | `bioskills/hi-c-analysis/`（9 个） |

## 模块 15 表型组与 GS（`scripts/phenomics/`）

| 场景 | 技能路径 |
|---|---|
| BLUP/遗传力 | `bioskills/experimental-design/randomization-blocking/` + `kdense-scientific-skills/skills/statsmodels/` |
| 基因组预测（GBLUP）与 ML | `kdense-scientific-skills/skills/scikit-learn/` + `bioskills/machine-learning/`（model-validation、omics-classifiers、biomarker-discovery、prediction-explanation） |
| AMMI/GGE 与功效 | `kdense-scientific-skills/skills/statistical-power/` |
| 图像表型 | `kdense-scientific-skills/skills/cellprofiler/`（PlantCV 工作流模板整理待做，见 roadmap） |
| 降维可视化 | `kdense-scientific-skills/skills/umap-learn/` |

## 模块 16 单细胞（`scripts/single_cell/`）

| 场景 | 技能路径 |
|---|---|
| 预处理/聚类/注释/整合 | `bioskills/single-cell/`（17 个：preprocessing、clustering、cell-annotation、markers-annotation、batch-integration、doublet-detection、data-io、differential-abundance 等）+ `kdense-scientific-skills/skills/anndata/`、`scanpy/` |
| 数据格式（h5ad） | `kdense-scientific-skills/skills/anndata/` |
| 批次整合（深度模型） | `kdense-scientific-skills/skills/scvi-tools/` |
| 轨迹/velocity/命运 | `bioskills/single-cell/trajectory-inference/` + `kdense-scientific-skills/skills/scvelo/` + `omicverse-skills/skills/`（single-cell-trajectory-inference、monocle2、sctour、via、cellrank-fate、rna-velocity、cytotrace2） |
| **基础模型（scGPT/Geneformer，RSTA 衔接）** | `omicverse-skills/skills/single-cell-foundation-model/` + `reference-label-transfer/`、`cross-modal-celltype-transfer/`、`popv-annotation/`、`cellvote-consensus/` |
| 细胞通信 | `sciagent-skills/skills/systems-biology-multiomics/cellchat-cell-communication/` + `bioskills/single-cell/cell-communication/` + `omicverse-skills/skills/`（single-cell-cellphonedb-communication、single-cell-liana-communication） |
| 调控子/SCENIC | `omicverse-skills/skills/single-cell-scenic/` + `bioskills/gene-regulatory-networks/scenic-regulons/` |
| 表型-细胞型关联（GWAS 衔接） | `omicverse-skills/skills/scissor-cell-selection/` |
| 剪接/代谢物通信/谱系 | `bioskills/single-cell/`（perturb-seq、lineage-tracing、metabolite-communication、cnv-inference、multimodal-integration） |

## 模块 17 空间组学（`scripts/spatial/`）

| 场景 | 技能路径 |
|---|---|
| 全流程（IO/预处理/domains/通信/反卷积） | `bioskills/spatial-transcriptomics/`（12 个：spatial-data-io、spatial-preprocessing、spatial-domains、spatial-neighbors、spatial-communication、spatial-deconvolution、spatial-statistics、spatial-visualization、image-analysis、high-resolution-binning 等） |
| 单细胞→空间映射 | `omicverse-skills/skills/single-to-spatial-mapping/` + `spatial-tutorials/` |

## 模块 18 进化/KA-KS/TF（`scripts/evolution/`、`scripts/tf/`）

| 场景 | 技能路径 |
|---|---|
| KA/KS 与正选择 | `bioskills/comparative-genomics/`（positive-selection、gene-family-evolution、ancestral-reconstruction、ortholog-inference） |
| 系统发育建树 | `kdense-scientific-skills/skills/phylogenetics/`（MAFFT+IQ-TREE+FastTree）+ `bioskills/alignment/`（multiple-alignment、alignment-trimming） |
| 树操作/可视化 | `kdense-scientific-skills/skills/etetoolkit/` + `bioskills/sequence-manipulation/codon-usage/` |
| 启动子 k-mer/motif | `bioskills/sequence-manipulation/motif-search/` + `bioskills/chip-seq/motif-analysis/` |
| GRN/调控网络 | `bioskills/gene-regulatory-networks/`（6 个：grn-inference、coexpression-networks、differential-networks、multiomics-grn、perturbation-simulation、scenic-regulons）+ `kdense-scientific-skills/skills/arboreto/`（GRNBoost2/GENIE3，v3.2 路径）+ `omicverse-skills/skills/bulk-wgcna-analysis/` |
| 序列基础 | `kdense-scientific-skills/skills/biopython/` |

## 跨模块层

| 场景 | 技能路径 |
|---|---|
| 出版级绘图规范（配色/多面板/统计标注） | `bioskills/data-visualization/`（19 个：color-palettes、multipanel-figures、statistical-annotation、heatmaps-clustering、volcano-and-ma-plots、circos-plots、upset-plots、sequence-logos、genome-tracks 等）+ `kdense-scientific-skills/skills/scientific-visualization/` |
| 分刊图规范（Nature/Cell/Science/Lancet/NEJM/PNAS/eLife） | `sciagent-skills/skills/scientific-writing/`（七刊 figure-guide）+ `nature-skills/skills/nature-figure/` + `cell-skills/skills/cell-cns-figure/` |
| 统计方法选择/报告 | `kdense-scientific-skills/skills/statistical-analysis/` + `sciagent-skills/skills/biostatistics/`（4 个） |
| 论文写作/润色/审稿回复/检索 | `nature-skills/skills/`（9 个）+ `cell-skills/skills/`（12 个）+ `sciagent-skills/skills/scientific-writing/`（manuscript-writing、peer-review-methodology、literature-review、pubmed/openalex/biorxiv 数据库） |
| 分析报告/表格输出 | `bioskills/reporting/`（rmarkdown-reports、quarto-reports、publication-tables、figure-export、jupyter-reports） |
| 数据库 API/批量下载 | `bioskills/database-access/`（15 个）+ `kdense-scientific-skills/skills/bioservices/`、`gget/` |
| VCF/BAM 程序化处理 | `kdense-scientific-skills/skills/pysam/` + `bioskills/variant-calling/`（13 个） |

---

## 维护说明

- 新增收编技能后，在本表对应模块行追加路径；删库时同步删行。
- 表中路径必须真实存在（用 `ls third_party/...` 核对），禁止凭记忆写路径。
- 上游库更新时：bioskills 已归档不再变；K-Dense/omicverse 活跃，升级后需重新核对。
