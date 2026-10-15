# NOTICE — 第三方组件溯源说明（omicverse-skills）

本目录包含从上游仓库复制的第三方 skill 源码，版权归原作者所有，按上游许可证使用。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/omicverse/omicverse-skills |
| 描述 | OmicVerse 官方 AI agent 技能库：单细胞/空间/_bulk/代谢/微生物组工作流执行指引（61 个技能，中英双语） |
| 获取时间 | 2026-10-15 |
| 获取版本 | commit `d9b0a82a14b5c8ac4190763eac7daf4f1fd2679f`（main HEAD，tarball） |
| 许可证 | MIT License, Copyright (c) 2026 Zehua Zeng（**LICENSE 原文已打开核验**，全文见本目录 `LICENSE`） |
| 上游主库许可说明 | ⚠️ 技能所驱动的 **omicverse Python 主包为 GPL-3.0**。本目录**仅收录 MIT 的技能文档**，不包含主包任何源码；主包由用户运行时从 PyPI/conda 安装，其许可随上游。同 bioservices 先例。 |

## 已包含的内容（全部 61 个技能）

| 分组 | 技能 | 与 RiceMultiOmics 的对应关系 |
|---|---|---|
| 单细胞核心 | single-cell-preprocessing / annotation / batch-integration / clustering-backends / differential-expression / differential-abundance / downstream-analysis | 模块 16 单细胞全流程 |
| 单细胞基础模型 | single-cell-foundation-model（scGPT/Geneformer/scFoundation/UCE/CellPLM）、reference-label-transfer、cross-modal-celltype-transfer、cellvote-consensus、cellmatch-ontology、popv-annotation、metatime-annotation | 模块 16 RSTA/scGPT 集成路线的直接参考 |
| 轨迹/命运 | single-cell-trajectory-inference / monocle2 / sctour / via / cellrank-fate / single-cellfate-analysis / rna-velocity / cytotrace2 / cnmf / lda-topic / scenic | 模块 16 轨迹与调控延伸 |
| 细胞通信 | single-cell-cellphonedb-communication、single-cell-liana-communication、scissor-cell-selection | 与 sciagent cellchat 互补（CellPhoneDB/LIANA 工具线）；scissor = 表型×scRNA 关联（GWAS 衔接） |
| 空间 | single-to-spatial-mapping、spatial-tutorials | 模块 17 空间组 |
| bulk/代谢 | bulk-wgcna-analysis、bulk-metabol-preprocessing / multivariate / pathway-multifactor / untargeted-lipidomics、bulk-stringdb-ppi、bulk-celltype-deconvolution、bulk-to-single-deconvolution、bulk-trajblend-interpolation | 代谢组模块 5-8（多因子通路）、WGCNA 模块 10、PPI |
| 微生物组 | microbiome-16s-amplicon-dada2 / da-comparison / meta-analysis / phylogeny、micro-metabol-paired | 根际微生物组候选扩展 |
| 数据/绘图 | data-io-loading、data-transform、data-stats-analysis、data-viz-plots、plotting-visualization、data-export-excel、data-export-pdf、datasets-loading | 通用数据操作与出版图 |
| 其他 | gsea-enrichment、tcga-preprocessing、bulk-deg/deseq2/combat、bulk-fastq-quantification、biocontext-knowledge | 富集互参；tcga 为人类数据（保留因与库内聚） |

## 明确未包含的内容

- 上游 `src/omicverse_skills/` 的 Python 包代码（`__init__.py`、加载器等）与 `pyproject.toml`、`logo/`——
  仅收录技能文档本体（各技能目录的 SKILL.md 与 references/scripts）。

## 使用说明

- 作为**方法论参考与流程规范**：需要时加载对应技能目录下的 `SKILL.md` 及其 `references/`。
- 这些技能描述的是 **omicverse 库（GPL-3.0）的 API 调用路径**；RiceMultiOmics 自研脚本不依赖 omicverse，
  引用时注意区分"方法思路"与"具体 API"。
- 二次分发本项目时，请保留本目录的 `LICENSE` 与 `NOTICE.md`。
