# NOTICE — 第三方组件溯源说明（K-Dense Scientific Agent Skills）

本目录包含从上游仓库复制的第三方 skill 源码，版权归原作者所有，按上游许可证使用。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/K-Dense-AI/scientific-agent-skills |
| 描述 | 177 个面向 AI agent 的科研技能（原 Claude Scientific Skills，41.5k+ stars，活跃维护） |
| 获取时间 | 2026-10-14 |
| 获取版本 | commit `a1043ef544e3a100c28e70a736b6cf32554d0677`（shallow clone） |
| 许可证 | MIT License, Copyright (c) 2025 K-Dense Inc.（全文见本目录 `LICENSE.md`） |
| 引用信息 | 见本目录 `CITATION.cff`（上游建议学术使用引用） |

## 已包含的内容（精选 28 个技能，2026-10-15 扩充）

首轮（5 个）+ 二轮（23 个），均取自同一上游 commit。

### 首轮收录（v3.1 期间）

| 技能目录 | 用途 | 与 RiceMultiOmics 的对应关系 |
|---|---|---|
| `skills/bulk-rnaseq/` | FASTQ→counts 流程规范（nf-core/rnaseq、Salmon/STAR/featureCounts、链特异性、参考基因组溯源） | 补强转录组模块的**上游**环节（我们默认从 count 矩阵开始） |
| `skills/pydeseq2/` | PyDESeq2 差异表达规范 | 与 `02_deseq2.R`（R 版 DESeq2）互参 |
| `skills/pathway-enrichment/` | 通路富集规范（过度代表分析与 GSEA 的最佳实践） | 与 `03_enrichment.R` / `04_gsea.R` 互参 |
| `skills/scientific-visualization/` | 出版级可视化规范（配色、字体、导出、多面板） | 与自研绘图主题（`config/nature_palette.json`、`plot_theme.*`）互参 |
| `skills/statistical-analysis/` | 统计分析方法选择与报告规范 | 全局方法论参考（样本量、检验选择、多重校正） |

### 二轮收录（2026-10-15，按水稻全组学互补性筛选）

| 技能目录 | 上游库许可 | 与 RiceMultiOmics 的对应关系 |
|---|---|---|
| `skills/arboreto/` | BSD-3 | **GRN 推断（GRNBoost2/GENIE3）**——共表达网络进阶的实现路径参考（含 Dask 兼容坑位说明） |
| `skills/bioservices/` | ⚠️ GPL-3（库） | UniProt/KEGG/Reactome API 接入——蛋白组 `protein_id→gene_id` 自动映射的参考 |
| `skills/statistical-power/` | MIT | 统计功效/样本量计算——GWAS/GS 实验设计参考 |
| `skills/anndata/` | BSD-3 | 单细胞数据格式（h5ad）——模块 16/17 数据结构参考 |
| `skills/scanpy/` | BSD-3 | 单细胞分析全流程——模块 16 补强 |
| `skills/scvi-tools/` | BSD-3 | 深度生成模型单细胞整合/批次校正——模块 16（scGPT 同类路线） |
| `skills/scvelo/` | BSD-3 | RNA velocity——单细胞轨迹延伸 |
| `skills/deeptools/` | BSD | ChIP/ATAC 信号矩阵与富集——模块 15 表观组 |
| `skills/matchms/` | Apache-2.0 | 质谱谱图匹配/相似度——代谢组模块质谱鉴定参考 |
| `skills/pyopenms/` | BSD-3 | 质谱数据处理——模块 14 蛋白组/代谢组参考 |
| `skills/biopython/` | Biopython（BSD 类） | 序列处理基础（KA/KS、启动子提取、FASTA/GenBank 解析） |
| `skills/etetoolkit/` | ⚠️ GPL-3+（库） | 系统发育树操作/可视化——KA/KS 进化模块参考 |
| `skills/phylogenetics/` | 上游可执行混合（MAFFT BSD / IQ-TREE GPL-2 / FastTree Apache/MIT） | MAFFT+IQ-TREE+FastTree 建树流程——进化分析参考 |
| `skills/pysam/` | MIT | BAM/BCF/VCF 程序化处理——模块 9 变异检测参考 |
| `skills/gget/` | BSD-2 | Enrichr/AlphaFold/BLAT 等快捷查询——富集与注释辅助 |
| `skills/networkx/` | BSD-3 | 网络分析——共表达/GRN 网络拓扑参考 |
| `skills/statsmodels/` | BSD-3 | 线性模型/混合模型——GWAS 统计基础参考 |
| `skills/scikit-learn/` | BSD-3 | 机器学习——GS 交叉验证、随机森林参考 |
| `skills/umap-learn/` | BSD-3 | 流形降维——组学数据可视化参考 |
| `skills/mageck/` | MIT | CRISPR 筛选分析（MAGeCK MLE/RRA）——基因编辑实验延伸 |
| `skills/primer-design/` | MIT | 引物设计（含 CRISPR/突变检测/qPCR 引物）——湿实验支撑 |
| `skills/experimental-design/` | MIT | 实验设计规范（重复、随机化、功效）——模块 12 参考 |
| `skills/cellprofiler/` | MIT | 图像表型分析流水线——模块 18 表型组参考 |

> ⚠️ **GPL 库说明**：`bioservices`、`etetoolkit` 及 `phylogenetics` 部分上游可执行文件为 GPL 系许可。
> 本目录**仅收录 K-Dense 的使用规范文档（MIT）**，不包含这些库的任何源码；
> 库本体由用户运行时从 PyPI/bioconda 安装，其许可随上游。

## 明确未包含的内容

- 上游其余 149 个技能（化学信息学、材料、天文、临床、量子计算、实验室自动化等与水稻多组学无关的领域）；
- 上游的 `docs/`、`tests/`、`scripts/`、`plugin.json` 等仓库级文件。

## 使用说明

- 作为**方法论参考与流程规范**：需要时加载对应技能目录下的 `SKILL.md` 及其 `references/`。
- 与 RiceMultiOmics 自研脚本相互独立：本包的分析脚本不依赖这些文件。
- 上游仍在活跃更新；如需最新版本请从上述仓库获取。二次分发本项目时，请保留本目录的 `LICENSE.md`、`CITATION.cff` 与 `NOTICE.md`。
