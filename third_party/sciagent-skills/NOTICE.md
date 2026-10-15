# NOTICE — 第三方组件溯源说明（SciAgent-Skills 精选子集）

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/jaechang-hits/SciAgent-Skills |
| 描述 | 203 个面向 AI agent 的生信/生命科学技能（BixBench-Verified-50 准确率 92%），12 个类别 |
| 获取时间 | 2026-10-14 |
| 获取版本 | commit `82c862c`（2026-09-29） |
| 许可证 | **CC BY 4.0**（Creative Commons Attribution 4.0 International, Copyright (c) 2024 jaechang-hits）——全文见本目录 `LICENSE` |
| 许可要点 | 允许复制、再分发、改编（含商业用途）；条件：**署名** + 注明是否修改 + 不得施加额外限制 |
| 上游附加声明 | 上游 LICENSE 文件载明：项目引用了 K-Dense-AI/claude-scientific-skills（MIT）与 Biomni（Apache-2.0）相关内容 |

## 已包含的内容（8 个类别，共 64 个技能）

| 类别目录 | 技能数 | 代表技能与用途 |
|---|---|---|
| `skills/systems-biology-multiomics/` | 11 | MOFA+ 多组学整合、KEGG/Reactome 通路、STRING PPI、组学分析指南 |
| `skills/scientific-writing/` | 24 | **期刊图规范（Nature/Cell/Science/Lancet/NEJM/PNAS/eLife 分刊指南）**、论文写作、文献检索、审稿方法 |
| `skills/data-visualization/` | 8 | 组学绘图、多面板图、显著性标注、seaborn/matplotlib/plotly 规范 |
| `skills/genomics-bioinformatics/rnaseq/` | 5 | DESeq2/PyDESeq2 差异表达、featureCounts、Salmon、gseapy 富集 |
| `skills/biostatistics/` | 4 | 统计分析、statsmodels、PyMC 贝叶斯、scikit-survival |
| `skills/genomics-bioinformatics/qc/` | 3 | fastp 预处理、MultiQC 报告、BUSCO 解读 |
| `skills/genomics-bioinformatics/variant/` | 6 | 补强基因组学模块（v2.0）：bcftools 操作、GATK、PLINK2-GWAS、SnpEff 注释、VCF 过滤、CNVkit |
| `skills/molecular-biology/` | 3 | 补强设计模块（v2.1）：sgRNA 设计指南、质粒注释（pLannotate）、RNA 二级结构（ViennaRNA） |

## 明确未包含的内容

- 上游其余约 149 个技能（单细胞、变异检测、结构生物学、药物发现、医学影像、实验室自动化等）；
- 上游仓库级文件（assets/、tests/、scripts/、templates/、legacy/ 等）。

## 改动说明（CC BY 4.0 要求注明改动）

本子集为**原样复制**：未修改任何技能文件内容，仅按与"水稻多组学 + 科研写作"的相关性
挑选了上述 6 个类别。二次分发本项目时请保留本目录的 `LICENSE` 与 `NOTICE.md`，
并保留对作者 jaechang-hits 的署名。
