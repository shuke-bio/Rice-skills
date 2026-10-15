# NOTICE — 第三方组件溯源说明（bioSkills）

本目录包含从上游仓库复制的第三方 skill 源码，版权归原作者所有，按上游许可证使用。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/GPTomics/bioSkills |
| 描述 | 562 个面向 AI agent 的生物信息学 SKILL.md（按领域组织） |
| 获取时间 | 2026-10-14 |
| 获取版本 | commit `d91ed3d563019e649dc854c56ccd62551359488a`（shallow clone） |
| 许可证 | MIT License, Copyright (c) 2026 GPTomics（全文见本目录 `LICENSE`） |
| 维护状态 | 上游仓库已归档（archived，不再更新）——内容为冻结快照，参考使用不受影响 |

## 已包含的内容（与 RiceMultiOmics 分析模块对口的 10 个类别，共 90 个技能）

| 类别目录 | 技能数 | 与 RiceMultiOmics 的对应关系 |
|---|---|---|
| `differential-expression/` | 6 | 补强转录组模块：批次校正、时间序列 DE、DE 结果与可视化 |
| `metabolomics/` | 9 | 补强代谢组模块：XCMS/MS-DIAL 上游预处理、代谢物注释、脂质组、同位素追踪、靶向分析 |
| `multi-omics-integration/` | 5 | 补强整合模块：数据协调、整合设计、MOFA/mixOmics 分析、相似性网络 |
| `pathway-analysis/` | 6 | 补强富集分析：GO/GSEA/KEGG/Reactome/WikiPathways |
| `data-visualization/` | 20 | 补强绘图：火山图/MA、UpSet、Circos、多面板、网络图、统计标注等 |
| `expression-matrix/` | 5 | 补强数据准备：计数矩阵导入、基因 ID 映射、元数据连接、归一化 |
| `variant-calling/` | 13 | 补强基因组学模块（v2.0）：GATK/DeepVariant 变异检测、VCF 规范化/统计/过滤、结构变异、注释 |
| `population-genetics/` | 7 | 补强基因组学模块（v2.0）：关联检验、连锁不平衡、PLINK、群体结构、选择统计 |
| `primer-design/` | 4 | 补强设计模块（v2.1）：qPCR 引物设计、引物特异性检查、引物验证基础 |
| `crispr-screens/` | 15 | 补强设计模块（v2.1）：CRISPR 筛选与 sgRNA 设计相关知识 |

## 二轮补充（2026-10-15，31 个类别共 264 个技能）

上游扩容后按水稻全组学缺口补收（同一 commit `d91ed3d`，已核验 MIT 原文）：

| 批次 | 类别（技能数） | 对应模块 |
|---|---|---|
| 单细胞/空间 | single-cell(17)、spatial-transcriptomics(12) | 模块 16/17 |
| 表观组 | chip-seq(12)、atac-seq(12)、methylation-analysis(10)、hi-c-analysis(9) | 模块 15 |
| 蛋白组/RNA 生物学 | proteomics(9)、alternative-splicing(9)、small-rna-seq(6)、epitranscriptomics(5) | 模块 14、miRNA/剪接 |
| 调控/进化 | gene-regulatory-networks(6)、comparative-genomics(13)、phylogenetics 未收（K-Dense 已覆盖） | v3.2 GRN、KA/KS |
| 基因组/变异 | genome-assembly(9)、genome-annotation(7)、genome-intervals(8)、genome-engineering(5)、copy-number(11)、phasing-imputation(4)、long-read-sequencing(9) | 泛基因组/QTL、CRISPR 编辑、群体遗传 |
| 因果/关联 | causal-genomics(11) | v3.2 共定位/MR 延伸 |
| 上游流程 | sequence-io(9)、sequence-manipulation(7)、database-access(15)、alignment(7)、alignment-files(10)、read-qc(7)、read-alignment(4)、rna-quantification(4) | FASTQ→count 上游补强 |
| 基建 | machine-learning(6)、reporting(6)、experimental-design(5) | GS/ML、报告、实验设计 |

> 注：首轮 10 类（90 技能）为早期快照，未随本轮同步更新；上游仓库已归档，无版本漂移风险。

## 明确未包含的内容

- 上游人类临床/免疫/化学类：clinical-databases、clinical-biostatistics、epidemiological-genomics、
  liquid-biopsy、immunoinformatics、tcr-bcr-analysis、chemoinformatics(20)、flow-cytometry、
  imaging-mass-cytometry、restriction-analysis、ribo-seq、clip-seq、rna-structure、
  systems-biology、metagenomics、microbiome、workflows(41)、resources、install 脚本等
  （单细胞通信已由 sciagent 的 cellchat 技能覆盖；宏基因组/根际微生物组为候选扩展，见 roadmap）。

## 使用说明

- 作为**方法论参考与流程规范**：需要时加载对应技能目录下的 `SKILL.md`（frontmatter 含 `name`/`description`/`tool_type` 等字段）。
- 与 RiceMultiOmics 自研脚本相互独立：本包的分析脚本不依赖这些文件；它们用于补充专家经验（版本兼容性、常见陷阱、参数选择）。
- 二次分发本项目时，请保留本目录的 `LICENSE` 与 `NOTICE.md`。
