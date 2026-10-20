# 场景卡：bulk RNA-seq（转录组）

## 1. 适用与前提

- 数据：count 矩阵（整数）+ 样本表（分组、批次、协变量）；
- 前提：有生物学重复（n≥2，推荐 n≥3）；参考基因组与注释版本一致；
- 目标：差异表达、标志物、通路/调控、与代谢组联合。

## 2. 三种设计对应三种模型

| 设计 | 设计公式（示例） | 说明 |
|---|---|---|
| 处理 vs 对照 | `~ batch + group` | 批次入模型，不然后果是把技术差异当处理效应 |
| 时间序列 | `~ group + time + group:time` | 关注交互项，而非各时间点单独检验 |
| 多环境/多地点 | `~ environment + group + environment:group` | 环境作为固定效应；环境数少时不做随机效应 |

## 3. 三档方案

### 最小可行（先看数据是否可信）
1. `01_qc_filter.py`：低表达过滤、样本相关性、PCA；
2. `02_deseq2.R`：仅在有重复时给差异表；n=1 只出描述性汇总；
3. `07_plot_transcriptome.R`：PCA/相关性/火山图初版。

**产出**：质控报告、PCA、差异表（若可用）。

### 标准（可解释的差异 + 通路）
1. `01_qc_filter.py` → `02_deseq2.R`（含批次）；
2. `03_enrichment.R`（GO/KEGG）→ `04_gsea.R`（排序基因集）；
3. `06_tf_analysis.R`（差异 TF 与靶基因）→ `07_plot_transcriptome.R`；
4. 敏感性分析：edgeR 与 DESeq2 取交集，报告一致比例。

**产出**：差异表、富集/GSEA 图、TF-靶基因表、出版级图。

### 进阶（网络与联合）
1. `05_wgcna.R`（n≥15 时建模块-性状关联）；
2. `integration/01_mofa.R` 或 `02_diablo.R`（与代谢组匹配样本）；
3. `integration/05_network_integration.R`（基因-代谢物网络）；
4. `qtl/02_colocalization.R`（有基因型时做 eQTL/mQTL 共定位）。

## 4. 文献依据（方法学）

- **DESeq2**：负二项 + 收缩估计，适合小 n、count 数据（Love et al., Genome Biol 2014, doi:10.1186/s13059-014-0550-8）；
- **edgeR**：精确检验/GLM，与 DESeq2 互补，用于敏感性分析（Robinson et al., Bioinformatics 2010, doi:10.1093/bioinformatics/btp616）；
- **limma-voom**：大样本或非 count 平台（Law et al., Genome Biol 2014, doi:10.1186/gb-2014-15-2-r29）；
- **WGCNA**：模块-性状关联（R 包实现；本目录不引用非白名单期刊的方法学原文，仅作为标准工具使用）；
- **GSEA**：排序基因集富集（Subramanian et al., PNAS 2005, doi:10.1073/pnas.0506580102）。

### 高分论文的设计启发（水稻）

- **单细胞/空间图谱类**（Nature 2025；Cell 2026）→ 先建参考图谱，再做标签迁移；bulk 数据可作为验证层；
- **多组学群体研究**（Genome Biol 2025，精米代谢组 + 高光谱 + mGWAS，doi:10.1186/s13059-025-03513-w）→ 表型代理 + 关联 + 共定位的三段式；
- **泛转录组/可变剪接**（Plant Cell 2024）→ 若关注剪接，需在差异之外单独做剪接分析，不能只看基因水平。

## 5. 失败模式与诊断

| 现象 | 可能原因 | 处理 |
|---|---|---|
| 差异基因上千、富集无结果 | 批次/质量差异主导 | 看 PCA 与批次交叉表，重跑含批次模型 |
| 两组完全分开但无差异 | 样本量不足或方差过大 | 报告效应量，不做强结论 |
| 候选基因与文献不符 | 注释版本不一致 | 统一 MSU/RAP-DB 并转换 ID |
| 富集结果随背景基因集变化 | 背景选择不当 | 用检测到的基因作背景，而非全基因组 |

## 6. 脚本映射

`scripts/transcriptomics/01–07` + `scripts/integration/01,02,05` + `scripts/qtl/02`。
