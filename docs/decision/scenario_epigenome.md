# 场景卡：表观基因组（WGBS / ATAC / ChIP）

## 1. 适用与前提

- WGBS：差异甲基化（DMR）；需报告覆盖度与亚硫酸盐转化率；
- ATAC/ChIP：差异可及/差异结合（DAR/DBR）；需报告 FRiP、TSS 富集、重复一致性；
- 前提：生物学重复 ≥2；材料与发育时期一致；批次记录完整。

## 2. 三档方案

### 最小可行
1. `epigenomics/01_methylation_analysis.R`（WGBS）或 `02_atac_chip_analysis.R`（ATAC/ChIP）；
2. 样本相关性/覆盖度检查；
3. 批次与技术协变量整理。

### 标准
1. DMR/DAR 检测（报告阈值与平滑参数）；
2. 区域注释（启动子/增强子/基因体）；
3. 与 RNA-seq 联合（若有）。

### 进阶
1. 三维基因组/染色质互作（Hi-C，若有）；
2. 与 GWAS/QTL 共定位；
3. 多组学整合（MOFA/DIABLO）。

## 3. 文献依据（白名单）

- 水稻单细胞染色质可及性图谱（Nature Plants 2025）——scATAC 与跨物种比较；
- 3D 表观基因组调控抽穗期（Genome Biology 2026）——Hi-C + 表观 + 表达的整合；
- 3D 基因组/表观/转录整合（NAR 2023）——多维修饰与节律基因表达；
- H3K9ac 标记的 3D 互作（Genome Biology 2024）——干旱下的增强子-启动子互作。

## 4. 失败模式

| 现象 | 原因 | 处理 |
|---|---|---|
| 覆盖度不足 | 测序深度低 | 补测或限定高覆盖区域 |
| 重复间相关性低 | 批次/材料差异 | 检查取材与建库批次 |
| DMR 全在重复区 | 比对歧义 | 过滤重复区，报告过滤规则 |
| ATAC 无 TSS 富集 | 建库失败 | 重做文库，不可强行分析 |

## 5. 脚本映射

`scripts/epigenomics/01_methylation_analysis.R`、`02_atac_chip_analysis.R`。
