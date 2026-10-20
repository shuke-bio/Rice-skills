# 场景卡：QTL 与 BSA

## 1. 适用与前提

- 家系设计：RIL / F2 / MAGIC / NAM / ROAM；**需要遗传图谱或高密度标记**；
- BSA-seq：双亲极端混池；需要足够混池个体与测序深度；
- 前提：表型可靠（区组/重复）；家系结构与系谱准确。

## 2. 先选设计

| 设计 | 建议规模 | 方法 |
|---|---|---|
| F2 / RIL | ≥150 个体/家系 | 连锁分析（LOD，置换检验） |
| MAGIC / NAM / ROAM | ≥300 家系 | 多亲本 QTL（两步法） |
| BSA-seq | 每个混池 ≥30 个体，≥30× | ΔSNP-index / G′ |
| 自然群体 | ≥200 | 关联分析（见 GWAS 场景卡） |

## 3. 三档方案

### 最小可行
1. `qtl/01_qtl_mapping.R`：遗传图谱/连锁分析；
2. 标记质控（若为芯片/重测序）：`gwas/01_plink_qc.sh`；
3. 表型分布与区组检查。

### 标准
1. `qtl/01_qtl_mapping.R`（LOD 阈值用置换检验）；
2. `qtl/04_qtl_region_snps.py`（区间 SNP 提取与注释）；
3. `qtl/05_bsa_scan.py`（BSA：ΔSNP-index/G′，模拟法置信线）。

### 进阶
1. `qtl/02_colocalization.R` / `03_coloc_bayesian.R`（与 GWAS/表达共定位）；
2. `gwas/07_haplotype.R`（候选基因单倍型）；
3. `design/01_sgrna_design.py`（编辑验证靶点）。

## 4. 文献依据（白名单）

- 多位点 QTN 模型（Molecular Plant 2022, doi:10.1016/j.molp.2022.02.012）——
  当 QTL 区间内候选位点多、连锁强时，可用多位点模型辅助缩小区间；
- 多亲本群体 QTL（Plant Communications 2025）——MAGIC/NAM/ROAM 的两步法；
- 水稻 QTL 精细定位与单倍型（JIPB 2024）——区间到候选基因的验证路径；
- 泛基因组 SV 关联（Nature Communications 2025）——结构变异可作为 QTL 区间的补充标记。

## 5. 失败模式

| 现象 | 原因 | 处理 |
|---|---|---|
| LOD 曲线全基因组都高 | 表型/家系错误 | 检查系谱、区组与表型分布 |
| 区间过大 | 标记密度低/重组少 | 增加标记或改用 BSA |
| BSA 两个混池无差异 | 混池个体不足/表型错误 | 增加个体与深度，重测表型 |
| 候选基因无法验证 | 区间仍有多个基因 | 单倍型 + 表达 + 编辑验证 |

## 6. 脚本映射

`scripts/qtl/01–05` + `scripts/gwas/01,07` + `scripts/design/01`。
