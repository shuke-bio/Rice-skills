# 场景卡：表型组与基因组预测（GS）

## 1. 适用与前提

- 数据：表型（多环境/多年份，最好有重复）+ 基因型；
- 前提：**样本 ID 严格对齐**；表型经过区组/环境校正（BLUP）；
- 目标：估计遗传力、评估预测准确率、排序候选品系。

## 2. 先做三件事

1. 表型与基因型样本取交集，记录交集前后 n；
2. 表型分布与离群检查（`phenomics/01_phenotype_blup.R`）；
3. 功效评估（`phenomics/05_power_analysis.R`）——n 太小就不要承诺预测精度。

## 3. 三档方案

### 最小可行
1. `phenomics/01_phenotype_blup.R`（BLUP/遗传力）；
2. `phenomics/05_power_analysis.R`（功效）；
3. 样本匹配表。

### 标准
1. `phenomics/03_gs_gblup.R`（VanRaden 亲缘 + 岭回归 + k 折 CV）；
2. 多环境：`phenomics/04_ammi_gge.R`（稳定性）；
3. 报告重复 CV 的均值±SD、预测-观测相关。

### 进阶
1. `phenomics/06_menet_gs.py`（MeNet 多环境预测，需 PyTorch）；
2. 多环境/多年份联合建模；
3. 与 GBLUP 基线对比（不做对比的深度模型结果不接受）。

## 4. 文献依据（白名单）

- MeNet：混合效应深度网络，水稻 12 性状×3 环境验证（Plant Communications 2025,
  doi:10.1016/j.xplc.2025.101620）；
- GBLUP 亲缘矩阵与交叉验证：方法实现见 `phenomics/03_gs_gblup.R` 与
  `docs/phenomics_usage.md`（不引用非白名单方法学原文）；
- AMMI/GGE 稳定性：实现见 `phenomics/04_ammi_gge.R`。

## 5. 失败模式

| 现象 | 原因 | 处理 |
|---|---|---|
| CV 准确率异常高（>0.9） | 数据泄漏/随机划分 | 按家系/环境分组划分 |
| 预测在新环境失效 | G×E 强 | 多环境联合模型，报告环境特有效应 |
| 遗传力估计接近 0 | 表型噪声/样本少 | 增加重复与环境，检查区组 |
| 深度模型不如 GBLUP | 样本不足/过拟合 | 以 GBLUP 为基线，缩小模型 |

## 6. 脚本映射

`scripts/phenomics/01,03,04,05,06` + `scripts/gwas/01`（基因型质控）。
