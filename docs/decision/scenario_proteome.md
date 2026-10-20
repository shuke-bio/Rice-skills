# 场景卡：蛋白质组

## 1. 适用与前提

- 数据：DIA/DDA 定量矩阵（LFQ/iBAQ/TMT）+ 样本表；
- 前提：报告缺失率与插补策略（MNAR 常见）；缺失 >50% 的蛋白不做差异；
- 目标：差异蛋白、功能富集、蛋白-基因映射、PTM、与转录组联合。

## 2. 三档方案

### 最小可行
1. `proteomics/01_preprocess.py`：缺失过滤/插补、归一化；
2. 样本相关性/PCA；
3. 缺失率与插补策略报告。

### 标准
1. `proteomics/02_diff_protein.R`（limma）；
2. `proteomics/03_plot_proteome.R`（火山/热图）；
3. `proteomics/04_uniprot_mapping.py`（蛋白→基因）；
4. `proteomics/05_interpro_domains.py`（结构域富集）。

### 进阶
1. 磷酸化等 PTM 分析（需 PTM 数据）；
2. 与转录组联合（蛋白-mRNA 一致性）；
3. 通路/网络整合。

## 3. 文献依据（白名单）

- 水稻蛋白质组全景与 m6A 的转录后调控（Nature Plants 2024）——蛋白组与表观/翻译调控联合；
- 质谱蛋白质组与 m6A（Nature Plants 2024）——蛋白丰度与修饰的整合设计；
- 蛋白组方法学与缺失值处理：实现见 `docs/proteomics_usage.md`（不引用非白名单原文）。

## 4. 失败模式

| 现象 | 原因 | 处理 |
|---|---|---|
| 差异蛋白很少 | 缺失率高/归一化不当 | 检查缺失机制，按组过滤 |
| 蛋白与 mRNA 不一致 | 转录后调控真实存在 | 作为结论而非错误 |
| 批次主导 PCA | 建库/上机批次 | 批次入模型或分层 |
| PTM 位点不可信 | 富集不足/定位概率低 | 提高阈值，报告 localization probability |

## 5. 脚本映射

`scripts/proteomics/01–05`。
