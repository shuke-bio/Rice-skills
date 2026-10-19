# QTN × 环境互作（QEI/QMI）：大规模气候-表型-基因组联合分析

> 方法学蒸馏笔记（自有语言重写）。对应脚本
> `scripts/gwas/09_multilocus_qtn.R --method multi_env`（IIIVmrMLM / Fast3VmrMLM
> 的多环境模块）。

## 1. 问题背景

基因型 × 环境互作（G×E）决定品种在不同年份、地点、栽培管理下的表现差异。
传统做法是分环境单独做 GWAS 再比较，功效低且无法区分「主效」与「互作效应」。
多环境联合分析把同一性状的各环境表型放进**一个混合模型**，同时估计：

- **QTN 主效应**（跨环境稳定的位点）
- **QEI**（QTN × Environment，位点效应随环境变化）
- 若管理措施（施肥/灌溉等）也编码为环境因子，则同时得到 **QMI**
  （QTN × Management）

## 2. 模型与识别策略

- 联合模型：表型向量按「性状 × 环境」展开，固定效应含环境主效应，随机效应
  含 QTN 主效应 + QEI + 全多基因背景（压缩方差组分）。显著性检验同样走
  「全基因组扫描 → 机器学习多基因座 → LRT」两步。
- IIIVmrMLM 的 `method="Multi_env"`：给 `n.en = c(2,2)` 指明每个性状的环境数；
  输出 `result.Q`（QTN 主效应）与 `result.QEI` 两个 sheet。
- Fast3VmrMLM 的 `Fast3VmrMLM_MEJA()`：同一思路的大数据实现，参数 `n_en`。
- 大规模气候数据（NSR 2026 论文）：把逐日气象因子聚合成环境协变量
  （温度/降水/辐射的均值与极端值等），与基因组一起做联合扫描，定位
  「气候响应型 QTN」，用于设计育种（breeding by design）——即挑选在目标
  气候条件下表达最优的等位组合。**论文验证数据：玉米、水稻、大豆的多环境
  联合分析（396 个已知基因 + 84 个已知 GEI 复现），另有 12 环境玉米数据
  解析开花期可塑性——水稻在验证物种之内。**

## 3. 数据准备要点

1. 表型宽表：`sample_id` + **性状×环境列**，同一性状的各环境列必须相邻
   （trait-major）；每列一个环境。
2. 环境标签建议写进列名（如 `trait1__E1`、`trait1__E2`），脚本会原样传给引擎。
3. 至少 2 个环境才能识别 QEI；环境数越多，QEI 功效越高但计算量线性增长。
4. 气象协变量：若要做「气候响应」分析，把环境级协变量（年均温、生长季降水、
   干旱指数等）整理成环境×样本的设计矩阵，作为额外协变量或环境分层变量。

## 4. 脚本示例

```bash
# 两个性状 × 各两个环境
Rscript scripts/gwas/09_multilocus_qtn.R --method multi_env \
    --genotype results/gwas/01_qc/qc_geno \
    --phenotype data/phenotype_wide.csv \
    --trait 1:2 --n-en 2,2 \
    --outdir results/gwas/09_multilocus_me

# --demo 自带合成多环境数据（两个性状 × 两个环境，含一个 QEI 位点）
Rscript scripts/gwas/09_multilocus_qtn.R --demo --method multi_env
```

输出中的 `qtn_results.csv` 会带 `result.Q` / `result.QEI` 的 sheet 标签；
显著 QEI 位点即「效应随环境变化」的候选，可进一步结合环境梯度做反应规范
（reaction norm）解读。

## 5. 结果解读

- **QTN 主效应**：跨环境稳定，适合广适性育种。
- **QEI**：只在特定环境显著，适合生态区专用品种或气候适应型设计。
- 强 QEI 且方向随环境翻转的位点，警惕「均值-方差」伪互作：检查各环境
  表型方差是否差异过大，必要时先做方差齐性变换。
- QEI 的显著性同样建议 Bonferroni（0.05 / 标记数）；多环境联合扫描的标记数
  按位点数计，不按「位点 × 环境」计。

## 6. 引用

1. Wang JT, Han XL, Zhao MM, Zhang HQ, Chen Y, Jiang QY, Zhang YM. A fast
   method for breeding by design via G × E interactions detected in large-scale
   climatic, phenomic and genomic data. **National Science Review** 2026;
   13(10): nwag095. doi:10.1093/nsr/nwag095
2. Li M, et al. A compressed variance component mixed model for detecting QTNs
   and QTN-by-environment and QTN-by-QTN interactions in GWAS.
   **Molecular Plant** 2022; 15(4): 630-650.
