# 标记成对上位性 GWAS（QTN×QTN，marker-pair epistasis）：Fast3VmrMLM

> 方法学蒸馏笔记（自有语言重写）。对应脚本
> `scripts/gwas/09_multilocus_qtn.R --method epistasis`。
> Fast3VmrMLM 为 GPL-3（包 DESCRIPTION 声明；仓库 LICENSE 文件写 BSD-3，
> 两者不一致，按保守原则当 GPL 处理），运行时安装，不打包。

## 1. 为什么上位性检测难

上位性（基因×基因互作）是「缺失遗传力」的重要来源。全基因组两两扫描的组合
数是 m(m-1)/2：1000 个标记约 50 万对，50 万个标记是 1.25×10¹¹ 对。既有方法
要么无法枚举（非穷举法忽略群体结构与多基因背景，假阳性高），要么只处理
二分类性状（BOOST/GBOOST），要么在作物群体上功效不足。

## 2. 方法框架（三步）

1. **全基因组两两扫描**：先用 PLINK `--epistasis` 对 aa / ad / da / dd 四种
   互作效应做穷举初筛（表型先经残差混合模型校正群体结构与多基因背景），
   保留最重要的约 100 万对。
2. **潜在关联对筛选**：条件期望 + GRAMMAR-gamma 近似的向量化 Wald 检验
   （并行）；对高 LD 区段的候选对做去共线性。
3. **加速经验贝叶斯 + LRT**：PCG 与 Woodbury 矩阵技术估计所有主效与互作效应，
   对非零效应对逐个似然比检验，输出显著（significant）与提示性（suggested）
   的 QQI（QTN-by-QTN interaction）。

**扩展模块**：`-Hap`（单倍型×单倍型互作）、`-NCII`（系谱推断基因型的育种群体）、
`-mQTL`（k-mer/SV/功能单倍型等多基因型标记的互作）。

## 3. 实测性能（论文报告）

- 模拟：功效 68.33%、FPR 0.0056‱、F1 = 0.70；优于 PLINK、FaST-LMM、REMMA、
  Next-Gen GWAS 与 3VmrMLM。
- 规模：50 万标记 × 1000 个体，7.23 h（30 线程）。
- 水稻真实数据（1439 / 3K / 18K）：58 个有实验证据的 GGIs + 219 个生信支持的
  候选 GGIs；平均每性状 0.26 / 0.91 / 5.45 h。
- 把 GGIs、基因、G×E 整合成遗传网络，18K 水稻中筛出 24 个产量关键基因、
  21 对优异基因单倍型组合。

## 4. 脚本用法

```bash
# IIIVmrMLM 路线（无需 PLINK，Windows 可用）
Rscript scripts/gwas/09_multilocus_qtn.R --method epistasis \
    --genotype results/gwas/01_qc/qc_geno \
    --phenotype data/phenotype.csv --trait 1 \
    --outdir results/gwas/09_epistasis

# Fast3VmrMLM 路线（大数据；需 PLINK 1.9 在 PATH）
Rscript scripts/gwas/09_multilocus_qtn.R --method epistasis \
    --engine Fast3VmrMLM ...
```

- 阈值经验：n < 4000 → `svpal=1e-5`、`svrad=1e5`；n > 4000 →
  `svpal = 0.05/[m(m-1)/2]`（m 为进入扫描的标记数）、`svrad=1e4`。
- 输出：`result.Q`（主效应 QTN）与 `result.QQI`（互作对）两个 sheet，
  脚本整理进 `qtn_results.csv`。
- `--demo`：200 个体 × 400 标记，IIIVmrMLM 实测约 6 min。

## 5. 解读与陷阱

- **组合爆炸**：标记数 > 10⁴ 时务必先做 LD 抽稀（脚本默认开启）。
- **假阳性**：两两扫描的检验数极大，必须用 Bonferroni 或置换阈值；报告时
  区分 significant 与 suggested。
- **生物学验证**：GGI 候选强烈依赖后续表达/单倍型/实验证据（论文中 58 个
  实验支持的互作是下限，219 个生信候选需验证）。
- **网络构建**：把 QQI 与 QTN、QEI 合并成遗传网络再挑 hub 基因，比单看
  显著对更稳健。

## 6. 引用

1. Wang JT, Zhao MM, Han XL, Chang XY, Ma XS, Jiang QY, Zhang YM. Fast3VmrMLM
   enables statistically powerful and computationally efficient large-scale
   marker-pair epistasis GWAS in sizable populations. **Plant Communications**
   2026 (accepted 2026-10-13). doi:10.1016/j.xplc.2026.102155
2. Li M, et al. A compressed variance component mixed model for detecting QTNs
   and QTN-by-environment and QTN-by-QTN interactions in GWAS.
   **Molecular Plant** 2022; 15(4): 630-650.
