# 二分类 / case-control 性状 GWAS：FastBiCmrMLM

> 方法学蒸馏笔记（自有语言重写）。对应脚本 `scripts/gwas/10_case_control_qtn.R`。
> 上游软件**无 LICENSE 文件**（GPL 风格研究软件，仅提供 Linux 二进制），
> 按运行时安装使用，本仓库不打包；使用前请自行评估并引用原文。

## 1. 为什么需要专门的方法

作物抗性性状（抗/感）、人类疾病（病例/对照）是**二分类**表型。把 0/1 当连续
性状拟合线性混合模型会：低估效应、放松方差假设、损失功效，且显性效应无法
被正确刻画。FastBiCmrMLM 用**混合逻辑回归**（mixed logistic model）代替线性
模型，同时保留 3VmrMLM 的压缩方差组分思想。

## 2. 模型要点

- 混合逻辑回归全模型：同时估计加性 + 显性效应及其多基因背景；把四个方差组分
  压缩为两个，降低求解复杂度。
- 三种计算实现按样本量分派（脚本 `--method auto` 自动选择）：

| 实现 | 适用样本量 | 说明 |
|---|---|---|
| `FastBiCmrMLM` | n < 1,000 | 主算法 |
| `FastBiCmrMLM-Time` | 1,000 ≤ n < 20,000 | 比主算法快数倍（WTCCC1 单疾病约 7 min） |
| `FastBiCmrMLM-RAM` | n ≥ 20,000 | 省内存；50 万 × 100 万标记约 13.5 h / 10 CPU / 141 GB |
| `FastBiCmrMLM-Hap` | 任意 | 用 LD 标记构建 bin 单倍型作标记（`c`、`NumHap`） |

- 与 GMMAT / SAIGE / fastGWA-GLMM 相比：在显性、小等位替代效应、稀有变异上
  检出更高；UK Biobank 规模可检出效应低至 0.03% 的变异。
- 也支持结构变异（前三种实现均可检测 SV）。

## 3. 脚本用法（10_case_control_qtn.R）

```bash
# Linux / WSL2（软件只有 Linux 版）
Rscript scripts/gwas/10_case_control_qtn.R \
    --genotype results/gwas/01_qc/qc_geno \
    --phenotype data/binary_trait.csv --trait blast_resistance \
    --method auto \
    --outdir results/gwas/10_case_control

# 连续表型按阈值二值化（<= 3 → 对照，> 3 → 病例）
Rscript scripts/gwas/10_case_control_qtn.R ... --cut-to-binary 3
```

- 输入：PLINK 前缀（或剂量 CSV + `--snp-info`）；表型 CSV 为 0/1 编码
  （`sample_id` + 性状列），多列时用 `--trait` 选一列。
- 参数默认：`--search-radius 1e6`（去共线性窗口 bp）、`--svpal 1e-5`、
  `--lod-threshold 3`。
- 输出：`qtn_results.csv`、`association_results.csv`、`engine_raw/`、
  `summary.txt`、`run_info.txt`。
- **Windows 明确不支持**：脚本在非 Linux 主机上给出 SKIP 与替代方案
  （数量性状走 `09_multilocus_qtn.R` 的 IIIVmrMLM）。

## 4. 安装（运行时）

```bash
# Linux / WSL2
wget https://github.com/YuanmingZhang65/FastBiCmrMLM/raw/main/FastBiCmrMLM_Linux.zip
unzip FastBiCmrMLM_Linux.zip
R -e 'install.packages("FastBiCmrMLM", repos = NULL)'
```

依赖：Rcpp/RcppArmadillo/data.table/RcppParallel/BH/bigmemory/dplyr/glmnet/SPAtest。

## 5. 选型建议

- 病例/对照比例悬殊（< 10%）时优先 `-Time`/`-RAM`，并用 `--svpal` 收紧阈值。
- bin 单倍型（`-Hap`）在稀有变异富集区段更有效，但解释成本更高。
- 数量性状不要用本方法：走 3VmrMLM 家族（09）；二分类性状用本方法。

## 6. 引用

Wang JT, Chang XY, Zhao Q, Zhang YM. FastBiCmrMLM: a fast and powerful
compressed variance component mixed logistic model for big genomic
case-control genome-wide association study. **Briefings in Bioinformatics**
2024; 25(4): bbae290. doi:10.1093/bib/bbae290
