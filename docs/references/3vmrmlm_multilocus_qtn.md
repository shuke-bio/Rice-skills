# 多位点 QTN 检测：3VmrMLM 模型族（Fast3VmrMLM / IIIVmrMLM）

> 方法学蒸馏笔记（自有语言重写，非论文原文复制）。对应脚本
> `scripts/gwas/09_multilocus_qtn.R`。上游软件 GPL-3，按运行时安装使用，
> 本仓库不打包其源码。

## 1. 一句话定位

与 `03_association.R` 的单位点扫描（GLM/MLM/FarmCPU/Blink，一次只检验一个
标记）不同，3VmrMLM 家族是**多位点联合模型**：在全基因组扫描后把候选标记
一起放进一个混合模型，用 EM 经验贝叶斯估计所有非零效应，再逐个做似然比检验。
它同时估计**加性 + 显性**效应、并**完整控制所有多基因背景**（加性、显性及
其交互背景），因此对微效位点与稀有等位基因的检出更敏感。

## 2. 模型要点（3VmrMLM）

- 压缩方差组分混合模型：把多个方差组分压缩为少数（一般 4 → 2），把「所有
  效应」压成一个效应向量、所有多基因背景压成一个背景向量，从而可在普通
  服务器上求解；QTN 检测时先估计 AA/Aa/aa 三种基因型效应，再一元方差分析
  分解为加性/显性。
- 两步框架（所有模块共用）：
  1. **全基因组扫描**：宽松阈值筛潜在关联标记（`svpal`）；
  2. **机器学习阶段**：EM 经验贝叶斯在多基因座模型中选择非零效应，似然比
     检验定显著位点。最终显著性建议用 Bonferroni（P = 0.05 / 标记数），
     `svpal` 只决定谁进入第二阶段，**不等于**最终阈值。
- Fast3VmrMLM 在此之上叠加七项大数据技术（表型校正最小二乘扫描、PCG 预条件
  共轭梯度、条件期望向量化 Wald 检验、GRAMMAR-gamma 近似、二进制基因型结构、
  Woodbury 矩阵降维、并行计算），把 50 万标记 × 1000 个体的扫描压到数小时级。

## 3. 引擎与模块对照

| 功能 | IIIVmrMLM（CRAN） | Fast3VmrMLM（GitHub） |
|---|---|---|
| 数量性状 QTN | `method="Single_env"` | `Fast3VmrMLM()` |
| 单倍型标记（bin/gene hap） | — | `genoType="Hap"`（`c_threshold`、`numofHaplotypes`） |
| 系谱推断基因型（NC II 杂交种） | — | `population="NCII"` + `filePedigree` |
| 多基因型分子标记（k-mer/SV/lncRNA） | — | `Fast3VmrMLM_SV()` / `_PGEN()` |
| QTN × 环境（QEI） | `method="Multi_env"` + `n.en=c(2,2)` | `Fast3VmrMLM_MEJA()`（见 `qei_environment_interaction.md`） |
| QTN × QTN（上位性） | `method="Epistasis"`（无需 PLINK） | `Fast3VmrMLM_QQI()`（需 PLINK 1.9，见 `epistasis_qqi.md`） |

## 4. 脚本用法（09_multilocus_qtn.R）

```bash
# 单环境多位点（PLINK 前缀 + 表型 CSV）
Rscript scripts/gwas/09_multilocus_qtn.R \
    --genotype results/gwas/01_qc/qc_geno \
    --phenotype data/phenotype.csv --trait 1:3 \
    --pop-structure results/gwas/01_qc/pca.csv --pop-type PC \
    --outdir results/gwas/09_multilocus

# 多环境（性状×环境宽表，trait 选性状、n-en 给每个性状的环境数）
Rscript scripts/gwas/09_multilocus_qtn.R --method multi_env \
    --genotype ... --phenotype ... --n-en 2,2

# 上位性
Rscript scripts/gwas/09_multilocus_qtn.R --method epistasis ...
```

- 输入基因型可直接是 PLINK 二进制（`--genotype` 指前缀），也可给 0/1/2 剂量
  CSV + `--snp-info`（脚本自动转 PLINK）。
- 参数默认值：`svpal=0.01`、`svrad=2e4`、`svmlod=3`。
- 输出：`qtn_results.csv`（trait/marker/chr/pos/lod/additive/dominance/pve/p）、
  `association_results.csv`（可直接喂 `05_plot_gwas.R`）、`engine_raw/`（引擎原始
  xlsx/csv）、`summary.txt`、`run_info.txt`。
- `--demo`：合成 200 个体 × 400 标记、植入 2 个 QTN（S120/S260），三模式均自检。
  实测单环境 15 s；epistasis 约 6 min（400 标记两两扫描）。

## 5. 阈值与参数经验

- `svpal`：n < 4000 用 0.01；n ≥ 4000 用 1e-5；想多进候选用 0.05，想收紧用 1e-5。
- `svrad`：去共线性窗口（bp），默认 2e4；物种连锁长度不同需调。
- 上位性模块：n < 4000 时 `svpal=1e-5`、`svrad=1e5`；n > 4000 时
  `svpal = 0.05/[m(m-1)/2]`、`svrad=1e4`。
- 群体结构：有 PCA 给 `--pop-structure` + `--pop-type PC`；有 ADMIXTURE Q 矩阵
  用 `Q`；不给则仅靠亲缘矩阵控制（引擎自行计算 K）。

## 6. 多亲本群体 QTL：mppQTL（11_multiparent_qtl.R）

MAGIC / NAM / ROAM 群体的多位点 QTL 定位，两步：单标记线性混合模型全扫描
（qtl2 的 founder 设计概率）→ group-lasso 选峰 → 多基因座模型 EM 经验贝叶斯 →
LRT 定显著 QTL。适合多亲本群体的微效 QTL；双亲群体仍用常规 QTL 流程。

```bash
Rscript scripts/gwas/11_multiparent_qtl.R --population magic \
    --genotype 4way_magic_genotype.csv --phenotype 4way_magic_phenotype.csv \
    --map 4way_magic_map.csv --cross-info 4way_magic_cross.csv \
    --crosstype riself4 --map-function haldane --cores 4 \
    --outdir results/gwas/11_multiparent
```

输入格式（原生 CSV，详见脚本头注释）：基因型行=标记、列=founder 等位字母 +
各系；表型行=Trait×Environment、列=系；`cross-info` 行=系、列=`founder_A..D`，
**每行必须是 1..n_founder 的排列**（qtl2 校验）；等位基因字母只认 A/T/C/G。

输出：每个 Trait×Env 一个 xlsx（`GWS` 全扫描 + `QTL` 显著位点两个 sheet）；
脚本整理为 `qtl_results.csv`（显著 QTL）与 `gws_scan.csv`（全扫描）。

## 7. 安装（运行时，全部 GPL-3，不打包）

```r
# IIIVmrMLM（CRAN，推荐；Windows/Linux 均可，自动匹配本机 R 版本）
install.packages("IIIVmrMLM")

# Fast3VmrMLM（GitHub；R 4.3 预编译二进制，R >= 4.4 可能加载失败）
# 见 https://github.com/YuanmingZhang65/Fast3VmrMLM 的 Installation 3.0
# 需要 Rcpp/RcppArmadillo/RcppParallel/data.table/openxlsx/BH/BEDMatrix/pgenlibr

# mppQTL（GitHub zip；依赖 qtl2）
install.packages("qtl2")
# download.file("https://github.com/YuanmingZhang65/IIIVmrMLM/raw/main/mppQTL_1.0.0.zip", "mppQTL.zip", mode="wb")
# install.packages("mppQTL.zip", repos=NULL)
```

版本实测（2026-10-19，R 4.6.0 / Windows）：IIIVmrMLM 与 mppQTL 正常加载并跑通；
Fast3VmrMLM 的 R 4.3 二进制在 R 4.6 下 DLL 入口点缺失（`LoadLibrary failure`），
脚本 `--engine auto` 会自动退回 IIIVmrMLM（方法等价），需要 Fast3VmrMLM 的
百万级标记吞吐时请用 R 4.3.x 环境。

## 8. 引用

1. Li M, Zhang YW, Zhang ZC, Xiang Y, Liu MH, Zhou YH, Zuo JF, Zhang HQ,
   Chen Y, Zhang YM. A compressed variance component mixed model for detecting
   QTNs and QTN-by-environment and QTN-by-QTN interactions in genome-wide
   association studies. **Molecular Plant** 2022; 15(4): 630-650.
2. Li M, Zhang YW, Xiang Y, Zhang YM. IIIVmrMLM: the R and C++ tools associated
   with 3VmrMLM. **Molecular Plant** 2022; 15(8): 1251-1253.
3. Wang JT, Chen Y, Shu GP, Zhao MM, Zheng A, Chang XY, Li GQ, Wang YB,
   Zhang YM. Fast3VmrMLM: a fast algorithm that integrates genome-wide scanning
   with machine learning to accelerate gene mining and breeding by design for
   polygenic traits in large-scale GWAS datasets. **Plant Communications**
   2025; 6(7): 101385. doi:10.1016/j.xplc.2025.101385
4. Li G, Zhou YH, Li HF, Zhang YM. A multi-locus linear mixed model methodology
   for detecting small-effect QTLs for quantitative traits in MAGIC, NAM, and
   ROAM populations. **Comput Biol Med** / **Computational and Structural
   Biotechnology Journal** 2023; 21: 2241-2252. doi:10.1016/j.csbj.2023.03.022
