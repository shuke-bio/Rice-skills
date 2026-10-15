# 表型组学模块使用指南（rice-phenomics，v3.0）

> Rice-skills v3.0 新增 **模块 16 · 表型组学**：表型数据清洗、BLUP 与遗传力估计、
> 光谱指数计算，以及图像表型与基因组选择（GS）的工具导航。
> 输出直接衔接模块 7（GWAS）与育种应用。

## 0. 流程总览

```
多点试验表型 ──01_phenotype_blup.R──▶ BLUP 表（GWAS/GS 的标准输入）+ 遗传力
冠层光谱 ──02_spectral_indices.R──▶ 植被指数表
图像数据 ──▶ PlantCV（外部工具，见 §4）
BLUP + 基因型 ──▶ 模块 7 GWAS / rrBLUP·GBLUP 基因组选择（见 §5）
```

## 1. 表型 QC 与 BLUP（`scripts/phenomics/01_phenotype_blup.R`）

```bash
Rscript scripts/phenomics/01_phenotype_blup.R \
    --pheno demo_data/phenomics/phenotype_long.csv \
    --env-col env --traits grain_yield,plant_height \
    --outdir results/phenomics/01_blup
```

**输入**：长表 `line_id, env, trait1, trait2, ...`（多年多点）；单环境表省略 `--env-col`。

**功能**：
- QC：MAD 法异常值标记（只标不删，输出 `outliers.csv` 供人工复核）；
- BLUP：环境调整 + 收缩加权（`blups.csv`）——**推荐作为 GWAS/GS 的表型输入**（消除环境与重复数不平衡）；
- 遗传力：`heritability.csv`（h2 = σ²g / (σ²g + σ²e/r̄)，line-mean 基础）；
- 图：各环境性状分布、性状相关性。

**方法说明**：方差组分用矩估计（ANOVA）；lme4 若安装会自动额外输出 REML 对照值。
严格分析（系谱、复杂设计）建议 lme4/ASReml，代码片段见脚本注释。

## 2. 光谱指数（`scripts/phenomics/02_spectral_indices.R`）

```bash
Rscript scripts/phenomics/02_spectral_indices.R \
    --spectra demo_data/phenomics/spectra.csv --outdir results/phenomics/02_indices
```

- 输入：宽表（sample_id + 波长列，列名含波长如 `X800`/`800nm`）；
- 自动定位波段（NIR 800 / Red 670 / Blue 470 / Red-edge 730 / Green 550 nm，±15 nm 容差，可 `--nir` 等覆盖）；
- 计算：**NDVI、EVI、SAVI、RVI、NDRE、WDRVI、GNDVI**；
- 输出：指数表 + 分布图 + 相关性图。指数可直接作为性状进 GWAS（替代/补充人工表型）。

## 3. 与 GWAS / 育种的衔接

```bash
# BLUP 表 -> 模块 7 GWAS（把 blups.csv 作为 phenotype.csv，列名对应 --trait）
Rscript scripts/gwas/03_association.R --genotype ... --snp-info ... \
    --phenotype results/phenomics/01_blup/blups.csv --trait grain_yield
```

## 4. 图像表型（外部工具指引）

| 工具 | 用途 | 获取 |
|---|---|---|
| **PlantCV**（MPL-2.0） | RGB/多光谱图像分割与形态特征（株高/叶面积/绿度） | `conda install -c conda-forge plantcv` · plantcv.danforthcenter.org |
| DeepLabCut | 器官关键点追踪（生长动态） | 运行时安装（PyPI） |
| CellProfiler | 批量图像特征提取 | cellprofiler.org |

典型 PlantCV 流程：`pcv.readimage → pcv.threshold.* → pcv.fill → pcv.analyze.size`，
输出形态学特征表 → 与本模块的 BLUP/GWAS 流程衔接（同表型长表格式）。

## 5. 基因组选择（GS，外部工具指引）

| 包 | 模型 | 获取 |
|---|---|---|
| **rrBLUP** | rrBLUP/GBLUP（R 包，经典） | `install.packages("rrBLUP")` |
| BGLR | 贝叶斯（BayesA/B/C） | `install.packages("BGLR")` |
| sommer | 混合模型（多性状/多环境） | `install.packages("sommer")` |

典型 GBLUP（RR-BLUP 等价）：
```r
# 基因型矩阵（marker x line，0/1/2）→ A.mat 构造亲缘 G → kin.blup
library(rrBLUP)
G <- A.mat(t(geno))
fit <- kin.blup(data = blups, geno = "line_id", pheno = "grain_yield", K = G)
# 预测值 fit$g；交叉验证报告预测准确率
```
> 本包模块 8 的 `vanraden_kinship()` 可直接生成 G 矩阵（HapMap/PLINK 输入），
> 与 rrBLUP 的 `A.mat` 等价（VanRaden 2008 方法 1）。

## 6. 水稻专用注意点

| 事项 | 说明 |
|---|---|
| 多点试验 | 至少 2–3 个环境才能分解 G×E；单点数据 h2 会偏高（假精度） |
| 遗传力水平 | 水稻产量类性状 h2 通常 0.3–0.6；株型/粒型 0.6–0.9；演示数据（yield 0.32 / height 0.70）即为典型量级 |
| 光谱采集 | 晴天 10:00–14:00、固定高度与视场角；盐/旱胁迫下 NDVI 下降、NDRE 更敏感 |
| BLUP vs 原始均值 | GWAS 用 BLUP 可减少环境噪声、提高检出；但需在方法中说明 BLUP 已用于表型标准化 |

> 演示数据：`python scripts/utils/make_demo_phenomics.py --outdir demo_data/phenomics`
> （60 品系 × 3 环境 × 2 性状 + 5 个注入异常值；40 小区光谱含盐处理效应）。

---

## v3.2 新增脚本

### 功效与样本量分析 `05_power_analysis.R`（零依赖）

实验设计阶段先回答"要多少样本/能检出多大效应"：

```bash
Rscript scripts/phenomics/05_power_analysis.R --demo        # 三段全部演示
Rscript scripts/phenomics/05_power_analysis.R \
    --ttest --gwas --gs --outdir results/phenomics/05_power
```

- **两群表型比较**（如处理 vs 对照）：给定 Cohen's d 与每组 n 算功效，或反推 80% 功效
  所需 n（与 R 内置 `power.t.test` 对拍一致）
- **GWAS 位点检测**：非中心 χ²（`ncp = n·PVE/(1−PVE)`），给定样本量/MAF 算检出功效，
  或反推 80% 功效所需 PVE；默认 Bonferroni 5e-8
- **GS 训练群大小**：Daetwyler 式 `ρ² = Nh²/(Nh²+Me)`，给定 h²/Me 算准确率或反推所需 N
- 输出 `power_ttest.csv / power_gwas.csv / gs_accuracy_curve.csv + power_report.txt`；
  `--demo` 全部锚点对拍通过（示例：d=0.8 需每组 26 样本达 80% 功效）
