# GAPIT GWAS 输入格式与模型选择——方法笔记

> 依据 GAPIT 官方用户手册、GitHub 仓库文档与各模型原始论文**用自己的语言重写**的方法笔记；
> 未复制任何第三方教程正文或源码。GAPIT 为 **GPL-2.0+**（DESCRIPTION 核验，2026-10），
> **运行时安装、不打包**（见 §7）。

## 1. 适用场景

- 需要用**社区经典 R 路线**复现或对比 GWAS 结果时：GAPIT 是作物 GWAS 领域引用最广的
  工具之一（Lipka et al. 2012 Bioinformatics 起，历经 v2/v3/v4 三代迭代），
  审稿人对其模型族（GLM/MLM/FarmCPU/BLINK）与输出格式有公认认知。
- 一份基因型要**同时跑多个模型**比较时：GAPIT 的 `Multiple_analysis=TRUE` 一次传多个
  model，统一输出各模型的 Manhattan/QQ 图与结果表，适合模型敏感性分析。
- 除了 GWAS 还想做**基因组选择（GS）**：同包内 gBLUP/cBLUP/sBLUP/GAGBLUP
  （Wang & Zhang 2026 MBE, GAPIT v4 主题即"GWAS 与 GP 整合"）。
- 数据是 **Numeric 矩阵 / HapMap / PLINK / VCF** 任意一种时的前处理决策（见 §3）。

## 2. 方法要点（语言重写）

**问题本质**：GWAS 的假阳性主要由群体结构和亲缘关系驱动——亲缘近的个体性状值天然接近，
与位点等位频率耦合后产生"假关联"；同时严格校正又会损失功效。各模型的差别本质是
"固定效应（Q/结构）+ 随机效应（K/亲缘）怎么进模型、以及算不算得动"。

GAPIT 实现的 GWAS 模型族（v3/v4 口径）：

| 模型 | 一句话原理 | 计算量 | 何时用 |
|---|---|---|---|
| GLM | 只含固定效应（Q 矩阵或标记协变量） | 低 | 小群体、结构弱、快速筛查 |
| MLM（P3D） | Q+K 混合线性模型，用 P3D 压缩方差组分估计（Zhang et al. 2010 Nat Genet） | 中 | **默认主力**：群体结构+亲缘同时校正 |
| CMLM | 对基因组分簇压缩 K 矩阵，把大群体拆成组内近似（Tang et al. 2016） | 中 | 大群体，MLM 太慢时的加速版 |
| SUPER | 逐位点迭代聚类，亲缘远者做回归、近者做聚类层级（Wang et al. 2014 PLoS One） | 高 | 结构复杂的自然群体 |
| MLMM | 逐步纳入多个背景 SNP 作为随机效应，替代单一 K（Segura et al. 2012） | 高 | 多位点遗传架构、MLM 过度校正时 |
| FarmCPU | 固定效应与随机效应**交替迭代**：先 GLM 筛协变量 SNP，再 MLM，循环直至收敛（Liu et al. 2016 PLoS Genet） | 中高 | **功效/效率均衡推荐**，大队列首选之一 |
| BLINK | 百万级个体与标记，压缩 K + 迭代（Huang et al. 2018 Gigascience） | 中高 | 超大队列（人/玉米级规模） |

GS 模型：gBLUP（默认）、cBLUP/sBLUP（Bayesian 压缩，Wang et al. 2018 Heredity）、
GAGBLUP——按预测目标选，GW/GS 共用一套 GD/GM 输入。

**实践建议（多模型并行的意义）**：不同模型对结构和亲缘的假设不同，单一模型的显著峰
可能是模型假象。发表级做法是 GLM + MLM + FarmCPU（或 BLINK）三者交集为主、
差集为候选补充，配合 QQ 图检验膨胀程度（见 §6）。

## 3. 输入格式详解

GAPIT 读三类信息：**Y**（表型）、**GD**（基因型数值）、**GM**（SNP 物理信息）。
表型与基因型的 taxa 顺序可以不同——GAPIT 自动取交集并对齐（官方口径），但仍建议
上游用同一 ID 体系，避免交集意外缩小。

| 格式 | GD 形态 | 适用 | 转换路径 |
|---|---|---|---|
| **Numeric（默认）** | 0/1/2 数值矩阵（多位点需编码为 0/1/2） | 本包 GWAS 主链的输出形态 | 本包 `01_plink_qc.sh` + `02_plink_to_matrix.R` 直接产 numeric |
| HapMap | HapMap 文本（单字符基因型符号） | TASSEL/VCF 生态常用交换格式 | 见 §4 |
| PLINK | .ped/.map 或二元 | Plink 生态 | `--file` 前缀族 |
| VCF | .vcf | 直接吃变异调用结果 | `--file` + 后缀选择 |

关键文件要求（官方手册口径重写）：

- **Y**：第一列 taxa，其余每列一个性状；缺失用 `NA`；多性状文件一次全跑。
- **GD（Numeric）**：行=SNP，列=taxa；编码 0/1/2（大群体可改 0/1）；首列为 SNP 名称。
- **GM**：行=SNP，第一列 SNP 名，后续列染色体、位置（bp）；与 GD 行序一致。
- **K（亲缘矩阵）**：可由 GAPIT 内部计算（VanRaden 等），也可外部提供（.csv）；
  需注意 K 的行序从第一行就开始是数据（无表头），外部 K 的 taxa 集合会被自动对齐。
- **CV（协变量）**：群体结构信息，通常取 PCA 前若干主成分（Q 矩阵）；行=taxa。

> HapMap 的前 11 列是元数据列（rs#、alleles、chrom、pos、strand、assembly#、
> center、protLSID、assayLSID、panelLSID、QCcode），第 12 列起才是 taxa 基因型；
> 用 R 读入时通常跳过前 11 列（`read.table` 后取 12:ncol）。Numeric 格式的缺失编码
> 用 `NA`；HapMap 的缺失字符是 `N`。两种格式混用是新手最常见的报错来源。

## 4. VCF → HapMap 前处理（Numeric 之外的第二条路）

当上游是变异检测流水线（GATK/bcftools）产出 VCF、而下游想用 GAPIT 或 TASSEL 时，
社区标准路径是 vcftools 过滤 + TASSEL 格式转换（TASSEL5 的 `run_pipeline.pl`）。
完整流程（各步独立可跳过）：

1. **定义 VCF 路径变量**（`-max-missing` 控制位点缺失率）。
2. **vcftools 质控过滤**（常见口径）：`--max-missing 0.8`（位点）、
   `--maf 0.05`（次等位基因频率）、`--recode` 输出过滤后 VCF。
3. **vcftools 转 PLINK 格式**（`--plink`），按需用 PLINK 进一步导出。
4. **TASSEL VCF→HapMap**：`run_pipeline.pl -Xmx20g -vcf input.vcf
   -export -exportType Hapmap`；TASSEL 要求样本名不含特殊字符，vcf 样本名需先规整。
5. **删除 HapMap 重复基因型列**：`sed` 去掉 `\-` 前的重复符号（HapMap 编码惯例）。

本包链路的对应关系：`scripts/variants/` 若有 VCF 产出的，`bcftools` 即可完成
Numeric 化，不需要再绕 HapMap；只有当对方工具（GAPIT/TASSEL 生态）强制要
HapMap 时才走 §4。**决策原则：链内保持 Numeric，出链才转 HapMap。**

## 5. 显著性判定与 QQ 解读

- **阈值双口径**：Bonferroni（0.05/标记数，严格、控制家系错误率）与 FDR
  （BH 校正后的观察/期望 p 值阈值）。本包 `03_association.R` 默认 Bonferroni，
  `--p-threshold` 可覆盖；建议两个口径都算，报告时给出判定依据。
- **QQ 图是模型体检**：观测 p 对期望 p 作图。早期偏上翘（bulge）= 群体结构/
  亲缘未校正干净，假阳性温床；中段贴合、仅尾部翘起 = 有真实信号。GAPIT 默认输出
  的 QQ 有时只显示部分 p 段，做模型诊断时建议用本包 `05_plot_gwas.R` 产全图对比。
- **Manhattan 读峰**：超阈值连续峰取 lead SNP（p 最小者）为峰代表，报告时附物理
  位置与候选基因交集（本包模块 7 后半链 `04_candidate_loci.R` 负责）。

## 6. 输出文件解读（GAPIT 工作目录）

GAPIT 把结果写进工作目录，常用产物（官方手册口径）：

| 文件/图 | 内容 |
|---|---|
| `GWAS.Results` 表 | 每 SNP 的 p 值等统计量（含 `P.value` 列） |
| Allelic effect 估计 | 加性效应值 |
| D 值与标准误 | 比值的方差组分诊断 |
| Pred | 表型预测值（GS 相关） |
| ROC | 有效/发现率评估表 |
| Manhattan / QQ | 每模型一套，多模型时按 model 后缀区分 |
| PCA 图 | 群体结构；`Inter.Plot=TRUE` 时给交互式 |
| 亲缘热图 / NJ 树 | K 矩阵与 `NJtree.group` 分组 |

本包对 GAPIT 输出的接法：`03_association.R --engine gapit` 会把
`GWAS` 结果表统一转成本包标准 `association_results.csv`（列：chrom/pos/p/…），
后续 `04_candidate_loci.R` / `05_plot_gwas.R` 无缝衔接，不依赖 GAPIT 的文件布局。

## 7. 与 Rice-skills 的衔接

| 场景 | 命令/路径 |
|---|---|
| 跑 GAPIT 引擎（默认 MLM） | `Rscript scripts/gwas/03_association.R --engine gapit …` |
| 指定模型（GLM/MLM/CMLM/SUPER/FarmCPU/MLMM/BLINK） | 加 `--model FarmCPU`（engine=gapit 时生效） |
| 多模型并行比较 | 对同一输入多次运行，`--outdir` 分开；或用 GAPIT 原生 `Multiple_analysis` |
| 显著位点→候选基因 | `04_candidate_loci.R`（下游一致） |
| 出版级图 | `05_plot_gwas.R`（Manhattan/QQ/总览，3 格式） |

**许可与安装**：GAPIT 为 **GPL-2.0+**（2026-10 核验 DESCRIPTION，无独立 LICENSE
文件；仓库未在任何位置声明更宽许可）。本包 MIT 主体**不打包、不内嵌其源码**；
用户自行安装属使用而非再分发：

```r
install.packages("remotes")
remotes::install_github("jiabowang/GAPIT")
```

引用口径：v4 → Wang & Zhang 2026 MBE (msag107)；v3 → Wang & Zhang 2021 GPB；
FarmCPU → Liu et al. 2016 PLoS Genet；BLINK → Huang et al. 2018 Gigascience；
SUPER → Wang et al. 2014 PLoS One；P3D → Zhang et al. 2010 Nat Genet；
MLMM → Segura et al. 2012 Nat Genet；
cBLUP/sBLUP → Wang et al. 2018 Heredity；GAPIT 初代 → Lipka et al. 2012 Bioinformatics。
