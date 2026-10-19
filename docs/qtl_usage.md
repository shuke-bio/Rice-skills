# QTL 映射模块使用指南（eQTL / mQTL / pQTL，v3.0）

> Rice-skills v3.0 新增 **模块 12 · QTL 映射**：把基因型变异与分子性状
> （基因表达 → eQTL；代谢物 → mQTL；蛋白 → pQTL）关联起来，支持 cis / trans 模式。
> 上游衔接：模块 10（变异检测）；下游衔接：模块 7（GWAS）、模块 11（泛基因组）。

## 0. 流程总览

```
基因型 (VCF→dosage) ──┐
                       ├── 01_qtl_mapping.R ──▶ cis/trans QTL 表 + Manhattan + 距离图
分子性状矩阵（表达/代谢/蛋白）──┘
```

**科学问题**：某变异是否调控某个基因的表达/代谢物积累？cis-eQTL 提供"变异→表达"的
因果证据链，是连接 GWAS 位点与候选基因的关键工具。

## 1. 输入格式

| 文件 | 格式 | 来源 |
|---|---|---|
| `--genotype` | SNP × 样本剂量矩阵（0/1/2） | 模块 10 输出的 VCF → PLINK → 模块 7 的 `02_plink_to_matrix.R` |
| `--snp-info` | snp_id, chr, pos | 同上 |
| `--features` | 特征 × 样本矩阵（**需预先标准化**：表达用 VST/TPM-log，代谢用 log+Pareto） | 模块 1 的 `02_deseq2.R` 输出的 VST 矩阵 / 模块 2 的预处理矩阵 |
| `--feature-info` | feature_id, chr, pos（cis 模式必需） | 基因坐标（GFF 提取）/代谢物可不做 cis |

## 2. 使用方法

```bash
# cis-eQTL（推荐首选：每个基因只测自身 ±1 Mb 的 SNP）
Rscript scripts/qtl/01_qtl_mapping.R \
    --genotype results/gwas/01_qc/matrix/genotype_matrix.csv \
    --snp-info demo_data/genomics/snp_info.csv \
    --features results/transcriptomics/02_deseq2/vst_matrix.csv \
    --feature-info gene_coords.csv \
    --mode cis --window 1000000 --outdir results/qtl/01_eqtl

# 全基因组 trans 扫描（内存大，chunk 控制）
Rscript scripts/qtl/01_qtl_mapping.R ... --mode both --trans-p 1e-5 --chunk 100

# mQTL（代谢物；无位置信息时只能做 trans）
Rscript scripts/qtl/01_qtl_mapping.R \
    --genotype ... --snp-info ... \
    --features results/metabolomics/01_preprocess/metabolites_normalized.csv \
    --mode trans --outdir results/qtl/01_mqtl
```

**关键参数**：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--mode` | cis | cis / trans / both |
| `--window` | 1e6 | cis 窗口（bp）；水稻 LD 衰减快（~100–200 kb 半衰减距离），常用 100 kb–1 Mb |
| `--min-maf` / `--max-missing` | 0.05 / 0.1 | 基因型 QC |
| `--trans-p` | 1e-5 | trans 每对保留阈值（BH 校正按全部测试数保守计算） |
| `--output-p` | 0.01 | cis 结果写入阈值（控制大文件） |
| `--chunk` | 100 | trans 分块大小（内存不够就调小） |

## 3. 输出

```
results/qtl/01_eqtl/
├── qtl_results.csv          # cis 对（p ≤ output-p）+ trans 保留对：beta/t/p/fdr/距离
├── qtl_significant.csv      # FDR < 0.05
├── qtl_top_per_feature.csv  # 每个基因/代谢物的最佳命中（cis 优先）
├── qtl_summary.txt          # 测试计数与概况
├── manhattan_qtl.*          # 每 SNP 最强信号的 Manhattan
├── cis_distance_plot.*      # cis 信号 vs 距离（经典 eQTL 图：峰常在基因附近）
└── run_info.txt
```

## 4. 水稻专用注意点

| 事项 | 说明 |
|---|---|
| 群体结构 | 结构化群体（籼/粳混合）中 eQTL 假阳性上升：先用亚群内分析，或在基因型 PCA 基础上扩展脚本加协变量 |
| LD 结构 | 水稻 LD 衰减快 → cis 窗口可用 100 kb–500 kb；显著 SNP 落在基因近端（<10 kb）是强 cis 证据 |
| 表达矩阵 | 必须取 normalized（如 VST）；不要直接输入 raw counts |
| 因果推断 | eQTL ≠ 因果：cis-eQTL 共定位（GWAS 位点与 eQTL 峰重叠）是更强证据；可用共定位方法（coloc）进一步分析 |
| 多组织 | 一个基因在不同组织可能有不同 eQTL；重要基因建议分组织分析 |

## 5. 与其它模块的衔接

| 任务 | 衔接方式 |
|---|---|
| GWAS 候选基因验证 | GWAS 显著位点（模块 7）× cis-eQTL（本模块）：位点调控哪个基因 = 优先候选 |
| 转录组差异基因 | 模块 1 的 DEG 中，哪些受 cis 变异调控（`qtl_top_per_feature.csv` 按基因过滤） |
| 代谢组关联 | mQTL 命中 × 模块 3 的基因-代谢物相关性网络 → 变异-基因-代谢物链 |
| Pan-GWAS | SV-eQTL 思路：把模块 11 的 PAV/SV 矩阵作为 `--genotype` 输入（二值编码） |
| 单倍型验证 | `qtl_significant.csv` 的 lead SNP → RiceVarMap2 在线单倍型-表型分析复核 |

> 引擎为自研实现（OLS + 保守 BH），无需 MatrixEQTL 等外部依赖；
> 大规模计算（上万基因 × 百万 SNP）建议用 MatrixEQTL/FastGWA 等专用工具，
> 本模块适合常规规模（≤ 数千特征）与教学/预分析。

---

## v3.2 新增脚本

### BSA / QTL-seq 完整链（05 → 04）

**1) 全基因组扫描 `05_bsa_scan.py`**（ΔSNP-index + G′ 滑窗，零外部依赖）

```bash
python scripts/qtl/05_bsa_scan.py \
    --vcf bulk.vcf --bulk-high HighPool --bulk-low LowPool \
    --window 2000000 --step 100000 --ci 0.99 \
    --outdir results/qtl/05_bsa_scan
```

- 输入：双混池 VCF（每池一个样本列，`AD` 字段承载池内等位计数；无 AD 时加 `--use-gt`）
- 统计：每 SNP `SNP-index = ALT/(REF+ALT)`；窗内均值差 ΔSNP-index 与 G′；显著性用
  正态模拟（`sd = sqrt(0.25*(1/n_high + 1/n_low))`）取 99% 分位作阈值
- 输出：`bsa_windows.csv`（每窗）、`bsa_peaks.csv`（候选区间）、`bsa_scan.png/pdf`（全基因组双轨图）
- `--demo` 自检：构建 3 条染色体并植入 Chr2 QTL → 精确回收、无 QTL 染色体零假阳性
- 与 QTLseqr 的关系：本脚本为轻量自研筛查；发表级分析可与 QTLseqr（R 包）互相印证，
  见 `docs/tools_navigation.md`

**2) 峰区间精细提取 `04_qtl_region_snps.py`**（把 05 找到的区间"放大看"）

```bash
python scripts/qtl/04_qtl_region_snps.py \
    --vcf filtered.vcf.gz --gff annotation.gff3.gz \
    --region Chr7:20000000-25000000 \
    --sample-a ParentA --sample-b ParentB --promoter 2000 \
    --outdir results/qtl/04_region_snps
```

- 功能：区间 VCF（双亲选择）→ 每 SNP 的 GT/深度/AD/突变频率/ΔAF/差异判定/携带归属 →
  GFF 五级注释（启动子/5'UTR/外显子/内含子/3'UTR/基因间区）→ 区间 VCF + **三表 Excel**
  （全部位点 / 仅差异位点 / 区间基因列表）
- 输出 Excel 需 `openpyxl`（缺失时自动降级三份 CSV）；`--demo` 内置示例 VCF+GFF 自检

**3) 贝叶斯共定位 `03_coloc_bayesian.R`**（02 筛查的定量升级）

```bash
Rscript scripts/qtl/03_coloc_bayesian.R \
    --gwas results/qtl/02_coloc/gwas_region.csv \
    --qtl  results/qtl/02_coloc/qtl_region.csv \
    --outdir results/qtl/03_coloc_bayes
```

- 自实现 coloc.abf（近似 Bayes 因子 + H0–H4 后验 + SNP 级共享因果后验），**零依赖**
- 区域 CSV 需 `snp_id, beta, se` 三列（两文件按 SNP 自动对齐）；`--demo` 共享因果
  PP4=0.994 / 不同因果 PP3>PP4
- 判定：PP4≥0.8 强共定位；0.5–0.8 提示；PP3>PP4 提示两个独立信号
