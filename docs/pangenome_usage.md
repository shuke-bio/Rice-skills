# 泛基因组分析模块使用指南（v3.0）

> Rice-skills v3.0 新增 **模块 11 · 泛基因组分析**：图形泛基因组构建 → 结构变异（SV）检测 →
> 存在/缺失变异（PAV）分析 → Pan-GWAS 关联。
> 上游衔接：模块 10（变异检测，提供参考比对基础）；下游衔接：模块 7（GWAS 引擎共用）。

## 0. 流程总览

```
多份组装 FASTA ──01_pangenome_graph.sh──▶ 图形泛基因组（GFA/VCF）
组装 FASTA ──02_sv_detect.sh──▶ 每样本 SV 表（SYRI）──▶ SV/PAV 矩阵
PAV 矩阵 ──03_pav_analysis.R──▶ 核心/可变/私有分类 + 结构 PCA + 性状关联
PAV/SV 矩阵 ──04_pan_gwas.R──▶ Pan-GWAS（Manhattan + QQ）
```

**适用场景**：多份高质量组装（≥10，推荐染色体级别）的泛基因组分析；材料少于 10 份或只有
重测序数据时，建议先做模块 7 的 SNP-GWAS，或只做 02+03 的 SV/PAV 检测。

**运行环境**：图构建与 SV 检测需 **Linux 服务器**（内存/磁盘消耗大）；03/04 为纯 R，本地可跑。

## 1. 依赖与工具获取

| 步骤 | 工具 | 安装/获取 |
|---|---|---|
| 图构建（快） | minigraph | `conda install -c bioconda minigraph` |
| 图构建（全） | minigraph-cactus | `conda install -c bioconda cactus` 或官方 Docker/Singularity |
| 图构建（无参考） | PGGB | `conda install -c bioconda pggb` |
| SV 检测 | minimap2 + SYRI + samtools | `conda install -c bioconda minimap2 syri samtools` |
| PAV 分析（第三方） | APAV | 官方仓库运行时获取：`git clone https://github.com/SJTU-CGM/APAV`（本包不打捆） |
| PAV/Pan-GWAS（本包自研） | 无额外依赖 | `03_pav_analysis.R` / `04_pan_gwas.R` 开箱即用 |

> 所有第三方工具均为**运行时由用户自行安装/获取**的独立程序（本包仅提供调用模板与文档），
> 请遵循各工具仓库的许可条款。

## 2. 步骤详解

### 2.1 图形泛基因组构建（`scripts/pangenome/01_pangenome_graph.sh`）

```bash
# 快速图（minigraph，适合先探索）
bash scripts/pangenome/01_pangenome_graph.sh \
    --ref IRGSP-1.0.fa --assemblies assemblies/ \
    --tool minigraph --outdir results/pangenome/01_graph --threads 32

# 标准图 + 可分型 VCF（minigraph-cactus，推荐）
bash scripts/pangenome/01_pangenome_graph.sh \
    --ref IRGSP-1.0.fa --assemblies assemblies/ \
    --tool minigraph-cactus --outdir results/pangenome/01_graph --threads 32
```

- 组装选择建议：3K-RGP 代表性材料 + 本实验室材料，覆盖籼/粳/野生稻多样性；
- minigraph-cactus 产出 `<name>.vcf.gz`（含基因型 SV）+ `.gbz`（图）——**SV 可以直接从 VCF 分型**，
  无需再做 2.2 的比对检测；
- 资源预估：minigraph ≈ 数小时 / 数十 GB 内存；minigraph-cactus ≈ 10–50 倍（准备大内存服务器）。

### 2.2 SV 检测（`scripts/pangenome/02_sv_detect.sh`）

```bash
# 单样本 / 批量（基于组装的 minimap2 + SYRI 路线）
bash scripts/pangenome/02_sv_detect.sh --ref IRGSP-1.0.fa \
    --queries assemblies/ --outdir results/pangenome/02_sv --threads 16
```

- 输出每样本 SV 表（`*.sv.syri.out`：INS/DEL/INV/TRA/CPG/CPL 等类型 + 坐标 + 大小）与汇总；
- `--preset asm10` 适合种内组装（约 1–5% 分歧）；精细抛光组装可用 asm5；
- **替代路线**：若已用 minigraph-cactus 建图，直接使用其 VCF 中的 SV 基因型，跳过本步骤。

### 2.3 PAV 分析（`scripts/pangenome/03_pav_analysis.R`）

```bash
Rscript scripts/pangenome/03_pav_analysis.R \
    --pav demo_data/pangenome/pav_matrix.csv \
    --meta demo_data/pangenome/marker_info.csv \
    --pheno demo_data/pangenome/phenotype.csv --trait grain_yield \
    --outdir results/pangenome/03_pav
```

**输入格式**（marker × sample 矩阵，第一列 marker_id）：

| 来源 | 说明 |
|---|---|
| APAV 输出 | PAV 表的覆盖度列（`--values cov --cov-threshold 0.5`） |
| minigraph-cactus VCF | `bcftools query -f '[%GT\t]\n'` 导出后转 0/1 |
| 自有数据 | 任何"基因/区域 × 材料"的存在矩阵（0/1 或覆盖度） |

**功能**：二值化 → 核心/可变/私有/缺失分类 → 统计摘要 → 频率分布图 → PAV-PCA 结构图 →
标记热图 → （可选）性状关联（连续性状 Wilcoxon / 分类性状 Fisher + BH 校正 + Manhattan）。

### 2.4 Pan-GWAS（`scripts/pangenome/04_pan_gwas.R`）

```bash
Rscript scripts/pangenome/04_pan_gwas.R \
    --pav demo_data/pangenome/pav_matrix.csv \
    --snp-info demo_data/pangenome/marker_info.csv \
    --phenotype demo_data/pangenome/phenotype.csv \
    --trait grain_yield --outdir results/pangenome/04_pan_gwas
```

- 与 SNP-GWAS 同一引擎（OLS/GLM + PCA 协变量），标记按存在频率（默认 0.05–0.95）与缺失率过滤；
- 输出全量结果表 + Bonferroni 显著表 + Manhattan/QQ 出版图；
- SV 剂量数据（0/1/2）同样支持（`--values auto` 保留原始编码）。

## 3. 输出结构

```
results/pangenome/
├── 01_graph/     # GFA / gbz / VCF + graph_log.txt
├── 02_sv/        # 每样本 SV 表 + sv_all.tsv + BAM
├── 03_pav/       # pav_class.csv + pav_stats.txt + 频率/PCA/热图/关联
└── 04_pan_gwas/  # pan_gwas_results.csv + 显著表 + Manhattan/QQ
```

## 4. 水稻专用注意点

| 事项 | 建议 |
|---|---|
| 参考骨干 | 以 IRGSP-1.0 为参考构建图，保证坐标与下游注释（RAP-DB）一致 |
| 材料选择 | 参考 3K-RGP（3,000+ 份重测序）与已发表染色体级组装（如 16 份代表性籼/粳/野生）；避免近缘重复材料 |
| 核心/可变基因 | 栽培稻核心基因占比约 60–80%（研究间差异大）；稀有/私有基因常富集于抗病基因家族（NBS-LRR） |
| SV 验证 | Pan-GWAS 命中位点建议用 PCR/长读长数据验证；已知案例：`Sub1A`、`OsGA20ox2`（sd1）等 |
| PAV 与表型 | 注意区分"存在/缺失"与"结构变化"——PAV 关联命中可能是连锁而非因果 |

## 5. 与其它模块的衔接

| 下游任务 | 衔接方式 |
|---|---|
| Pan-GWAS → 候选基因 | 显著标记（`pan_gwas_significant.csv`）× 标记区间基因（对照 `marker_info.csv` 坐标 + GFF） |
| PAV 分类 → 候选基因 | `pav_class.csv` 中 private/稀有标记的基因列表 |
| 与 SNP-GWAS 互补 | 同一性状分别跑模块 7 与 Pan-GWAS，比较命中区域（SV 常捕获 SNP 遗漏的因果变异） |
| 单倍型/表达验证 | 候选基因 → 模块 1（转录组）/模块 3（联合分析）验证表达与代谢关联 |
| 数据库复核 | RiceVarMap2 / RFGB / Rice SNP-Seek 在线单倍型分析（见 `docs/rice_databases.md`） |
