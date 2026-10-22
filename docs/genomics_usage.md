# 基因组学模块使用指南（GWAS + 群体遗传）

> Rice-skills v2.0 新增两大模块：**模块 7 · GWAS 全链**与**模块 8 · 群体遗传**。
> 全部脚本与主包共用同一套 Nature 主题（`config/nature_palette.json`）、工具层
> （`scripts/utils/plot_theme.R`、`genomics_utils.R`）与溯源规范（`run_info.txt`）。

## 0. 输入数据格式（两种路线）

**路线 A · 剂量矩阵（推荐，零外部依赖，全流程可跑）**

| 文件 | 格式 |
|---|---|
| `genotype.csv` | 第一列 `snp_id`，其余列=样本；剂量 0/1/2，缺失留空 |
| `snp_info.csv` | `snp_id, chr, pos`（pos 单位 bp） |
| `phenotype.csv` | 第一列 `sample_id`，其余列为性状；做 Fst 时需一列分组（如 `group`） |
| `annotation.gff3`（可选） | 候选基因注释（MSU/RAP-DB 官方 GFF3，或跳过） |

**路线 B · VCF（标准生信路线，需要 PLINK）**

```bash
bash scripts/gwas/01_plink_qc.sh --vcf genotypes.vcf.gz --outdir results/gwas/01_qc
Rscript scripts/gwas/02_plink_to_matrix.R \
    --raw results/gwas/01_qc/rice.raw --map results/gwas/01_qc/rice.map \
    --prune results/gwas/01_qc/rice.prune.in --outdir results/gwas/01_qc/matrix
```

未装 PLINK 时脚本会给出安装指引（conda 或官方二进制）；路线 A 不需要任何二进制。

**快速体验**：`python scripts/utils/make_demo_genomics.py --outdir demo_data/genomics`
生成 300 样本 × 3000 SNP 的示例数据（含亚群结构、LD 块、4 个因果位点与注释）。

## 1. 模块 7 · GWAS 全链

```bash
# ① 关联分析（三引擎可选：auto / simple / rmvp / gapit）
Rscript scripts/gwas/03_association.R \
    --genotype demo_data/genomics/genotype.csv \
    --snp-info demo_data/genomics/snp_info.csv \
    --phenotype demo_data/genomics/phenotype.csv \
    --trait salt_tolerance_index --npc 5 --outdir results/gwas/03_association

# ② 显著位点 → 候选区间 → 候选基因
Rscript scripts/gwas/04_candidate_loci.R \
    --assoc results/gwas/03_association/association_results.csv \
    --snp-info demo_data/genomics/snp_info.csv \
    --gff demo_data/genomics/annotation.gff3 \
    --flank 100000 --outdir results/gwas/04_candidates

# ③ 出版级图（Manhattan + QQ + 组合总览，各 3 格式）
Rscript scripts/gwas/05_plot_gwas.R \
    --assoc results/gwas/03_association/association_results.csv \
    --snp-info demo_data/genomics/snp_info.csv --outdir results/gwas/05_figures
```

**引擎说明**：

| 引擎 | 模型 | 依赖 | 适用 |
|---|---|---|---|
| `simple` | 逐 SNP OLS（GLM 族）+ PCA 协变量 | 无（内置） | 探索性分析、小数据、快速验证 |
| `rmvp` | GLM/MLM/FarmCPU | `install.packages("rMVP")`（Apache-2.0） | **发表级主推**（内存高效、多模型） |
| `gapit` | GLM/MLM/CMLM/SUPER/MLMM/FarmCPU/Blink（`--model` 选择，默认 MLM） | GitHub 安装 GAPIT（GPL-2.0+，运行时安装不打包） | 社区经典路线 |

显著性阈值：默认 Bonferroni 0.05/n；`--p-threshold` 可覆盖。

## 2. 模块 8 · 群体遗传

```bash
# ① PCA 群体结构 + VanRaden 亲缘矩阵
Rscript scripts/popgen/01_pca_structure.R \
    --genotype demo_data/genomics/genotype.csv --snp-info demo_data/genomics/snp_info.csv \
    --phenotype demo_data/genomics/phenotype.csv --group-col group \
    --outdir results/popgen/01_pca

# ② LD 衰减（自动标出 r²=0.2 的半衰减距离，指导候选窗口大小）
Rscript scripts/popgen/02_ld_decay.R \
    --genotype demo_data/genomics/genotype.csv --snp-info demo_data/genomics/snp_info.csv \
    --max-dist 500000 --outdir results/popgen/02_ld

# ③ 选择信号扫描：Hudson Fst（两组）+ 核苷酸多样性 π
Rscript scripts/popgen/03_selection_stats.R \
    --genotype demo_data/genomics/genotype.csv --snp-info demo_data/genomics/snp_info.csv \
    --phenotype demo_data/genomics/phenotype.csv --group-col group \
    --group1 indica --group2 japonica --window 100000 --outdir results/popgen/03_selection

# ④ 系统发育树（IBS 距离 + NJ/UPGMA，导出 Newick）
Rscript scripts/popgen/06_phylogeny.R \
    --genotype demo_data/genomics/genotype.csv --snp-info demo_data/genomics/snp_info.csv \
    --phenotype demo_data/genomics/phenotype.csv --group-col group \
    --outdir results/popgen/06_phylogeny

# ⑤ ADMIXTURE（Linux/WSL 上运行；无二进制则按提示跳过）
bash scripts/popgen/04_admixture.sh --bed plink/rice.bed --k-min 2 --k-max 8 \
    --outdir results/popgen/04_admixture
# ⑥ 祖先成分堆积图 + CV 误差曲线（对任意 .Q 矩阵均可绘图）
Rscript scripts/popgen/05_plot_admixture.R \
    --q-dir results/popgen/04_admixture --cv results/popgen/04_admixture/cv_errors.txt \
    --outdir results/popgen/05_admixture_figures
```

**依赖与降级**：核心计算全部为自实现（base R + data.table/ape），
无 ADMIXTURE/PLINK/vcftools 二进制也能完成 PCA、LD、Fst、π、系统树与分析级 GWAS；
ADMIXTURE 与 PLINK 仅在上游（.Q 矩阵、VCF 转换）需要，脚本均提供安装指引与替代方案。

## 3. 输出说明

| 目录 | 关键产物 |
|---|---|
| `03_association/` | `association_results.csv`（全量）、`association_significant.csv`、`pca_scores.csv` |
| `04_candidates/` | `gwas_loci.csv`（lead SNP/区间）、`candidate_genes.csv`（与 GFF 注释重叠） |
| `05_figures/` | `manhattan.*`、`qq.*`、`fig_gwas_overview.*` |
| `popgen/01_pca/` | `pca_scores.csv`、`kinship_matrix.csv`、`pca_scatter.*`、`fig_structure_overview.*` |
| `popgen/02_ld/` | `ld_decay.csv`、`ld_decay_curve.*` |
| `popgen/03_selection/` | `fst_per_snp.csv`、`fst_windows.csv`、`pi_windows.csv`、`fst_genome_scatter.*` |
| `popgen/05_admix/` | `admixture_k<K>.*`、`cv_curve.*` |
| `popgen/06_phylogeny/` | `tree.nwk`、`phylogeny.*` |

所有目录含 `run_info.txt`（参数 + 包版本 + 种子），图均为 PDF+TIFF+PNG（300 dpi）。

## 4. 与转录组/代谢组模块的衔接

- GWAS 候选区间（`gwas_loci.csv`）可与转录组差异表达、WGCNA 模块交叉比对，
  得到"变异 × 表达"的高置信候选（eQTL 思路的手工版）；
- 群体结构 PCA 得分可直接作为主包分析的分层变量；
- 候选基因列表可回到 `06_tf_analysis.R` / `03_enrichment.R` 做通路层面的解读。

---

## v3.2 新增脚本

### 2D-GWAS 成对性状关联 `08_2d_association.R`（零依赖）

找"只在某个背景性状的极端群体中才显现"的关联位点（标准 GWAS 会漏掉的
条件特异/上位性效应；方法思路整理自 PATOWAS，见 `docs/references/patowas-2d-association.md`）：

```bash
Rscript scripts/gwas/08_2d_association.R \
    --geno demo_data/genomics/genotype_matrix.csv \
    --pheno demo_data/genomics/phenotypes.csv \
    --traits traitA,traitB --tail 0.1 --p-cutoff 0.01 \
    --outdir results/gwas/08_2d
```

- 逻辑：按性状 A 取上下极端组（`--tail`，默认 10%）→ 组内分别做 B ~ SNP 检验 →
  只在单侧极端组显著的位点记为"A 条件特异"，方向 +A / -A
- 输出 `2d_results.csv`（每 SNP 双侧 p 值、Δlog10p、方向）+ `2d_scatter.png`（高 × 低组散点）
- `--demo` 自检：植入的上位性 SNP 以 +A 方向正确回收；命中后建议回到标准 GWAS
  （模块 7）与共定位（02/03）在条件性状上做验证

## v3.3 新增脚本：3VmrMLM 家族（多位点 / QEI / 上位性 / case-control / 多亲本）

> 方法笔记：`docs/references/3vmrmlm_multilocus_qtn.md`、
> `fastbicmrmlm_case_control.md`、`qei_environment_interaction.md`、
> `epistasis_qqi.md`。全部引擎为 **GPL-3 或未声明许可，运行时安装，不打包**。
> 实测（2026-10-19，R 4.6.0 / Windows）：IIIVmrMLM、mppQTL 可正常加载；
> Fast3VmrMLM 二进制按 R 4.3 构建，R ≥ 4.4 下 DLL 加载失败（脚本自动退回
> IIIVmrMLM）；FastBiCmrMLM 仅提供 Linux 版。

### 09 多位点 QTN / QEI / 上位性 `09_multilocus_qtn.R`

与 `03_association.R` 的单位点扫描互补：多位点联合模型 + 全多基因背景控制，
对微效/稀有位点更敏感；支持 QTN×环境与 QTN×QTN 互作。

```bash
# 单环境多位点（直接吃 PLINK 前缀，接 01_plink_qc.sh 的产出）
Rscript scripts/gwas/09_multilocus_qtn.R \
    --genotype results/gwas/01_qc/qc_geno \
    --phenotype data/phenotype.csv --trait 1:3 \
    --pop-structure results/gwas/01_qc/pca.csv --pop-type PC \
    --outdir results/gwas/09_multilocus

# 多环境（QEI；性状×环境宽表，trait 选性状、n-en 给环境数）
Rscript scripts/gwas/09_multilocus_qtn.R --method multi_env \
    --genotype ... --phenotype data/phenotype_wide.csv --n-en 2,2

# 上位性（QQI）
Rscript scripts/gwas/09_multilocus_qtn.R --method epistasis \
    --genotype ... --phenotype ... --trait 1

# 自检（三模式都可）
Rscript scripts/gwas/09_multilocus_qtn.R --demo
Rscript scripts/gwas/09_multilocus_qtn.R --demo --method multi_env
Rscript scripts/gwas/09_multilocus_qtn.R --demo --method epistasis
```

- 引擎：`--engine auto`（默认，Fast3VmrMLM 可用则用，否则 IIIVmrMLM）/
  `IIIVmrMLM` / `Fast3VmrMLM`
- 关键参数：`--svpal`（扫描阈值，n<4000 用 0.01，n≥4000 用 1e-5）、
  `--svrad`（去共线性窗口 bp）、`--svmlod`（提示性 LOD）
- 输出 `qtn_results.csv`（trait/marker/chr/pos/lod/additive/dominance/pve/p）与
  `association_results.csv`（可直接接 `05_plot_gwas.R`）
- demo 实测：单环境 15 s、multi_env 22 s、epistasis 约 6 min（200×400）

### 10 二分类 / case-control `10_case_control_qtn.R`（Linux/WSL）

```bash
Rscript scripts/gwas/10_case_control_qtn.R \
    --genotype results/gwas/01_qc/qc_geno \
    --phenotype data/binary_trait.csv --trait blast_resistance \
    --method auto --outdir results/gwas/10_case_control
```

- 0/1 表型；连续表型用 `--cut-to-binary <阈值>` 二值化
- 算法按样本量自动分派：<1k 主算法 / 1k–20k `-Time` / ≥20k `-RAM`（省内存）
- **Windows 不支持**：脚本在非 Linux 主机输出 SKIP + 替代路径（数量性状走 09）
- 实测：Windows 端 skip 分支验证通过；Linux 实跑请在 Linux/WSL 环境执行 `--demo`

### 11 多亲本群体 QTL `11_multiparent_qtl.R`

```bash
Rscript scripts/gwas/11_multiparent_qtl.R --population magic \
    --genotype 4way_magic_genotype.csv --phenotype 4way_magic_phenotype.csv \
    --map 4way_magic_map.csv --cross-info 4way_magic_cross.csv \
    --crosstype riself4 --map-function haldane --cores 4 \
    --outdir results/gwas/11_multiparent

Rscript scripts/gwas/11_multiparent_qtl.R --demo   # 150 系 × 240 标记合成 MAGIC
```

- 适用 MAGIC/NAM/ROAM；`--crosstype` 如 `riself4/8/16`、`risib4/8/16`、`magic19`、`nam5/26`、`roam14`
- 基因型 CSV 行=标记、列=`founder_*` + 各系；`--cross-info` 每行必须是
  1..n_founder 的**排列**（qtl2 校验）；等位字母只认 A/T/C/G
- 输出 `qtl_results.csv`（显著 QTL）+ `gws_scan.csv`（全扫描）+ 每性状×环境 xlsx
- demo 实测：32 s，植入 QTL 精确命中 D2M30
