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
生成 300 样本 × 3000 SNP 的合成数据（含亚群结构、LD 块、4 个因果位点与注释）。

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
| `gapit` | MLM/BLINK/FarmCPU | GitHub 安装 GAPIT | 社区经典路线 |

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

所有目录含 `run_info.txt`（参数 + 包版本 + 种子），图均为 PDF+SVG+TIFF（300 dpi）。

## 4. 与转录组/代谢组模块的衔接

- GWAS 候选区间（`gwas_loci.csv`）可与转录组差异表达、WGCNA 模块交叉比对，
  得到"变异 × 表达"的高置信候选（eQTL 思路的手工版）；
- 群体结构 PCA 得分可直接作为主包分析的分层变量；
- 候选基因列表可回到 `06_tf_analysis.R` / `03_enrichment.R` 做通路层面的解读。
