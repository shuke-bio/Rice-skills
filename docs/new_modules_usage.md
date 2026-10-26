# 新增模块使用说明（v3.4–v3.8）

> 本文覆盖 P0/P1/P2 阶段新增的脚本。所有脚本均支持 `--help`；带 `--demo` 的脚本可先用合成数据验证环境。
> 详细方法边界见各脚本头部注释；本文只列输入、输出与典型命令。

## 转录组

### `transcriptomics/00_fastq_counts.sh`
FASTQ → counts（fastp + HISAT2 + featureCounts），默认只打印计划，`--run` 执行。

```bash
bash scripts/transcriptomics/00_fastq_counts.sh \
  --r1 R1.fq.gz --r2 R2.fq.gz --reference genome.fa --annotation genes.gtf \
  --sample S01 --outdir results/transcriptomics/00_counts --run
```
输出：`S01_counts.txt`、BAM、fastp 报告。需要 fastp/HISAT2/samtools/featureCounts。

### `transcriptomics/08_tf_prioritize.py`
PURE-like TF 排序：共表达 + motif + ChIP/DAP 证据。

```bash
python scripts/transcriptomics/08_tf_prioritize.py --expression expr.csv \
  --deg deg.csv --tfs tfs.csv --motif motif.csv --binding binding.csv
```
输出：`tf_ranking.csv`（含各分量 z-score 与综合分）。

### `transcriptomics/09_dtu.py`
isoform 水平 DTU（卡方 + BH）。

```bash
python scripts/transcriptomics/09_dtu.py --counts iso_counts.csv \
  --info iso_info.csv --groups groups.csv
```
输出：`dtu_results.csv`、`isoform_proportions.csv`。仅支持两组。

### `transcriptomics/10_ssgsea.R`
ssGSEA/GSVA 通路活性；无 GSVA 时退化为 set mean z-score。

```bash
Rscript scripts/transcriptomics/10_ssgsea.R --expression expr.csv --gmt sets.gmt
```
输出：`ssgsea_scores.csv`。

## 代谢组

### `metabolomics/08_mgwas.R`
代谢物 GWAS（MLM）与 mQTL（标记回归 + LOD）。

```bash
Rscript scripts/metabolomics/08_mgwas.R --metabolites met.csv \
  --genotype geno.csv --snp-info snp_info.csv --mode gwas
```
输出：`mgwas_results.csv` 或 `mqtl_results.csv`、`significant.csv`。

### `metabolomics/09_annotate.py`
合并 SIRIUS 与 GNPS 注释。

```bash
python scripts/metabolomics/09_annotate.py --features features.csv \
  --sirius sirius.csv --gnps gnps.csv
```
输出：`annotated_features.csv`。只合并不运行注释工具。

## GWAS / QTL

### `gwas/12_kmer_gwas.py`
k-mer 剂量 GWAS（KMERIA 简化版）。

```bash
python scripts/gwas/12_kmer_gwas.py --kmer-matrix kmers.csv \
  --phenotype phe.csv --genome genome.fa
```
输出：`kmer_gwas_results.csv`、显著表（含基因组定位）。未实现图泛基因组回贴。

### `qtl/06_twas.R`
weighted-sum TWAS。

```bash
Rscript scripts/qtl/06_twas.R --eqtl eqtl_weights.csv --gwas gwas_summary.csv
```
输出：`twas_results.csv`。非完整 FUSION，不含 LD 参考。

### `qtl/07_mr.R`
两样本 MR：IVW、MR-Egger、加权中位数、异质性 Q。

```bash
Rscript scripts/qtl/07_mr.R --exposure eqtl.csv --outcome gwas.csv
```
输出：`mr_results.csv`、`mr_snps.csv`。无 MR-PRESSO/Steiger。

### `qtl/08_finemap_susie.R`
SuSiE（需 susieR）或近似 ABF fine-mapping。

```bash
Rscript scripts/qtl/08_finemap_susie.R --z z.csv --ld ld.csv --n 500 --method susie
```
输出：`finemap_results.csv`、`credible_sets.csv`。ABF 回退不考虑 LD。

## 变异 / 群体遗传

### `variants/05_consensus_call.sh`
bcftools + FreeBayes + GATK 共识变异（≥2 支持）。

```bash
bash scripts/variants/05_consensus_call.sh --reference ref.fa --reads1 R1.fq --reads2 R2.fq --sample S01
```
输出：`isec/0000.vcf`。默认只打印计划，`--run` 执行。

### `variants/06_sv_cnv.sh`
Sniffles2/cuteSV/Delly SV/CNV 计划。

```bash
bash scripts/variants/06_sv_cnv.sh --bam sample.bam --reference ref.fa --mode ont
```
输出：`tool_plan.txt`，`--run` 后生成 VCF。

### `popgen/08_selection_scan.sh`
selscan iHS/CLR + norm + Dsuite。

```bash
bash scripts/popgen/08_selection_scan.sh --vcf phased.vcf.gz --map chr.map \
  --sets sets.txt --tree tree.nwk
```
输出：iHS/CLR 结果、Dsuite D/Fbranch 计划。

## 泛基因组

### `pangenome/07_t2t_qc.sh`
T2T 组装 QC（N50/GC/gap；可选 BUSCO/Merqury/RepeatMasker/Liftoff）。

```bash
bash scripts/pangenome/07_t2t_qc.sh --assembly asm.fa --run-tools \
  --busco-lineage poales_odb10
```
输出：`qc_summary.tsv`、`tool_plan.txt`、工具日志。

### `pangenome/08_graph_genotype.sh`
vg Giraffe + PanGenie 图基因型计划/执行。

```bash
bash scripts/pangenome/08_graph_genotype.sh --graph rice.gfa --reads sample.fq --run
```
输出：`sample.vcf`（vg）。

### `pangenome/09_graph_visualize.sh`
odgi 图可视化。

```bash
bash scripts/pangenome/09_graph_visualize.sh --graph rice.gfa --run
```
输出：`graph.png`。

## 蛋白 / 表观 / 表型

### `proteomics/06_ptm_ppi.py`
PTM Fisher 富集 + STRING PPI 边/hub。

```bash
python scripts/proteomics/06_ptm_ppi.py --proteins ptm.csv --string string.tsv
```
输出：`ptm_enrichment.csv`、`ppi_edges.csv`、`ppi_hubs.csv`。

### `epigenomics/03_hic_analysis.sh`
Hi-C balance/zoomify/TAD/loops 计划。

```bash
bash scripts/epigenomics/03_hic_analysis.sh --matrix sample.cool --resolution 10000
```

### `epigenomics/04_scatac.py`
scATAC TF-IDF/LSI/KMeans。

```bash
python scripts/epigenomics/04_scatac.py --matrix scatac.csv --n-clusters 8
```

### `phenomics/09_dl_gs.py`
DL-GS 基线：Ridge/ElasticNet/MLP/CNN。

```bash
python scripts/phenomics/09_dl_gs.py --genotype geno.csv --phenotype phe.csv --method all
```

### `phenomics/10_temporal_gwas.py`
逐时间点 + 生长速率 GWAS。

```bash
python scripts/phenomics/10_temporal_gwas.py --genotype geno.csv --phenotype long.csv
```

## 单细胞 / 空间

### `single_cell/02_integration_annotation.py`
PCA/批次中心化 + KMeans + marker 注释。

```bash
python scripts/single_cell/02_integration_annotation.py --matrix expr.csv \
  --metadata meta.csv --markers markers.csv
```

### `single_cell/03_cell_communication.py`
配体-受体打分 + 置换检验。

```bash
python scripts/single_cell/03_cell_communication.py --matrix expr.csv --labels labels.csv
```

### `spatial/02_deconvolution.py`
NNLS 反卷积。

```bash
python scripts/spatial/02_deconvolution.py --sc sc.csv --labels labels.csv --spatial spots.csv
```

## 设计

### `design/04_offtarget.py`
SpCas9 全基因组脱靶（≤50 Mb 内置）。

```bash
python scripts/design/04_offtarget.py --sgrna sgrna.csv --genome region.fa
```

### `design/05_base_prime_edit.py`
CBE/ABE 与 pegRNA 设计。

```bash
python scripts/design/05_base_prime_edit.py --sequence ATG... --position 10 --edit C2T
```

## 工程工具

### `utils/export_workflow.py`
导出 Nextflow/Snakemake/Dockerfile。

### `utils/benchmark_suite.py`
演示回归：`--quick` 3 个，`--all` 6 个；CI 使用 `--all`。

### `utils/hpc_submit.py`
SLURM/PBS 提交脚本生成。

### `utils/build_skill_import.py`
生成 description ≤1024 的导入版 `SKILL.md`（不改主文件）。

### `integration/08_causal_network.py`
基因-代谢物关联网络，可合并外部 MR 边。

## v3.9.0 新增

### `transcriptomics/11_batch_timecourse.R`
ComBat（sva 可选）+ 时间序列模型：每基因 `lm(y ~ group*time)` vs `lm(y ~ group)`，F 检验 + BH。

```bash
Rscript scripts/transcriptomics/11_batch_timecourse.R --counts counts.csv --metadata meta.csv
```
输出：`timecourse_results.csv`（含 beta_group/beta_time/beta_interaction）、`batch_corrected_log2cpm.csv`。

### `qtl/06_twas.R --ld`
新增 LD 矩阵输入，TWAS 分母改为 `sqrt(w' R w)`；不给 `--ld` 时等价于独立 SNP。

```bash
Rscript scripts/qtl/06_twas.R --eqtl weights.csv --gwas gwas.csv --ld ld.csv
```

### `gwas/13_burden_test.py`
rare-variant gene-based burden：MAF ≤ 阈值、每基因 collapse、PCA 校正 OLS + BH。

```bash
python scripts/gwas/13_burden_test.py --genotype geno.csv --snp-info info.csv --phenotype phe.csv
```
输出：`burden_results.csv`。

### `qtl/10_qtl_epistasis.R`
成对 QTL 上位性：`y ~ g1*g2 + PCs`，交互项 p + BH。

```bash
Rscript scripts/qtl/10_qtl_epistasis.R --genotype geno.csv --phenotype phe.csv --pairs pairs.csv
```
输出：`epistasis_results.csv`。

## v3.10.0 新增

### `metabolomics/10_sirius.sh`
SIRIUS formula/structure 计划/执行；GNPS 生成上传说明（网页流程，无 CLI）。

```bash
bash scripts/metabolomics/10_sirius.sh --mzml data/mzml --run
```

### `metabolomics/11_gemma_mgwas.sh`
GEMMA kinship + AI-REML LMM。

```bash
bash scripts/metabolomics/11_gemma_mgwas.sh --bfile geno --pheno met1.txt --run
```

### `qtl/11_bsa_ed_sl.py`
BSA/QTL-seq 的 SNP index、ΔSNP-index、ED、SL（Fisher），支持窗口平滑。

```bash
python scripts/qtl/11_bsa_ed_sl.py --counts bsa_counts.csv --window 1000000
```

### `proteomics/07_alphafold.py`
AlphaFold DB 查询（pLDDT/PDB/PAE），可下载 PDB。

```bash
python scripts/proteomics/07_alphafold.py --uniprot ids.csv --download
```

### `phenomics/11_multitrait_gs.R`
逐性状 GBLUP + GEBV 相关（不是联合 MTM）。

```bash
Rscript scripts/phenomics/11_multitrait_gs.R --genotype geno.csv --phenotype phe.csv --traits trait1,trait2
```

## v3.11.0 新增

### `integration/09_mofa_interpret.R`
MOFA 因子解释：载荷平方和方差、top 特征、因子-性状相关。

### `integration/10_diablo_tune.R`
DIABLO 调参：生成 `diablo_plan.R`；有 mixOmics 且 `--run` 时执行。

### `single_cell/04_scvi_cellchat.py`
scVI 整合（有 scvi-tools 且 `--run` 时执行）+ CellChat R 计划。

### `spatial/03_cell2location.py`
cell2location/COMMOT 计划生成器。

### `variants/07_deepvariant.sh` / `08_phasing.sh`
DeepVariant 与 WhatsHap/SHAPEIT 计划/执行封装。

## v3.12.0 新增

### `qtl/12_fusion_weights.R`
FUSION-like cis-eQTL 权重（ridge 闭式解，无需 glmnet）；输出 `fusion_weights.csv` / `fusion_models.csv`。

### `qtl/13_cim.R`
复合区间作图：单标记选 cofactor，再逐点 `y ~ SNP + cofactors + PCs`，输出 LOD。

### `pangenome/10_phg.sh`
PHG（Practical Haplotype Graph）流程计划/执行。

### `epigenomics/05_archr_signac.R`
ArchR / Signac 可执行计划生成。

### `design/06_edit_efficiency.py`
碱基/先导编辑 guide 启发式排序（GC/homopolymer/窗口），并生成 BE-Hive/PRIDICT 复核计划。

### `gwas/14_mtag.py`
MTAG 输入校验与命令计划（SNP, Z, N + LD score）。

### `qtl/07_mr.R --presso`
有 MRPRESSO 包时运行 MR-PRESSO，否则写 `mr_presso_plan.txt`。
