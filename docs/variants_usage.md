# 变异检测与注释模块使用指南（v3.0）

> Rice-skills v3.0 新增 **模块 10 · 变异检测与注释**：从测序原始数据（FASTQ）到
> 分析就绪变异集（过滤 + 功能注释）的完整上游流程。
> 下游衔接：模块 7（GWAS）、模块 8（群体遗传）、模块 11（泛基因组）。

## 0. 流程总览

```
FASTQ ──01_bwa_align.sh──▶ BAM ──02_gatk_call.sh──▶ raw VCF
     ──03_vcf_filter.sh──▶ filtered VCF ──04_snpeff_annotate.sh──▶ 注释 VCF + 影响表
     ──▶ scripts/gwas/01_plink_qc.sh（GWAS / 群体遗传）
```

**适用场景**：自有重测序数据（WGS / GBS / 外显子捕获）；如果你从测序公司直接拿到
**VCF 或 PLINK 文件**，跳过 01–02，直接从 `03_vcf_filter.sh`（或模块 7 的 PLINK 质控）接入即可。

**运行环境**：比对与变异检测建议在 **Linux 服务器 / WSL2** 上运行（内存与磁盘开销大）；
Windows 原生环境可用 conda 安装工具，但性能与兼容性以 WSL2 为佳。

## 1. 依赖安装

```bash
conda install -c bioconda bwa samtools bcftools gatk4 snpeff
```

| 工具 | 用途 | 备注 |
|---|---|---|
| bwa | 短读长比对（BWA-MEM） | — |
| samtools | BAM 排序/去重/统计 | — |
| gatk4 | HaplotypeCaller 变异检测、硬过滤 | conda 包自带 Java 17 |
| bcftools | 轻量变异检测/过滤/统计（替代引擎） | 内存占用低 |
| snpeff | 变异功能注释 | 首次运行自动下载水稻数据库 |

> 本包**不捆绑**上述二进制文件——脚本运行时调用你自己安装的工具；未安装时脚本会
> 给出对应安装指引后退出。

## 2. 步骤详解

### 2.1 比对（`scripts/variants/01_bwa_align.sh`）

```bash
bash scripts/variants/01_bwa_align.sh \
    --fastq1 sample_R1.fq.gz --fastq2 sample_R2.fq.gz \
    --ref IRGSP-1.0_genome.fa --sample IR24 \
    --outdir results/variants/01_align --threads 16
```

- 流程：bwa index（仅第一次，自动复用）→ bwa mem → fixmate → markdup → index → QC；
- 自动添加 `@RG` 读组（GATK 必需）；输出 `<sample>.sorted.bam` + flagstat/idxstats/stats 质控报告；
- **参考基因组**：水稻推荐 IRGSP-1.0（Nipponbare，RAP-DB）或 MSU7；同一项目所有样本共用同一参考。

### 2.2 变异检测（`scripts/variants/02_gatk_call.sh`）

```bash
# 单样本
bash scripts/variants/02_gatk_call.sh --bam sample.sorted.bam \
    --ref IRGSP-1.0_genome.fa --sample IR24 --outdir results/variants/02_call

# 多样本联合分型（推荐：一起 call 减少批次差异）
bash scripts/variants/02_gatk_call.sh --bamdir results/variants/01_align \
    --ref IRGSP-1.0_genome.fa --outdir results/variants/02_call --java-mem 16g

# 轻量引擎（小数据/内存受限）
bash scripts/variants/02_gatk_call.sh --engine bcftools --bam s.sorted.bam \
    --ref IRGSP-1.0_genome.fa --outdir results/variants/02_call
```

- GATK 模式：逐样本 HaplotypeCaller（GVCF）→ GenomicsDBImport → GenotypeGVCFs，输出 `cohort.vcf.gz`；
- bcftools 模式：mpileup + call，适合小规模数据快速出结果；
- 大规模（>50 样本）自动按 `--intervals` 分级并行（按染色体拆分后分别 call 再合并）。

### 2.3 过滤（`scripts/variants/03_vcf_filter.sh`）

```bash
bash scripts/variants/03_vcf_filter.sh --vcf cohort.vcf.gz \
    --ref IRGSP-1.0_genome.fa --outdir results/variants/03_filter
```

默认阈值（GATK 硬过滤，已按水稻调整）：`QUAL<30 | QD<2 | FS>60 | MQ<40 | SOR>3 |
MQRankSum<-12.5 | ReadPosRankSum<-8`，可用 `--min-qual` 等参数覆盖。

- GATK 模式先软标注（`flagged.vcf.gz` 保留审计轨迹）再硬剔除；bcftools 模式等价表达式；
- 输出 `filtered.vcf.gz` —— 这是进入 GWAS/群体遗传的标准输入。

### 2.4 注释（`scripts/variants/04_snpeff_annotate.sh`）

```bash
bash scripts/variants/04_snpeff_annotate.sh \
    --vcf results/variants/03_filter/filtered.vcf.gz \
    --outdir results/variants/04_annot
```

- 默认库 `Oryza_sativa`（Ensembl，首次自动下载）；自定义注释（如 RAP-DB）用
  `--build-custom --gff ... --fasta ... --db-name ...`；
- 输出：注释 VCF（ANN 字段）、**扁平影响表** `ann_table.tsv`（CHROM/POS/REF/ALT/基因/效应/影响/HGVS）、
  `snpeff_summary.html` 统计报告；
- 常用下游筛选：
  ```bash
  grep -P "\tHIGH\t" results/variants/04_annot/ann_table.tsv   # 高影响位点（终止/移码/剪接）
  grep "Os01g0100100" results/variants/04_annot/ann_table.tsv  # 按基因查
  ```

## 3. 输出结构

```
results/variants/
├── 01_align/     # BAM + flagstat/idxstats/stats + align_log.txt
├── 02_call/      # cohort.vcf.gz (+gVCF) + call_log.txt
├── 03_filter/    # filtered.vcf.gz + filter_stats.txt
└── 04_annot/     # annotated.vcf.gz + ann_table.tsv + snpeff_summary.html
```

## 4. 水稻专用注意点

| 事项 | 建议 |
|---|---|
| 参考基因组 | IRGSP-1.0（RAP-DB）或 MSU7；基因 ID 体系与后续注释保持一致 |
| 测序深度 | GWAS 群体 ≥10×；候选基因验证 ≥30×；GBS 按酶切方案评估覆盖均匀性 |
| 群体规模 | 联合分型（GenomicsDBImport）比逐样本 call 后合并更稳，建议 ≥20 样本时启用 |
| HWE 过滤 | 水稻自花授粉，偏分离普遍——`01_plink_qc.sh` 的 HWE 阈值不必过严（默认 1e-6 已宽松） |
| 注释库 | SnpEff `Oryza_sativa` 基于 Ensembl 注释；与 RAP-DB/MSU 基因号需用内置映射表转换（`assets/annotations/rapdb_to_msu_map.csv`） |

## 5. 与其它模块的衔接

| 下游任务 | 输入与命令 |
|---|---|
| GWAS（模块 7） | `bash scripts/gwas/01_plink_qc.sh --vcf filtered.vcf.gz ...` |
| 群体遗传（模块 8） | 同上得到 PLINK fileset → PCA/LD/Fst/ADMIXTURE |
| 泛基因组（模块 11） | 以图形泛基因组为参考时，把 `--ref` 换成泛基因组主路径序列 |
| 候选基因注释 | `ann_table.tsv` × GWAS 候选区间（`scripts/gwas/04_candidate_loci.R`） |
| 单倍型分析 | filtered VCF → 按基因区间提取基因型（RFGB/RiceVarMap2 在线复核） |
