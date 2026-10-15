# 表观基因组学模块使用指南（rice-epigenomics，v3.0）

> Rice-skills v3.0 新增 **模块 15 · 表观基因组学**：DNA 甲基化（WGBS）差异甲基化区域（DMR）
> 与染色质可及性/结合（ATAC-seq / ChIP-seq）差异区域（DAR/DHR）分析。
> **分工**：上游比对与 peak calling 用标准外部工具（本文给命令模板）；下游统计分析与出版图由本包脚本完成。

## 0. 流程总览

```
WGBS:  FASTQ ─Bismark─▶ 甲基化调用（CX_report）─区域聚合─▶ 区域甲基化矩阵
       ─01_methylation_analysis.R─▶ DMR 表 + 火山/热图
ATAC/ChIP: FASTQ ─bowtie2─▶ BAM ─MACS2─▶ peak ─计数矩阵─▶ 02_atac_chip_analysis.R
       ─▶ DAR/DHR 表 + PCA/火山/热图
```

## 1. 上游流程模板（外部工具，按需自行安装）

### 1.1 WGBS（DNA 甲基化）

```bash
# 比对（Bismark；水稻参考需先 bismark_genome_preparation）
bismark --genome /ref/IRGSP-1.0 -1 R1.fq.gz -2 R2.fq.gz -o bam/
# 去重 + 甲基化提取
deduplicate_bismark -p --bam bam/*_bismark_bt2_pe.bam
bismark_methylation_extractor --CX --cytosine_report \
    --genome_folder /ref/IRGSP-1.0 --bedGraph sample.dedup.bam
# → sample.CX_report.txt（chr, pos, strand, count_meth, count_unmeth, context, 三核苷酸）
```

**区域聚合**（本包需要的输入）：把 CX_report 按基因/转座子/自定义窗口聚合：
```bash
# 示例：按 200 bp 窗口聚合 CG 甲基化（awk/或 Python）
# 输出格式：region_id, chr, start, end, context, sampleA, sampleB, ...
```
> 也可直接用公开甲基化数据（RiceENCODE 等）的区域水平表。

### 1.2 ATAC-seq / ChIP-seq

```bash
# 比对（Bowtie2）+ 过滤
bowtie2 -x /ref/IRGSP-1.0 -1 R1.fq.gz -2 R2.fq.gz | samtools sort -o sample.bam
samtools index sample.bam
# 小鼠/人黑名单在水稻不存在；用重复区域黑名单（可选自定义）
# peak calling
macs2 callpeak -t sample.bam -c input.bam -f BAMPE -g 3.7e8 -n sample --outdir peaks/
# ATAC 加 --nomodel --shift -100 --extsize 200；ChIP 用 --broad（组蛋白修饰）
# 计数矩阵（featureCounts/多样本合并 counts）
featureCounts -T 8 -a merged_peaks.saf -F SAF -o peak_counts.txt *.bam
```

## 2. 下游分析（本包脚本）

### 2.1 DMR 分析（`scripts/epigenomics/01_methylation_analysis.R`）

```bash
Rscript scripts/epigenomics/01_methylation_analysis.R \
    --matrix demo_data/epigenome/methylation_regions.csv \
    --metadata demo_data/epigenome/metadata.csv \
    --ref control --alt salt --delta 0.1 --fdr 0.05 \
    --outdir results/epigenome/01_dmr
```

| 参数 | 默认 | 说明 |
|---|---|---|
| `--matrix` | — | 区域 × 样本矩阵（region_id, chr, start, end [,context] + 甲基化比例 0-1 或 0-100 自动识别） |
| `--delta` | 0.1 | 组间甲基化差阈值（分数尺度，即 10%） |
| `--test` | t | t / wilcox |

输出：`dmr_all.csv`（delta/p/fdr）、`dmr_significant.csv`（hyper/hypo）、甲基化分布箱线、
火山图、top-DMR 热图。

### 2.2 DAR / DHR 分析（`scripts/epigenomics/02_atac_chip_analysis.R`）

```bash
Rscript scripts/epigenomics/02_atac_chip_analysis.R \
    --matrix demo_data/epigenome/peak_counts.csv \
    --metadata demo_data/epigenome/metadata.csv \
    --ref control --alt salt --lf2fc 1 --fdr 0.05 \
    --outdir results/epigenome/02_dar
```

- CPM 过滤（默认为样本数一半的样本中 CPM≥1）→ CPM+log2 → PCA → limma（装了自动用）或 Wilcoxon+BH；
- 输出：`dar_all.csv`、`dar_significant.csv`（open/up 与 closed/down）、PCA、火山、top-peak 热图。

## 3. 水稻专用注意点

| 事项 | 说明 |
|---|---|
| 甲基化上下文 | 水稻 CG（~40-60%）、CHG、CHH 三类各有生物学含义；建议**分 context 分别聚合与检验**（矩阵里加 context 列），CHH 在生殖/胁迫中动态更强 |
| 重复元素 | 水稻转座子（TRE 家族）甲基化是研究热点；区域注释建议覆盖基因 + TE 两类 |
| 黑名单 | 水稻无通用黑名单；着丝粒/重复区信号异常时应自行过滤 |
| ATAC 参数 | 水稻基因组 3.7e8（MACS2 -g 3.7e8）；peak 宽度通常 200-800 bp |
| 组蛋白修饰 | H3K4me3/H3K27ac（激活）、H3K27me3（Polycomb）等；broad peaks 用 `--broad` |
| 关联转录组 | DMR/DAR → 邻近基因 → 与模块 1 的 DEG 交集：`bedtools intersect` 或坐标匹配 |

## 4. 与其它模块的衔接

| 任务 | 衔接方式 |
|---|---|
| DMR/DAR → 基因 | 区域坐标 × GFF 注释求交（`scripts/gwas/04_candidate_loci.R` 的坐标匹配思路） |
| 与转录组联合 | DAR/DMR 邻近基因 × 模块 1 DEG → "表观调控候选"列表 |
| 与 GWAS 联合 | DAR 区域 × GWAS 候选区间 → 变异可能通过染色质状态起作用的证据 |
| 公共数据对照 | RiceENCODE（ricencode.cn）查询组织特异性 ATAC/甲基化对照 |
| 上游工具文档 | Bismark / Bowtie2 / MACS2 / featureCounts 官方文档 |

> 演示数据：`python scripts/utils/make_demo_epigenome.py --outdir demo_data/epigenome`
> （3000 区域 + 4000 peak × 20 样本；真值 DMR/DAR 各含 ground_truth 表；
> 演示检出：DMR 194/200、DAR 248/250，均 0 假阳性）。
