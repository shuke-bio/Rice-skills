# 小 RNA 与可变剪接分析指南（v3.1）

> 水稻 **miRNA / siRNA** 与 **可变剪接（AS）** 分析的流程模板与外部工具指引。
> 这两类分析的工具链较重（比对/折叠/定量），本包以"输入准备 + 质量控制 + 与已有模块衔接"
> 的方式集成；核心统计可回流入模块 1/3 的差异分析与联合分析。

## 1. miRNA / 小 RNA 分析

### 1.1 流程总览

```
小 RNA 测序 FASTQ（18-30 nt 富集）
  → cutadapt 去接头/低质量
  → miRDeep2（已知 miRNA 定量 + 新 miRNA 预测；Bowtie 比对水稻基因组）
  → miRNA 表达矩阵（miRNA × 样本）
  → 差异 miRNA（DESeq2 —— 用本包模块 1 的 02_deseq2.R，输入换成 miRNA 矩阵）
  → 靶基因预测（psRNATarget 在线 / TargetFinder；水稻 miRNA 参考 miRBase/miR-Rice）
  → miRNA-靶基因调控对 → 与模块 1 DEG 求交 → 调控网络（模块 3 的 05_network_integration.R）
```

### 1.2 关键命令模板（外部工具，按需安装）

```bash
# 接头去除（小 RNA 接头按建库试剂盒定）
cutadapt -a TGGAATTCTCGGGTGCCAAGG -q 20 -m 18 -M 30 -o trimmed.fq raw.fq

# miRDeep2（需要：水稻基因组 FASTA、比对好的 BAM 或直接 fasta、miRNA 前体参考）
mapper.pl trimmed.fq -c -e -h -i -j -l 18 -m -p rice_bowtie_index -s collapser.csv
miRDeep2.pl trimmed.fq rice_genome.fa mature_miRBase.fa hairpin.fa collapser.csv ... 2>report.log

# 定量表整理为"miRNA × 样本"矩阵后：
Rscript scripts/transcriptomics/02_deseq2.R --counts mirna_matrix.csv \
    --metadata metadata.csv --contrast salt_vs_control --outdir results/mirna/02_de
```

### 1.3 水稻专用资源

| 资源 | 用途 |
|---|---|
| **miRBase**（mirbase.org，osa- 前缀） | 已知水稻 mature/前体 miRNA 参考 |
| **miR-Rice / PmiREN**（pmiren.com） | 水稻专门 miRNA 数据库（比 miRBase 更全的水稻新发现） |
| **psRNATarget**（plantgrn.noble.org/psRNATarget） | 植物 miRNA 靶基因预测（在线） |
| **TAS / phasiRNA 说明** | 21-nt phasiRNA 与 TAS 基因座在水稻生殖发育中重要——小 RNA 长度分布图先看 21/24 nt 峰 |

### 1.4 注意点

- 小 RNA 文库链特异性强：miRDeep2 结果按链读取；不同建库（DpnII 适配）会影响比对方向；
- 24 nt siRNA 峰常来自转座子/重复区甲基化（联动模块 15 的甲基化数据解读）；
- miRNA 定量的归一化推荐 DESeq2 尺寸因子或 spike-in（外参）；TPM 用于展示。

## 2. 可变剪接（Alternative Splicing）

### 2.1 流程总览

```
RNA-seq FASTQ（建议 ≥ 100 bp 双端，最好有长读长互补）
  → HISAT2 比对（--dta）→ StringTie 定量 + 组装
  → rMATS（成对样本的 AS 差异：SE/A5SS/A3SS/MXE/RI）或 SUPPA2（基于 TPM 的 ΔPSI）
  → 差异 AS 事件表（event, PSI, FDR）
  → 事件注释（suppa2 cluster / RMATS2simplerep）→ 受影响基因 → 与 DEG 交叉
  → 蛋白域破坏检查（AS 事件是否造成移码/提前终止 → 功能影响）
```

### 2.2 关键命令模板（外部工具）

```bash
hisat2 -p 16 --dta -x rice_index -1 R1.fq.gz -2 R2.fq.gz | samtools sort -o sample.bam
stringtie sample.bam -G annotation.gtf -o sample.gtf -A sample_gene_abundance.tsv

# rMATS（两组比较；需要 BAM 列表与 GTF）
rmats.py --b1 salt_b1.txt --b2 ctrl_b1.txt --gtf IRGSP_annotation.gtf \
    -od rmats_out --tstat 4 -t paired --nthread 8

# SUPPA2（轻量，基于 transcript TPM）
suppa.py generateEvents -i transcripts.gtf -o events -f SE A5 A3 MX RI
suppa.py psiPerEvent -i transcripts.gtf -e events.ioe -m tpm.tsv -o psi
suppa.py diffSplice -m tpm.tsv -i events.ioe -p psi_* -e cond1_rep1.psicon,cond1_rep2.psicon -e cond2... -o diff
```

### 2.3 水稻专用资源与注意点

| 事项 | 说明 |
|---|---|
| 参考转录本 | IRGSP/RAP-DB 的 GTF；水稻 AS 比例约 30-40% 基因有 AS（低于动物） |
| 长读长 | Iso-Seq/ONT 直接转录本定相，RSTA 类图谱项目已产出全长转录本参考 |
| 差异 AS vs 差异表达 | 分开判读：AS 影响（ΔPSI）与总表达量变化（DE）是两种调控 |
| 功能影响 | SE/RI 移码 → NMD 或截短蛋白；与蛋白组（模块 14）的肽段覆盖交叉验证 |
| 盐/旱胁迫 | 水稻胁迫下剪接因子（SR 蛋白）自身常差异表达——联动模块 1 的 TF 分析 |

## 3. 与本包的衔接

| 任务 | 衔接方式 |
|---|---|
| 差异 miRNA → 靶基因 | 靶基因列表 → 模块 1 的 `03_enrichment.R`（GO/KEGG） |
| miRNA-DEG 反负相关 | miRNA 表达 × 靶基因表达 Spearman（负相关筛调控对）→ 模块 3 网络 |
| 差异 AS 基因 | 与 DEG 交集 → "转录与剪接双重调控"基因清单 |
| 空间/单细胞中的剪接 | Scanpy 的 spliced/unspliced 层 → scVelo（模块 17 已引） |

---

## v3.2 新增：本地 miRNA 靶基因预测 `02_mirna_targets.py`

零依赖的 psRNATarget 式打分（前文流程模板的"本地预筛"环节；正式分析仍建议
psRNATarget 在线版 + 降解组验证）：

```bash
python scripts/evolution/02_mirna_targets.py \
    --mirna mirna.fasta --transcripts cds.fasta \
    --max-score 4.0 --outdir results/evolution/02_mirna_targets
```

- 打分规则：21 nt 窗口反平行配对；种子区（miRNA 第 2–8 位）错配重罚、3' 区轻罚、
  G:U 计半配对；**cleavage 级** = 种子完美 + 第 10–11 位完美配对，其余为 inhibition 级
- 输出 `mirna_targets.csv`（miRNA/靶标/位置/分数/模式/错配数）；`--demo` 自检
  （完美位点 score=0 判 cleavage、种子错配判 inhibition、无位点靶标不出现在结果）
- 与 miRDeep2 流程衔接：miRDeep2 得到成熟 miRNA 序列 → 本脚本全 CDS 扫描 → 候选靶
  基因列表 → 模块 1 富集 / 表达负相关分析
