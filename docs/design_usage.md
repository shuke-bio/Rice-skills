# 设计模块使用指南（引物设计 + CRISPR 靶点设计，v2.1）

> Rice-skills v2.1 新增 **模块 9 · 设计与验证**：本地离线完成
> **sgRNA 靶点设计 → 双靶点配对 → pYLCRISPR/Cas9 构建引物生成 → 通用引物设计**。
> 配套实验操作参考见 [lab_protocols.md](lab_protocols.md)。

## 0. 输入

| 文件 | 说明 |
|---|---|
| `gene.fa` | 目标基因序列（FASTA；单条记录，50-5000 bp 均可） |
| `gene.gff3`（可选） | 与 FASTA 坐标对齐的 GFF3（CDS/exon 记录）——用于"优先 CDS"筛选与 qPCR 跨内含子设计 |

演示数据：`python scripts/utils/make_demo_genomics.py --outdir demo_data/genomics`
会生成 `demo_gene.fa` + `demo_gene.gff3`（3 外显子演示基因）。

## 1. sgRNA 靶点设计（`scripts/design/01_sgrna_design.py`）

```bash
python scripts/design/01_sgrna_design.py \
    --fasta demo_data/genomics/demo_gene.fa \
    --gff demo_data/genomics/demo_gene.gff3 \
    --outdir results/design/01_sgrna --n-pairs 3
```

双链扫描 SpCas9 靶点（20 nt + NGG），按标准方案（Ma et al. 2015）的七条原则逐条筛选：

| # | 原则 | 实现 |
|---|---|---|
| 1 | 编辑区优先 CDS | 有 GFF 时：非 CDS 靶点直接剔除 |
| 2 | GC 含量 40–75%（偏好 50–75%） | `--min-gc/--max-gc`（默认 40–80 判定，评分向 60% 收敛） |
| 3 | 禁连续 ≥4 个 T（Pol III 终止信号） | 直接剔除 |
| 4 | 二级结构：与 sgRNA 骨架连续配对 ≥7 bp 剔除 | G-U 摆动配对模拟（等价于 mfold 检查的规则版） |
| 5 | 靶点选基因中上游 | 评分项：CDS 5' 段权重最高 |
| 6 | 避免 BsaI 位点干扰构建 | 靶点+PAM 内含 GGTCTC/GAGACC 时标注警告 |
| 7 | 优先 On-score 高 / off-target 少 | 透明评分 0–100 + 序列内近匹配统计（≤3 错配计数） |

**输出**：

| 文件 | 内容 |
|---|---|
| `sgrna_candidates.csv` | 全部候选：序列、链、坐标、GC、连续T、骨架配对、是否在 CDS、序列内脱靶数、评分、通过与否+原因 |
| `sgrna_recommended_pairs.csv` | **推荐的 3 组双靶点**（默认间隔 ≥100 bp、两者均通过、GC 差 ≤15%） |
| `construction_primers.csv` | 每组双靶点的 **gRT1+/U6aT1-/gRT2+/U6bT2-**（含"G 开头取 N19、否则 N20"规则）+ SP-L1/SP-R + 全套通用引物 |
| `sgrna_map.{pdf,svg,tiff}` | 基因结构图 + 靶点位置箭头（正/负链分色）+ 推荐配对区域 |

**重要提醒**：
- 全基因组脱靶评估（植物中通常非瓶颈，但建议做）请用 **CRISPR-P 2.0**（http://crispr.hzau.edu.cn/CRISPR2/）或
  **Cas-OFFinder**（http://www.rgenome.net/cas-offinder/）——本脚本只做序列内特异性检查；
- 二级结构深度预测可用 mfold / OligoEvaluator 复核（与原则 4 的规则化检查互补）；
- 构建引物的接头序列转录自官方方案文档，**下单前请与最新版官方文件核对一次**。

## 2. 引物设计（`scripts/design/02_primer_design.py`，primer3 引擎）

依赖：`pip install primer3-py`（GPL-2.0，仅作为外部引擎调用）。

```bash
# ① qPCR 跨内含子引物（要求给定 GFF 的 exon 记录）
python scripts/design/02_primer_design.py --mode qpcr \
    --fasta demo_data/genomics/demo_gene.fa --gff demo_data/genomics/demo_gene.gff3 \
    --outdir results/design/02_primers_qpcr

# ② CRISPR 靶点区验证引物（扩增靶点 ±100 bp，供 Sanger 验证）
python scripts/design/02_primer_design.py --mode validation \
    --fasta demo_data/genomics/demo_gene.fa \
    --target-start 300 --target-end 330 --flank 100 \
    --outdir results/design/02_primers_val

# ③ 通用 PCR / 克隆引物
python scripts/design/02_primer_design.py --mode generic \
    --fasta demo_data/genomics/demo_gene.fa --target-start 1 --target-end 600 \
    --outdir results/design/02_primers_generic
```

- **qPCR 模式**在拼接后的转录本上设计，并**逐个内含子结合点定向设计**（扩增子必跨结合点），
  基因组 DNA 污染无法扩增；输出同时给出转录本与基因组两套坐标；
- 每对引物附**完整质量指标**：Tm、GC、发夹 Tm、同源二聚体 Tm、异源二聚体 Tm、端稳定性、Primer3 penalty；
- 质量判读：|Tm(F)-Tm(R)| ≤ 2 °C；hairpin/homodimer Tm 低于退火温度即可；异源二聚体 Tm 越低越好。

## 3. 输出结构

```
results/design/
├── 01_sgrna/              # 靶点表 + 推荐双靶点 + 构建引物 + 靶点分布图
├── 02_primers_qpcr/       # qPCR 引物 + 覆盖图
├── 02_primers_val/        # 验证引物 + 覆盖图
└── 02_primers_generic/    # 通用引物 + 覆盖图
```
所有图 PDF+SVG+TIFF，所有目录含 `run_info.txt` 溯源。

## 4. 免责与引用

- 评分模型为**透明启发式**（非 Doench/Shen 等训练模型），用于优先级排序；最终靶点建议
  在 skl.scau.edu.cn / CRISPR-P 等平台复核后再下单；
- 构建体系引用：Ma X, Zhang Q, Zhu Q, et al. *A robust CRISPR/Cas9 system for convenient
  high-efficiency multiplex genome editing in monocot and dicot plants.* Mol Plant, 2015;
  实验流程参考 Ma X, Liu Y-G. Curr Protoc Mol Biol, 2016.

---

## v3.2 新增脚本

### qPCR 相对定量 `03_qpcr_analysis.py`（设计-验证闭环）

引物设计后的湿实验数据落地：把 qPCR 仪导出的 **Ct/Cp 值**变成出版级相对表达图：

```bash
python scripts/design/03_qpcr_analysis.py \
    --ct-table ct_values.csv --ref-gene OsUBQ5 --calibrator Control \
    --outdir results/design/03_qpcr

# 附转录组交叉验证（qPCR vs RNA-seq 相关性）
python scripts/design/03_qpcr_analysis.py --ct-table ct_values.csv \
    --ref-gene OsUBQ5 --calibrator Control \
    --rna-table rnaseq_values.csv --rna-format fpkm \
    --outdir results/design/03_qpcr
```

**输入格式（长表，从仪器导出后整理一次即可）**：

```csv
gene,group,bio_rep,ct
GeneA,Control,1,18.20
GeneA,Control,1,18.15     # 同一 bio_rep 多行 = 技术重复，自动取均值
GeneA,Treat,1,16.10
RefGene,Control,1,22.00   # 内参基因以普通 gene 行出现（--ref-gene 指定名称）
```

**计算**：技术重复均值 → ΔCt（同组同重复）→ ΔΔCt（校准组）→ **2^-ΔΔCt**；
校准组均值≈1（生物个体差异保留在 SD 中）。

**输出**：
- `ddct_table.csv`（逐重复明细：Ct/ΔCt/ΔΔCt/相对表达）
- `rel_expression_summary.csv`（gene × group：mean/SD/SEM/n）
- `qpcr_<gene>.pdf/png`（每基因条形图，mean±SD，Nature 配色）
- `qpcr_all_genes_grid.pdf`（≤24 基因总览）
- 加 `--rna-table` 时：`correlation_qpcr_rnaseq.csv` + 散点图（Pearson/Spearman）

**注意事项**：
- "Undetermined"/空值自动剔除并告警；缺内参 Ct 的样本自动跳过
- `--rna-format counts` 时自动 log2(x+1) 变换后再相关
- 组间显著性检验建议在 2^-ΔΔCt 值上用 Welch t 或 Mann-Whitney
  （样本量估计见 `phenomics/05_power_analysis.R`）
- `--demo` 用内置示例 Ct 值（已知 4×/2×/0.5× 效应）自检，数值精确回收
