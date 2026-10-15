# 蛋白质组学模块使用指南（rice-proteomics，v3.0）

> Rice-skills v3.0 新增 **模块 14 · 蛋白质组学**：DIA/DDA 定量数据（LFQ / iBAQ / TMT 比率）
> 的预处理、质控、差异蛋白分析（limma）、出版级火山图/热图，并与转录组、GWAS 模块衔接。
> 与转录组+代谢组形成"转录-蛋白-代谢"三角验证体系。

## 0. 流程总览

```
蛋白定量矩阵 ──01_preprocess.py──▶ 过滤+归一化矩阵 ──02_diff_protein.R──▶ 差异蛋白表
    ──03_plot_proteome.R──▶ 火山图 + 热图
    ──▶ 富集分析（复用模块 1 的 03_enrichment.R，蛋白→基因 ID 映射后）
```

## 1. 输入格式

| 文件 | 格式 | 说明 |
|---|---|---|
| 蛋白定量矩阵 | `protein_id` + 样本列（CSV） | LFQ/iBAQ 强度或 log 值；**行=蛋白，列=样本** |
| 样本元数据 | `sample_id, group [, batch, ...]` | 分组比较必需 |
| 蛋白-基因映射 | `protein_id, gene_id`（可选） | UniProt→RAP-DB/MSU 映射，供富集分析衔接 |

> 从测序/质谱公司拿到的 DIA-NN / MaxQuant / Spectronaut 输出：导出
> "Protein.Group / Protein.Ids + 各样本强度"为宽表即可（保留蛋白 ID 为第一列）。

## 2. 步骤详解

### 2.1 预处理（`scripts/proteomics/01_preprocess.py`）

```bash
python scripts/proteomics/01_preprocess.py \
    --abundance demo_data/proteomics/proteins.csv \
    --metadata demo_data/proteomics/metadata.csv \
    --outdir results/proteomics/01_preprocess --min-valid 0.7 --impute halfmin
```

| 参数 | 默认 | 说明 |
|---|---|---|
| `--min-valid` | 0.7 | 蛋白在 ≥70% 样本中检出才保留（蛋白组惯例） |
| `--log2` | auto | 强度 >50 自动 log2；可强令 always/never |
| `--impute` | halfmin | halfmin（log 空间最低值再低 1 单位）/knn/none |
| `--batch-col` | — | 批次列名 → 逐批次中位数中心化 |

输出 QC：缺失比例条形图、PCA（分组着色）、样本相关性热图、组内 CV 箱线图。

### 2.2 差异蛋白（`scripts/proteomics/02_diff_protein.R`）

```bash
Rscript scripts/proteomics/02_diff_protein.R \
    --matrix results/proteomics/01_preprocess/proteins_filtered.csv \
    --metadata demo_data/proteomics/metadata.csv \
    --ref control --alt salt --lf2fc 1 --fdr 0.05 \
    --outdir results/proteomics/02_diff
# 加批次校正：--batch-col batch（limma 引擎）
```

- 引擎：**limma**（装了则自动用，moderated t，推荐）；未装则回退 Welch t + BH（可用但功效低）；
- 输出全表（log2FC/stat/p/fdr/组均值）+ 显著表（双重阈值）+ 摘要。

### 2.3 出版图（`scripts/proteomics/03_plot_proteome.R`）

```bash
Rscript scripts/proteomics/03_plot_proteome.R \
    --diff results/proteomics/02_diff/protein_diff_all.csv \
    --matrix results/proteomics/01_preprocess/proteins_filtered.csv \
    --metadata demo_data/proteomics/metadata.csv --top 50 \
    --outdir results/proteomics/03_plots
```

- 火山图（上/下调着色 + top 蛋白标注：ggrepel）；
- Top N 差异蛋白 z-score 热图（样本按组排列）；
- 均为 PDF/SVG/TIFF 三格式。

### 2.4 功能富集与网络（复用现有模块）

```bash
# 差异蛋白 gene_id 列 → 模块 1 的富集脚本（GO/KEGG）
# 04_ppi 步骤：蛋白列表 → STRING (https://string-db.org) 提交 / STRINGdb 包
```

## 3. 水稻专用注意点

| 事项 | 说明 |
|---|---|
| 蛋白 ID 体系 | UniProt（Q...）→ RAP-DB（Os...）/ MSU（LOC_Os...）：用 `protein_gene_map.csv` 建立映射；水稻 UniProt 参考 proteome UP000059680 |
| 过滤阈值 | 70% 是常用起点；检测深度低的样品组可放宽到 50% 并注明 |
| 缺失机制 | 质谱缺失多为 MNAR（低丰度→缺失）→ halfmin 填充合理；KNN 用于随机缺失 |
| TMT 数据 | 比率数据先 log2 再归一化；注意 TMT 批次（--batch-col） |
| 与转录组对照 | mRNA-蛋白一致性常仅 0.3–0.5（翻译调控真实存在），不一致的基因往往是调控关键 |

## 4. 与其它模块的衔接

| 任务 | 衔接方式 |
|---|---|
| 磷酸化（如有） | 位点定量表同流程处理；激酶推断用 KSEA（文档级指导） |
| 富集分析 | `protein_gene_map.csv` 映射后 → 模块 1 `03_enrichment.R` / `04_gsea.R` |
| 与转录组联合 | 同组设计下 mRNA vs 蛋白 log2FC 散点（四象限分类：一致上调/下调/翻译调控） |
| 与 GWAS 衔接 | 差异蛋白对应基因 × GWAS 候选区间 → 功能候选优先级 |
| 多组学整合 | 蛋白矩阵可作为模块 3 的第三组学接入 MOFA/DIABLO（输入格式同代谢组） |

> 演示数据：`python scripts/utils/make_demo_proteomics.py --outdir demo_data/proteomics`
> （1200 蛋白 × 24 样本，80 个真差异蛋白含 ground_truth.csv；演示检出 69/80、0 假阳性）。

---

## v3.2 新增脚本

### UniProt 映射自动化 `04_uniprot_mapping.py`

从 UniProt 参考 proteome（UP000059680，日本晴）自动生成 `protein_gene_map.csv`
（蛋白 → MSU/RAP-DB 基因 ID），供本模块与模块 1 富集分析衔接：

```bash
python scripts/proteomics/04_uniprot_mapping.py --outdir results/proteomics/04_mapping
# 网络受限（TLS 拦截）：--insecure；或浏览器下载 TSV 后 --tsv-file 传入（完全离线）
python scripts/proteomics/04_uniprot_mapping.py --tsv-file uniprot_up000059680.tsv --outdir ...
```

- 解析规则：Gene Names 同时含 MSU（`LOC_Os04g55560`）与 RAP-DB（`Os04g0649100`）两种
  风格时 **MSU 优先**；版本号自动剥离；输出含 `id_style` 列与 `unmatched.csv` 诊断
- 首次下载后自动缓存（`data/annotations/`），之后离线可用；`--demo` 合成数据自检

### 差异蛋白结构域富集 `05_interpro_domains.py`

InterProScan TSV（`-T` 输出或 EBI 在线版）解析 + 结构域过度表示分析（蛋白版 GO 富集）：

```bash
python scripts/proteomics/05_interpro_domains.py \
    --interpro results/proteomics/interproscan.tsv \
    --de results/proteomics/02_diffpro/de_proteins.csv \
    --background results/proteomics/01_preprocess/proteins_filtered.csv \
    --outdir results/proteomics/05_domains
```

- 超几何检验 + BH FDR；`--level signature` 可切到 PFAM 签名级视图
- **强烈建议传 `--background`**（实测蛋白全集）——质谱缺失非随机，用全域蛋白做背景会偏
- 输出 `domain_enrichment.csv`（含 fold、p_adj）+ top 结构域条形图；`--demo` 自检
