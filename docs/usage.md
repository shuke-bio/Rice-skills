# 使用手册（Usage）

## 1. 准备输入数据

三个 CSV 文件（UTF-8、首行表头、逗号分隔）：

| 文件 | 首列 | 必需列 | 说明 |
|---|---|---|---|
| `counts.csv` | `gene_id` | 样本列（整数计数） | 行=基因，列=样本；**不要**预先归一化 |
| `metabolites.csv` | `metabolite_id` | 样本列；建议 `metabolite_name`/`class`/`kegg_id` | 缺失值留空即可 |
| `metadata.csv` | `sample_id` | `group` | 可选 `batch`、`variety` 及 WGCNA 性状列 |

自检清单：
1. 样本 ID 与矩阵列名**完全一致**（大小写敏感、无多余空格）；
2. 每组至少 3 个生物学重复（少于 2 个脚本会警告）；
3. count 为整数，代谢物丰度为数值；
4. `group` 列中对照组排第一（DESeq2 默认以首个水平为参照，也可用 `--contrast` 显式指定）。

没有数据？先生成示例数据试跑：

```bash
python scripts/utils/make_demo_data.py --outdir demo_data
```

### 1.1 测序公司/代谢组平台数据的适配（推荐先跑）

不同公司（诺禾/华大/百迈客/迈维/欧易等）交付的表格列名、ID 体系、分隔符各不相同，
用适配工具一步转成标准格式，并生成校验报告：

```bash
# 转录组（自动识别 Gene ID 列、MSU/RAP-DB/Ensembl ID 体系、非整数与 TPM 风险）
python scripts/utils/adapt_company_data.py --mode transcriptome \
    --input 公司表达量矩阵.xlsx --metadata metadata.csv --outdir adapted/
#   若 ID 是 RAP-DB/Ensembl，可加 --id-map 转换表.csv（两列 from,to）
#   若确实是"以浮点存储的原始 count"，可加 --round-counts
#   注意：DESeq2 需要原始 count；公司给的 TPM/FPKM 不能直接用于差异表达，需索取 count 矩阵

# 代谢组（支持中英文表头：名称/分类/KEGG ID → metabolite_name/class/kegg_id）
python scripts/utils/adapt_company_data.py --mode metabolome \
    --input 公司代谢组丰度表.csv --outdir adapted/
```

输出：`adapted/counts.csv` 或 `adapted/metabolites.csv` + `adapted/adaptation_report.txt`
（含行列统计、ID 体系检测、缺失值比例、样本与元数据匹配情况、修复建议）。

支持的输入：`.csv` / `.tsv` / `.txt`（自动识别分隔符）/ `.xlsx`（需 `pip install openpyxl`）。
Excel 含合并单元格、脚注行时请先清理或另存为 CSV。

### 1.2 水稻注释资源（ID 转换 / TF 列表 / GSEA 基因集）——已内置，开箱即用

仓库内置官方数据快照 `assets/annotations/`，**无需任何下载**即可使用：

| 文件 | 用途 |
|---|---|
| `rapdb_to_msu_map.csv` / `msu_to_rapdb_map.csv` | ID 转换，直接用于 `adapt_company_data.py --id-map` |
| `rice_TF_PlantTFDB_msu.csv`（1,862 TF / 56 家族） | 直接用于 `06_tf_analysis.R --tf-annotation` |
| `rice_TF_families.gmt`（53 个基因集） | 直接用于 `04_gsea.R --gmt` |

需要更新到官方最新版本时：

```bash
python scripts/utils/fetch_rice_annotations.py --force   # 自动从 RAP-DB / PlantTFDB 官方源刷新
#   网络拦截 TLS（报错含 "self-signed certificate"）→ 加 --insecure
#   RAP-DB 返回 403（反爬）或无法联网 → 浏览器下载官方文件后传本地路径：
#   --rapdb-url RAP-MSU_2018-03-29.txt.gz --tf-url Osj_TF_list.txt.gz
#   不加 --force 时优先使用内置快照（离线模式），缓存写入 data/annotations/
```

## 2. 一条命令跑全流程

```bash
bash run_pipeline.sh \
    --counts demo_data/counts.csv \
    --metadata demo_data/metadata.csv \
    --metabolites demo_data/metabolites.csv \
    --outdir results --stages all
```

常用参数：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--stages` | all | `all` / `transcriptome` / `metabolome` / `integration` |
| `--design` | `~ group` | DESeq2 设计公式，如 `~ batch + group` |
| `--contrast` | 全部两两 | 形如 `salt_vs_control`（分子_分母） |
| `--gmt` | 无 | GSEA 基因集（`.gmt`）；不提供则跳过 GSEA |
| `--tf-annotation` | 无 | PlantTFDB/iTAK 注释 csv；不提供则跳过 TF 分析 |
| `--traits` | 无 | WGCNA 性状列，逗号分隔，如 `tiller_number,Na_content` |
| `--padj` / `--log2fc` | 0.05 / 1 | 显著性阈值 |
| `--vip` | 1 | OPLS-DA VIP 阈值 |
| `--cor-cutoff` | 0.8 | 基因-代谢物相关性阈值 |
| `--python` / `--rscript` | 自动探测 | 指定解释器路径 |
| `--seed` | 42 | 随机种子 |

`run_pipeline.sh` 会把每一步的输出自动串给下一步（如 QCd 后的 count 矩阵 → DESeq2），
可选步骤（GSEA/TF/MOFA 等）在缺少输入或依赖时自动跳过并在
`results/pipeline_log.txt` 中记录 `SKIP`。

## 3. 分模块、分步运行

### 模块 1：转录组

```bash
# 01 质控过滤
python scripts/transcriptomics/01_qc_filter.py \
    --counts counts.csv --metadata metadata.csv \
    --outdir results/transcriptome/01_qc \
    --min-counts 10 --min-samples 3

# 02 差异表达（DESeq2）
Rscript scripts/transcriptomics/02_deseq2.R \
    --counts results/transcriptome/01_qc/counts_filtered.csv \
    --metadata metadata.csv --outdir results/transcriptome/02_de \
    --design "~ group" --contrast salt_vs_control --padj 0.05 --log2fc 1

# 03 GO/KEGG 富集
Rscript scripts/transcriptomics/03_enrichment.R \
    --de-sig results/transcriptome/02_de/de_sig_salt_vs_control.csv \
    --de-all results/transcriptome/02_de/de_all_salt_vs_control.csv \
    --outdir results/transcriptome/03_enrichment --org osa --ontologies "BP,CC,MF"

# 04 GSEA（需自备 GMT）
Rscript scripts/transcriptomics/04_gsea.R \
    --de-all results/transcriptome/02_de/de_all_salt_vs_control.csv \
    --gmt rice_kegg.gmt --outdir results/transcriptome/04_gsea --n-top 20

# 05 WGCNA
Rscript scripts/transcriptomics/05_wgcna.R \
    --counts results/transcriptome/01_qc/counts_filtered.csv \
    --metadata metadata.csv --outdir results/transcriptome/05_wgcna \
    --traits tiller_number,Na_content --top-n 10000 --min-module-size 30

# 06 转录因子（可选，需注释文件；--ppi 追加 STRING 蛋白互作网络）
Rscript scripts/transcriptomics/06_tf_analysis.R \
    --de-all results/transcriptome/02_de/de_all_salt_vs_control.csv \
    --tf-annotation planttfdb_rice.csv --outdir results/transcriptome/06_tf

# 07 出版图汇编
Rscript scripts/transcriptomics/07_plot_transcriptome.R \
    --results-dir results/transcriptome --outdir results/transcriptome/07_figures
```

### 模块 2：代谢组

```bash
# 01 预处理
python scripts/metabolomics/01_preprocess.py \
    --abundance metabolites.csv --metadata metadata.csv \
    --outdir results/metabolome/01_preprocess \
    --impute half-min --normalize auto

# 02 OPLS-DA + VIP + 差异代谢物
Rscript scripts/metabolomics/02_oplsda.R \
    --data results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --metadata metadata.csv --outdir results/metabolome/02_oplsda \
    --vip 1 --pvalue 0.05 --permutations 200

# 03 代谢通路富集（需 kegg_id 列）
Rscript scripts/metabolomics/03_pathway.R \
    --diff-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
    --outdir results/metabolome/03_pathway --org osa

# 04 化学分类统计
Rscript scripts/metabolomics/04_metabolite_classification.R \
    --diff-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
    --all-metabolites results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --outdir results/metabolome/04_classification

# 05 生物标志物（ROC + 随机森林）
Rscript scripts/metabolomics/05_biomarker.R \
    --data results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --metadata metadata.csv --outdir results/metabolome/05_biomarker --n-top 10

# 06 代谢物相关性网络
Rscript scripts/metabolomics/06_metabolite_network.R \
    --data results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --outdir results/metabolome/06_network --cor-cutoff 0.8 --method spearman

# 07 出版图汇编
Rscript scripts/metabolomics/07_plot_metabolome.R \
    --results-dir results/metabolome --outdir results/metabolome/07_figures
```

### 模块 3：多组学整合（转录组 + 代谢组）

```bash
# 01 MOFA 潜在因子
Rscript scripts/integration/01_mofa.R \
    --rna results/transcriptome/02_de/normalized_counts_vst.csv \
    --metabolites results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --metadata metadata.csv --outdir results/integration/01_mofa --n-factors 10

# 02 DIABLO（或 --method o2pls）
Rscript scripts/integration/02_diablo.R \
    --rna results/transcriptome/02_de/normalized_counts_vst.csv \
    --metabolites results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --metadata metadata.csv --outdir results/integration/02_diablo \
    --ncomp 2 --keepx "10,10"

# 03 基因-代谢物相关性
Rscript scripts/integration/03_correlation.R \
    --de-genes results/transcriptome/02_de/de_sig_salt_vs_control.csv \
    --de-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
    --rna results/transcriptome/02_de/normalized_counts_vst.csv \
    --metabolites results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --outdir results/integration/03_correlation --cor-cutoff 0.8 --method spearman

# 04 联合通路分析
Rscript scripts/integration/04_joint_pathway.R \
    --de-genes results/transcriptome/02_de/de_sig_salt_vs_control.csv \
    --de-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
    --outdir results/integration/04_joint_pathway --org osa

# 05 多组学网络（需 WGCNA 的 gene_modules.csv）
Rscript scripts/integration/05_network_integration.R \
    --gene-modules results/transcriptome/05_wgcna/gene_modules.csv \
    --rna results/transcriptome/02_de/normalized_counts_vst.csv \
    --metabolites results/metabolome/01_preprocess/metabolome_preprocessed.csv \
    --metadata metadata.csv --outdir results/integration/05_network

# 06 整合出版图
Rscript scripts/integration/06_plot_integration.R \
    --results-dir results/integration \
    --transcriptome-dir results/transcriptome \
    --metabolome-dir results/metabolome \
    --outdir results/integration/06_figures
```

## 4. 参数速查表

| 参数 | 适用脚本 | 默认 | 调优建议 |
|---|---|---|---|
| `--padj` | 02/03/04/06、代谢 02/03 | 0.05 | 宽松筛选可用 0.1 |
| `--log2fc` | 02/06 | 1 | 弱效应可降到 0.5 |
| `--vip` | 代谢 02 | 1 | 严格可用 1.5 |
| `--cor-cutoff` | 整合 03/05 | 0.8 | 样本少时抬高到 0.85 |
| `--top-n` | 05 | 10000 | 内存不足降到 5000 |
| `--min-module-size` | 05 | 30 | 模块过碎可提到 50 |
| `--n-factors` | 整合 01 | 10 | 先看方差解释曲线再定 |
| `--keepx` | 整合 02 | 10,10 | 用 perf 结果调（见输出表） |
| `--impute` | 代谢 01 | half-min | 缺失>20% 试 `knn` |
| `--normalize` | 代谢 01 | auto | `auto` = sum→log2→pareto |
| `--block` | 整合 03 | 500 | 显存/内存紧张时减小 |

## 5. 输出与图表

- 每张图同时输出 **PDF（矢量可编辑）+ SVG + TIFF（300 dpi）**；单栏 89 mm、宽图 120/180 mm。
- 每个结果目录都有 `run_info.txt`（运行时间、全部参数、关键包版本、随机种子）。
- 完整输出清单见 [output_description.md](output_description.md)。

## 6. 常见问题（FAQ）

**Q1 只想做两组差异表达，不做后面步骤？**
只跑 01 + 02 即可，或 `bash run_pipeline.sh --stages transcriptome` 后忽略后续目录（可选步骤会自动跳过）。

**Q2 多组实验怎么做比较？**
`--contrast` 用逗号分隔多个（`salt_vs_control,drought_vs_control`），不写则自动做所有两两组合。

**Q3 有批次效应？**
`metadata.csv` 加 `batch` 列，分析时 `--design "~ batch + group"`。

**Q4 水稻没有 MSigDB，GSEA 怎么办？**
用 KEGG osa 或 PlantTFDB 自建 GMT（脚本头部注释给出构建方法），或跳过 GSEA 用 ORA 富集（03）替代。

**Q5 结果能复现吗？**
所有脚本固定 `seed=42` 并记录在 `run_info.txt`；但 DESeq2/WGCNA 结果会随包版本变化，升级依赖后建议重跑比对。

**Q6 想换物种？**
改三处：注释包（`org.Osativa.eg.db`）、KEGG 代码（`--org`，如 `zma`/`ath`）、TF 注释来源。流程不变。

---

## v3.2 新增脚本（模块 3 扩展）

### TF 调控网络推断 `07_grn_boost.py`

从表达矩阵推断**有方向**的 TF → 靶基因调控边（GRNBoost2 算法内核的 sklearn 实现，
无需安装 Java/arboreto；WGCNA 给共表达模块，本脚本给调控边排序）：

```bash
python scripts/integration/07_grn_boost.py \
    --expr demo_data/transcriptomics/counts_norm.csv \
    --tfs assets/annotations/rice_tfs_planttfdb.csv \
    --top-edges 5 --outdir results/integration/07_grn
```

- 对每个靶基因以候选 TF 为特征拟合梯度提升回归（250 树、subsample 0.632），
  特征重要性即调控分数；`--all-genes` 可放开全部基因作候选
- 输出 `grn_edges.csv`（全量排序边）+ `top_edges.csv`（每靶基因 top-K）；`--demo` 自检
- 重要程度是关联分数、非直接结合证据——与 SCENIC/WGCNA、motif 分析（tf 模块）交叉验证
