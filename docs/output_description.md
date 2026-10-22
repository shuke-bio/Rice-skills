# 输出结果说明（Output description）

所有结果写入 `--outdir`（默认 `results/`），按「模块/步骤」分目录。
**每张图同时输出 PDF + TIFF + PNG(300 dpi) 三种格式**，下文用 `*.{pdf,tiff,png}` 表示。
每个目录均含 `run_info.txt`：运行时间、全部参数、关键包版本、随机种子（默认 42）。

```
results/
├── pipeline_log.txt                        # 全流程日志（含 PASS/FAIL/SKIP，由 run_pipeline.sh 生成）
├── transcriptome/
│   ├── 01_qc/
│   ├── 02_de/
│   ├── 03_enrichment/
│   ├── 04_gsea/
│   ├── 05_wgcna/
│   ├── 06_tf/
│   └── 07_figures/
├── metabolome/
│   ├── 01_preprocess/
│   ├── 02_oplsda/
│   ├── 03_pathway/
│   ├── 04_classification/
│   ├── 05_biomarker/
│   ├── 06_network/
│   └── 07_figures/
└── integration/
    ├── 01_mofa/
    ├── 02_diablo/
    ├── 03_correlation/
    ├── 04_joint_pathway/
    ├── 05_network/
    └── 06_figures/
```

---

## 模块 1：转录组（results/transcriptome/）

### 01_qc/ — 数据质控与过滤
| 文件 | 说明 |
|---|---|
| `counts_filtered.csv` | 过滤低表达基因后的 count 矩阵（下游 DESeq2/WGCNA 的输入） |
| `qc_metrics.csv` | 每样本指标：文库大小、检出基因数、测序深度相关统计 |
| `qc_report.txt` | 文本报告：过滤前后基因数、离群样本提示（与中位相关性偏低者） |
| `library_size.*` | 每样本文库大小条形图（检查测序深度一致性） |
| `detected_genes.*` | 每样本检出基因数图（检查样本质量一致性） |
| `sample_pca.*` | PCA 散点图：颜色=分组，形状=批次；判断组间分离与批次效应 |
| `sample_correlation_heatmap.*` | 样本间相关性热图（高变基因）；识别离群/错标样本 |
| `run_info.txt` | 运行溯源 |

### 02_de/ — 差异表达
| 文件 | 说明 |
|---|---|
| `normalized_counts_vst.csv` | VST 归一化表达矩阵（下游 WGCNA/整合分析输入） |
| `de_all_<contrast>.csv` | 全部基因：`gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj` |
| `de_sig_<contrast>.csv` | 显著差异基因（padj ≤ 阈值且 \|log2FC\| ≥ 阈值），含 `significant=UP/DOWN` 列 |
| `de_summary.csv` | 每个比较的上调/下调/不显著基因计数 |
| `dispersion_estimates.*` | DESeq2 离散度诊断图（评估模型拟合） |
| `run_info.txt` | 运行溯源（含 design 公式与 contrast） |

### 03_enrichment/ — GO/KEGG 富集
| 文件 | 说明 |
|---|---|
| `go_<ONT>_<UP|DOWN>.csv` | GO 富集结果（BP/CC/MF × 上调/下调），含 Description、GeneRatio、BgRatio、p、padj、geneID |
| `kegg_<UP|DOWN>.csv` | KEGG 通路富集结果（organism = `osa`） |
| `enrichment_overview.csv` | 各类 top 10 汇总，便于快速制表 |
| `dotplot_go.*` / `dotplot_kegg.*` | GO/KEGG 点图（点大小=基因数，颜色=padj） |
| `run_info.txt` | 运行溯源 |

### 04_gsea/ — 基因集富集分析
| 文件 | 说明 |
|---|---|
| `gsea_results.csv` | 全部基因集：ID、NES、pvalue、padj、set size、core enrichment 基因 |
| `gsea_nes_bubble.*` | NES 气泡图（正负按语义色着色） |
| `gsea_enrichplot_<set>.*` | top N 通路的富集曲线（running ES） |
| `run_info.txt` | 运行溯源（含 GMT 来源） |

### 05_wgcna/ — 共表达网络
| 文件 | 说明 |
|---|---|
| `gene_modules.csv` | 基因 → 模块归属（下游整合分析 05 的输入） |
| `module_trait_correlation.csv` | 模块-性状相关系数与 p 值（长表） |
| `module_trait_heatmap.*` | 模块-性状关联热图（标注 r 与 p） |
| `module_dendrogram.*` | 基因树状图 + 模块颜色条 |
| `soft_threshold_selection.*` | 软阈值筛选图（scale-free 拟合指数） |
| `hub_genes.csv` | 各模块 hub 基因（按 kME 排序） |
| `module_top_hub.csv` | 每模块的 top hub 基因 |
| `module_enrichment.csv` | 主要模块的 GO 富集结果 |
| `run_info.txt` | 运行溯源 |

### 06_tf/ — 转录因子（可选）
| 文件 | 说明 |
|---|---|
| `de_TF.csv` | 差异表达转录因子（含家族注释与显著性） |
| `tf_family_stats.csv` | TF 家族统计（up/down 数量、富集 p 值） |
| `tf_family_top20.*` | TF 家族条形图 |
| `ppi_edges.csv` / `ppi_nodes.csv` | 蛋白互作网络边/点表（`--ppi` 时生成） |
| `ppi_network.*` | PPI 网络图 |
| `run_info.txt` | 运行溯源 |

### 07_figures/ — 出版图汇编
| 文件 | 说明 |
|---|---|
| `fig1_qc_overview.*` | QC 总览：文库大小 + 检出基因 + PCA + 相关性热图 |
| `fig2_de_overview.*` | 差异表达总览：火山图（标注 top5 基因）+ MA 图 + DE 计数 |
| `fig3_enrichment_overview.*` | 富集总览：GO 点图 + KEGG 点图 |
| `fig4_wgcna_overview.*` | WGCNA 总览：模块-性状热图 + 模块大小 + hub 基因表达热图 |
| `fig5_tf_overview.*` | TF 家族总览（有注释文件时） |

---

## 模块 2：代谢组（results/metabolome/）

### 01_preprocess/
| 文件 | 说明 |
|---|---|
| `metabolome_preprocessed.csv` | 预处理后矩阵（保留 `metabolite_id` 与注释列，下游统一输入） |
| `missing_value_summary.csv` | 每代谢物/每样本缺失统计 |
| `pca_scores.csv` | 归一化后 PCA 得分（样本 × PC） |
| `missing_value_barplot.*` | 缺失值分布图 |
| `pca_before.*` / `pca_after.*` | 归一化前后 PCA 对比 |
| `boxplot_normalization.*` | 归一化前后箱线图（评估离散度改善） |
| `run_info.txt` | 运行溯源（含填充与归一化方法） |

### 02_oplsda/
| 文件 | 说明 |
|---|---|
| `oplsda_scores.csv` | OPLS-DA 样本得分（含分组） |
| `oplsda_loadings.csv` | 各代谢物载荷（p(1)、p(corr) 等） |
| `vip_scores.csv` | 全部代谢物 VIP 值 |
| `diff_metabolites.csv` | **差异代谢物**（VIP ≥ 阈值且 p < 阈值）：VIP、p、padj、log2FC、significant |
| `oplsda_score_plot.*` | OPLS-DA 得分图（组间分离程度） |
| `s_plot.*` | S-plot：筛选贡献最大的代谢物 |
| `vip_barplot.*` | VIP 前 30 条形图 |
| `permutation_test.*` | 置换检验（R2Y/Q2 分布，判断模型是否过拟合） |
| `pca_scores.*` / `pca_scores.csv` | 无监督 PCA 基线 |
| `run_info.txt` | 运行溯源 |

### 03_pathway/
| 文件 | 说明 |
|---|---|
| `pathway_enrichment.csv` | 代谢通路富集：pathway_id、name、hits、fold enrichment、p、padj、compound_ids |
| `pathway_bubble.*` | 通路气泡图 |

### 04_classification/
| 文件 | 说明 |
|---|---|
| `class_counts.csv` | 各化学类别（黄酮/生物碱/脂质…）的差异数量、背景数量、占比与富集 p 值 |
| `class_barplot.*` | 化学类别组成图 |
| `class_upset.*` | 类别交集 UpSet 图（UpSetR 可用时） |

### 05_biomarker/
| 文件 | 说明 |
|---|---|
| `biomarker_roc.csv` | 候选代谢物 AUC 与 95% CI、p 值 |
| `roc_coordinates.csv` | ROC 曲线坐标（便于自定义绘图） |
| `rf_importance.csv` | 随机森林重要性（MeanDecreaseAccuracy / Gini） |
| `biomarker_panel.csv` | ROC(AUC 高) 与 RF(重要性高) 的**交集标志物组合** |
| `roc_curves.*` | ROC 曲线面板 |
| `rf_importance_plot.*` | 重要性排序条形图 |

### 06_network/
| 文件 | 说明 |
|---|---|
| `metabolite_edges.csv` | 显著相关的代谢物对（from, to, r, padj） |
| `metabolite_nodes.csv` | 节点属性：degree、betweenness、是否枢纽 |
| `network_plot.*` | 代谢物相关性网络图 |

### 07_figures/
`fig1_preprocess_overview.*`（PCA + 箱线 + 缺失值）、`fig2_oplsda_overview.*`（得分 + S-plot + VIP）、
`fig3_biomarker_overview.*`（ROC + RF 重要性）、`fig4_network_overview.*`（代谢网络）。

---

## 模块 3：多组学整合（results/integration/）

### 01_mofa/
| 文件 | 说明 |
|---|---|
| `mofa_model.hdf5` | MOFA 模型文件（可重复加载分析） |
| `variance_explained.csv` / `variance_explained.*` | 每个因子在各组学解释的方差比例 |
| `factor_scores.csv` | 样本 × 因子得分 |
| `factor_scores_heatmap.*` | 因子得分热图（按分组注释） |
| `factor_loadings_rna.csv` / `factor_loadings_metabolites.csv` | 各因子 top 特征载荷 |
| `top_features_per_factor.csv` | 每因子最重要特征汇总 |

### 02_diablo/
| 文件 | 说明 |
|---|---|
| `diablo_scores.csv` | 样本联合得分（各成分） |
| `diablo_loadings.csv` | 各组学特征载荷 |
| `diablo_perf_error_rates.csv` | 交叉验证错误率（用于选择 keepX/成分数） |
| `diablo_perf.*` | 性能评估图 |
| `diablo_selected_features.csv` | 被选中的关键特征（基因/代谢物） |
| `sample_projection.*` | 样本投影图 |
| `circos_loadings.*` | 载荷相关性 Circos 图（RNA ↔ 代谢物） |
| `o2pls_scores.csv` / `o2pls_loadings.*` | `--method o2pls` 时的等价输出 |

### 03_correlation/
| 文件 | 说明 |
|---|---|
| `gene_metabolite_pairs_sig.csv` | 显著基因-代谢物对：r、p、padj、方向（\|r\| ≥ 阈值） |
| `correlation_heatmap.*` | 显著对相关矩阵热图 |
| `correlation_network_edges.csv` / `correlation_network_nodes.csv` | 二分网络边表/点表（供 05/06 使用） |

### 04_joint_pathway/
| 文件 | 说明 |
|---|---|
| `gene_kegg_enrichment.csv` | 差异基因 KEGG 富集结果 |
| `metabolite_pathway_enrichment.csv` | 差异代谢物通路富集结果 |
| `joint_pathway_table.csv` | **联合通路表**：pathway、基因侧 padj、代谢物侧 padj、category（both / transcriptome_only / metabolome_only） |
| `joint_pathway_dotplot.*` | 联合通路点图（按 category 着色） |
| `pathway_links.csv` | 节点-通路链接表（供桑基图使用） |

### 05_network/
| 文件 | 说明 |
|---|---|
| `module_eigengenes.csv` | WGCNA 模块特征基因矩阵 |
| `module_metabolite_cor.csv` | 模块 × 代谢物相关性长表 |
| `module_metabolite_heatmap.*` | 模块-代谢物相关热图 |
| `multiomics_network_edges.csv` / `multiomics_network_nodes.csv` | 三层网络（基因-模块-代谢物）边/点表 |
| `hub_table.csv` | 枢纽节点表：degree、kME、所属模块 |
| `network_plot.*` | 多组学网络图 |

### 06_figures/
`fig1a_variance_explained.*` + `fig1b_factor_scores.*`（MOFA 整合总览）、
`fig2a_correlation_heatmap.*` + `fig2b_correlation_network.*`（相关性总览）、
`fig3_joint_pathway*`（联合通路）、`fig4_sankey*`（基因→通路→代谢物）、
`fig5_circos*`（按通路分区的 Circos）、`fig6a/b/c_upset*`（通路集合关系）。
缺失依赖包时对应面板跳过，日志会提示。

---

## 图表规范

- **配色**：`config/nature_palette.json`，低饱和色盲友好；上调红 `#CE6A6B`、下调蓝 `#5B84AE`、不显著灰 `#C9C4BE`、高亮琥珀 `#E8A33D`。
- **字体**：Arial / Helvetica（PDF 中文本可编辑）；正文 8 pt、标题 9 pt 加粗、刻度与图例 7 pt；线宽 0.5 pt。
- **尺寸**：单栏 89 mm，1.5 栏 120 mm，双栏 180 mm；TIFF 300 dpi。
- **后端**：R 脚本用 ggplot2/ComplexHeatmap，Python 脚本用 matplotlib，两端共用同一调色板，风格一致。

## 可复现性

每个 `run_info.txt` 记录：脚本名、运行时间、R/Python 版本、全部参数、关键包版本、随机种子。
撰写论文方法部分时可直接引用这些文件。
