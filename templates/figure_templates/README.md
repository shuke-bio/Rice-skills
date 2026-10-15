# templates/figure_templates/ — 出版级绘图模板库

> 本目录提供 Rice-skills 全部出版图表的**可复制起点模板**：每个模板函数
> 返回一张符合 Nature 绘图规范的图（配色来自 `config/nature_palette.json`，
> Arial/Helvetica 8 pt，PDF + SVG + TIFF 300 dpi，单栏 89 mm）。
> **`scripts/` 下的成品绘图脚本即由这些模板生成**——写新图时先在这里找模板，
> 复制改写，而不是从零开始。

## 两个后端如何选

模板提供 R（ggplot2/pheatmap）与 Python（matplotlib/seaborn）两个实现，
**视觉完全一致**（共用同一调色板与主题）。选择规则：

| 场景 | 推荐后端 | 模板文件 |
|---|---|---|
| 分析主栈是 R（DESeq2/WGCNA/OPLS-DA/MOFA 的结果直接成图） | R | `ggplot2_templates.R` |
| 质控/预处理脚本是 Python，或需要嵌入 Python 流水线 | Python | `matplotlib_templates.py` |
| 批量 QC 图（PCA、箱线、缺失值） | Python（脚本已在 Python） | `matplotlib_templates.py` |
| 需要注释热图且环境有 pheatmap / ComplexHeatmap | R | `ggplot2_templates.R` |

两个后端的函数一一对应，但 Python 版未提供 `template_bar()`（R 独有）；
Python 版热图用 `seaborn.heatmap` + 手工注释色条，R 版用 `pheatmap`。

## 模板清单

| 模板函数 | 图类型 | R 保存方式 | Python 保存方式 |
|---|---|---|---|
| `template_volcano()` | 火山图（语义色 + 阈值线 + top5 基因标注） | `save_nature_figure()` | `save_nature_figure()` |
| `template_ma()` | MA 图（log10 baseMean vs log2FC） | 同上 | 同上 |
| `template_pca()` | PCA 得分图（95% 置信椭圆，可加 shape） | 同上 | 同上 |
| `template_dotplot()` | 富集点图（GeneRatio + Count + -log10 padj，按 Category 分面） | 同上 | 同上 |
| `template_bubble()` | 富集气泡图（NES/折叠富集 + set size + padj） | 同上 | 同上 |
| `template_heatmap()` | 热图（diverging 相关矩阵 / sequential 表达矩阵，groups5 注释色） | `save_pheatmap_nature()` | `save_nature_figure()` |
| `template_bar()` | 通用水平条形图 | `save_nature_figure()` | —（无 Python 版） |

## 如何使用（R 后端）

模板文件**不主动加载主题**，调用前必须先 source 共享主题（仓库根目录执行）：

```r
source("scripts/utils/plot_theme.R")                        # 主题与调色板
source("templates/figure_templates/ggplot2_templates.R")   # 模板库

# 复制改写：把模板函数体搬进 scripts/transcriptomics/07_plot_transcriptome.R，
# 数据来源换成流水线输出的 CSV，保存调用保持不变：
de <- read.csv("results/transcriptome/02_de/de_all_salt_vs_control.csv")
p <- template_volcano(de)
save_nature_figure(p, "results/transcriptome/07_figures/fig_volcano",
                   width_mm = 89, height_mm = 80)
```

注意：

- **颜色一律来自主题 API**（`rice_palette()` / `nature_semantic()` /
  `diverging_palette()` / `sequential_palette()` / `scale_*_nature()`），
  模板内没有硬编码调色板外的颜色；改写时也请保持这一纪律。
- `template_heatmap()` 返回 pheatmap 对象，必须用 `save_pheatmap_nature()`
  保存（`save_nature_figure()` 只接受 ggplot/ComplexHeatmap 对象）。
- `significant` 列缺省时可由阈值自动推导（volcano/MA）。

## 如何使用（Python 后端）

```bash
# 方式 1（推荐）：把模板复制到主题旁边使用
cp templates/figure_templates/matplotlib_templates.py scripts/utils/
cd scripts/utils && python matplotlib_templates.py   # 自带合成数据演示

# 方式 2：在仓库根目录直接运行（自动把 scripts/utils 加入 sys.path）
python templates/figure_templates/matplotlib_templates.py
```

```python
# 方式 3：在自己的脚本中导入
import sys
sys.path.insert(0, "scripts/utils")
from matplotlib_templates import template_volcano
from plot_theme import save_nature_figure

fig, ax = template_volcano(de_df)          # 返回 (fig, ax)
save_nature_figure(fig, "results/fig_volcano", width_mm=89, height_mm=80)
```

`python matplotlib_templates.py` 会以 numpy 合成数据（seed=42）演示全部
模板并输出到**系统临时目录**（不污染仓库），可直接运行验证环境。

## 输入数据列名约定

模板假定流水线脚本输出的列名如下（与 `SKILL.md` 第 2 节、
`scripts/utils/data_io.py` 的加载约定一致）：

| 模板 | 数据框 | 必需列 | 可选列 / 说明 |
|---|---|---|---|
| volcano / MA | 差异表达表 | `gene_id`, `log2FoldChange` | `padj`（`significant` 缺省时用于推导）, `baseMean`（MA 必需）, `significant` |
| PCA | 得分表 | `sample_id`, `PC1`, `PC2` | 元数据按 `sample_id` 索引；`color_col`/`shape_col` 指定上色/形状列 |
| dotplot | 富集表 | `Description`, `p.adjust` | `GeneRatio`（或 `Count`）, `Count`, `Category`（分面） |
| bubble | 富集/GSEA 表 | `Description` + p 值列 | `NES`（或折叠富集列）, `set_size`/`Count`（点大小） |
| heatmap | 矩阵 + 注释表 | 行/列名 | 注释表按矩阵列名索引（对称相关矩阵同时注释行） |
| bar | 任意表 | 值列 + 类别列 | `fill` 分组列可选 |

`significant` 列的取值约定：`Up` / `DOWN` / `NS`（大小写与同义词均被容忍，
内部归一化为三级语义色：上调红 `#CE6A6B`、下调蓝 `#5B84AE`、不显著灰
`#C9C4BE`）。

## 绘图规范要点（改写时不可破坏）

- 语义色固定：up=红 / down=蓝 / not_significant=灰 / highlight=琥珀
  `#E8A33D`；阈值线用 `reference_line` 灰色虚线。
- 字体 Arial/Helvetica，正文 8 pt、标题 9 pt 加粗、刻度/图例 7 pt，线宽 0.5 pt。
- 单栏 89 mm / 1.5 栏 120 mm / 双栏 180 mm；每张图同时输出
  PDF + SVG + TIFF（300 dpi）。
- 每张图对应明确科学结论；图注科学含义写在生成脚本注释中（英文注释）。
