# 空间组学模块使用指南（rice-spatial，v3.0）

> Rice-skills v3.0 新增 **模块 18 · 空间组学**：Stereo-seq / Visium 空间数据的
> 域识别、空间可变基因（SVG）与空间图分析，并对接 **RSTA（水稻时空图谱）** 与
> **RICE scGPT 基础模型**的本地运行。
> 上游衔接：测序公司/平台输出（SAW / Space Ranger）；下游衔接：单细胞模块（模块 17）。

## 0. 流程总览

```
Stereo-seq（SAW 输出 bin 矩阵）/ Visium（Space Ranger）
    ──01_spatial_pipeline.py──▶ 表达域 + 空间散点 + Moran's I SVG + 邻域富集
    ├─▶ 细胞类型解卷积（Cell2location / Tangram / RCTD，外部工具）
    ├─▶ RSTA 参考对照 / RICE scGPT 映射（§4）
    └─▶ 空间通讯（CellChat / Squidpy，外部工具）
```

## 1. 环境安装

```bash
conda create -n spatial python=3.11 -y && conda activate spatial
pip install scanpy squidpy          # 核心（squidpy 提供空间图/自相关统计）
pip install cell2location           # 解卷积（可选，需 GPU 更佳）
```

## 2. 数据准备

| 平台 | 上游 | 输入形态 |
|---|---|---|
| **Stereo-seq** | 官方 **SAW** 管线（gef → bin 矩阵；`github.com/STOmics/SAW`，运行时自行获取） | bin 矩阵 CSV/h5ad + 坐标（barcode, x, y） |
| **Visium** | Space Ranger | `filtered_feature_bc_matrix` + `spatial/` 目录（h5ad 直接含 obsm['spatial']） |

## 3. 主流程（`scripts/spatial/01_spatial_pipeline.py`）

```bash
python scripts/spatial/01_spatial_pipeline.py \
    --input demo_data/spatial/counts.csv --coords demo_data/spatial/coords.csv \
    --genes GeneD1_1,GeneD2_1,GeneD3_1 \
    --resolution 0.5 --outdir results/spatial/01_pipeline
```

| 参数 | 默认 | 说明 |
|---|---|---|
| `--input` / `--coords` | — / — | 矩阵 + 坐标；或直接 `.h5ad`（含 obsm['spatial']） |
| `--genes` | — | 需要在组织图上展示的基因 |
| `--resolution` / `--n-hvg` | 0.5 / 2000 | 域识别粒度 |
| `--n-svg` | 2000 | Moran's I 检验基因数（squidpy） |

**输出**：`domain_map`（表达域空间图）、`gene_maps`、`spatial_domains.csv`、
`svg_results.csv`（Moran's I 排序）、`neighbourhood_enrichment`（squidpy）、`processed.h5ad`。

**方法说明**：域识别采用"表达聚类 + 空间可视化"的基线方案（Leiden on PCA）；
更高级的空间域方法（BayesSpace、STAGATE、GraphST）与解卷积（Cell2location、Tangram）
为外部工具，可基于 `processed.h5ad` 继续。

## 4. RSTA 集成（重点）

**RSTA（Rice Spatio-Temporal Atlas，南方科技大学等，*Cell* 论文）**提供：
- 40,155 注释基因；**788,030 单核 / 47 组织**；**43 个组织的空间基因表达**；**202 细胞类型 marker**；
- 配套 **RICE scGPT 单细胞基础模型**（本包所集成的核心资产）。

> ⚠️ **获取方式**：RSTA 代码仓库（`github.com/RiceSpatioTemporalAtlas/RSTA`）**未附许可证文件**，
> 本包按合规规则**不复制其任何代码**——以下为"用户运行时自行获取与运行"的操作指引。

### 4.1 本地运行 RICE scGPT（操作指引）

```bash
# 1) 获取官方代码（用户自行执行）
git clone https://github.com/RiceSpatioTemporalAtlas/RSTA.git
cd RSTA

# 2) 阅读官方说明与模型目录
#    README.md
#    01.Rice-scGPT_model/    ← 模型与教程（Tutorial_Annotation.ipynb 等）
#    00.Stereo-seq_analysis/ ← 时空图谱分析代码
#    04.Common_scripts/      ← 通用脚本

# 3) 按 01.Rice-scGPT_model/ 内 notebook 的说明准备 Python 环境与模型权重
#    （环境依赖与权重获取方式以官方 notebook/文档为准）
```

**典型用途**（对应官方 Tutorial）：

| 任务 | 说明 |
|---|---|
| 细胞类型注释 | 对缺乏 marker 的细胞群给出注释建议（与 RSTA 202 类对齐） |
| 跨数据集整合 | 把自有 scRNA-seq 映射进 RSTA 参考空间，便于比较与复用 |
| 参考映射验证 | 本包模块 17 的 marker 打分结果 → scGPT 复核 |

### 4.2 时空图谱工作流（整理要点）

官方 `00.Stereo-seq_analysis` 的分析组织方式（要点 distillation，供自有数据参考）：

1. **空间 bin 矩阵质控**（低表达 bin 过滤、总计数分布）；
2. **表达域识别 + 空间可视化**（对应本包 §3 流程）；
3. **空间 marker 与梯度分析**（本包 `svg_results.csv` + `gene_maps` 覆盖此步）；
4. **单细胞—空间联合**（snRNA-seq 与空间数据的标签转移 / 解卷积）；
5. **发育轨迹的空间映射**（沿组织轴的空间拟时序）。

> 本包脚本覆盖第 2–3 步的基线实现；第 1/4/5 步与 RSTA 的专用管线配合使用。

## 5. 下游分析（外部工具）

| 分析 | 工具 | 说明 |
|---|---|---|
| 细胞类型解卷积 | **Cell2location**（GPU）/ Tangram / RCTD | 参考 scRNA-seq → 每 bin 细胞类型比例 |
| 空间通讯 | CellChat / Squidpy ligrec | 配体-受体空间共表达 |
| 与单细胞整合 | Harmony / scGPT 嵌入 | 空间与单细胞的统一表征 |
| 空间 TWAS 类分析 | 自定义 | 空间域差异基因 × GWAS 区间 |

## 6. 水稻专用注意点

| 事项 | 说明 |
|---|---|
| 组织形态 | 水稻叶片维管束/叶肉/表皮带状结构 → bin 尺寸影响域边界（Stereo-seq 常见 bin20-bin50） |
| 叶绿体 | 绿色组织叶绿体基因占比高，注意归一化与过滤策略 |
| 参考图谱 | 优先对照 RSTA 对应组织的空间模式；胚乳背腹轴等已知空间分工可作阳性对照 |
| 空间与单细胞 | 先 scRNA-seq 注释，再做解卷积；纯空间数据的域注释参考 RSTA 202 细胞类型 |

> 演示数据：`python scripts/utils/make_demo_spatial.py --outdir demo_data/spatial`
> （50×50 网格 2500 bins / 4 个象限域 / 含空间梯度基因；需 scanpy+squidpy 跑主流程）。
