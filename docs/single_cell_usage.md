# 单细胞组学模块使用指南（rice-single-cell，v3.0）

> Rice-skills v3.0 新增 **模块 17 · 单细胞组学**：水稻 scRNA-seq 的质控、聚类、
> marker 注释、轨迹分析流程，并对接 **RSTA 水稻时空图谱资源与 RICE scGPT 基础模型**。
> 上游衔接：测序公司返回的 Cell Ranger 输出；下游衔接：空间组学（模块 18）。

## 0. 流程总览

```
Cell Ranger 输出 / 表达矩阵 ──01_scanpy_pipeline.py──▶ 聚类 + marker + UMAP + 轨迹
    ├─▶ marker 注释（本包 marker_sets / RSTA 202 细胞类型 markers）
    ├─▶ RICE scGPT 辅助注释与整合（本地运行，见 §3）
    └─▶ 细胞通讯 CellChat / 调控网络 SCENIC（外部工具，见 §4）
```

## 1. 环境安装

```bash
# 推荐独立的 conda 环境（scanpy 生态依赖较新）
conda create -n sc python=3.11 -y && conda activate sc
pip install scanpy leidenalg igraph
pip install harmony-pytorch            # 多样本批次整合（可选）
pip install scvelo                     # RNA velocity（可选）
```

## 2. 主流程（`scripts/single_cell/01_scanpy_pipeline.py`）

```bash
# 基本用法（Cell Ranger 目录或 CSV/h5ad）
python scripts/single_cell/01_scanpy_pipeline.py \
    --input filtered_feature_bc_matrix/ \
    --outdir results/single_cell/01_pipeline \
    --min-genes 200 --max-mito 10 --resolution 0.8

# 多样本批次整合 + 参考 marker 打分
python scripts/single_cell/01_scanpy_pipeline.py --input counts.csv \
    --metadata sample_meta.csv --batch-key sample \
    --markers marker_sets.csv \
    --outdir results/single_cell/01_pipeline
```

| 参数 | 默认 | 说明 |
|---|---|---|
| `--input` | — | Cell Ranger 目录 / `.h5ad` / counts CSV（cells × genes 自动识别方向） |
| `--min-genes` / `--max-mito` / `--min-cells` | 200 / 10 / 3 | 细胞与基因过滤 |
| `--n-hvg` / `--n-pcs` / `--resolution` | 2000 / 30 / 0.8 | 特征选择与聚类粒度 |
| `--batch-key` | — | metadata 列名 → Harmony 整合（未装 harmony 自动回退） |
| `--markers` | — | 细胞类型 marker 集（第一列类型名 + 后续列为基因）→ score_genes 打分与自动注释 |

**输出**：QC 图、`umap_clusters`、`marker_genes.csv`（Wilcoxon 每簇 top 基因）、
`cluster_counts.csv`、`cell_scores.csv` + `cell_type_assignment.csv`（marker 打分）、
`paga_trajectory`、`processed.h5ad`（供下游工具直接读取）。

## 3. RSTA 资源与 RICE scGPT（重点）

**RSTA（Rice Spatio-Temporal Atlas）**：覆盖水稻"种子到种子"全生命周期的时空图谱
（78.8 万单核 / 47 组织 / 202 细胞类型 / 43 组织空间数据），配套 **RICE scGPT 单细胞基础模型**。

> ⚠️ RSTA 仓库（github.com/RiceSpatioTemporalAtlas/RSTA）**未附许可证文件**——本包不复制其代码，
> 采用"**用户运行时自行获取**"模式：按官方仓库/论文提供的入口下载与使用。

### 3.1 获取与本地运行（参考流程）

```bash
# 1) 克隆官方仓库（由用户自行执行；模型权重按仓库 Release/官网说明下载）
git clone https://github.com/RiceSpatioTemporalAtlas/RSTA.git
# 2) 按 01.Rice-scGPT_model/ 目录内 notebook 说明准备环境与模型权重
# 3) 典型用途（参考该目录 Tutorial_Annotation.ipynb 的调用方式）：
#    - 细胞类型注释（对缺乏 marker 的细胞群）
#    - 跨数据集整合（统一嵌入空间）
```

### 3.2 与本包流程的衔接

| 场景 | 做法 |
|---|---|
| 自动注释 | 本包 `--markers` 打分给出初步注释 → RICE scGPT 复核/补充无 marker 类型 |
| 参考映射 | 把自己的数据映射到 RSTA 参考图谱（按 scGPT 微调管线） |
| marker 参考 | 用 RSTA 的 202 细胞类型 marker 构建 `marker_sets.csv` 后直接喂给本包脚本 |

## 4. 下游分析（外部工具）

| 分析 | 工具 | 说明 |
|---|---|---|
| 细胞通讯 | **CellChat**（R）/ NicheNet | 基于配体-受体库（水稻可用拟南芥同源注释） |
| 拟时序/RNA velocity | scVelo / Monocle3 | `processed.h5ad` 直接输入 |
| 调控网络 | SCENIC / pySCENIC | TF regulon（PlantTFDB 水稻 TF 列表见 `assets/annotations/`） |
| 差异表达（伪批量） | 本包模块 1 | 按细胞类型聚合成伪批量 counts → DESeq2 流程 |

## 5. 水稻专用注意点

| 事项 | 说明 |
|---|---|
| 注释来自参考 | 水稻单细胞注释依赖参考图谱（RSTA、各组织论文）；无 marker 时优先 scGPT 辅助 |
| 原生质体影响 | 原生质体制备会诱导应激基因（如热激、氧化胁迫）——QC 时结合空液滴/应激分数排查 |
| 叶绿体/线粒体基因 | 叶绿体基因高表达细胞需过滤（如 rbcL/psbA 家族）；视组织调整 |
| 聚类精度 | 植物细胞类型差异常小于动物——resolution 0.6–1.2 间调参，结合 marker 稳定性 |
| 与转录组衔接 | 单细胞伪批量（per cell type）→ 模块 1 差异分析 → 与组织水平 RNA-seq 对照 |

> 演示数据：`python scripts/utils/make_demo_single_cell.py --outdir demo_data/single_cell`
> （1500 细胞 × 1000 基因 / 3 类 + marker 集 + 真值标签；需安装 scanpy 后跑主流程）。
