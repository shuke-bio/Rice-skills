# examples/ — 输入数据格式模板

> 本目录是 **数据格式模板**，不含任何真实实验数据。把三个文件的结构复制为自己的
> 数据即可直接进入分析流程；输入格式的权威定义见 [SKILL.md 第 2 节](../SKILL.md)。

## 目录内容

```
examples/
├── README.md                        # 本文件
├── metadata.csv                     # 样本元数据模板（6 个虚拟样本）
├── transcriptome_data/
│   └── README.md                    # 转录组 count 矩阵格式说明（含示例小表）
└── metabolome_data/
    └── README.md                    # 代谢组丰度矩阵格式说明（含示例小表）
```

## 三个标准输入（要点速览）

所有输入均为 **UTF-8 编码、逗号分隔（CSV）、首行表头**。详细列定义与常见错误见
各子目录 README；此处只列最容易出错的点。

| 文件 | 结构 | 硬性要求 |
|---|---|---|
| `counts.csv`（转录组） | 行 = 基因，列 = 样本 | 整数计数；**不允许缺失值**；基因 ID 唯一 |
| `metabolites.csv`（代谢组） | 行 = 代谢物，列 = 样本 | 允许空单元格（缺失值）；推荐 `metabolite_name`/`class`/`kegg_id` 注释列 |
| `metadata.csv`（元数据） | 首列 `sample_id` | **必含 `group` 列**；样本 ID 与两个矩阵的列名**完全一致**（大小写敏感、无空格） |

样本 ID 在三个文件之间必须严格一致，这是全流程自动对齐的前提；矩阵的样本列
顺序无需与元数据一致（脚本内部按 ID 对齐）。缺失值规则：count 矩阵不允许，
代谢组矩阵允许（由 `scripts/metabolomics/01_preprocess.py` 填充）。

## 用示例数据快速跑通流程

仓库不附真实数据，但提供确定性示例数据生成器（seed=42，可复现）：

```bash
# 在仓库根目录执行
python scripts/utils/make_demo_data.py --outdir demo_data
```

生成 `demo_data/` 下的三个文件（`counts.csv` 2000 基因、`metabolites.csv` 300
代谢物含缺失值、`metadata.csv` 6 样本 control/salt 两组、带批次效应与部分
差异基因/代谢物）。随后即可运行：

```bash
python scripts/transcriptomics/01_qc_filter.py \
    --counts demo_data/counts.csv --metadata demo_data/metadata.csv
python scripts/metabolomics/01_preprocess.py \
    --abundance demo_data/metabolites.csv --metadata demo_data/metadata.csv
```

或一条命令跑完整流程：

```bash
bash run_pipeline.sh --counts demo_data/counts.csv \
    --metadata demo_data/metadata.csv \
    --metabolites demo_data/metabolites.csv \
    --outdir results --stages all
```

> 注意：示例数据的样本 ID 由组名前 3 个字母生成（`control → CON_1`、
> `salt → SAL_1`），与本模板的 `CTRL_/SALT_` 不同——这只是命名习惯差异，
> 格式完全一致。示例数据仅供验证流程可运行，**不得用于任何生物学结论**。

## 换成自己的数据：检查清单

提交数据前逐项确认（任何一项不满足，脚本都会给出明确报错或结果不可信）：

- [ ] **样本 ID 一致性**：`metadata.csv` 的 `sample_id` 与 count/代谢物矩阵的
      列名逐一对应（大小写敏感、无空格、无多余后缀）；建议复制本模板的
      `metadata.csv` 改写而非从零新建。
- [ ] **group 列**：元数据必含 `group` 列，取值简洁明确（如 `control`/`salt`）；
      多组比较时每组名称固定，后续对比（如 `salt_vs_control`）以组名为准。
- [ ] **生物学重复 ≥ 3**：每组至少 3 个独立生物学重复（2 重复无法估计离散度，
      差异表达统计不可靠）；检查是否有技术重复被误当成生物学重复。
- [ ] **count 为整数**：count 矩阵全部为非负整数，无小数、无科学计数法文本、
      无千分位逗号；表头除基因 ID 列外全部是样本 ID。
- [ ] **代谢物缺失值处理**：缺失值留空单元格即可（不要填 0 或 "NA" 之外的
      文本）；预处理脚本会自动做缺失值汇总与填充（半最小值/KNN）及
      归一化（总和/log/Pareto），**无需手工归一化**；低丰度代谢物缺失多属
      正常现象，脚本会报告缺失比例。
- [ ] **文件编码与分隔符**：均为 UTF-8 + 逗号；Excel 另存为 CSV 时选择
      "CSV UTF-8"，避免 GBK 编码和分号分隔。
- [ ] **基因/代谢物 ID 唯一**：无重复行 ID；代谢物 ID 与 `kegg_id` 的对应
      关系正确（`kegg_id` 用于通路富集，错误的映射会静默影响富集结果）。

## 相关文档

- [SKILL.md 第 2 节](../SKILL.md) — 输入数据格式权威定义
- `transcriptome_data/README.md` — count 矩阵详解与 DESeq2 衔接说明
- `metabolome_data/README.md` — 代谢组矩阵详解与注释列用途
- `../docs/usage.md` — 各脚本命令行参数
