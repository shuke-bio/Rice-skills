# examples/metabolome_data/ — 代谢组丰度矩阵格式说明

> 本目录只提供格式说明与示例，不含真实数据。权威定义见
> [SKILL.md 第 2 节](../../SKILL.md)。

## 文件结构

- **行 = 代谢物**，每行一个代谢物。
- **第一列**必须是 `metabolite_id`（代谢物唯一 ID，如 `M0001`）。
- **其余列**：先是可选的注释列，然后是样本列（样本 ID 必须与
  `metadata.csv` 的 `sample_id` 完全一致，顺序不限）。
- 表头必须存在；文件为 UTF-8、逗号分隔。

**必需**：`metabolite_id` + 各样本列。
**推荐注释列**（按此顺序放在 ID 列之后）：

| 注释列 | 用途 |
|---|---|
| `metabolite_name` | 代谢物名称，仅用于图表展示与结果可读性 |
| `class` | 化学分类（黄酮/生物碱/脂质/氨基酸等），**供 `04_metabolite_classification.R` 做差异代谢物的分类组成统计** |
| `kegg_id` | KEGG compound ID（如 `C12271`），**供 `03_pathway.R` / `04_joint_pathway.R` 做通路富集**；无对应 KEGG 条目时留空 |

## 示例小表（4 代谢物 × 6 样本）

`metabolites.csv`：

| metabolite_id | metabolite_name | class | kegg_id | CTRL_1 | CTRL_2 | CTRL_3 | SALT_1 | SALT_2 | SALT_3 |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| M0001 | Flavonoid-3-glucoside | Flavonoids | C12271 | 15200 | 14100 | 15900 | 38900 | 41200 | 37600 |
| M0002 | Proline-betaine | Alkaloids | C00719 | 860 | 910 | 795 | 2540 | 2380 | 2610 |
| M0003 | Linolenic-acid | Lipids | C06427 | 3210 | 2980 |  | 1520 | 1400 | 1610 |
| M0004 | Unknown-compound-4 | Terpenoids |  | 45 | 52 | 38 | 44 | 51 | 40 |

对应的 CSV 原文（前两行）：

```csv
metabolite_id,metabolite_name,class,kegg_id,CTRL_1,CTRL_2,CTRL_3,SALT_1,SALT_2,SALT_3
M0001,Flavonoid-3-glucoside,Flavonoids,C12271,15200,14100,15900,38900,41200,37600
```

## 缺失值

- 代谢组矩阵**允许缺失值**，用**空单元格**表示（如示例中 `M0003` 的
  `CTRL_3`）。不要填 `0`、`NA`、`-` 或 `null` 等文本——空单元格才能被
  正确识别为缺失。
- 低丰度代谢物在部分样本中低于检测限而缺失，属正常现象；脚本会先报告整体
  缺失比例，再按半最小值或 KNN 填充（`01_preprocess.py`）。
- 若某代谢物缺失比例过高（如 > 50%），建议自行删除该行。

## 归一化

**无需手工归一化**。`scripts/metabolomics/01_preprocess.py` 依次完成：
缺失值汇总与填充 → 总和归一化 / log 变换 / Pareto 缩放（方法可调）→
异常样本检查（归一化前后 PCA 与箱线图）。提交**原始仪器响应峰面积或强度**
即可，不要提前 log 或标准化——否则归一化步骤会被重复施加。

## 与其他步骤的衔接

- `class` 列 → `04_metabolite_classification.R`：按化学分类统计差异代谢物
  分布（分类组成图）。
- `kegg_id` 列 → `03_pathway.R`（代谢通路富集）与 `04_joint_pathway.R`
  （基因–代谢物联合通路）：没有 `kegg_id` 的通路富集步骤会被跳过并提示。
- 预处理输出 → `02_oplsda.R`（OPLS-DA / VIP）、`05_biomarker.R`
  （ROC / 随机森林）、`06_metabolite_network.R`（相关性网络）。
