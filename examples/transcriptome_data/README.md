# examples/transcriptome_data/ — 转录组 count 矩阵格式说明

> 本目录只提供格式说明与示例，不含真实数据。权威定义见
> [SKILL.md 第 2 节](../../SKILL.md)。

## 文件结构

- **行 = 基因**，每行一个基因，行数 = 基因数（几千到几万）。
- **列 = 样本**，每列一个样本，列数 = 样本数。
- **第一列**是基因 ID（列名习惯写作 `gene_id`，脚本会自动识别首列为 ID 列）。
- **其余列**是样本 ID，必须与 `metadata.csv` 的 `sample_id` 完全一致
  （大小写敏感、无空格）；列顺序无需与元数据一致，脚本内部按 ID 对齐。
- 表头（第一行）必须存在；文件为 UTF-8、逗号分隔。

## 示例小表（5 基因 × 6 样本）

`counts.csv`：

| gene_id | CTRL_1 | CTRL_2 | CTRL_3 | SALT_1 | SALT_2 | SALT_3 |
|---|---:|---:|---:|---:|---:|---:|
| LOC_Os01g01010 | 1520 | 1380 | 1602 | 4210 | 3980 | 4150 |
| LOC_Os01g01019 | 88 | 102 | 76 | 231 | 205 | 198 |
| LOC_Os02g02100 | 0 | 1 | 0 | 0 | 0 | 2 |
| LOC_Os03g02870 | 4521 | 4890 | 4312 | 398 | 445 | 376 |
| LOC_Os05g01120 | 12 | 9 | 15 | 9 | 11 | 8 |

对应的 CSV 原文（前两行）：

```csv
gene_id,CTRL_1,CTRL_2,CTRL_3,SALT_1,SALT_2,SALT_3
LOC_Os01g01010,1520,1380,1602,4210,3980,4150
```

基因 ID 使用水稻 MSU/RAP 命名（`LOC_Os` 格式）；如使用其他 ID 体系
（如 RAP ID、基因符号），全流程保持统一即可，富集步骤需要 ID 能映射到
`org.Osativa.eg.db`。

## 数值要求

- **整数计数**：全部为非负整数（0 允许）。不要小数、科学计数法、千分位
  逗号或引号。
- **不允许缺失值**：count 矩阵不能有空单元格。原始数据中未检出的基因写
  `0`，不要留空。
- 低表达基因（如示例中 `LOC_Os02g02100`）无需手动删除，QC 步骤会按
  `--min-counts 10 --min-samples 3` 过滤。

## 常见错误清单

| 错误 | 现象 / 后果 | 正确做法 |
|---|---|---|
| 重复基因 ID | 脚本报 "duplicate gene_id" 并中止 | 去重（合并或删除重复行）后再提交 |
| 非整数值（小数、`1.2e3`、带引号数字） | 读取报错或被当作文本 | 全部改为非负整数 |
| 样本名与元数据不一致（`ctrl_1` vs `CTRL_1`、多余空格、`.1` 后缀） | 样本对齐失败或样本被丢弃 | 严格复制 `metadata.csv` 的 `sample_id` |
| 逗号/制表符混乱（Excel 存成 TSV、分号分隔） | 解析失败或只有一列 | 统一 UTF-8 + 逗号分隔；Excel 另存为选 "CSV UTF-8" |
| 表头缺失或首列无列名 | 首行被当成数据 | 保留表头；首列列名可为空但建议写 `gene_id` |
| count 矩阵含 TPM/FPKM 等连续值 | 差异表达模型不适用 | DESeq2 必须用原始整数 count |

## 与 DESeq2 输入的衔接

`scripts/transcriptomics/02_deseq2.R` 直接消费本格式（QC 过滤后的
`counts_filtered.csv` 亦为本格式），内部处理：

1. 以首列为行名（基因 ID），其余列为样本计数矩阵；
2. 按 `metadata.csv` 的 `sample_id` 对齐样本顺序；
3. 构建 `DESeqDataSet`：设计公式默认为 `~ group`，传 `--design "~ batch + group"`
   可校正批次；分组由 `--contrast salt_vs_control` 指定（对照组在前）；
4. 输出 `de_all_*.csv`（含 `log2FoldChange`/`pvalue`/`padj`/`significant` 列，
   可直接用于绘图模板 `template_volcano()` / `template_ma()`）与
   VST 归一化矩阵（用于热图）。

因此：**不要在提交前做任何归一化**（CPM/TPM/FPKM 均不适合 DESeq2），
也不要预先做差异筛选——归一化与检验由脚本完成。
