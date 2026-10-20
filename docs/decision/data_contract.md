# 数据契约（Data Contract）

> 目的：在跑任何脚本之前，用 10 个问题确定"这批数据能回答什么、不能回答什么"。
> 配套：`scripts/plan/01_make_plan.py`（读 JSON → 输出三档方案）。

## 一、10 个必答问题

| # | 字段 | 可选值 | 为什么重要 |
|---|---|---|---|
| 1 | `assay` | bulk_rna / metabolome / resequencing / rna_metabolome | 决定可用方法与统计模型 |
| 2 | `goal` | differential_expression / marker_discovery / gwas / prediction / integration / qtl | 决定主线与验收标准 |
| 3 | `design` | treatment_control / timecourse / multi_environment / natural_population / ril_magic_nam / case_control / paired | 决定随机效应与协变量 |
| 4 | `n_per_group` | 整数（GWAS 为群体大小） | 决定功效与能否做推断 |
| 5 | `replicates` | true / false | **无重复不能做推断统计** |
| 6 | `factors` | 处理、时间、地点、基因型等 | 设计公式的唯一来源 |
| 7 | `batch` | true / false | 批次与技术混杂会制造假差异 |
| 8 | `reference` | IRGSP-1.0 / 其他组装 + 注释版本 | 注释错版本会让候选基因全错 |
| 9 | `phenotype` | continuous / binary / none | 二分类需逻辑回归与更大样本 |
| 10 | `environments` | 整数 | ≥2 时才谈 G×E 与稳定性 |

## 二、跑之前的五项体检（不通过就不要往下走）

1. **样本匹配**：多组学按样本 ID 取交集，记录交集前后 n；
2. **批次结构**：批次 × 分组交叉表——若完全混杂，效应不可分离；
3. **重复与离群**：PCA/相关矩阵找离群样本，记录剔除理由；
4. **数据分布**：counts 是否整数、代谢物是否相对丰度（组成型）、蛋白是否 LFQ；
5. **注释版本**：基因 ID 与参考基因组必须一致（MSU/RAP-DB 需转换）。

## 三、三条红线

| 红线 | 后果 | 允许的替代 |
|---|---|---|
| 无生物学重复（n=1/组） | P 值不可解释 | 描述性分析 + 独立验证 |
| 批次与处理完全混杂 | 无法区分处理效应 | 仅报告相关性，不写因果 |
| GWAS 未校正群体结构 | 假阳性大量增加 | 先 PCA/亲缘校正，或改用家系设计 |

## 四、JSON 模板

```json
{
  "assay": "bulk_rna",
  "goal": "differential_expression",
  "design": "treatment_control",
  "n_per_group": 3,
  "replicates": true,
  "factors": ["treatment"],
  "batch": false,
  "reference": "IRGSP-1.0",
  "phenotype": "none",
  "environments": 1
}
```

运行：

```bash
python scripts/plan/01_make_plan.py --input my_data.json --outdir results/plan
python scripts/plan/01_make_plan.py --demo     # 4 个内置场景，含 n=1 红线
```
