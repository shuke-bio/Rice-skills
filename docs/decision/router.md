# 路由表：用户一句话 → 该问什么 → 用哪张卡

> 目标：用户说"我有批数据想分析"时，**先补字段、再出方案**，不要直接跑脚本。

## 1. 意图识别

| 用户说法 | 先补问 | 场景卡 | 生成器路由 |
|---|---|---|---|
| "转录组差异表达 / 找差异基因" | 分组、每组 n、是否有重复、批次、注释版本 | `scenario_rnaseq.md` | `bulk_rna` + `differential_expression` |
| "GWAS 找位点 / 关联分析" | 群体大小、表型类型、群体结构、参考基因组、环境数 | `scenario_gwas.md` | `resequencing` + `gwas` |
| "QTL / BSA / 混池测序" | 家系类型与规模、表型、是否有亲本 | `scenario_qtl.md` | `goal=qtl` 或 `design=ril_magic_nam` |
| "基因组选择 / 预测品种" | 表型环境数、基因型、交叉验证设计 | `scenario_gs.md` | `goal=prediction` |
| "单细胞 / 空间" | 样本数、组织、批次、是否有参考图谱 | `scenario_singlecell.md` | `single_cell` / `spatial` |
| "代谢组差异 / 标志物" | 靶向或非靶向、样本量、缺失率 | `scenario_metabolome.md` | `metabolome` + `marker_discovery` |
| "转录组 + 代谢组联合" | 样本能否匹配、n、批次 | `scenario_multiomics.md` | `rna_metabolome` + `integration` |
| "甲基化 / ATAC / ChIP" | 组学类型、重复、覆盖度 | `scenario_epigenome.md` | `wgbs` / `atac` |
| "蛋白组差异" | 定量方式、缺失率、重复 | `scenario_proteome.md` | `proteome` |
| "泛基因组 / SV / PAV" | 组装数、注释是否统一、参考基因组 | `scenario_pangenome.md` | `pangenome` |

## 2. 标准流程（agent 应固定执行）

1. **补齐字段**：按 `data_contract.md` 问 10 项；能从上下文推断的直接填，不要重复问；
2. **数据体检**（有文件时）：
   `python scripts/plan/02_check_data.py --matrix ... --metadata ... --outdir results/qc`
3. **生成方案**：
   `python scripts/plan/01_make_plan.py --input data.json --outdir results/plan`
4. **先看阻断项**：有 `blockers` 就停下，向用户说明原因并给替代方案；
5. **让用户选档**：最小可行 / 标准 / 进阶；
6. **执行**：先 `--demo` 验证环境，再传真实参数（见 `docs/*_usage.md`）。

## 3. 不该做的事

- 用户只说"帮我分析"就直接跑脚本——缺少设计字段，结果不可解释；
- 有阻断项仍输出 P 值或"显著"结论；
- 把场景卡的示例数值当成用户数据的结论；
- 引用非白名单期刊作为方法依据（见 `README.md` 的准入标准）。
