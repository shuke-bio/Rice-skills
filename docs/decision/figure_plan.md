# 分析方案与图形方案确认（Plan & Figure Approval）

> 硬性规则：**在用户明确确认"分析方案 + 图形方案"之前，不执行任何分析脚本**。
> 本文件定义确认单包含什么、图形规格、以及确认与执行的记录方式。

## 1. 确认单必须包含

1. **数据契约摘要**：assay / goal / design / n / 重复 / 批次；
2. **阻断项**（若有）：先解决，不允许带阻断执行；
3. **风险提示**：样本量、批次混杂、注释版本等；
4. **分析步骤**：逐条列出脚本与目的；
5. **图形方案**：每张图的名称、类型、规格、生成脚本；
6. **预期产出**：表格、图、报告清单。

生成：

```bash
python scripts/plan/01_make_plan.py --input data.json --outdir results/plan
python scripts/plan/05_run_plan.py --plan results/plan/plan.json --tier minimal
```

## 2. 图形规格（默认 CNS 风格）

| 项 | 默认值 |
|---|---|
| 宽度 | 单栏 85 mm / 双栏 170 mm |
| 分辨率 | ≥300 dpi（位图）；矢量优先 |
| 格式 | PDF + SVG + TIFF |
| 字体 | Arial / Helvetica，7–9 pt |
| 配色 | `config/nature_palette.json`（低饱和、色盲友好） |
| 坐标轴 | 完整刻度、单位、样本量标注 |
| 统计标注 | 检验方法、P/FDR、效应量（按需） |

图名与类型由方案自动给出；具体样式由各脚本内置的绘图规范控制，无需另调。

## 3. 确认与执行记录

| 文件 | 作用 |
|---|---|
| `approval.json` | 记录方案哈希、档位、图形数、确认时间 |
| `run.log` | 执行日志（跳过/失败原因逐条记录） |

```bash
# 用户确认后
python scripts/plan/05_run_plan.py --plan results/plan/plan.json --tier minimal --approve
# 再执行（会先做依赖体检；缺依赖则停止）
python scripts/plan/05_run_plan.py --plan results/plan/plan.json --tier minimal --execute
```

保护机制：

- 没有 `approval.json` → 拒绝执行；
- 方案文件被修改（哈希不一致）或档位不一致 → 拒绝执行；
- 有阻断项 → 拒绝执行；
- 依赖未就绪 → 停止并指向 `deps_report.md`；
- 不支持 `--demo` 的脚本 → 跳过并提示需真实参数，不伪造执行。

## 4. 为什么必须先确认

- 分析方案决定"结论能不能成立"（如无重复不能做推断）；
- 图形方案决定"结果能不能发表"（尺寸/字体/配色不合规会返工）；
- 先确认再执行，避免跑了几小时才发现设计或图表方向不对。
