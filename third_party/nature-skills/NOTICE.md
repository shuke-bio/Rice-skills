# NOTICE — 第三方组件溯源说明

本目录包含从上游仓库复制的第三方 skill 源码，版权归原作者所有，按上游许可证使用。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/xiaolai/nature-skills |
| 获取时间 | 2026-10-14（两次：v2.1.1 时 shallow clone commit `c1ba255bedf3eefb9c336b7ae20dfe6973c157be`；v3.1 时经镜像获取 `main` 分支 tarball，无 git 元数据） |
| 许可证 | MIT License, Copyright (c) 2026 Yuan Yizhe（全文见本目录 `LICENSE`） |

## 已包含的内容（按需选择，仅复制源码与参考资料）

| 目录 | 上游对应 | 用途 |
|---|---|---|
| `skills/nature-figure/` | `skills/nature-figure/` | Nature 级图表规范与 matplotlib 模板（SKILL.md、README.md、references/、evals/） |
| `skills/nature-polishing/` | `skills/nature-polishing/` | 英文论文润色规范（25 条规则、12 步工作流） |
| `skills/nature-writing/` | `skills/nature-writing/` | 论文各章节写作指南（摘要/引言/方法/结果/讨论/结论等） |
| `skills/nature-academic-search/` | 同名 | 学术检索（PubMed/MEDLINE 查询、DOI 提取、RIS/BibTeX 导出） |
| `skills/nature-citation/` | 同名 | 引用格式转换与规范 |
| `skills/nature-data/` | 同名 | 数据呈现规范 |
| `skills/nature-paper2ppt/` | 同名 | 论文转 Nature 风格中文汇报 PPTX |
| `skills/nature-reader/` | 同名 | 全文中英对照精读辅助（图表感知） |
| `skills/nature-response/` | 同名 | 审稿意见逐条回复信（response to reviewers）工作流 |

## 明确未包含的内容

- **`assets/` 目录（约 30 MB）**：包括 `assets/figures4papers/`（来自已发表论文的演示脚本与预览图，第三方素材）
  与 `assets/chart-atlas/`（图集预览 PNG）。出于仓库体积与第三方素材许可审慎考虑未复制；
  如需原始演示资产，请从上游仓库获取。
- 上游仓库的 `.claude-plugin/`、CI 配置等仓库基础设施未包含。

## 使用说明

- 这些文件以"参考规范与模板库"形式提供：RiceMultiOmics 自身的分析脚本使用自研绘图主题
  （`config/nature_palette.json` + `scripts/utils/plot_theme.R` / `plot_theme.py`），与上游代码相互独立。
- 需要 Nature 级图表规范、英文润色、论文写作、文献检索、审稿回复或汇报 PPT 辅助时，
  加载对应的 `SKILL.md` 及其 `references/`。
- 如二次分发本项目，请保留本目录的 `LICENSE` 与 `NOTICE.md`。
