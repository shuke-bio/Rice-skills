# NOTICE — 第三方组件溯源说明（ltc 艺术/自然风离散配色）

本目录收录 **ltc** R 包（CRAN: ltc）的许可证、描述与调色板数据。
**这是数据组件，不是技能库**（无 SKILL.md），与技能库相互独立。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/loukesio/ltc-color-palettes |
| 描述 | ltc：Artistic and Nature-Inspired Color Palettes——艺术与自然灵感的离散配色，每套配色带背景故事 |
| 获取时间 | 2026-10-18 |
| 获取版本 | commit `6e69caa`（v0.4.0.9000，shallow clone，默认分支 main） |
| 许可证 | **MIT License, Copyright (c) 2021 ltc authors**（LICENSE.md 原文已核验，全文见本目录；R 惯例 `MIT + file LICENSE`） |
| 维护状态 | 活跃 |

## 已包含的内容

| 文件 | 说明 |
|---|---|
| `@palettes.csv` | 32 套配色 hex 数据，字段 `id,p_sour,p_name,c_numb,colors`；从上游 `R/ltc_functions.R` 的 `palettes` 列表逐字转录（含大写归一化），未改动色值 |
| `LICENSE.md` | 上游 MIT 许可证原文（未改动） |
| `DESCRIPTION` | 上游包元数据（作者 Loukas Theodosiou 等，未改动） |
| `README.md` | 上游说明（未改动） |

## 明确未包含的内容

- 上游 R 包**源码其余部分**（`R/` 其余文件、`data-raw/` 画廊生成脚本、Docker 相关文件）：
  仅取许可证、包描述、说明与调色板数据本身（hex 值为事实数据）。
- `.git/`：仅取文件快照，不带版本历史。

## 使用说明

- **查询式使用**：读 `@palettes.csv` 取 hex，配合 ggplot2 `scale_color_manual()` /
  `scale_fill_manual()`；选型方法论见 `docs/references/plot_color_palettes.md`。
- **完整功能**（backstory 文档、色盲模拟对比）：`install.packages("ltc")`。
- 二次分发本项目时，请保留本目录的 `LICENSE.md` 与 `NOTICE.md`。
