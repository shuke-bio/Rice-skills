# NOTICE — 第三方组件溯源说明（ggpubfigs 色盲友好配色 + 出版主题）

本目录收录 **ggpubfigs** R 包的许可证、描述与调色板数据。
**这是数据组件，不是技能库**（无 SKILL.md），与技能库相互独立。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/JLSteenwyk/ggpubfigs |
| 描述 | 面向发表级插图的 ggplot2 简易主题 + 色盲友好配色（bright_seven、tol_eight、tableau_ten、glasbey_twelve 等 15 套） |
| 获取时间 | 2026-10-18 |
| 获取版本 | commit `7cf2193`（v1.1.4，shallow clone，默认分支 master） |
| 许可证 | **MIT License, Copyright (c) 2020 Jacob L. Steenwyk**（LICENSE.md 原文已核验，全文见本目录） |
| 维护状态 | 活跃 |

## 已包含的内容

| 文件 | 说明 |
|---|---|
| `@palettes.csv` | `friendly_pals` 15 套配色 hex 数据，字段 `id,p_sour,p_name,c_numb,colors`；从上游 `R/colors.R` 逐字转录，未改动色值。上游注明这些配色整理自 Paul Tol/SRON、davidmathlogic、IBM 等公开来源（事实数据） |
| `LICENSE.md` | 上游 MIT 许可证原文（未改动） |
| `DESCRIPTION` | 上游包元数据（作者 Jacob L. Steenwyk，未改动） |
| `README.md` | 上游说明（未改动） |

## 明确未包含的内容

- 上游 ggplot2 **主题函数源码**（`R/theme_*.R`）与文档站点：仅取许可证、包描述、
  说明与调色板数据本身；主题属代码表达，需要时按 §使用说明 安装上游包。
- `.git/`：仅取文件快照，不带版本历史。

## 使用说明

- **查询式使用**：读 `@palettes.csv` 取 hex；分类变量配色优先用 `tol_eight`（Paul Tol，
  色盲安全）、`tableau_ten`（Tableau 10）、`glasbey_twelve`（多类别）等；
  选型方法论见 `docs/references/plot_color_palettes.md`。
- **完整功能**（`ggpubfigs::npg()`、`theme_grid()` 等 ggplot2 扩展）：`install.packages("ggpubfigs")`。
- 二次分发本项目时，请保留本目录的 `LICENSE.md` 与 `NOTICE.md`。
