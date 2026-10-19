# NOTICE — 第三方组件溯源说明（colors 科研绘图配色推荐器）

本目录包含从上游仓库复制的第三方程序源码与数据，版权归原作者所有，按上游 MIT 许可证使用。
**这是工具组件，不是技能库**（无 SKILL.md），与 bioskills 等技能库相互独立。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/lcpmgh/colors |
| 描述 | 科研绘图配色推荐器：R/Shiny 交互程序，含 229 套论文插图配色方案（顶刊风格 + R 包内置），支持双栏对比、透明度调节、RGB/HEX 信息查看 |
| 获取时间 | 2026-10-18 |
| 获取版本 | commit `36b2cc64994d1c463dc3c382c605d7cc63399162`（v4.0，shallow clone） |
| 许可证 | **MIT License, Copyright (c) 2025 Liang Chen**（LICENSE 原文已核验，全文见本目录 `LICENSE`） |
| 维护状态 | 活跃（末次提交 2026-06-26） |

## 已包含的内容

| 文件 | 说明 |
|---|---|
| `app.R` | Shiny 主程序（v4.0 双栏对比布局），未经改动 |
| `@palettes.csv` | 配色数据库：229 套方案，字段 `id,p_sour,p_name,c_numb,colors`；来源构成：知乎精选 177、ggsci 28、grDevices 16、RColorBrewer 8；每套 2–51 个颜色。**已知上游笔误**：第 124 行（id 123，知乎来源 `nameless` 6 色套）末色为 `#T3ABCF`——上游原文即如此（非本库转录错误），为保持逐字转录**未改动**；取用该套时请手动把第 6 色改为合法 hex（如 `#53ABCF`） |
| `README.md` | 上游说明（含配色来源文章链接与推荐工具） |
| `LICENSE` | 上游 MIT 许可证原文（未改动） |
| `preview.png`、`www/` | 上游预览图与 app 图标资源 |

## 明确未包含的内容

- `articles/`：上游备份的 8 篇知乎文章（.mhtml）与 1 个第三方网页（微生信 colors.html）。
  **这些是第三方版权内容，MIT 许可证不覆盖，不复制、不分发**；其配色已以上游整理后的
  hex 值形式体现在 `@palettes.csv` 中（事实数据），溯源链接见 `docs/references/plot_color_palettes.md`。
- `oldversion/v3.0.1`：已被 v4.0 取代的历史版本。
- `.git/`：仅取文件快照，不带版本历史。

## 使用说明

- **本地运行**：`shiny::runApp("third_party/colors")`（R 依赖见 `docs/references/plot_color_palettes.md` §4）。
- **查询式使用**：直接读 `@palettes.csv` 取 hex，配合 ggplot2 的 `scale_color_manual()` /
  `scale_fill_manual()` 用于本库各出图脚本；方法论与选型速查见
  `docs/references/plot_color_palettes.md`（自有语言蒸馏，非上游代码）。
- 二次分发本项目时，请保留本目录的 `LICENSE` 与 `NOTICE.md`。
