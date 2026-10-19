# NOTICE — 第三方组件溯源说明（colormaps 连续科学色图）

本目录收录**连续型科学色图**（sequential / diverging / cyclic / rainbow 感知均匀色图，共 133 套）的
蒸馏数据与许可证原文。与 `colors/`（离散配色 + Shiny 工具）互补：colors 管"分类变量配色"，
本目录管"连续变量色图"（热图、等值线、密度图、梯度图）。
**这是数据组件，不是技能库**（无 SKILL.md）。

**收录准则（2026-10-18 用户审定）**：只收**仍在活跃维护**的上游（末次提交 ≤1 年内），
以许可证原文核验为前提。长期停更库不收录。

## 上游信息与许可核验

所有许可证均**打开原文核验**（LICENSE / pyproject），非仅凭 GitHub 标签：

| 家族 | 上游 | 获取版本 | 许可证（原文已核验） | 末次更新 | 规模 |
|---|---|---|---|---|---|
| cmasher | https://github.com/1313e/CMasher（默认分支 master） | commit `4360a56`（v1.9.4） | **BSD-3-Clause, © 2019-2021 Ellert van der Velden**（`LICENSES/CMasher-BSD-3-Clause.txt`） | 2026-10-09 | 53 套 |
| colorcet | https://github.com/holoviz/colorcet | commit `076eb84` | **CC BY 4.0, © Peter Kovesi**（`LICENSES/colorcet-CC-BY-4.0.txt`）；再分发须署名 | 2026-10-12 | 80 套 |

## 已包含的内容

| 文件 | 说明 |
|---|---|
| `@colormaps.csv` | 133 套连续色图的 hex 数据，字段 `c_sour,c_name,c_type,n_stops,colors`（`;` 分隔 hex 列表，个数 = `n_stops`）。色阶数如实转录上游文件（116 套 256 级、cmasher 循环/发散类 510–511 级、cmasher `neutral` 2560 级），未做重采样。`c_type` 取值为 sequential / diverging / cyclic / rainbow，取自上游官方分类：cmasher 用其 `_known_cmap_types.py` 类型注册表；colorcet 按其官方 CSV 类名（cyclic-/diverging-/rainbow_/linear-） |
| `LICENSES/` | 上游许可证原文（未改动） |
| `README.md` | 数据 Schema 与快速用法 |

## 明确未包含的内容

- 各上游**源代码、构建文件、测试、文档站点、示例图、CVD 诊断图**：仅取许可证与色图数据本身
  （hex 值为事实数据；颜色数值本身不受版权保护，受保护的是代码与表达）。
- CMasher 各色图的 `_viscm.png`/`.png` 预览图与 `app_data.txt.gz`。
- colorcet 的 `_s25.csv`（25 级降采样重复文件）与 `assets/Glasbey`（256 色分类图，属离散分类用途，
  与 colors 组件定位重叠）。
- `.git/`：仅取文件快照，不带版本历史。

## 调研后未收录的上游（长期停更，按收录准则剔除）

| 家族 | 断更时间 | 说明 |
|---|---|---|
| crameri（Scientific colour maps，batlow/vik/broc 等） | Zenodo v8.0.0 发布于 2023-06-14（3 年无新版） | 上游 GitHub 仓库与作者账户均已 404、官方网站本机不可达、项目转向私有资助。若需该家族，用 `remotes::install_github("thomasp85/scico")` 或 `pip install cmcrameri` 运行时安装，引用 Crameri et al. 2020 (Nat Commun 11, 5444) |
| cmocean（balance/thermal/phase 等） | 2024-04-05（约 2.5 年） | `pip install cmocean`；引用 Thyng et al. 2016 |
| viridis（viridis/magma/plasma 等） | 2024-01-28（约 2.8 年） | matplotlib/seaborn 内置零成本；R 用 `install.packages("viridis")` |
| scico（Crameri 家族 R 移植） | 2023-08-14（3 年） | `install.packages("scico")` |

## 使用说明

- **查询式使用（推荐）**：读 `@colormaps.csv` 取 hex 数组，用 ggplot2 `scale_fill_gradientn(colours=)` /
  matplotlib `LinearSegmentedColormap.from_list()` 出图；方法笔记见
  `docs/references/plot_color_palettes.md` §8（含离散/连续选型决策与代码）。
- **完整功能**（交互浏览器、其他格式）：`pip install cmasher` / `pip install colorcet`。
- **引用口径**：CMasher 请引 van der Velden (2021, JOSS)；colorcet 请引 Kovesi (2015)。
- 二次分发本项目时，请保留本目录全部 `LICENSES/` 与 `NOTICE.md`（colorcet 的 CC BY 4.0 要求署名）。
