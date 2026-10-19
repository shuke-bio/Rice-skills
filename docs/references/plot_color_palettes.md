# 论文插图配色方法论（colors 蒸馏笔记）

> **定位**：本文是方法笔记（自有语言撰写），不复制上游代码。配套组件：
> - `third_party/colors/`（lcpmgh/colors，v4.0，MIT，Copyright (c) 2025 Liang Chen，
>   溯源见该目录 NOTICE.md）；配色数据源：`third_party/colors/@palettes.csv`（229 套离散）
> - `third_party/colormaps/`（CMasher BSD-3-Clause + colorcet CC BY 4.0）；
>   `@colormaps.csv`（133 张连续色图，见 §8）
> - `third_party/ltc/`（MIT，32 套自然风离散）、`third_party/ggpubfigs/`（MIT，
>   15 套色盲友好离散）
>
> 如需交互式挑选，本地运行 Shiny app；日常更快的方式是按下文直接查表。

## 1. 适用场景

- 为本库所有 ggplot2 出版图（模块默认 300 dpi CNS 级输出）选择**离散分组色**
-  VOLCANO/MA/热图/散点/箱线/富集条形图等需要 3–10 个可区分颜色时的取色
-  投稿前统一全文图表色系（同一项目固定 1–2 套 palette，避免一栏一色）
-  Quick-look 的连续型数据（表达量、pH、海拔）选 sequential/diverging 色系

## 2. 配色方案的三个类型（先分类，再选色）

| 类型 | 用途 | 判据 | 示例 |
|---|---|---|---|
| **离散 categorical** | 分组、样品、通路、物种 | 色相彼此拉开、明度接近、数量=组数 | NPG、Lancet、JCO、D3 |
| **连续 sequential** | 数值大小（0→max） | 单色相由浅入深或感知均匀 | viridis 系、 Blues/YlOrRd |
| **发散 diverging** | 有中点意义的数值（log2FC±、差异±） | 两端对比色经中性色过渡 | RdBu、PiYG、BrBG |

顶刊插图的"高级感"主要来自三点：**低饱和、明度一致、灰白底**。§2–§7 只覆盖离散型
（论文中使用频率最高）；连续/发散色图的方法论与查表见 **§8**（`colormaps` 组件）。

## 3. 配色数据库概况（@palettes.csv，229 套）

- **字段**：`id, p_sour, p_name, c_numb, colors`（`colors` 为逗号分隔 hex 列表）
- **来源构成**：知乎顶刊配色文章精选 177 套；R 包内置 52 套（ggsci 28、grDevices 16、
  RColorBrewer 8）
- **色数分布**：每套 2–51 色，主力区间 4–10 色（离散分组最常用段）
- **筛选方法**：按 `c_numb == 实际组数` 过滤，再目测挑一套；组数 >10 时优先合并类目
  或换连续型编码，不要硬上 10 色以上离散 palette

### 3.1 常用期刊风格 palette（摘自 CSV，ggsci 来源）

| 名称 | 色数 | Hex（按使用顺序） |
|---|---|---|
| NPG（Nature Publishing Group） | 10 | `#E64B35 #4DBBD5 #00A087 #3C5488 #F39B7F #8491B4 #91D1C2 #DC0000 #7E6148 #B09C85` |
| Lancet | 9 | `#00468B #ED0000 #42B540 #0099B4 #925E9F #FDAF91 #AD002A #ADB6B6 #1B1919` |
| NEJM | 8 | `#BC3C29 #0072B5 #E18727 #20854E #7876B1 #6F99AD #FFDC91 #EE4C97` |
| JAMA | 7 | `#374E55 #DF8F44 #00A1D5 #B24745 #79AF97 #6A6599 #80796B` |
| JCO（J Clin Oncol） | 10 | `#0073C2 #EFC000 #868686 #CD534C #7AA6DC #003C67 #8F7700 #3B3B3B #A73030 #4A6990` |
| AAAS/Science | 10 | `#3B4992 #EE0000 #008B45 #631879 #008280 #BB0021 #5F559B #A20056 #808180 #1B1919` |
| Tron（浅色底电影感） | 7 | `#FF410D #6EE2FF #F7C530 #95CC5E #D0DFE6 #F79D1E #748AA6` |

> 使用惯例：按组数取前 N 色；同页图表保持同一套 palette 与同一顺序（图例顺序 =
> hex 顺序），全文统一。

## 4. 三条使用路径

### 4.1 查表法（推荐日常默认，零依赖）

```r
# 从 CSV 取一套 palette 直接用于 ggplot2
pal <- "third_party/colors/@palettes.csv"
df  <- read.csv(pal, stringsAsFactors = FALSE)
hex <- strsplit(df$colors[df$p_name == "pal_npg"][1], ",", fixed = TRUE)[[1]]
hex <- trimws(hex)                       # 去空格
ggplot(df2, aes(x, y, color = group)) +
  geom_point(size = 2) +
  scale_color_manual(values = head(hex, length(unique(df2$group)))) +
  theme_classic(base_size = 12)          # 顶刊默认白底无网格
```

填充色用 `scale_fill_manual(values = hex)`；`alpha` 统一 0.8–0.9（出版图忌全透明）。

### 4.2 本地 Shiny app（交互挑选、双栏对比）

```r
# 一次性安装依赖（R >= 4.1）
install.packages(c("shiny", "ggplot2", "stringr", "magrittr", "data.table",
                   "shinyWidgets", "colourpicker", "DT"))
# 运行（在本仓库根目录）
shiny::runApp("third_party/colors")
```

功能：按色数筛选 → 双栏对比 → 查看 RGB/HEX 与透明度效果。挑中后记下名字，回 §4.1 取用。

### 4.3 程序化生成（连续/发散、或需更多控制时）

- **连续/发散色图（本包默认）**：查 `third_party/colormaps/@colormaps.csv`
  （133 张感知均匀色图，选型与代码见 §8）
- `ggsci::scale_color_npg()` 等——与 §3.1 同源，适合快速调用
- `rcartocolor`（CARTO 系列）、`scico`（Scientific colour maps，感知均匀，含色盲友好变体）
- `RColorBrewer` 的 `scale_*_brewer(palette = "Set2"/"Dark2", direction = -1)`

## 5. 选色原则（论文场景）

1. **组数 ≤ 8 用离散 palette；>8 合并类目**，或改为"高亮一组 + 其余灰阶"的叙事配色
2. **相邻组避开红-绿相邻**（最常见色盲障碍）；NPG/JCO 前几色已规避，直接沿用最稳
3. **灰度打印可辨**：选中的 palette 转灰度后明度应有层次（Lancet/JAMA 表现好；
   可以先 `colorspace::desaturate()` 自查）
4. **不靠颜色单独承载信息**：热图等务必配数值标注或形状辅助
5. **全文一致**：一套 palette 走到底；子图需强调时用加粗/边框/星标，而非另换色系
6. **白底黑字灰网格**（或 classictheme 去网格）是顶刊默认，palette 只服务数据元素

## 6. 与 Rice-skills 出图模块的衔接

- 各模块 ggplot2 脚本（转录组/代谢组/联合/GWAS/群体遗传/表型组等）均支持在
  `scale_color_manual()/scale_fill_manual()` 处注入 hex；建议在自己的项目级主题文件
  （或 .Rprofile 片段）固定一套 palette 常量全文复用——本包 R 侧共用
  `scripts/utils/plot_theme.R` + `config/nature_palette.json`，Python 侧对应
  `scripts/utils/plot_theme.py`
- 需要 CNS 级排版时配合：白底、Arial/Helvetica 或按要求字体（见默认出图规范）、
  300 dpi、色块图例置于轴外侧

## 7. 许可与引用

- 组件：MIT（Copyright (c) 2025 Liang Chen），保留 `third_party/colors/LICENSE` 与 NOTICE.md
- 知乎来源配色的原作者链接见上游 `README.md`（本库不复制其文章正文）
- 若论文方法学部分需致谢工具，建议引用：lcpmgh/colors (https://github.com/lcpmgh/colors) 及
  上述知乎文章/`ggsci`（Nan Xiao, https://github.com/nanxstats/ggsci）

## 8. 连续色图方法论（colormaps 组件）

### 8.1 五个类型：先分类，再选色

连续标量（表达量、log2FC、相关、角度）上色与离散分组是两套逻辑，核心是先判断数据是
**单向、双向有中心、还是循环**：

| 类型 | 数据特征 | 典型场景 | 关键要求 |
|---|---|---|---|
| **sequential** | 单向递增（0→max） | 表达量/丰度/覆盖度/深度的热图、密度图 | 感知均匀，浅端可读 |
| **diverging** | 双向、有中点意义（0 或均值） | log2FC 火山图/热图、相关系数、偏差 | **中心浅、两端深**，两臂明度对称 |
| **cyclic** | 循环量，无起点终点 | 风向/角度、相位、24 h 节律、圆周树图 | **首尾同色**，环绕处无色带断崖 |
| **rainbow** | 连续但为视觉风格服务 | 一般数量编码**应避免**（jet 式视觉畸变）；仅用于非量化自然变量或艺术化呈现 | 慎用，见 8.4 |
| **special** | 特殊物理量/地形专用 | 海水层化、地形高程、黑洞径向速度图等 | 本组件不收录，按 §8.5 运行时安装 |

本组件（`third_party/colormaps/@colormaps.csv`，133 张 = CMasher 53 + colorcet 80）
只覆盖前四类；色阶数如实转录上游文件（116 张 256 级，cmasher 循环/发散类 510–511 级，
`neutral` 2560 级，未重采样），全部经上游感知均匀性与色盲友好（CVD）诊断。

### 8.2 类型怎么判：数据驱动 + 上游官方分类交叉校验

对上游给出的每张色图，判定流程（本库建库时实际使用）：

1. **首尾色距**：计算 hex 首尾两点在 RGB 空间的归一化距离。**< 0.05 判 cyclic**
   ——循环色图首尾必须同色（如 colorcet `cyclic_*`、cmasher `seasons`/`copper`）。
   阈值不能放宽：深色发散色图两端也可能接近（实测 cmasher 类发散图曾被误判），
   这类必须靠第 3 步官方分类兜底。
2. **亮度剖面**：沿 256 阶计算感知亮度（如 L\*）曲线。
   - 呈 **U 型或 ∩ 型**（中段与两端明显不同，偏移 > 0.06）→ **diverging**
   - 单调（或近单调）递增/递减 → **sequential 候选**
3. **上游官方分类交叉校验**（最高优先）：CMasher 的 `_known_cmap_types.py`、
   colorcet 官方 CSV 的文件名前缀（`diverging_`/`cyclic_`/`rainbow_`/`linear_`/
   `glasbey_` 等）、crameri 官方 PDF 版本史。文件名带 `rainbow_` 且官方归为虹谱的
   才标 rainbow（实测 colorcet 5 张）。
> 结论：**启发式只能当初筛**，最终类型以官方分类为准（本库 CSV 的 `c_type` 列已是
> 二者校正后的结果）。

### 8.3 快速选型表

| 数据 | 用（均可在 CSV 中按名检索） |
|---|---|
| 单向递增（表达量、丰度、深度） | cmasher：rainforest、ember、apple、jungle、ocean；colorcet：linear_* 系列 |
| 双向有中心（logFC、相关、距零偏差） | cmasher：redshift、iceburn、watermelon、holly；colorcet：diverging_* 系列 |
| 循环量（角度、相位、节律） | cmasher：copper、emergency、seasons；colorcet：cyclic_* 系列 |
| 黑洞中心发散（径向速度图等） | cmasher：iceburn、redshift（black-centered diverging） |
| >30 类分类变量 | 回 `colors/@palettes.csv`（离散 229 套），**勿**把连续色图切成 30 段当分类色 |

### 8.4 三条铁律（论文场景）

1. **发散数据必须配 diverging 色图**。用 sequential 表达 ±log2FC 会让"显著下调"
   看起来比"显著上调"浅——方向感全丢。中心务必压在白色/浅灰，两臂同色系镜像
   （红-蓝、绿-紫均可，避开红-绿相邻组合）。
2. **不要用 rainbow/jet 类表达数量**。色相跳跃会被读成不存在的分组边界，
   且感知亮度不均匀（同样的 Δ 在黄段和在蓝段视觉差几倍）。只有在对"连续但非量化"
   的变量做风格化呈现时才考虑（本组件仅收录 colorcet 官方 5 张 rainbow）。
3. **热图色图与数值区间匹配**。发散色图的中点要对齐实际中点（如 log2FC 的 0），
   而不是数据均值；离群值可设 `limits` 截断，避免全图挤进两端。

### 8.5 用法：查询式取 hex（本包默认，零依赖）

CSV 字段：`c_sour, c_name, c_type, n_stops, colors`（`;` 分隔 hex 列表，阶数见 `n_stops`）。

```r
# R：读 CSV 取 hex，喂给 ggplot2 / ComplexHeatmap
cm   <- read.csv("third_party/colormaps/@colormaps.csv")
hexs <- strsplit(cm$colors[cm$c_name == "redshift"], ";")[[1]]
# volcano/heatmap 用发散色图（中心对齐 0）：
ggplot(df, aes(x, y)) + geom_tile(aes(fill = z)) +
  scale_fill_gradientn(colours = hexs, limits = c(-4, 4), oob = scales::squish) +
  theme_classic(base_size = 12)
```

```python
# Python：pandas + matplotlib
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
cm = pd.read_csv("third_party/colormaps/@colormaps.csv")
hexs = cm.loc[cm.c_name == "rainforest", "colors"].iloc[0].split(";")
cmap = LinearSegmentedColormap.from_list("rainforest", hexs)
# 交互功能（在线预览、其他软件格式）：pip install cmasher / pip install colorcet
```

连续色图同时适用于 **fill 与 color**：热图/密度用 `scale_fill_gradientn()`，
轨迹线用 `scale_color_gradientn()`；每个模块出图脚本的 `plot_theme.R` /
`plot_theme.py` 已开放 hex 注入点。

### 8.6 许可与引用

- 组件：CMasher 各图 BSD-3-Clause（© 2019-2021 Ellert van der Velden，
  1313e/CMasher）；colorcet 各图 CC BY 4.0（© Peter Kovesi，
  holoviz/colorcet）。LICENSE 原文见 `third_party/colormaps/LICENSES/`，
  溯源与剔除范围见 `third_party/colormaps/NOTICE.md`
- **停更家族不收录**（2026-10-18 收录准则，末次提交 ≤ 约 1 年）：crameri（2023-06，
  Zenodo 存档）/ cmocean（2024-04）/ scico（2023-08）/ viridis（2024-01）**未克隆入包**，
  仅在该组件 NOTICE.md 的「调研后未收录的上游」表留档并附运行时安装命令
- 致谢建议：CMasher（van der Velden, https://github.com/1313e/CMasher）与
  colorcet（Kovesi, https://github.com/holoviz/colorcet）；
  如需 special 类色图，安装 `cmocean`/`crameri` 后按其文档引用
