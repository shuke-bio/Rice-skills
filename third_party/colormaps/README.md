# colormaps — 连续科学色图数据组件

133 套**感知均匀连续色图**（sequential / diverging / cyclic / rainbow）的
hex 数据，来自 2 个活跃维护的上游家族（cmasher / colorcet）。
色阶数按上游 `_hex.txt` 如实转录（见 `n_stops` 列）：116 套 256 级，
cmasher 循环/发散类 510–511 级，cmasher `neutral` 2560 级。
许可证原文见 `LICENSES/`，溯源明细与收录准则见 `NOTICE.md`。

## 数据 Schema（@colormaps.csv）

| 列 | 含义 |
|---|---|
| `c_sour` | 家族：cmasher / colorcet |
| `c_name` | 色图名（与上游一致；colorcet 用官方 CSV 基名，如 `diverging_bwr_40-95_c42`） |
| `c_type` | sequential / diverging / cyclic / rainbow |
| `n_stops` | 色阶数（如实转录上游文件行数：116 套 256、cmasher 循环/发散 510–511、`neutral` 2560） |
| `colors` | hex 列表，`;` 分隔，个数 = `n_stops` |

## 快速选型

| 数据 | 用 |
|---|---|
| 单向递增（表达量、丰度、深度） | cmasher：rainforest、ember、apple、jungle；colorcet：linear_* 系列 |
| 双向有中心（logFC、相关、距零偏差） | cmasher：redshift、iceburn、watermelon、holly；colorcet：diverging_* 系列 |
| 循环量（角度、相位） | cmasher：copper、emergency、seasons；colorcet：cyclic_* 系列 |
| 黑洞中心发散（径向速度图等） | cmasher：iceburn、redshift（black-centered diverging） |
| >30 类分类变量 | 见 `colors/@palettes.csv`（离散配色 229 套），勿用连续色图切段 |

原则：**发散数据必须用 diverging 色图**（中心浅、两端深）；避免 rainbow/jet 类表达
数量（视觉畸变）；全部 133 套均经上游感知均匀性与色盲友好（CVD）诊断。

## 用法

查询式（本包默认，无需装上游）：

```r
# R：读 CSV 取 hex，喂给 ggplot2
cm <- read.csv("third_party/colormaps/@colormaps.csv")
hexs <- strsplit(cm$colors[cm$c_name == "redshift"], ";")[[1]]
ggplot(df, aes(x, y)) + geom_tile(aes(fill = z)) +
  scale_fill_gradientn(colours = hexs)
```

```python
# Python：pandas + matplotlib
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
cm = pd.read_csv("third_party/colormaps/@colormaps.csv")
hexs = cm.loc[cm.c_name == "rainforest", "colors"].iloc[0].split(";")
cmap = LinearSegmentedColormap.from_list("rainforest", hexs)
```

完整功能（交互选色、其他软件格式）安装上游：`pip install cmasher` / `pip install colorcet`。
停更家族（crameri/cmocean/viridis）的运行时安装命令见 `NOTICE.md` §未收录上游。
选型方法论与期刊配色查表见 `docs/references/plot_color_palettes.md`。
