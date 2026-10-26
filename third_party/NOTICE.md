# NOTICE — third_party 第三方组件总索引

本目录收录 6 个精选上游技能库 + 4 个工具/数据组件（colors 配色推荐器、colormaps 连续色图数据、
ltc 与 ggpubfigs 离散配色数据）。**每个库的详细溯源（获取版本、commit、收录范围、
明确未包含内容）见各子目录的 `NOTICE.md`；各库许可证原文见其 `LICENSE*` 文件。**

| 子目录 | 上游 | 许可（LICENSE 原文已核验） | 规模 | 上游状态 |
|---|---|---|---|---|
| `bioskills/` | github.com/GPTomics/bioSkills | MIT | 41 类 354 技能 | 已归档快照 d91ed3d |
| `kdense-scientific-skills/` | github.com/K-Dense-AI/scientific-agent-skills | MIT | 28 技能 | 活跃（a1043ef） |
| `omicverse-skills/` | github.com/omicverse/omicverse-skills | MIT 文档（驱动主包 GPL-3.0，运行时安装） | 61 技能 | 活跃（d9b0a82） |
| `sciagent-skills/` | github.com/jaechang-hits/SciAgent-Skills | CC BY 4.0 | 64 技能 | 已归档快照 |
| `nature-skills/` | github.com/xiaolai/nature-skills | MIT | 22 技能 | 已归档快照 |
| `cell-skills/` | github.com/yrui-cmd/Cell + Cell_CNS_Figure | MIT ×2（LICENSE-Cell / LICENSE-Cell_CNS_Figure） | 12 技能 | 已归档快照 |
| `colors/` | github.com/lcpmgh/colors | MIT | 229 套配色 + Shiny app（**工具组件，非技能库**） | 活跃（36b2cc6） |
| `colormaps/` | github.com/1313e/CMasher + github.com/holoviz/colorcet | BSD-3-Clause / CC BY 4.0 | 133 张科学连续色图 hex 数据 + 快速选型 README（**数据组件**） | 双库均活跃（2026-10） |
| `ltc/` | github.com/loukesio/ltc-color-palettes | MIT | 32 套自然风调色板 hex 数据（**数据组件**） | 活跃（2026-08-30） |
| `ggpubfigs/` | github.com/JLSteenwyk/ggpubfigs | MIT | 15 套色盲友好调色板 hex 数据（**数据组件**） | 活跃（2026-05-12） |

**汇总**：541 个 SKILL.md + 4 个工具/数据组件（colors、colormaps、ltc、ggpubfigs），约 30 MB。
Rice-skills 自研部分（scripts/、docs/、config/、templates/、tests/、examples/、assets/）不依赖这些文件；
它们用于按 `docs/skill_routing.md` 路由表提供专家经验参考；四个配色组件共同服务于科研绘图取色速查
（离散配色见 `colors/` 与 `ltc/ggpubfigs/`，连续色图见 `colormaps/`，方法论见
`docs/references/plot_color_palettes.md`）。
**收录准则**（2026-10-18 起执行）：上游末次提交不超过约 1 年且许可证原文可核验；长期停更的上游
（如 crameri / cmocean / scico / viridis）不克隆入包，仅在 `colormaps/NOTICE.md` 的
「调研后未收录的上游」表中留档并附运行时安装命令。
二次分发时请保留本目录全部 `LICENSE*` 与各库 `NOTICE.md`。

## 隐私脱敏说明

为保护个人隐私，本目录上游文档中的个人邮箱与本地用户路径已做脱敏处理
（邮箱替换为“已脱敏”，`/Users/<name>/` 替换为 `/Users/demo/`），不影响任何功能；
原始信息请见各上游仓库。

## 运行时外部工具（未打包）

本仓库不复制、不打包下列工具的源码或二进制，仅在脚本中生成命令或调用本机安装版本：
fastp、HISAT2、samtools、featureCounts、PLINK、rMVP、GAPIT、Cas-OFFinder、CRISPOR、
vg、PanGenie、minigraph-cactus、odgi、cooler、HiCExplorer、Sniffles2、cuteSV、Delly、
Manta、SURVIVOR、selscan、Dsuite、SIRIUS、GNPS、susieR、rrBLUP、BGLR、GSVA、
harmonypy、celltypist、scrublet、cell2location、COMMOT、BE-Hive、PRIDICT 等。

完整清单、官方入口与许可核验状态见 `docs/runtime_dependencies.md`。分发 Docker 镜像或
商业部署前，必须逐一核对上游 LICENSE 原文并把结论补入本文件；GPL 系工具只通过命令行
调用，不并入本仓库 MIT 主体。
