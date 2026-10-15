# NOTICE — 第三方组件溯源说明（Cell 科研工作流技能集）

本目录包含从上游仓库复制的第三方 skill 源码，版权归原作者所有，按上游许可证使用。

## 上游信息（两个来源）

| 项 | 来源 1 | 来源 2 |
|---|---|---|
| 仓库 | https://github.com/yrui-cmd/Cell | https://github.com/yrui-cmd/Cell_CNS_Figure |
| 获取时间 | 2026-10-14 | 2026-10-14 |
| 获取版本 | commit `9ee2480`（2026-10-12） | commit `6ea19c4`（2026-10-12） |
| 许可证 | MIT License, Copyright (c) 2026 yutuur (木纹) | MIT License, Copyright (c) 2026 yrui-cmd |
| LICENSE 文件 | 本目录 `LICENSE-Cell` | 本目录 `LICENSE-Cell_CNS_Figure` |

## 已包含的内容（12 个技能，248 个文件）

| 技能目录 | 用途 | 来源 |
|---|---|---|
| `skills/cell/` | 总入口路由（按科研任务分流到成员技能） | 来源 1 |
| `skills/cell-brainstorm/` | 选题头脑风暴（证据驱动的科研选题） | 来源 1 |
| `skills/cell-plan/` | 研究计划（按题目/投稿目标安排路线） | 来源 1 |
| `skills/cell-review/` | 顶刊综述 | 来源 1 |
| `skills/cell-manuscript-editing/` | 英文论文深度润色 | 来源 1 |
| `skills/cell-reviewer-response/` | 审稿回复（point-by-point） | 来源 1 |
| `skills/cell-submission/` | 投稿材料整理 | 来源 1 |
| `skills/cell-visualization-code/` | Python/MATLAB 科研可视化代码 | 来源 1 |
| `skills/cell-data-figure/` | 数据图复现（本地脚本，无外部依赖） | 来源 1 |
| `skills/cell_su7/` | 科研图转矢量（PPT / Adobe Illustrator 流程） | 来源 1 |
| `skills/cell-ppt-edited/` | PowerPoint 原生对象编辑 | 来源 1 |
| `skills/cell-cns-figure/` | CNS 风格科研图重建（独立仓库版，含 docs/INSTALL.md） | 来源 2 |

另收录：`README-Cell.md`（上游总览文档）。

## ⚠️ 外部依赖提示（使用前必读）

- **`skills/cell-cns-figure/`**：依赖**外部付费服务**（额度制：入口 1 约 20 额度 / 入口 2 约 45 额度），
  需按该技能 `docs/INSTALL.md` 配置第三方 API；会产生真实费用，提交任务前请自行确认。
- **`skills/cell_su7/`**：矢量转换流程需要 Windows + PowerPoint 或 Adobe Illustrator
  （据上游说明 macOS 路径仅有模拟测试）。
- **`skills/cell-ppt-edited/`**：需要 Windows 桌面版 Microsoft PowerPoint（不支持 WPS / macOS）。
- 各技能目录内的 `requirements.txt` 需按需单独安装（如 cell-cns-figure 需要
  Pillow / svgpathtools / cssselect2 / tinycss2）。

## 使用说明

- 与既有集成互不冲突：`third_party/nature-skills`（Nature 绘图/润色/写作）、
  `third_party/bioskills`、`third_party/kdense-scientific-skills` 各自独立目录，按需加载。
- 二次分发本项目时，请保留本目录的两个 LICENSE 文件与 `NOTICE.md`。
