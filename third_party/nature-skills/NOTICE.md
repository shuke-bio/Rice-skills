# NOTICE — 第三方组件溯源说明

本目录包含从上游仓库复制的第三方 skill 源码，版权归原作者所有，按上游许可证使用。

## 上游信息

| 项 | 内容 |
|---|---|
| 仓库 | https://github.com/Yuan1z0825/nature-skills |
| 获取时间 | 2026-10-16，shallow clone commit `7c1cdd2f4c484d4038ea36f1b7397f61237a4f90`（2026-10-15 23:31 +0800） |
| 许可证 | **Apache License 2.0**（全文见本目录 `LICENSE`） |
| 历史版本 | 本包 v3.1–v3.2 曾收录同一技能族的 v1 版 9 个技能（当时经镜像 xiaolai/nature-skills 获取，MIT License, Copyright (c) 2026 Yuan Yizhe）；本次已整体替换为作者主仓最新版（v2 静态/动态分层架构，仓库级许可证升级为 Apache-2.0） |

## 已包含的内容（22 个技能，仅复制源码与参考资料，排除演示资产）

| 目录 | 用途 | 许可 |
|---|---|---|
| `skills/nature-writing/` | 论文各章节写作（摘要/引言/方法/结果/讨论/结论） | Apache-2.0 |
| `skills/nature-polishing/` | 英文论文润色（25 条规则、12 步工作流） | Apache-2.0 |
| `skills/nature-academic-search/` | 多源文献检索、MeSH 策略、引用文件管理（RIS/BibTeX） | Apache-2.0 |
| `skills/nature-reader/` | 论文全文中英对照精读辅助（图表感知） | Apache-2.0 |
| `skills/nature-paper-card/` | 单篇论文结构化精读卡片（01–16 节，来源可约束） | Apache-2.0 |
| `skills/nature-paper-xray/` | 论文深读"X 光"剖析（还原作者真实思路与证据边界） | Apache-2.0 |
| `skills/nature-citation/` | Nature/CNS 严格引用规范 | Apache-2.0 |
| `skills/nature-ref-verifier/` | 参考文献逐条多源交叉验证（投稿前查引用错误） | Apache-2.0 |
| `skills/nature-statistics/` | 统计报告规范审查（伪重复/多重比较/显著性滥用） | Apache-2.0 |
| `skills/nature-reviewer/` | Nature 风格预审稿模拟（审稿人视角压力测试） | Apache-2.0 |
| `skills/nature-response/` | 审稿意见逐条回复信工作流 | Apache-2.0 |
| `skills/nature-figure/` | Nature 级图表规范与 matplotlib/R 模板 | Apache-2.0 |
| `skills/nature-data/` | 数据可用性声明（Data Availability）规范 | Apache-2.0 |
| `skills/nature-downloader/` | 论文全文合法下载（出版社 API / OA / CNKI / 机构浏览器会话） | 自带 MIT LICENSE（Copyright (c) 2026 baihe26），与 Apache-2.0 兼容 |
| `skills/nature-image2ppt/` | 图片/截图/扫描 PDF 还原为对象级可编辑 PPTX | 自带 MIT LICENSE，与 Apache-2.0 兼容 |
| `skills/nature-literature-pipeline/` | 文献监测管线（多源检索→六维评分→精读→推送→归档，支持 cron） | SKILL.md 声明 MIT（Jiahao8595） |
| `skills/nature-paper2ppt/` | 论文一键转 Nature 风格中文汇报 PPTX | Apache-2.0 |
| `skills/nature-paper-trans/` | 英文论文 PDF 逐页翻译为简体中文图片型 PDF | Apache-2.0 |
| `skills/nature-experiment-log/` | 标准化实验日志记录（图片/语音/文字→YAML Markdown） | SKILL.md 声明 MIT |
| `skills/nature-proposal-writer/` | 研究计划/开题报告/项目申请书写作 | SKILL.md 声明 MIT |
| `skills/nature-paper-to-patent/` | 论文→发明专利技术交底书/专利草案 | Apache-2.0 |
| `skills/nature-shared/` | 上游内部共享参考包（被 writing/polishing/response/reader/paper2ppt 依赖） | Apache-2.0 |

## 明确未包含的内容

- **所有技能的 `assets/` 目录**：演示图片、示例 PDF、图集预览等（上游约 4.6 MB + nature-figure 内约 34 MB + nature-paper-trans 内约 1.2 MB），
  含已发表论文的演示素材，出于仓库体积与第三方素材许可审慎考虑未复制；需要时请从上游仓库获取。
- **`tests/` 以外的大型 `downloads/`、`examples/`、`__pycache__` 等运行产物**未包含。
- **`skills/nature-social-science/`**：社会科学方向选题/写作工作台，与本包水稻科研定位无关，未收录（ Apache-2.0，需要者可自行从上游获取）。
- 上游仓库的 `.claude-plugin/`、CI 配置、`index.html`、`docs/`、`scripts/` 等仓库基础设施未包含。

## 使用说明

- 这些文件以"参考规范与模板库"形式提供：Rice-skills 自身的分析脚本使用自研绘图主题
  （`config/nature_palette.json` + `scripts/utils/plot_theme.R` / `plot_theme.py`），与上游代码相互独立。
- 完整流程路由（分析→绘图→结果总结→文献→写作→投稿/汇报）见根目录 `SKILL.md` §5。
- 二次分发本项目时请保留本目录的 `LICENSE` 与 `NOTICE.md`；对含自带 MIT 声明的技能，
  其上游许可条款同样继续有效。
