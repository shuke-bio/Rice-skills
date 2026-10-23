# 技能市场上架文案（可直接复制）

## 基本信息

| 字段 | 内容 |
|---|---|
| 技能 ID（`name`） | `rice-skills` |
| 中文标题（若平台支持） | 水稻组学分析（Rice-skills） |
| 副标题 | 从数据到发表图，一套可复现的水稻组学分析技能 |
| 分类 | 数据分析 / 科研 |
| 标签 | 水稻、组学、GWAS、转录组、代谢组、育种、生物信息、科研绘图 |

## 短描述（≤60 字，列表卡片用）

> 水稻转录组、代谢组、GWAS、QTL、泛基因组、表型组、单细胞与空间分析；先出方案与图形方案，确认后执行。

## 长描述（详情页用）

**水稻全组学 + 育种分析，一条链路跑通。**

- **18 个模块、91 个 R/Python 脚本**：转录组、代谢组、多组学整合、变异检测、GWAS、群体遗传、泛基因组、QTL/BSA、蛋白组、表观组、表型组/GS、单细胞与空间、引物与 CRISPR 设计。
- **先方案，后执行**：读取数据契约后先给「阻断项 + 三档分析方案 + 图形方案」，你确认后才运行，避免白跑。
- **依赖透明**：缺 R/Python 包先列清单，经你同意才安装；系统工具只给指引。
- **出版级出图**：统一调色板，PDF/TIFF/PNG、300 dpi，单栏/双栏尺寸，跑完即得。
- **可复现**：固定随机种子，每次输出 `run_info.txt`（参数、包版本）。

**适用场景**

- 水稻转录组差异表达、富集、WGCNA、TF 调控
- 代谢组 OPLS-DA、通路、标志物
- GWAS（GLM/MLM/FarmCPU、多位点 QTN/QEI/上位性、case-control、多亲本 QTL、Meta-GWAS）
- QTL/BSA、泛基因组 SV/PAV、群体遗传
- 表型组与基因组选择、单细胞与空间转录组
- qPCR 引物、CRISPR sgRNA、qPCR 相对定量

**使用示例**

> 这批水稻转录组数据做差异分析，出火山图和富集图。

> 300 份材料、两个环境的株高，先做 MLM，再做多位点 QTN。

> 给 demo_gene.fa 设计 CRISPR sgRNA 和 qPCR 引物。

**说明**

- 示例数据为程序生成的演示数据，仅用于流程验证，不作生物学结论；
- 需要本机具备 R（≥4.3）与 Python（≥3.10），或按提示安装依赖；
- 技能执行前会先给方案，确认后才动数据。

## SKILL.md 里的 description（GitHub 原版，1203 字符）

```yaml
description: 水稻（Oryza sativa）多组学分析 Skill（全组学平台）。覆盖转录组、代谢组及联合分析（数据质控、DESeq2/edgeR、GO/KEGG/Reactome 富集、GSEA、WGCNA、转录因子与 PPI；OPLS-DA/VIP、代谢通路、ROC/随机森林、代谢网络；MOFA/DIABLO/O2PLS）；基因组学（变异检测 BWA→GATK→SnpEff、GWAS 三引擎、多位点 QTN·QEI·上位性（3VmrMLM/Fast3VmrMLM）、case-control QTN（FastBiCmrMLM，Linux）、多亲本 QTL（mppQTL）、Meta-GWAS、群体遗传、泛基因组 SV/PAV/Pan-GWAS、eQTL/mQTL QTL 映射、BSA/QTL-seq 扫描与区间提取）；蛋白质组（limma 差异蛋白、火山/热图）；表观基因组（WGBS/DMR、ATAC/ChIP 差异区域）；表型组（BLUP/遗传力、光谱指数、GS 指引）；单细胞与空间组学（Scanpy/Squidpy 流程、RSTA/RICE scGPT 集成指引）；qPCR 引物设计；CRISPR sgRNA 设计（含序列内脱靶统计；全基因组脱靶需外部工具）；qPCR 相对定量（2^-ΔΔCt）与转录组交叉验证；配套分子实验协议参考。**分析即出图：每个分析步骤自动生成 CNS 级出版图表**（低饱和色盲友好配色、Arial/Helvetica 7-9 pt、PDF/TIFF/PNG ≥300 dpi；绘图内置于脚本 + config/nature_palette.json，跑完即得，无需另行调用绘图技能）。**默认自动链**：跑分析 → 自动出图 → 自动产出结果解读总结（无需用户开口）；文献检索/下载、论文写作/润色、统计审查/预审稿/审稿回复、论文转 PPT 等（third_party 内置 22 个 Apache-2.0 nature-skills）为**按需加载**，用户明确要求时才触发（两级路由见 SKILL.md §5）。当用户提到水稻转录组、水稻代谢组、水稻GWAS、多位点GWAS、上位性、QTN×环境、QEI、case-control GWAS、MAGIC/NAM/ROAM、群体结构/群体遗传、泛基因组、eQTL、蛋白质组、表观组/甲基化/ATAC、表型/BLUP、单细胞/scRNA、空间转录组/Stereo-seq、引物设计、CRISPR 靶点、qPCR 定量、多组学联合分析、WGCNA、GSEA、DESeq2、OPLS-DA、MOFA、DIABLO、差异表达、富集分析、CNS 级绘图、发表图/论文配图、画图、写论文、论文润色、英文润色、查文献、文献检索、文献调研、下载论文、精读论文、总结结果、结果解读、审稿回复、论文转PPT、开题报告等关键词时触发。
```

## 写作要点

1. **第一句写"是什么 + 给谁"**，不要以"本技能……"开头；
2. **能力用动词开头**（做差异分析、出图、定位 QTL），不要堆形容词；
3. **触发词放最后**，覆盖用户真实说法（水稻、转录组、GWAS、火山图……）；
4. **不写无法验证的承诺**（"最强""一键发 Nature"）；
5. **说明边界**（示例数据、需装 R/Python、先方案后执行），减少差评。
