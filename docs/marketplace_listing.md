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

## SKILL.md 里的 description（922 字符，保留全部触发词）

```yaml
description: 水稻（Oryza sativa）多组学与育种分析技能：转录组（QC、DESeq2/edgeR、GO/KEGG/Reactome、GSEA、WGCNA、TF/PPI）、代谢组（OPLS-DA/VIP、通路、ROC/随机森林、网络）、多组学（MOFA、DIABLO、O2PLS、联合通路、Circos、GRNBoost2）、变异检测（BWA→GATK→SnpEff）、GWAS（GLM/MLM/FarmCPU、3VmrMLM/Fast3VmrMLM 多位点 QTN·QEI·上位性、FastBiCmrMLM 二分类、mppQTL 多亲本、Meta-GWAS、共定位 PP4）、群体遗传（PCA/亲缘、LD、Fst/π、ADMIXTURE、XP-EHH）、泛基因组（SV/PAV、Pan-GWAS、eQTL/mQTL/pQTL）、QTL/BSA、蛋白组（limma、UniProt、InterPro）、表观组（WGBS/DMR、ATAC/ChIP）、表型组（BLUP/遗传力、光谱、GBLUP、AMMI/GGE）、单细胞/空间（Scanpy、Squidpy、RSTA/scGPT）、qPCR 引物、CRISPR sgRNA（序列内脱靶）、qPCR 2^-ΔΔCt、实验协议。分析即出 CNS 级图（PDF/TIFF/PNG、300 dpi、共享调色板）；默认先给分析与图形方案，确认后执行；缺依赖先列清单，同意后安装。触发词：水稻转录组、水稻代谢组、水稻GWAS、多位点GWAS、上位性、QTN×环境、QEI、case-control GWAS、MAGIC/NAM/ROAM、群体结构/群体遗传、泛基因组、eQTL、蛋白质组、表观组/甲基化/ATAC、表型/BLUP、单细胞/scRNA、空间转录组/Stereo-seq、引物设计、CRISPR 靶点、qPCR 定量、多组学联合分析、WGCNA、GSEA、DESeq2、OPLS-DA、MOFA、DIABLO、差异表达、富集分析、CNS 级绘图、发表图/论文配图、画图、写论文、论文润色、英文润色、查文献、文献检索、文献调研、下载论文、精读论文、总结结果、结果解读、审稿回复、论文转PPT、开题报告。
```

## 写作要点

1. **第一句写"是什么 + 给谁"**，不要以"本技能……"开头；
2. **能力用动词开头**（做差异分析、出图、定位 QTL），不要堆形容词；
3. **触发词放最后**，覆盖用户真实说法（水稻、转录组、GWAS、火山图……）；
4. **不写无法验证的承诺**（"最强""一键发 Nature"）；
5. **说明边界**（示例数据、需装 R/Python、先方案后执行），减少差评。
