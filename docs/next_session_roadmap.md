# 下阶段路线图（新会话从这里继续）

> **交接快照（2026-10-14 深夜）**：v3.1.0 完成，commit `93e5100`，59 分析脚本 + 17 工具，
> 自检 39/39 通过（38 项整轮 + XP-EHH 修复后等价补验）。桌面 ZIP 仍是 v2.1.1 旧版（未重打包，遵用户指示）。
> 本文档承接用户在会话结束前提出的三个方向。

## 1. 数据库补充（✅ 本轮已完成大部分，见 docs/rice_databases.md）

用户点名的 7 个站点已补入 `rice_databases.md`：
Ensembl Plants FTP、InterPro（蛋白结构域）、ROAD（水稻 KEGG 面板）、RGAP R7 下载页、
UniProt（UP000059680 参考 proteome）、PlantPAN 4.0（启动子/motif）、agriGO v2（水稻 GO 富集）。

**遗留小项（✅ 2026-10-15 全部完成）**：
- [x] UniProt UP000059680 映射自动化 → `scripts/proteomics/04_uniprot_mapping.py`（下载/缓存/离线兜底，
      MSU 优先双风格解析，输出模块 1 富集衔接的 protein_gene_map.csv；真实下载待网络良好时跑一次缓存）
- [x] InterProScan TSV 解析 → `scripts/proteomics/05_interpro_domains.py`（含差异蛋白结构域富集）

## 2. GitHub 水稻 skills 调研（✅ 已完成 2026-10-15：K-Dense 二轮 + bioSkills 扩容 + omicverse-skills 新收）

**成果总览**（第三方技能总数 185 → 444）：
1. **K-Dense 二轮 +23**（5→28）：单细胞 anndata/scanpy/scvi-tools/scvelo、GRN arboreto/networkx、
   数据库 bioservices/gget、质谱 matchms/pyopenms、表观 deeptools、进化 biopython/etetoolkit/phylogenetics、
   变异 pysam、统计 statistical-power/statsmodels/scikit-learn/umap-learn、
   实验 mageck/primer-design/experimental-design/cellprofiler。commit `4df4059`。
2. **bioSkills 扩容 +264**（90→354，41 类）：单细胞(17)/空间(12)/表观 4 类(43)/蛋白组(9)/
   剪接+smallRNA+表转录(20)/GRN(6)/因果(11)/比较基因组(13)/长读(9)/基因组组装注释区间编辑(29)/
   CNV+相位(15)/序列比对QC定量(65)/数据库(15)/ML+报告+实验设计(17)。
   上游 commit `d91ed3d`（2026-08-15，仓库已归档，无漂移）。
3. **omicverse-skills 新收 61**（新库 third_party/omicverse-skills/，MIT 文档，上游 commit `d9b0a82`）：
   单细胞全流程 + **基础模型 scGPT/Geneformer/scFoundation**（模块 16 RSTA 路线直接参考）+
   CellPhoneDB/LIANA/Scissor + 空间映射 + bulk WGCNA + 代谢组多因子通路 + 微生物组。
   ⚠️ 驱动主包 omicverse 为 GPL-3.0——仅收 MIT 文档，主包运行时安装，NOTICE 已注明。

**许可核验**（全部打开原文验证）：K-Dense MIT（2025 K-Dense Inc.）✅、bioSkills MIT（2026 GPTomics）✅、
omicverse-skills MIT（2026 Zehua Zeng）✅；bioservices/etetoolkit/omicverse 主包等 GPL 库一律
"仅收 MIT 文档 + 库运行时安装 + NOTICE 注明"模式。
**未收但记录在案**：PlantCV（MPL-2.0 不在白名单→整理工作流模板路线，未做）、
bioSkills 的 metagenomics/microbiome/workflows(41)/chemoinformatics/临床类（候选扩展）、
K-Dense 的 qiime2/scikit-bio/esm/nmrglue 等（候选）。
**上游提取物**：`%TEMP%\kdense_ok`（K-Dense a1043ef）、`%TEMP%\bioskills-up`（d91ed3d）、
`%TEMP%\ov-skills-x`（d9b0a82）——临时目录可能被清理，需要时重新下载。

**许可合规红线（不变）**：只收 MIT/BSD/Apache-2.0/CC BY 4.0（需保留 LICENSE + 写 NOTICE）；
GPL 系独立目录；无 LICENSE 的（如 RSTA、APAV、ChromSkills）一律"运行时获取"模式，不复制源码。
网络提示：github 直连不稳时用 **ghfast.top** 镜像 tarball（`https://ghfast.top/https://github.com/<owner>/<repo>/archive/refs/heads/main.tar.gz`）。

**搜索关键词**：`rice bioinformatics skill`、`plant genomics claude skill`、`GWAS skill`、
`bioinformatics skills collection`、`agronomy LLM skill`、`breeding analysis skill`。

**已知候选金矿（按优先级）**：
1. **K-Dense（github.com/K-Dense-AI/claude-scientific-skills，MIT）**：177 个技能我们只收了 5 个。
   重点翻：proteomics、scanpy/anndata、gsea、pathway、biomart/ensembl 接入、PDB/结构类——
   逐个查 `skills/` 目录，挑与水稻全组学互补的收编。
2. **omics 类**：omicsverse（已在 third_party 有空间相关）、Scanpy 生态官方教程技能化、
   CellChat/NicheNet 的 R 流程技能。
3. **PlantCV**（MPL-2.0）：图像表型——目前只有文档指引，可考虑把其工作流模板做成技能目录。
4. **Galaxy/EF tools 的配方技能**：无许可风险的方法整理。

**收编流程**（每个库）：克隆/镜像 → 逐 SKILL.md 审阅 → LICENSE 核验（打开原文！不能只信 GitHub 标签）
→ 复制 + NOTICE 溯源（仓库/commit/收录范围/改动）→ README 计数更新 → 自检。

## 3. 候选分析清单（✅ 2026-10-15 八项全部实现并 demo 自检通过）

| 分析 | 状态 | 脚本 |
|---|---|---|
| Bayesian coloc | ✅ | `scripts/qtl/03_coloc_bayesian.R`（coloc.abf 自实现零依赖，PP4=0.994 demo） |
| motif → TF 匹配 | ✅ | `scripts/tf/02_motif_tf_match.py`（MEME PWM 库，衔接 01 k-mer） |
| 差异蛋白结构域富集 | ✅ | `scripts/proteomics/05_interpro_domains.py`（超几何 + BH） |
| UniProt 映射自动化 | ✅ | `scripts/proteomics/04_uniprot_mapping.py` |
| 共表达网络进阶 | ✅ | `scripts/integration/07_grn_boost.py`（GRNBoost2 式，sklearn） |
| miRNA 靶基因预测本地化 | ✅ | `scripts/evolution/02_mirna_targets.py`（psRNATarget 式打分） |
| 表型功效分析 | ✅ | `scripts/phenomics/05_power_analysis.R`（t/GWAS/GS 三段，与 power.t.test 对拍） |
| 成对性状关联 2D-GWAS | ✅ | `scripts/gwas/08_2d_association.R`（PATOWAS 方法整理） |

PlantCV 整理 → `docs/references/plantcv_workflows.md`（MPL-2.0 方法笔记，未复制源码）。
各分析使用时先查 `docs/skill_routing.md` 对应行的第三方技能经验（版本坑/参数选择）。

## 4. 工程待办（不变）

- [ ] **装 16 个缺失 R 包**（DESeq2、edgeR、fgsea、MOFA2、mixOmics、ropls、ComplexHeatmap、
      MetaboAnalystR、pROC、randomForest、UpSetR、ggalluvial、circlize、ggraph、patchwork、
      org.Osativa.eg.db）→ 跑模块 1-3 的 R 全流程验证
- [ ] **重打包 ZIP**（用户要求时一次性做：当前桌面 ZIP 是 v2.1.1 旧版）
- [ ] **GitHub 推送**（建 repo + push，需用户 GitHub 账号）
- [ ] 真实数据试跑（用户的水稻 GWAS/组学数据）

## 5. 已验证模块一览（v3.2 终态）

模块 1-18 全部就位；v3.1 增补：单倍型（07_haplotype.R）、共定位筛查（02_colocalization.R）、
GBLUP（03_gs_gblup.R）、AMMI/GGE（04_ammi_gge.R）、XP-EHH（07_xpehh.R）、KA/KS（01_ka_ks.py）、
启动子 k-mer（01_promoter_kmer.py）、miRNA/剪接文档；nature-skills 9 项；
第三方技能库 8 个（bioskills 90 / SciAgent 64 / Cell 12 / K-Dense **28**（2026-10-15 扩充）/ nature-skills 9 / 等）。
