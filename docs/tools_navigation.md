# 水稻生信工具导航（GitHub 开源生态盘点）

> **更新**：2026-10-18。本文件是**导航与许可评估**，不是代码收录——表中工具均通过
> 官方仓库获取，本包**不复制**其中任何源代码（已收录的 MIT / BSD-3 / CC BY 组件除外，见 §7）。收录/引用任何工具前请按项目红线
> **打开其 LICENSE 原文核验**（GitHub 标签不可全信）。
> 与 `rice_databases.md`（数据库导航）互补：本文件面向**分析工具与流程**。

## 0. 一句话结论

- **TBtools 不能收**：仅"非商用免费"、无开源许可证文件——学术上可自行下载使用，
  但**不得打包进任何仓库分发**（商业使用还需作者授权）。
- 许可干净的（MIT/Apache）如 rMVP、RicePilaf、Et-GWAS、hifiasm 可用于参考与集成。
- **colors**（lcpmgh/colors，MIT 已核验）：科研绘图配色推荐器，已按 `third_party/NOTICE.md`
  流程收录为 `third_party/colors/`（见 §7）。
- **colormaps / ltc / ggpubfigs**（2026-10-18 新增，均 MIT/BSD-3/CC BY 原文已核验）：
  连续色图数据（CMasher + colorcet 共 133 张）与两套离散配色数据（ltc 32 套、ggpubfigs 15 套），
  已按同一流程收录为 `third_party/colormaps/`、`third_party/ltc/`、`third_party/ggpubfigs/`（见 §7）。
- 大多数专业工具（QTLseqr、GAPIT、easyQTLseq、deepBSA 等）**未声明许可证**——
  使用没问题，但“收录”前必须联系作者或打开 LICENSE 核验。
- 本包已内置等效功能清单见第 6 节——**先查已内置再找外部工具**。

## 1. BSA / QTL-seq（混池测序定位）

| 工具 | 用途 | 仓库 | 许可（GitHub 标注） | 与本包关系 |
|---|---|---|---|---|
| **QTLseqr** ★85 | NGS-BSA 全基因组扫描（ΔSNP-index / G′ 两种统计） | `bmansfeld/QTLseqr` | **无 LICENSE 文件**（使用前核验 DESCRIPTION） | 互补：QTLseqr 全基因组扫峰 → 本包 `qtl/04_qtl_region_snps.py` 对峰区间做精细 SNP 提取+基因注释 |
| **Et-GWAS** | 极端表型 GWAS（3000 水稻基因组，低/高/对照三池） | `IRRI-South-Asia-Hub/Et-GWAS` | **MIT** ✅ | 互补：稀有变异关联路线；许可干净可放心学习代码 |
| easyQTLseq | QTL-seq 一体化（R Shiny） | `yzeng2017/easyQTLseq` | 未声明 | 备选 GUI 方案 |
| deepBSA | BSA 深度学习 + 传统方法 | `biollh/deepBSA` | 未声明 | 方法参考 |

## 2. GWAS 引擎

| 工具 | 用途 | 仓库 | 许可 | 与本包关系 |
|---|---|---|---|---|
| **rMVP** ★341 | 高性能 GWAS（C 内核，百万级标记） | `xiaolei-lab/rMVP` | **Apache-2.0** ✅ | 本包 `gwas/` 文档中为可选引擎之一 |
| **GAPIT** ★244 | 经典 GWAS R 包（GLM/MLM/CMLM/SUPER/MLMM/FarmCPU/Blink + GS） | `jiabowang/GAPIT` | **GPL-2.0+**（DESCRIPTION 核验，2026-10） | 本包 `gwas/03_association.R --engine gapit --model` 可调用；MIT 主体不打包，运行时安装 |
| **IIIVmrMLM** | 3VmrMLM 多位点 QTN / QEI / 上位性（R+C++） | `YuanmingZhang65/IIIVmrMLM`（CRAN） | **GPL-3**（DESCRIPTION 核验，2026-10） | 本包 `gwas/09_multilocus_qtn.R` 默认引擎；运行时安装，不打包 |
| **Fast3VmrMLM** | 大数据 3VmrMLM（百万级标记；Hap/NCII/mQTL/MEJA/QQI 模块） | `YuanmingZhang65/Fast3VmrMLM` | **GPL-3**（DESCRIPTION 声明；仓库 LICENSE 文件写 BSD-3，**两者不一致，按 GPL 保守处理**） | `gwas/09_multilocus_qtn.R --engine Fast3VmrMLM`；预编译二进制按 R 4.3，R≥4.4 可能 DLL 加载失败 |
| **FastBiCmrMLM** | 二分类/case-control 性状多位点 GWAS（混合逻辑回归） | `YuanmingZhang65/FastBiCmrMLM` | **无 LICENSE 文件**（仓库无许可声明，使用需自行评估） | `gwas/10_case_control_qtn.R`；**仅 Linux**，运行时安装 |
| **mppQTL** | 多亲本群体（MAGIC/NAM/ROAM）多位点 QTL | `YuanmingZhang65/mppQTL`（随 IIIVmrMLM 仓库分发） | **GPL-3**（DESCRIPTION 核验，2026-10） | `gwas/11_multiparent_qtl.R`；依赖 qtl2，运行时安装 |
| TASSEL | 关联+连锁分析（Java GUI） | `tassel` | 按项目网页说明 | 备选 |
| **RiceG2G** | GWAS/QTL 位点**因果基因优先级**（多证据打分） | `xhhuanglab/RiceG2G` | 未声明 | 互补：接在本包 GWAS/共定位下游做候选排序 |
| **RicePilaf** | Post-GWAS/QTL 仪表板（泛基因组/共表达/调控/表观/富集整合） | `bioinfodlsu/rice-pilaf` | **MIT** ✅ | 互补：功能解读层；论文 GigaScience 2024 |
| geneNR | GWAS/QTL 候选基因窗口检索（R，CRAN） | CRAN / `Nirmalaruban/geneNR` | 查 CRAN 页面 | 轻量替代 |

## 3. 变异检测与注释基础设施（事实标准）

| 工具 | 用途 | 仓库 | 许可 | 与本包关系 |
|---|---|---|---|---|
| BWA / samtools / bcftools | 比对 / BAM 处理 / VCF 操作 | `samtools/bcftools` 等 | bcftools: 标注 NOASSERTION（**用前打开 LICENSE 核验**，htslib 系为宽松许可） | 本包 `variants/01-03` 调用 |
| **SnpEff** | 变异功能注释 | `pcingola/SnpEff` | 查仓库（论文 MIT-ish 需核验） | 本包 `variants/04_snpeff_annotate.sh` 调用 |
| **hifiasm** ★816 | PacBio HiFi 组装（水稻泛基因组常用） | `chhylp123/hifiasm` | **MIT** ✅ | 与 `pangenome/` 上游衔接 |
| GATK | 金标准变异检测 | broadinstitute/gatk | 宽松（BSD 系，需核验） | 本包 `variants/02_gatk_call.sh` 调用 |

## 4. 综合工具箱（TBtools 类，GUI）

| 工具 | 说明 | 许可红线 |
|---|---|---|
| **TBtools** | 130+ 功能 GUI 工具箱（序列提取/GO/热图/Circos/QTL-seq 图…），湿实验首选 | ⚠️ **仅非商用免费，无开源许可证**——只能自行下载使用（学术免费），**不可打包/镜像/再分发**；收费培训等商业场景需作者授权 |
| Rice Galaxy（IRRI） | 网页版组学分析平台（含 3K 基因组数据） | GPL 系（工具各自许可） |
| IRRI Bioinformatics 工具集 | SNPseek / BRS-utilities / Crop Galaxy 等 | 逐仓库核验 |

## 5. 数据库与在线平台（详见 rice_databases.md）

- **SNP-Seek**（IRRI，3K Rice Genomes 的 SNP/表型查询与下载）
- **RicePilaf 在线版**（post-GWAS 交互分析 demo 站）
- RAP-DB / RiceXPro / PlantTFDB / agriGO / ROAD / PlantPAN — 已收录于 `rice_databases.md`

## 6. 先查本包已内置（避免重复造轮子）

| 需求 | 本包位置 |
|---|---|
| GWAS 关联（GLM/可选 rMVP/GAPIT 路线） | `scripts/gwas/03_association.R` + `01-02`（PLINK QC） |
| BSA/QTL-seq 峰区间精细 SNP 提取 + 基因区域注释 | `scripts/qtl/04_qtl_region_snps.py`（**v3.2 新增**） |
| 共定位（筛查 + Bayesian PP4） | `scripts/qtl/02_colocalization.R` + `03_coloc_bayesian.R` |
| 变异检测链（BWA→GATK→过滤→SnpEff） | `scripts/variants/01-04` |
| 功能富集（GO/KEGG/GSEA） | `scripts/transcriptomics/03-04` |
| 出版图（Manhattan/QQ/热图/Circos 等） | 全模块统一 Nature 调色板（`config/nature_palette.json`）；离散分组配色速查 `third_party/colors`（229 套，MIT）、`third_party/ltc`（32 套自然风，MIT）、`third_party/ggpubfigs`（15 套色盲友好，MIT）；连续色图速查 `third_party/colormaps`（133 张，BSD-3 / CC BY 4.0）——均见 §7 |

## 7. 出版图配色组件（colors / colormaps / ltc / ggpubfigs，均已收录）

| 工具 | 用途 | 仓库 | 许可 | 与本包关系 |
|---|---|---|---|---|
| **colors** ★818 | 科研绘图配色推荐器（R/Shiny 双栏对比，229 套配色：顶刊风格 177 + ggsci/grDevices/RColorBrewer 52，每套 2–51 色） | `lcpmgh/colors` | **MIT**（LICENSE 原文已核验，Copyright (c) 2025 Liang Chen） | **已收录** `third_party/colors/`（v4.0, commit `36b2cc6`；剔除知乎抓取文章 `articles/` 与旧版本 `oldversion/`）；方法笔记 `docs/references/plot_color_palettes.md`（选色原则 + ggplot2 取用代码） |
| **CMasher** ★400 | 科学连续色图（23 基础 + 30 混合，科研配色重头） | `1313e/CMasher` | **BSD-3-Clause**（LICENSE 原文已核验，© 2019-2021 Ellert van der Velden） | **已收录** `third_party/colormaps/@colormaps.csv`（含 colorcet 共 133 张；末次提交 2026-10-09，活跃） |
| **colorcet** ★210 | 感知均匀连续色图（widely used by HoloViews/PyViz） | `holoviz/colorcet` | **CC BY 4.0**（LICENSE 原文已核验，© Peter Kovesi；再分发须署名） | **已收录** `third_party/colormaps/@colormaps.csv`（80 张；末次提交 2026-10-12，活跃） |
| **ltc** ★303 | 自然风调色板（基于 Books/Chart 等视觉来源，每套含 rgb hex） | `loukesio/ltc-color-palettes` | **MIT**（LICENSE 原文已核验，© 2021 ltc authors） | **已收录** `third_party/ltc/@palettes.csv`（32 套；末次提交 2026-08-30，活跃） |
| **ggpubfigs** ★132 | 色盲友好出版配色（NPG/Lancet/NEJM/JAMA/JCO/AAAS/Tol/Viridis 等 15 套） | `JLSteenwyk/ggpubfigs` | **MIT**（LICENSE 原文已核验，Copyright (c) 2020 Jacob L. Steenwyk） | **已收录** `third_party/ggpubfigs/@palettes.csv`（15 套；末次提交 2026-05-12，活跃） |

用法：离散配色查询式直接读各组件 `@palettes.csv` 配 `scale_color_manual()`；连续色图查询式读
`third_party/colormaps/@colormaps.csv`（字段 `c_sour/c_name/c_type/n_stops/colors`，
选型表与 R/Python 取用代码见 `third_party/colormaps/README.md`）；交互挑选本地
`shiny::runApp("third_party/colors")`（依赖 shiny/ggplot2/DT/colourpicker 等 8 个 CRAN 包）。

> 长期停更的上游（crameri / cmocean / scico / viridis，末次提交 2023–2024）按 2026-10-18
> 收录准则**不入包**，仅在 `third_party/colormaps/NOTICE.md` 的「调研后未收录的上游」表留档，
> 附运行时安装命令。

## 8. 获取方式约定

```bash
# 推荐的运行时获取（不进仓库）：
git clone https://github.com/<owner>/<tool>.git   # 到工作区自行使用
# 或：conda install -c bioconda <tool> / install.packages("<pkg>") / pip install <pkg>
```

仅当某工具许可为 MIT/BSD/Apache/CC-BY 且确有必要时，才考虑按 `third_party/NOTICE.md`
的流程收录（保留 LICENSE + 写 NOTICE 溯源）；GPL 系只能作独立目录组件；
**无许可证工具一律不收录**。
