<div align="center">

<img src="assets/cover.jpg" alt="Rice-skills" width="100%">

**English** | [简体中文](README.zh-CN.md)

<h1>Rice-skills</h1>

<h3>水稻（<i>Oryza sativa</i>）核心组学 + 育种分析技能包 + CNS 级绘图</h3>

<b>Transcriptome · Metabolome · Integration · GWAS · Causal inference · Population genetics ·<br>
Pangenome · QTL · Proteomics · Epigenomics · Phenomics · Single-cell · Spatial · Editing · Workflows</b>

<br>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Wiki](https://img.shields.io/badge/wiki-Rice--skills-blue)](https://github.com/shuke-bio/Rice-skills/wiki)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23586529.svg)](https://doi.org/10.5281/zenodo.23586529)
[![Download ZIP](https://img.shields.io/badge/download-ZIP-success.svg)](https://github.com/shuke-bio/Rice-skills/archive/refs/tags/v3.12.0.zip)
![Version](https://img.shields.io/badge/version-3.12.0-brightgreen)
![Modules](https://img.shields.io/badge/modules-18-purple)
![Scripts](https://img.shields.io/badge/scripts-130-teal)
![Skills](https://img.shields.io/badge/curated%20skills-541-orange)
![R >= 4.3](https://img.shields.io/badge/R-%3E%3D4.3-276DC3)
![Python >= 3.10](https://img.shields.io/badge/Python-%3E%3D3.10-3776AB)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)

</div>

---

## 📑 Contents

- [What you can do](#-what-you-can-do)
- [Requirements](#-requirements)
- [Install](#-install)
- [Quick start — just talk to your AI](#-quick-start--just-talk-to-your-ai)
- [What's inside](#-whats-inside)
- [Curated third-party skills (541)](#-curated-third-party-skills-541-skillmd--4-tooldata-components)
- [Testing](#-testing)
- [Troubleshooting](#-troubleshooting)
- [License](#-license)
- [Documentation](#-documentation)

## 🎯 What you can do

**Transcriptome**
- FASTQ → counts (fastp → HISAT2 → featureCounts) and count matrix → differential expression (QC → DESeq2 / edgeR)
- ComBat batch correction and time-course models (`group × time` F-test)
- Volcano / heatmap / batch-aware models (batch as a covariate) / time-course contrasts via design
- GO / KEGG / Reactome enrichment, GSEA and ssGSEA/GSVA pathway activity
- Isoform-level differential transcript usage (DTU)
- WGCNA co-expression modules and hub genes
- Transcription-factor identification (PlantTFDB) and TF–target networks; PPI via external STRING import

**Metabolome**
- Preprocessing, missing-value handling, normalization and QC
- OPLS-DA with VIP-based marker selection
- KEGG metabolite pathway enrichment
- Biomarker discovery (ROC + random forest) and metabolite–trait networks
- SIRIUS formula/structure execution + GNPS upload guide; annotation merge
- GEMMA AI-REML mGWAS (kinship + LMM)

**Multi-omics integration**
- MOFA latent factors; DIABLO discriminant integration (O2PLS-style PLS is an approximation — see script note)
- Gene–metabolite correlation networks, joint KEGG pathways and Circos visualization
- GRNBoost2-style TF regulatory network inference
- MOFA factor interpretation (variance/loadings/trait correlation); DIABLO tuning

**Variant calling & annotation**
- Full chain: BWA alignment → GATK variant calling → VCF filtering → SnpEff annotation
- Consensus calling (bcftools + FreeBayes + GATK), DeepVariant, phasing, SV/CNV

**GWAS & genomics**
- PLINK QC → rMVP / GAPIT / built-in GLM association tests
- Publication-grade Manhattan / QQ plots and significant-locus annotation
- Meta-GWAS · haplotype (EM) analysis · 2D pairwise-trait GWAS
- Bayesian colocalization (PP4) · power and sample-size analysis
- **Multi-locus QTNs / QTN×environment (QEI) / epistasis (QTN×QTN)** — 3VmrMLM / Fast3VmrMLM
- **Binary (case-control) multi-locus GWAS** — FastBiCmrMLM (Linux / WSL)
- **Multi-parent population QTL** — MAGIC / NAM / ROAM (mppQTL)
- **TWAS with LD**, two-sample MR (IVW/Egger/weighted median + MR-PRESSO), SuSiE/ABF fine-mapping
- mGWAS/mQTL, k-mer GWAS, rare-variant burden, MTAG plan, FUSION weights, CIM, QTL epistasis, BSA ED/SL

**Population genetics**
- PCA structure + kinship · LD decay · Fst / π · ADMIXTURE · phylogeny
- XP-EHH, iHS/CLR and Dsuite introgression plans

**Pangenome & QTL**
- SV / PAV detection and Pan-GWAS — three routes: assembly (minimap2 + SYRI),
  graph (minigraph-cactus VCF genotyping) and read-based (HitSV short-read /
  long-read / hybrid calling with genotyping)
- eQTL / mQTL / pQTL mapping and GWAS×QTL colocalization
- T2T assembly QC, vg/PanGenie graph genotyping, PHG, odgi graph visualization
- Genome-wide BSA / QTL-seq scan (ΔSNP-index + G′, simulation-based confidence)
- QTL/BSA interval SNP extraction with gene-region annotation (promoter/UTR/exon/intron) — region VCF + 3-sheet Excel

**Proteomics**
- Preprocessing and differential abundance; UniProt protein→gene mapping automation
- InterProScan domain enrichment; PTM enrichment, STRING PPI and AlphaFold DB query

**Epigenomics**
- WGBS methylation DMR analysis; ATAC / ChIP differential regions
- Hi-C TAD/loops plan, scATAC (TF-IDF/LSI), ArchR/Signac plans

**Phenomics & breeding**
- BLUP / heritability; spectral indices (NDVI …)
- GBLUP genomic selection with cross-validation; AMMI / GGE biplots
- DL-GS baseline (Ridge/ElasticNet/MLP/CNN), multi-trait GBLUP, temporal GWAS

**Single-cell & spatial**
- Scanpy pipeline (clustering / marker genes / PAGA) plus integration/annotation and ligand–receptor communication
- scVI/CellChat plans (optional dependencies)
- Squidpy spatial domains and Moran's I; RSTA / scGPT foundation-model integration guide
- NNLS deconvolution and cell2location/COMMOT plans

**Evolution & regulation**
- KA/KS gene-family evolution
- Promoter k-mer enrichment; motif→TF matching; miRNA target prediction
- Alternative-splicing analysis playbook (rMATS / SUPPA2) and small-RNA pipeline (miRDeep2)

**Design & lab**
- qPCR primer design; CRISPR sgRNA design with in-sequence off-target statistics
- Base/prime editing design, edit-efficiency ranking, genome-wide off-target (≤50 Mb built-in)
- **qPCR relative quantification** (2^-ddCt, per-gene bar plots with SD) and
  qPCR↔RNA-seq correlation validation
- Lab protocols: knockout constructs, DNA extraction, edit detection, LUC/GUS assays …

**Automatic by default (no extra instruction needed)**
- CNS-grade figures — PDF / TIFF / PNG at 300 dpi with one shared palette,
  auto-generated by every analysis script (plotting skills auto-routed by chart type)
- Result interpretation & summary (key findings, statistical wording, conclusion boundaries)

**On demand (only when you ask — never force-routed)**
- Literature search (multi-source), lawful full-text download, structured paper reading
  - **Built-in sources**: PubMed (E-utilities + MeSH), CrossRef (DOI / citation checks), arXiv, OpenAlex, bioRxiv; Scopus / ScienceDirect require local `pybliometrics` configuration
  - Semantic Scholar / medRxiv are secondary leads only and are **not** automated sources
  - Last-resort routes: Google Scholar · Web of Science · Scopus · CNKI (institutional access or manual; scraped sources are flagged as unstable when used)
  - Full text only via lawful channels: publisher APIs (Elsevier / Springer Nature / IEEE, your own key), open access (PMC / Unpaywall), CNKI institutional sessions — no Sci-Hub-style grey routes
- Manuscript drafting & polishing (Nature-style), citation formatting & reference verification
- Statistical-reporting review, pre-submission reviewer simulation, point-by-point response letters
- Paper-to-PPT, paper translation, proposals, patent drafts, experiment logs
- Literature monitoring pipelines (daily/weekly tracking of a research topic)

## 🖥 Requirements

| Component | Version | Needed for |
|---|---|---|
| R | ≥ 4.3 (4.6 tested) | all `*.R` scripts (base stats; Bioconductor packages installed on demand) |
| Python | ≥ 3.10 (3.14 tested) | all `*.py` scripts (numpy / pandas / scipy / scikit-learn / matplotlib) |
| Conda | any recent | recommended for one-command environment setup |
| External tools | PLINK, fastp, HISAT2, GATK, InterProScan, … | optional — only the modules that wrap them; core analysis and every `--demo` run work without them |

No GPU required. Everything runs on a laptop-class machine; the heaviest
single-cell / pangenome steps scale with your data, not with the package.

## 📦 Install

**Rice-skills** is a plain folder of skills + scripts — pick any of these routes:

**A. skills CLI (Codex / Claude Code / OpenClaw …)**

```bash
npx skills add shuke-bio/Rice-skills --list          # browse available skills
npx skills add shuke-bio/Rice-skills --skill rice-skills --yes --copy
# or install as a project-level skill with an AGENTS.md binding:
# npx rice-skills install --project --write-agents
```

**B. Git clone**

```bash
git clone https://github.com/shuke-bio/Rice-skills.git
cd Rice-skills
conda env create -f environment.yml && conda activate rice-skills   # optional, for running scripts
```

**C. Point your agent at it** — paste the repo link into Claude / Codex / any
agentic CLI and say *"install the rice skills from this repo"*; or hit the
**download ZIP** badge at the top for the current source snapshot and unzip it
into your agent's skills directory.

Then run the pipeline exactly like any local tool — see [Quick start](#-quick-start--just-talk-to-your-ai).

## 💬 Quick start — just talk to your AI

Once installed, these prompts work out of the box:

| What you want | Say this to your agent |
|---|---|
| Run a full transcriptome + metabolome pipeline | `用仓库自带示例数据跑一遍 run_pipeline.sh 全流程，看每个模块的输出结构。` |
| Do a GWAS on my own data | `我有基因型CSV、SNP信息和表型CSV，按 docs/genomics_usage.md 帮我跑完整GWAS流程。` |
| Design sgRNA + primers for my gene | `给 demo_gene.fa 设计 CRISPR sgRNA 和 qPCR 引物，输出到 results/design/。` |
| Test Bayesian colocalization | `运行 scripts/qtl/03_coloc_bayesian.R --demo 给我看 PP4 输出格式。` |
| Locate QTLs from two bulk pools (BSA) | `用 05_bsa_scan.py 对这两个混池 VCF 做全基因组扫描，输出候选区间和 ΔSNP-index 曲线。` |
| Zoom into a QTL interval | `用 04_qtl_region_snps.py 提取 Chr7:20-25Mb 区间的双亲差异 SNP，标注基因区域，出 Excel。` |
| Estimate my GWAS sample size | `我的 GWAS 大约 300 份材料，用 05_power_analysis.R 估计能检出多大的效应。` |
| Predict miRNA targets | `用 02_mirna_targets.py 把这批 miRNA 对水稻 CDS 做靶基因预测。` |
| Match enriched motifs to TFs | `01_promoter_kmer.py 找到的 top k-mer，用 02_motif_tf_match.py 匹配已知 TF。` |
| Build a regulatory network | `用 07_grn_boost.py 从表达矩阵推断 TF 调控网络，top 边输出。` |
| Summarize my results | `根据 results/ 里的差异表达和富集结果写一份结果解读，指出关键发现和结论边界。` |
| Find literature on a gene | `查一下 OsGATA 近期高影响力文献，按 SKILL.md §5.3 走检索→下载→精读。` |
| Polish my manuscript | `把这份中文稿的 Discussion 按 nature-polishing 规范润色成 Nature 风格英文。` |
| Prepare for submission | `投稿前帮我校验参考文献、审查统计报告、做一次预审稿模拟。` |

Not sure which skill to use? Just describe the task — the SKILL.md trigger list
(rice transcriptome / GWAS / pangenome / CRISPR / 多组学联合分析 …) routes
automatically.

## 🧭 What's inside

<details open>
<summary><b>Module map (18 modules, 130 R/Python scripts + 23 shell)</b></summary>

- **Module 1 · Transcriptome**: FASTQ→counts → QC → DESeq2/edgeR DE → GO/KEGG/Reactome → GSEA/ssGSEA →
  WGCNA → TF ranking (PlantTFDB/PURE-like) → PPI; ComBat + time-course, DTU
- **Module 2 · Metabolome**: imputation/normalization → OPLS-DA + VIP → KEGG →
  biomarker (ROC + RF) → metabolite network; mGWAS/mQTL, SIRIUS, GEMMA AI-REML
- **Module 3 · Integration**: MOFA (+ factor interpretation) → DIABLO/O2PLS (+ tuning) → joint KEGG → Circos
- **Module 7 · GWAS**: PLINK QC → rMVP/GAPIT/GLM → Manhattan/QQ → loci annotation ·
  multi-locus QTN / QEI / epistasis · case-control · multi-parent QTL ·
  TWAS (LD-aware) · MR (+MR-PRESSO) · fine-mapping · mGWAS · k-mer GWAS · burden · MTAG plan
- **Module 11 · Pangenome**: graph/SV/PAV/Pan-GWAS · Oryza genus workflow · SINE annotation ·
  T2T QC · vg/PanGenie genotyping · PHG · odgi visualization
- **Module 16 · Phenomics**: BLUP/GBLUP/AMMI · multi-environment GP · RGB phenotyping ·
  chloroplast counting · DL-GS baseline · multi-trait GBLUP · temporal GWAS
- **Module 8 · Population genetics**: PCA + VanRaden kinship → LD decay →
  Fst/π → ADMIXTURE → phylogeny → XP-EHH
- **Modules 9–18**: CRISPR/base/prime editing + off-target · consensus/DeepVariant/phasing ·
  pangenome SV/PAV · QTL/BSA/ED-SL · proteomics (PTM/PPI/AlphaFold) · epigenomics
  (WGBS/ATAC/ChIP/Hi-C/scATAC) · phenomics · single-cell (Scanpy/scVI/CellChat) ·
  spatial (Squidpy/NNLS/cell2location) · KA/KS evolution
- **v3.2 additions**: Bayesian colocalization (self-implemented, dependency-free) ·
  2D-GWAS · GRNBoost2-style networks · power analysis · UniProt mapping ·
  InterProScan enrichment · miRNA targets · motif→TF matching
- **v3.3.1 additions**: multi-locus QTN / QTN×environment (QEI) / marker-pair
  epistasis (3VmrMLM family) · binary case-control GWAS (FastBiCmrMLM) ·
  multi-parent population QTL (mppQTL) — all runtime-installed, method notes in
  `docs/references/`
- **v3.3.2 additions**: literature-driven wrappers — MeNet (multi-environment GS) ·
  ELMERF (salt-stress RGB phenotyping) · Chloro-Count · AnnoSINE · Oryza genus
  pangenome stage runner — plus a 70-method distilled index graded S/A1/A2/A3/B
  (`docs/references/literature_methods_2022_2026.md`)
- **v3.3.3 additions**: decision layer — data contract, 10 scenario cards, three-tier
  plan generator with blockers (`scripts/plan/01_make_plan.py`), data QC
  (`02_check_data.py`), 22 validated method cards and 10 English scenario cards;
  citations restricted to the A1–A3 journal whitelist (`docs/decision/`)
- **v3.3.4 additions**: dependency preflight (`scripts/plan/04_check_deps.py`) — scans the
  scripts in a chosen plan tier for R/Python/tool dependencies, prints exact install
  commands, installs only after explicit consent, re-checks and logs; system tools are
  never auto-installed
- **v3.3.5 additions**: plan & figure approval gate — each tier lists the figures
  (type, 85/170 mm, 300 dpi, PDF+TIFF+PNG, script); `scripts/plan/05_run_plan.py`
  prints the confirmation sheet and refuses to run until approved (`approval.json`
  with a plan hash), then checks dependencies before executing demo-capable scripts
- **v3.4–v3.12 additions**: off-target/TWAS/MR/mGWAS/k-mer/burden · T2T/graph genotype/PHG ·
  single-cell/spatial upgrades · editing design · workflow export/HPC/CI ·
  now 130 R/Python + 23 shell
- **v3.3.6 change**: default figure formats are now **PDF / TIFF / PNG** (previously
  PDF/SVG/TIFF) across the shared R and Python helpers, 19 plotting scripts and the
  configuration; SVG remains available only on explicit request
- **v3.4.0 additions**: genome-wide CRISPR off-target search (`design/04_offtarget.py`),
  TWAS (`qtl/06_twas.R`) and two-sample Mendelian randomization (`qtl/07_mr.R`); all with
  synthetic demos and explicit method boundaries
- **v3.4.1 additions**: metabolite GWAS/mQTL (`metabolomics/08_mgwas.R`) and
  k-mer dosage GWAS (`gwas/12_kmer_gwas.py`, KMERIA-style with quantile dosage, PCA,
  kinship and eigen-transformed MLM; significant k-mers map back to the genome)
- **v3.4.2 additions**: T2T assembly QC (`pangenome/07_t2t_qc.sh`), graph pangenome
  genotyping plan/runner (`pangenome/08_graph_genotype.sh`, vg Giraffe + PanGenie) and a
  DL genomic-selection baseline (`phenomics/09_dl_gs.py`, Ridge/ElasticNet/MLP/CNN CV)
- **v3.5.0 additions**: PURE-like TF prioritization (`transcriptomics/08_tf_prioritize.py`),
  consensus variant calling (`variants/05_consensus_call.sh`) and fine-mapping
  (`qtl/08_finemap_susie.R`, SuSiE when installed, approximate ABF fallback)
- **v3.6.0 additions**: single-cell integration/annotation (`single_cell/02_integration_annotation.py`),
  cell-cell communication (`single_cell/03_cell_communication.py`), spatial NNLS deconvolution
  (`spatial/02_deconvolution.py`), SV/CNV plan/runner (`variants/06_sv_cnv.sh`), iHS/CLR+Dsuite
  (`popgen/08_selection_scan.sh`), metabolite annotation merge (`metabolomics/09_annotate.py`)
  and DTU (`transcriptomics/09_dtu.py`)
- **v3.7.0 additions**: base/prime editing design (`design/05_base_prime_edit.py`), Hi-C plan
  (`epigenomics/03_hic_analysis.sh`), scATAC (`epigenomics/04_scatac.py`), workflow export
  (`utils/export_workflow.py`), benchmark suite (`utils/benchmark_suite.py`), HPC submit
  (`utils/hpc_submit.py`), causal network (`integration/08_causal_network.py`) and temporal
  GWAS (`phenomics/10_temporal_gwas.py`)
- **v3.7.1 quick wins**: sgRNA `--genome` genome-wide off-target, GS `--method rrblup|bayesb`,
  Scanpy `--celltypist-model`/`--doublet`, spatial `--deconv-*`, and README script-count checks
  in `release_check.py`
- **v3.8.0 additions**: FASTQ→counts (`transcriptomics/00_fastq_counts.sh`), ssGSEA
  (`transcriptomics/10_ssgsea.R`), PTM/PPI (`proteomics/06_ptm_ppi.py`), import builder
  (`utils/build_skill_import.py`), CI workflow, and runtime-dependency/licence inventory
- **v3.9.0 additions**: ComBat + time-course model (`transcriptomics/11_batch_timecourse.R`),
  LD-aware TWAS (`qtl/06_twas.R --ld`), rare-variant burden (`gwas/13_burden_test.py`) and
  QTL epistasis (`qtl/10_qtl_epistasis.R`)
- **v3.10.0 additions**: SIRIUS/GEMMA execution wrappers, BSA ED/SL, AlphaFold DB query
  and multi-trait GBLUP summary
- **v3.11.0 additions**: MOFA factor interpretation, DIABLO tuning plan, scVI/CellChat,
  cell2location/COMMOT, DeepVariant and phasing wrappers; plan generator now includes
  the new modules
- **v3.12.0 additions**: FUSION weight training, CIM, PHG, ArchR/Signac plans, edit-efficiency
  heuristic ranking and MTAG plan; MR now supports `--presso`
- **P2 also includes** pangenome graph visualization (`pangenome/09_graph_visualize.sh`, odgi)
- **Molecular lab companion**: rice CRISPR knockout reference + molecular protocols
</details>

<details open>
<summary><b>Repository structure</b></summary>

```
Rice-skills/
├── SKILL.md                        # Skill manifest (frontmatter + full usage)
├── environment.yml                 # Conda env (name: rice-skills)
├── run_pipeline.sh                 # One-command transcriptome+metabolome run
├── config/nature_palette.json      # Shared Nature-style palette (R + Python)
├── assets/                         # Cover image + bundled annotation snapshots
├── demo_data/                      # Example demo data for every module
├── scripts/                        # 130 R/Python scripts + 23 shell wrappers, 18 modules
├── templates/                      # ggplot2/matplotlib figure templates + Rmd report
├── examples/                       # Input format specs + metadata templates
├── third_party/                    # 541 curated skills + palettes, 6 skill libraries + 4 tool/data components (NOTICE.md)
├── docs/                           # Per-module usage guides & lab protocols (Chinese)
└── tests/                          # Pipeline self-check + utility smoke tests
```

</details>

## 🧠 Curated third-party skills (541 SKILL.md + 4 tool/data components)

| Directory | Upstream | License | Size |
|---|---|---|---|
| `bioskills/` | [GPTomics/bioSkills](https://github.com/GPTomics/bioSkills) | MIT | 354 skills |
| `omicverse-skills/` | [omicverse/omicverse-skills](https://github.com/omicverse/omicverse-skills) | MIT docs | 61 skills |
| `sciagent-skills/` | [jaechang-hits/SciAgent-Skills](https://github.com/jaechang-hits/SciAgent-Skills) | CC BY 4.0 | 64 skills |
| `kdense-scientific-skills/` | [K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills) | MIT | 28 skills |
| `cell-skills/` | [yrui-cmd/Cell](https://github.com/yrui-cmd/Cell) + CNS_Figure | MIT | 12 skills |
| `nature-skills/` | [Yuan1z0825/nature-skills](https://github.com/Yuan1z0825/nature-skills) | Apache-2.0 | 22 skills |
| `colors/` | [lcpmgh/colors](https://github.com/lcpmgh/colors) | MIT | 229 palettes + Shiny app |
| `colormaps/` | [1313e/CMasher](https://github.com/1313e/CMasher) + [holoviz/colorcet](https://github.com/holoviz/colorcet) | BSD-3-Clause / CC BY 4.0 | 133 scientific colormaps (hex data) |
| `ltc/` | [loukesio/ltc-color-palettes](https://github.com/loukesio/ltc-color-palettes) | MIT | 32 palettes (hex data) |
| `ggpubfigs/` | [JLSteenwyk/ggpubfigs](https://github.com/JLSteenwyk/ggpubfigs) | MIT | 15 colorblind-friendly palettes (hex data) |

Skills are **methodology references only** — the pipeline never depends on them.
Load the relevant `SKILL.md` for expert guidance via the routing table in
[`docs/skill_routing.md`](docs/skill_routing.md).
The four color components are plotting aids, not skills: `colors/` holds 229
publication-grade discrete palettes plus a Shiny picker; `ltc/` (32) and
`ggpubfigs/` (15) add discrete palettes; `colormaps/` holds 133 scientific
continuous colormaps. Selection guides and R/Python usage are in
`third_party/colormaps/README.md` and the distilled method note
[`docs/references/plot_color_palettes.md`](docs/references/plot_color_palettes.md).

## 🧪 Testing

```bash
bash tests/test_pipeline.sh    # end-to-end smoke test
Rscript tests/test_utils.R     # R utility layer
python tests/test_utils.py     # Python utility layer
```

## 🛠 Troubleshooting

| Symptom | Fix |
|---|---|
| Download fails with `self-signed certificate` / TLS errors | Your network intercepts TLS. Use the `--insecure` flag (annotation fetch, UniProt mapping), or download the file in a browser and pass it via `--tsv-file` / `--local`. |
| A Bioconductor R package is missing | Install on demand: `BiocManager::install("DESeq2")` — scripts report exactly which package they need. The bundled annotation snapshot (`assets/annotations/`) keeps core pipelines working offline. |
| Not sure whether a script works on your machine | Run its self-test first: most scripts accept `--demo` (e.g. `python scripts/integration/07_grn_boost.py --demo`) — a demo run on built-in example data, zero input needed. |
| Figures look different from the paper style | Ensure `config/nature_palette.json` is present and scripts are run from the repository root (the theme path is repo-relative). |
| Everything is slow on first run | First runs compile/load heavy stacks (tidyverse, Bioconductor); subsequent runs are much faster. Pin `--seed` for identical results. |

## 📄 License

MIT — see [LICENSE](LICENSE). Third-party components keep their own licenses;
see [`third_party/NOTICE.md`](third_party/NOTICE.md) and per-library files.

**Cite this work.** Fu, Liang. *Rice-skills* (2026). Concept DOI (all versions): [10.5281/zenodo.23586529](https://doi.org/10.5281/zenodo.23586529). Machine-readable metadata: [`CITATION.cff`](CITATION.cff).

## 📬 Contact

**Questions or issues?** Open a [GitHub issue](https://github.com/shuke-bio/Rice-skills/issues) first — public and searchable, so the answer helps other users too. For anything that should not be public, open an issue without sensitive details and we will follow up.

## 📚 Documentation

Per-module guides and lab protocols live in [`docs/`](docs/) (Chinese):
[installation](docs/installation.md) · [pipeline usage](docs/usage.md) ·
[GWAS & popgen](docs/genomics_usage.md) · [design](docs/design_usage.md) ·
[pangenome](docs/pangenome_usage.md) · [QTL](docs/qtl_usage.md) ·
[skill routing](docs/skill_routing.md) · [tool navigation](docs/tools_navigation.md) · [promoter/TF usage](docs/tf_promoter_usage.md).

<div align="right">

[English](README.md) | [简体中文](README.zh-CN.md)

</div>
