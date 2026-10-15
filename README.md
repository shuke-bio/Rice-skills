<div align="center">

<img src="assets/cover.jpg" alt="Rice-skills" width="100%">

**English** | [简体中文](README.zh-CN.md)

<h1>Rice-skills</h1>

<h3>水稻（<i>Oryza sativa</i>）核心组学 + 育种分析技能包 + CNS 级绘图</h3>

<b>Transcriptome · Metabolome · Integration · GWAS · Population genetics · Pangenome · QTL ·<br>
Proteomics · Epigenomics · Phenomics · Single-cell · Spatial · Evolution · CRISPR design</b>

<br>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Version](https://img.shields.io/badge/version-3.2.0-brightgreen)
![Modules](https://img.shields.io/badge/modules-18-purple)
![Scripts](https://img.shields.io/badge/scripts-76-teal)
![Skills](https://img.shields.io/badge/curated%20skills-528-orange)
![R >= 4.3](https://img.shields.io/badge/R-%3E%3D4.3-276DC3)
![Python >= 3.10](https://img.shields.io/badge/Python-%3E%3D3.10-3776AB)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)

</div>

---

## 📑 Contents

- [Requirements](#-requirements)
- [Install](#-install)
- [Quick start — just talk to your AI](#-quick-start--just-talk-to-your-ai)
- [What's inside](#-whats-inside)
- [Curated third-party skills (528)](#-curated-third-party-skills-528)
- [Testing](#-testing)
- [Troubleshooting](#-troubleshooting)
- [License](#-license)
- [Documentation](#-documentation)

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
```

**B. Git clone**

```bash
git clone https://github.com/shuke-bio/Rice-skills.git
cd Rice-skills
conda env create -f environment.yml && conda activate rice-skills   # optional, for running scripts
```

**C. Point your agent at it** — paste the repo link into Claude / Codex / any
agentic CLI and say *"install the rice skills from this repo"*; or download the
ZIP from Releases and unzip into your agent's skills directory.

Then run the pipeline exactly like any local tool — see [Quick start](#-quick-start--just-talk-to-your-ai).

## 💬 Quick start — just talk to your AI

Once installed, these prompts work out of the box:

| What you want | Say this to your agent |
|---|---|
| Run a full transcriptome + metabolome pipeline | `用合成数据跑一遍 run_pipeline.sh 全流程，看每个模块的输出结构。` |
| Do a GWAS on my own data | `我有基因型CSV、SNP信息和表型CSV，按 docs/genomics_usage.md 帮我跑完整GWAS流程。` |
| Design sgRNA + primers for my gene | `给 demo_gene.fa 设计 CRISPR sgRNA 和 qPCR 引物，输出到 results/design/。` |
| Test Bayesian colocalization | `运行 scripts/qtl/03_coloc_bayesian.R --demo 给我看 PP4 输出格式。` |
| Estimate my GWAS sample size | `我的 GWAS 大约 300 份材料，用 05_power_analysis.R 估计能检出多大的效应。` |
| Predict miRNA targets | `用 02_mirna_targets.py 把这批 miRNA 对水稻 CDS 做靶基因预测。` |
| Match enriched motifs to TFs | `01_promoter_kmer.py 找到的 top k-mer，用 02_motif_tf_match.py 匹配已知 TF。` |
| Build a regulatory network | `用 07_grn_boost.py 从表达矩阵推断 TF 调控网络，top 边输出。` |

Not sure which skill to use? Just describe the task — the SKILL.md trigger list
(rice transcriptome / GWAS / pangenome / CRISPR / 多组学联合分析 …) routes
automatically.

## 🧭 What's inside

<details open>
<summary><b>Module map (18 modules, 76 scripts)</b></summary>

- **Module 1 · Transcriptome**: QC → DESeq2/edgeR DE → GO/KEGG/Reactome → GSEA →
  WGCNA → TF (PlantTFDB) → PPI
- **Module 2 · Metabolome**: imputation/normalization → OPLS-DA + VIP → KEGG →
  biomarker (ROC + RF) → metabolite network
- **Module 3 · Integration**: MOFA → DIABLO/O2PLS → joint KEGG → Circos
- **Module 7 · GWAS**: PLINK QC → rMVP/GAPIT/GLM → Manhattan/QQ → loci annotation
- **Module 8 · Population genetics**: PCA + VanRaden kinship → LD decay →
  Fst/π → ADMIXTURE → phylogeny → XP-EHH
- **Modules 9–18**: CRISPR/primer design · variant calling (BWA→GATK→SnpEff) ·
  pangenome SV/PAV · QTL (eQTL/mQTL) · proteomics · epigenomics (WGBS/ATAC/ChIP) ·
  phenomics (BLUP/GBLUP/AMMI) · single-cell (Scanpy) · spatial (Squidpy) ·
  KA/KS evolution
- **v3.2 additions**: Bayesian colocalization (self-implemented, dependency-free) ·
  2D-GWAS · GRNBoost2-style networks · power analysis · UniProt mapping ·
  InterProScan enrichment · miRNA targets · motif→TF matching
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
├── demo_data/                      # Synthetic demo data for every module
├── scripts/                        # 76 R/Python analysis scripts (18 modules)
├── templates/                      # ggplot2/matplotlib figure templates + Rmd report
├── examples/                       # Input format specs + metadata templates
├── third_party/                    # 528 curated skills, 6 upstream libraries (NOTICE.md)
├── docs/                           # Per-module usage guides & lab protocols (Chinese)
└── tests/                          # Pipeline self-check + utility smoke tests
```

</details>

## 🧠 Curated third-party skills (528 SKILL.md)

| Directory | Upstream | License | Size |
|---|---|---|---|
| `bioskills/` | [GPTomics/bioSkills](https://github.com/GPTomics/bioSkills) | MIT | 354 skills |
| `omicverse-skills/` | [omicverse/omicverse-skills](https://github.com/omicverse/omicverse-skills) | MIT docs | 61 skills |
| `sciagent-skills/` | [jaechang-hits/SciAgent-Skills](https://github.com/jaechang-hits/SciAgent-Skills) | CC BY 4.0 | 64 skills |
| `kdense-scientific-skills/` | [K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills) | MIT | 28 skills |
| `cell-skills/` | [yrui-cmd/Cell](https://github.com/yrui-cmd/Cell) + CNS_Figure | MIT | 12 skills |
| `nature-skills/` | [xiaolai/nature-skills](https://github.com/xiaolai/nature-skills) | MIT | 9 skills |

Skills are **methodology references only** — the pipeline never depends on them.
Load the relevant `SKILL.md` for expert guidance via the routing table in
[`docs/skill_routing.md`](docs/skill_routing.md).

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
| Not sure whether a script works on your machine | Run its self-test first: most scripts accept `--demo` (e.g. `python scripts/integration/07_grn_boost.py --demo`) — a synthetic run with no inputs needed. |
| Figures look different from the paper style | Ensure `config/nature_palette.json` is present and scripts are run from the repository root (the theme path is repo-relative). |
| Everything is slow on first run | First runs compile/load heavy stacks (tidyverse, Bioconductor); subsequent runs are much faster. Pin `--seed` for identical results. |

## 📄 License

MIT — see [LICENSE](LICENSE). Third-party components keep their own licenses;
see [`third_party/NOTICE.md`](third_party/NOTICE.md) and per-library files.

## 📚 Documentation

Per-module guides and lab protocols live in [`docs/`](docs/) (Chinese):
[installation](docs/installation.md) · [pipeline usage](docs/usage.md) ·
[GWAS & popgen](docs/genomics_usage.md) · [design](docs/design_usage.md) ·
[pangenome](docs/pangenome_usage.md) · [QTL](docs/qtl_usage.md) ·
[skill routing](docs/skill_routing.md).

<div align="right">

[English](README.md) | [简体中文](README.zh-CN.md)

</div>
