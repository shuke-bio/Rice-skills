# Changelog

All notable changes to Rice-skills are documented here.
Versioning follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [3.3.6] — 2026-10-23

### Changed
- **Default figure formats switched to PDF / TIFF / PNG** (was PDF / SVG / TIFF):
  shared helpers `scripts/utils/plot_theme.R` / `plot_theme.py`, `config/nature_palette.json`,
  19 plotting scripts with inline save loops, plus docs and SKILL wording.
  PNG is written at 300 dpi; SVG is no longer emitted by default (still supported when
  requested explicitly via the R helper's `formats` argument).
- Restored TIFF/PNG completeness in `design/03_qpcr_analysis.py` and
  `qtl/05_bsa_scan.py` (previously PDF/PNG only).

## [3.3.5] — 2026-10-23

### Added
- **Plan & figure approval gate** — every tier now carries a figure plan (name, type,
  85/170 mm, 300 dpi, PDF+SVG+TIFF, generating script); `scripts/plan/05_run_plan.py`
  prints the confirmation sheet and refuses to execute until the user approves
  (`--approve` writes `approval.json` with a plan hash; `--execute` verifies it,
  runs the dependency preflight, then executes demo-capable scripts and skips the rest
  with an explicit reason).
- `docs/decision/figure_plan.md`: confirmation sheet contents, figure standards and
  approval/execution records.

## [3.3.4] — 2026-10-23

### Added
- **Dependency preflight with consent gating** — `scripts/plan/04_check_deps.py` scans the
  scripts referenced by a plan (or given directly), checks Python imports, R packages and
  external tools, and reports what is missing with exact install commands.
  Default is check-only; installation requires `--install --yes` (non-interactive refusal
  exits with code 2), re-checks afterwards and logs to `results/deps/install.log`.
  System tools (PLINK/BLAST/TRF/Snakemake/…) are never auto-installed.
- `docs/decision/dependencies.md`: flow, safety boundaries, what is/isn't installed,
  detection rules and FAQ.

## [3.3.3] — 2026-10-21

### Added
- **Scenario cards for the remaining high-frequency designs**: single-cell & spatial,
  metabolome, QTL/BSA, phenomics/GS, pangenome, epigenome and proteome; each card lists
  prerequisites, three tiers, failure modes and whitelist-only sources.
- `scripts/plan/02_check_data.py`: pre-analysis health check (missingness, library size,
  constant features, PCA outliers, batch×group confounding, multi-omics overlap);
  `--demo` plants four defect classes and asserts they are detected.
- `docs/decision/router.md` + `template.json`: intent-to-question routing and a ready
  data contract.
- `docs/decision/cards/`: 22 machine-readable method cards, validated by
  `scripts/plan/03_validate_cards.py` (required fields, script existence, journal whitelist,
  DOI format); includes ELMERF, AnnoSINE, Chloro-Count, GBLUP, AMMI/GGE,
  FastBiCmrMLM and the Oryza pangenome workflow.
- `docs/decision/en/`: English versions of all 10 scenario cards, same structure and
  whitelist-only citation policy.

### Added
- **Decision layer for analysis planning** — `scripts/plan/01_make_plan.py` turns a
  10-field data contract (JSON) into blockers + three tiers (minimal / standard /
  advanced), each mapped to real scripts with method citations; `--demo` covers four
  scenarios (RNA-seq n=3, RNA-seq n=1 guardrail, GWAS 300 accessions, multi-omics n=60).
- `docs/decision/`: data contract, scenario cards for bulk RNA-seq, GWAS and
  transcriptome×metabolome integration; citation policy restricted to the A1–A3
  journal whitelist (non-whitelisted method papers are named as tools only).

## [3.3.2] — 2026-10-19

### Added
- **Literature-driven method landings (S/A1/A2/A3/B grading)** — five new wrappers plus a
  full distilled index covering 70 screened papers (2022–2026):
  - `scripts/pangenome/05_oryza_pangenome_pipeline.sh` — genus-scale Oryza pan-genome
    Snakemake stage runner (MIT upstream, not vendored)
  - `scripts/pangenome/06_annosine.sh` — SINE annotation (AnnoSINE, MIT; needs blastn/TRF/node)
  - `scripts/phenomics/06_menet_gs.py` — multi-environment genomic prediction (MeNet;
    academic-use licence, runtime only)
  - `scripts/phenomics/07_elmerf_salt.py` — rice seedling salt-stress RGB phenotyping
    (ELMERF, Apache-2.0)
  - `scripts/phenomics/08_chlorocount.py` — chloroplast/bundle-sheath quantification
    (Chloro-Count, BSD-3-Clause)
  - All five expose `--demo`: demo fixture + dependency report; heavy upstream
    dependencies (torch/mmseg/blastn/TRF/snakemake/weights) yield an explicit SKIP with
    install commands instead of a silent failure.
- `docs/references/literature_methods_2022_2026.md` — distilled index of 70 screened
  papers graded S (6) / A1 (4) / A2 (39) / A3 (5) / B (16), with inputs/outputs, licence
  status and landing route for each.

## [3.3.1] — 2026-10-19

### Added
- **3VmrMLM-family GWAS scripts (three new wrappers)** — closing the multi-locus /
  interaction / binary-trait / multi-parent gaps left by the single-locus engines
  in `03_association.R`:
  - `scripts/gwas/09_multilocus_qtn.R` — multi-locus QTN / QTN×environment (QEI) /
    marker-pair epistasis (QQI) via IIIVmrMLM (CRAN) or Fast3VmrMLM (big-data);
    accepts PLINK prefixes or dosage CSV; `--demo` self-test in all three modes
    (measured 15 s / 22 s / ~6 min on 200 x 400 simulated data).
  - `scripts/gwas/10_case_control_qtn.R` — binary / case-control traits via
    FastBiCmrMLM mixed logistic model (algorithm auto-selected by sample size;
    Linux/WSL only — non-Linux hosts get an explicit SKIP with the alternative path).
  - `scripts/gwas/11_multiparent_qtl.R` — MAGIC / NAM / ROAM multi-parent QTL via
    mppQTL (two-stage multi-locus model); `--demo` simulated 4-way MAGIC self-test
    (150 lines x 240 markers, ~32 s).
  - `scripts/utils/vmrml_utils.R` — shared helpers: PLINK binary writer (verified
    round-trip + independent decoder), demo-data generators, engine result
    normalisation, runtime-install guidance.
- **Four distilled method notes** (`docs/references/`): `3vmrmlm_multilocus_qtn.md`
  (model, engines, Hap/NCII/mQTL modules, mppQTL section, install guide, citations),
  `fastbicmrmlm_case_control.md`, `qei_environment_interaction.md`,
  `epistasis_qqi.md`. Papers cited by DOI; no text copied.
- `docs/genomics_usage.md` v3.3 section with runnable examples; `tools_navigation.md`
  §2 now lists IIIVmrMLM / Fast3VmrMLM / FastBiCmrMLM / mppQTL with verified
  license status (GPL-3 / no LICENSE file — runtime install only, never vendored);
  SKILL.md module table, description and trigger words updated.

### Notes
- Measured on R 4.6.0 / Windows (2026-10-19): IIIVmrMLM and mppQTL load and run;
  Fast3VmrMLM's R 4.3-built DLL fails to load on R 4.6 (auto-fallback to
  IIIVmrMLM); FastBiCmrMLM has no Windows build (Linux/WSL only).

## [3.3.0] — 2026-10-19

### Added
- **Publication color-palette component (colors)** — vendored
  `third_party/colors/` from **lcpmgh/colors** (v4.0, commit `36b2cc6`,
  MIT, Copyright (c) 2025 Liang Chen; LICENSE原文已核验): Shiny 配色推荐器
  + 229 套论文插图配色（`@palettes.csv`：知乎顶刊精选 177 + ggsci/grDevices/
  RColorBrewer 52，2–51 色/套）。上游 `articles/`（知乎抓取的第三方文章）
  与 `oldversion/` **明确剔除**（MIT 不覆盖），详见 `third_party/colors/NOTICE.md`。
- **Scientific continuous-colormap component (colormaps)** — vendored
  `third_party/colormaps/@colormaps.csv` (**133 张科学连续色图 hex 数据**):
  CMasher 53 张（BSD-3-Clause, © 2019-2021 Ellert van der Velden,
  `1313e/CMasher`, 默认分支 `master`）+ colorcet 80 张（CC BY 4.0,
  © Peter Kovesi, `holoviz/colorcet`）。两上游均活跃
  （末次提交 2026-10-09 / 2026-10-12），LICENSE 原文已核验。字段
  `c_sour/c_name/c_type/n_stops/colors`（`;` 分隔 hex），类型分布
  sequential 78 / diverging 31 / cyclic 19 / rainbow 5；附快速选型表与
  R/Python 取用代码（`third_party/colormaps/README.md`）。**剔除决策**
  （2026-10-18 收录准则：末次提交 ≤ 约 1 年）：crameri / cmocean / scico /
  viridis 长期停更（2023–2024）**不入包**，仅在
  `third_party/colormaps/NOTICE.md` 的「调研后未收录的上游」表留档并附
  运行时安装命令。
- **Discrete palette data components (ltc + ggpubfigs)** — vendored
  `third_party/ltc/@palettes.csv`（**32 套**自然风调色板，MIT,
  `loukesio/ltc-color-palettes`, © 2021 ltc authors，末次提交 2026-08-30）与
  `third_party/ggpubfigs/@palettes.csv`（**15 套**色盲友好出版配色，MIT,
  `JLSteenwyk/ggpubfigs`, © 2020 Jacob L. Steenwyk，末次提交 2026-05-12）。
  二者均由上游源码（R 列表）转录为 hex 数据，LICENSE 原文已核验，
  各附 LICENSE + NOTICE + README。
- `docs/references/plot_color_palettes.md` — 蒸馏方法笔记（自有语言）：配色三类型
  （离散/连续/发散）、NPG/Lancet/NEJM/JAMA/JCO/AAAS 期刊 palette hex 直取、
  六条选色原则（≤8 组、避红绿相邻、灰度可辨、全文一致等）、ggplot2
  `scale_color_manual()` 取用代码、本地 Shiny 运行与安装命令；
  **§8 连续色图五类型方法论**（sequential/diverging/cyclic/rainbow/special
  判定：cyclic=首尾色距 < 0.05、diverging=U/∩型亮度剖面，其余为 sequential；
  上游官方分类交叉校验）+ colormaps 组件选型与 R/Python 取用代码。
- SKILL.md §4/§5.1、docs/tools_navigation.md §0/§6/§7 挂载 colors 配色速查入口；
  `third_party/NOTICE.md` 总索引、README 双语同步四个配色组件
  （colors / colormaps / ltc / ggpubfigs）规模与收录准则。
- **Read-based SV detection route (module 11 · pangenome)** —
  `scripts/pangenome/02_sv_detect.sh --mode read-based` wraps **HitSV**
  (Li G et al., *Nat Biotechnol* 2026, doi:10.1038/s41587-026-03332-9): pure
  short-read (`--seq srs`), pure long-read (`--seq lrs`, `--platform
  HIFI|ONT_Q20|ERR_PRONE|ASM`) and hybrid (`--seq hybrid`) calling + genotyping,
  outputs a genotyped VCF that feeds `03_pav_analysis.R` / `04_pan_gwas.R`.
  Ships a `--demo` end-to-end self-test (300 kb reference, simulated
  ~10× 10 kb HiFi-like reads carrying a known 500 bp DEL + 1 kb INS,
  DEL-overlap verification). **HitSV is GPL-3.0 and deliberately NOT bundled** —
  the script prints runtime-install instructions (GitHub Releases static
  binary) and records the provenance in `run_info.txt`, keeping the MIT
  licensing model of this repo intact.
- `docs/references/hitsv_sv_calling.md` — distilled method note (own wording,
  no paper text/figures copied): SV-detection route decision guide
  (assembly / graph / SRS / LRS / hybrid / ultra-low-depth), HitSV method
  essentials (repeat- & signature-density-aware breakpoint recognition,
  haplotype-resolved local assembly), rice-relevant preset guidance
  (`-P asm10` intraspecies default), integration paths, reporting checklist,
  and references (paper, tool repo, CC BY 4.0 preprint, call-result repos).
- **GAPIT multi-model support (module 7 · GWAS)** — `scripts/gwas/03_association.R`
  gains `--model` (engine `gapit` only): GLM / MLM / CMLM / SUPER / MLMM /
  FarmCPU / Blink, case-insensitive with canonical normalization; default MLM
  keeps existing behavior; the chosen model is recorded in the engine label,
  `qc_summary.txt` and `run_info.txt`. Additions are distilled method notes only
  — **no GAPIT source bundled** (GPL-2.0+ verified via DESCRIPTION 2026-10;
  runtime install `remotes::install_github("jiabowang/GAPIT")`, instructions
  printed by the script).
- `docs/references/gapit_gwas_formats.md` — distilled method note (own wording,
  no third-party tutorial text or source copied): GAPIT input formats
  (Y/GD/GM/K/CV; Numeric vs HapMap vs PLINK vs VCF), VCF→HapMap preprocessing
  (vcftools + TASSEL route), the seven-model decision table (GLM/MLM/CMLM/
  SUPER/MLMM/FarmCPU/Blink with complexity and use-case guidance), output-file
  interpretation, and citation list (GAPIT v1–v4, FarmCPU, BLINK, SUPER, P3D,
  cBLUP/sBLUP papers).

### Changed
- **Literature source coverage documented** — both READMEs and SKILL.md §5.3 now state
  exactly which databases the literature-search chain covers: T1 APIs PubMed
  (E-utilities/MeSH), CrossRef, arXiv; T2 Semantic Scholar / bioRxiv / medRxiv /
  OpenAlex; last-resort Google Scholar / Web of Science / Scopus / CNKI (institutional
  access or manual, flagged as unstable). Full-text download stays on lawful channels
  only (publisher APIs / OA·PMC·Unpaywall / CNKI sessions; no Sci-Hub-style routes).
- **Colormap stop-count claims corrected** — `@colormaps.csv` stop counts are transcribed
  verbatim from upstream `_hex.txt` files (116 maps at 256 stops, cmasher cyclic/diverging
  at 510–511, cmasher `neutral` at 2560); README/NOTICE/method-note wording updated to
  match (no resampling). A known upstream hex typo in `colors/@palettes.csv`
  (`#T3ABCF`, palette id 123) is now documented in `third_party/colors/NOTICE.md`
  instead of being silently kept.
- **Public-facing wording sweep** — bundled demo data is now called
  "示例数据 / 演示数据" everywhere in README (both languages), SKILL.md,
  CITATION.cff, CONTRIBUTING.md and docs/ (previously used "simulated data" wording).
  The demo data itself is unchanged (still program-generated placeholders,
  seed=42, never for biological conclusions); only the shop-window wording
  changed so first-time visitors do not read the package as fake-data-only.
- **`third_party/nature-skills/` upgraded to the author's latest release** — 22 skills
  (was 9), now under Apache-2.0 (repo-level LICENSE, upstream
  [Yuan1z0825/nature-skills](https://github.com/Yuan1z0825/nature-skills) commit
  `7c1cdd2`); demo `assets/` and the off-topic `nature-social-science` skill excluded
  (full provenance in `third_party/nature-skills/NOTICE.md`). New capabilities:
  paper full-text download, structured paper cards / x-ray deep reading, reference
  cross-verification, statistical-reporting review, pre-submission reviewer
  simulation, literature monitoring pipelines, paper translation, proposals,
  patent drafts, experiment logs.
- **Two-tier skill routing (SKILL.md §5)** — **default route (automatic, no user
  instruction needed)**: CNS-grade plotting (§5.1, 40 skills routed by chart type)
  plus result interpretation & summary (§5.2). **On-demand (loaded only when the
  user explicitly asks — never force-routed)**: literature search/download/reading
  (§5.3), manuscript writing / polishing / submission (§5.4), presentation
  conversion (§5.5). §5.6 documents that live web query depends on the host
  agent's tools rather than third_party skills. SKILL.md §3 now shows the default
  chain and the on-demand chain separately.
- Curated third-party skill count 528 → 541; `docs/skill_routing.md` and both
  READMEs (capability lists, prompt tables, license notes) updated accordingly.

## [3.2.0] — 2026-10-15

### Added
- **Bayesian colocalization** (`scripts/qtl/03_coloc_bayesian.R`) — self-implemented
  coloc.abf (approximate Bayes factors, H0–H4 posteriors, per-SNP shared-causal
  posterior), zero extra dependencies.
- **2D-GWAS** (`scripts/gwas/08_2d_association.R`) — pairwise trait-conditioned
  association for epistasis-like background effects.
- **GRN inference** (`scripts/integration/07_grn_boost.py`) — GRNBoost2-style TF
  regulatory network via scikit-learn.
- **Power & sample-size analysis** (`scripts/phenomics/05_power_analysis.R`) —
  t-test / GWAS (non-central χ²) / GS training-population sizing; cross-checked
  against R's `power.test`.
- **UniProt mapping automation** (`scripts/proteomics/04_uniprot_mapping.py`) —
  proteome UP000059680 → protein→gene map, cache + offline fallbacks.
- **InterProScan domain enrichment** (`scripts/proteomics/05_interpro_domains.py`).
- **miRNA target prediction** (`scripts/evolution/02_mirna_targets.py`) —
  psRNATarget-style seed/pairing scoring, dependency-free.
- **motif→TF matching** (`scripts/tf/02_motif_tf_match.py`) — MEME PWM libraries,
  log-odds matching, TF ranking.
- **QTL/BSA interval SNP extraction** (`scripts/qtl/04_qtl_region_snps.py`) —
  two-parent genotype/ΔAF comparison, promoter/UTR/exon/intron gene-region
  annotation, region VCF + 3-sheet Excel output.
- **Genome-wide BSA / QTL-seq scan** (`scripts/qtl/05_bsa_scan.py`) —
  ΔSNP-index + G′ sliding windows with simulation-based confidence lines,
  candidate-region merging, genome-wide two-track figure; completes the
  BSA chain upstream of 04.
- Tool navigation (`docs/tools_navigation.md`), promoter/TF usage guide
  (`docs/tf_promoter_usage.md`), v3.2 smoke tests in `tests/test_pipeline.sh`.
- `docs/skill_routing.md` — 528 curated skills routed to analysis modules.
- Bilingual README (English + Chinese), cover image, workflow prompt tables.
- CITATION.cff, CONTRIBUTING.md, this changelog.

## [3.1.0] — 2026-10

### Added
- Haplotype EM analysis (`07_haplotype.R`), GWAS×QTL colocalization screen
  (`02_colocalization.R`), GBLUP genomic selection with CV (`03_gs_gblup.R`),
  AMMI/GGE biplots (`04_ammi_gge.R`), XP-EHH selection scans (`07_xpehh.R`),
  KA/KS evolution (`01_ka_ks.py`), promoter k-mer enrichment (`01_promoter_kmer.py`).
- Small-RNA / splicing playbook; molecular lab protocols; PlantCV workflow note.

## [3.0.0] — 2026-09

### Added
- Full-omics expansion: variant calling, pangenome, QTL (eQTL/mQTL/pQTL),
  proteomics, epigenomics (WGBS/ATAC/ChIP), phenomics (BLUP/spectral/GS/AMMI),
  single-cell (Scanpy), spatial (Squidpy + RSTA/scGPT integration guide),
  rice database resource map.

## [2.1.0] — 2026-08

### Added
- primer & CRISPR target design module; Meta-GWAS; lab protocol companions.

## [2.0.0] — 2026-07

### Added
- Genomics expansion: variant detection chain (BWA→GATK→SnpEff), GWAS pipeline
  (PLINK QC → rMVP/GAPIT/GLM → Manhattan/QQ → loci), population genetics
  (PCA/structure, LD decay, Fst/π, ADMIXTURE, phylogeny).

## [1.0.0] — 2026-06

### Added
- Initial release: rice transcriptome + metabolome + integration pipelines
  (QC → DE → enrichment → WGCNA; OPLS-DA → pathway → biomarker; MOFA/DIABLO),
  dual R/Python backends, shared Nature-style figure specification, simulated
  demo data and end-to-end self-tests.

## [Unreleased]

### Added
- **qPCR relative quantification** (`scripts/design/03_qpcr_analysis.py`) —
  2^-ddCt from long-format Ct tables (technical-replicate averaging, calibrator
  group = 1), per-gene bar plots with SD, all-gene grid, and optional
  qPCR↔RNA-seq Pearson/Spearman cross-check with scatter.
