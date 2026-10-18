# Changelog

All notable changes to Rice-skills are documented here.
Versioning follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **Read-based SV detection route (module 11 · pangenome)** —
  `scripts/pangenome/02_sv_detect.sh --mode read-based` wraps **HitSV**
  (Li G et al., *Nat Biotechnol* 2026, doi:10.1038/s41587-026-03332-9): pure
  short-read (`--seq srs`), pure long-read (`--seq lrs`, `--platform
  HIFI|ONT_Q20|ERR_PRONE|ASM`) and hybrid (`--seq hybrid`) calling + genotyping,
  outputs a genotyped VCF that feeds `03_pav_analysis.R` / `04_pan_gwas.R`.
  Ships a `--demo` synthetic end-to-end self-test (300 kb reference, simulated
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
  dual R/Python backends, shared Nature-style figure specification, synthetic
  demo data and end-to-end self-tests.

## [Unreleased]

### Added
- **qPCR relative quantification** (`scripts/design/03_qpcr_analysis.py`) —
  2^-ddCt from long-format Ct tables (technical-replicate averaging, calibrator
  group = 1), per-gene bar plots with SD, all-gene grid, and optional
  qPCR↔RNA-seq Pearson/Spearman cross-check with scatter.
