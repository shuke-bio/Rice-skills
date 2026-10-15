# Changelog

All notable changes to Rice-skills are documented here.
Versioning follows [Semantic Versioning](https://semver.org/).

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
