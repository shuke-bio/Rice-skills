# Scenario: Multi-omics integration (transcriptome x metabolome)

## 1. Scope and prerequisites

- Input: >=2 omics layers from the same samples (typically RNA-seq + metabolome).
- Requirements: sample matching by ID; n>=20 before factor models.
- Goals: cross-layer associations, modules, candidate gene-metabolite pairs, mQTL/eQTL colocalization.

## 2. Can these layers be integrated?

| Check | Alternative if it fails |
|---|---|
| Sample IDs match one-to-one | Independent analyses per layer |
| Same collection batch | Model batch, or stratify |
| n >= 20 | Correlation/co-expression only; no MOFA/DIABLO |
| Relative abundances | Apply CLR/log-ratio before correlation |

## 3. Three tiers

### Minimal
1. `metabolomics/01_preprocess.py`: missing-value handling, normalization, QC.
2. `transcriptomics/01_qc_filter.py`: expression QC.
3. Intersect sample IDs; report n before and after.

### Standard
1. `metabolomics/02_oplsda.R` (discrimination, VIP).
2. `integration/02_diablo.R` (supervised; requires labels).
3. `integration/01_mofa.R` (unsupervised factors; n>=20).
4. Gene-metabolite Spearman correlation with BH correction.

### Advanced
1. `integration/05_network_integration.R` (module-metabolite network).
2. `qtl/02_colocalization.R` / `03_coloc_bayesian.R` (mQTL/eQTL colocalization).
3. `gwas/07_haplotype.R` (candidate haplotypes).

## 4. Sources (whitelist only)

- Rice polished-grain metabolome + hyperspectral + mGWAS (Genome Biology 2025, doi:10.1186/s13059-025-03513-w).
- DIABLO: supervised multi-omics discrimination (Bioinformatics 2019, doi:10.1093/bioinformatics/btz158).
- MOFA: standard factor model; non-whitelisted method paper is not cited as a design source.

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| Factors reflect batch only | Batch unmodelled | Model/stratify by batch |
| Strong but implausible correlations | Compositional data | CLR transform; inspect drivers |
| Low variance explained | Layers truly unrelated | Report the negative result |
| High PP4 without functional support | LD/power artefacts | Add expression/haplotype/edit validation |

## 6. Script map

`scripts/metabolomics/01-07`, `scripts/transcriptomics/01-07`,
`scripts/integration/01,02,05`, `scripts/qtl/02,03`.
