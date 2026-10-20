# Scenario: Proteomics

## 1. Scope and prerequisites

- Input: DIA/DDA quantification matrix (LFQ/iBAQ/TMT) + sample sheet.
- Requirements: report missingness and imputation (MNAR is common); proteins with >50% missing are not tested.
- Goals: differential proteins, enrichment, protein-gene mapping, PTMs, integration with RNA-seq.

## 2. Three tiers

### Minimal
1. `proteomics/01_preprocess.py`: missing filter/imputation, normalization.
2. Sample correlation/PCA.
3. Report missingness and imputation strategy.

### Standard
1. `proteomics/02_diff_protein.R` (limma).
2. `proteomics/03_plot_proteome.R` (volcano/heatmap).
3. `proteomics/04_uniprot_mapping.py` (protein-to-gene).
4. `proteomics/05_interpro_domains.py` (domain enrichment).

### Advanced
1. PTM analysis (phosphoproteomics when available).
2. Integration with transcriptome (protein-mRNA concordance).
3. Pathway/network integration.

## 3. Sources (whitelist only)

- Rice proteome landscape and m6A-mediated post-transcriptional regulation (Nature Plants 2024).
- Mass-spectrometry proteome and m6A (Nature Plants 2024).
- Missing-value handling: see `docs/proteomics_usage.md` (no non-whitelisted citation).

## 4. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| Few differential proteins | High missingness/normalization | Inspect missingness mechanism; filter by group |
| Protein and mRNA disagree | Real post-transcriptional regulation | Report as a finding, not an error |
| Batch dominates PCA | Library/run batch | Model batch or stratify |
| Unreliable PTM sites | Low enrichment/localization | Raise thresholds; report localization probability |

## 5. Script map

`scripts/proteomics/01-05`.
