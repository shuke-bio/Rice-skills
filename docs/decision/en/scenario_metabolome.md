# Scenario: Metabolomics (stand-alone)

## 1. Scope and prerequisites

- Input: targeted or untargeted metabolite abundance matrix + sample sheet.
- Requirements: report missingness and imputation; untargeted data are semi-quantitative.
- Goals: differential metabolites, pathways, biomarkers, networks, or integration.

## 2. Data type first

| Type | Allowed | Not allowed |
|---|---|---|
| Targeted (absolute) | Group comparisons, pathways, biomarkers | - |
| Untargeted (relative) | Relative changes, pathways, networks | Absolute concentrations, cross-batch comparison |
| Compositional | Correlation/discrimination after CLR | Raw proportions in regression |

## 3. Three tiers

### Minimal
1. `metabolomics/01_preprocess.py`: missing filter/imputation, normalization, QC.
2. PCA and sample correlation for outliers.
3. Report detected metabolites and missingness.

### Standard
1. `metabolomics/02_oplsda.R`: OPLS-DA, VIP, permutation test.
2. `metabolomics/03_pathway.R`: pathway enrichment.
3. `metabolomics/05_biomarker.R`: ROC/random forest with independent validation.

### Advanced
1. `metabolomics/06_metabolite_network.R`: co-expression network.
2. `metabolomics/04_metabolite_classification.R`: annotation classes.
3. `qtl/02_colocalization.R`: mQTL/colocalization when genotypes exist.

## 4. Sources (whitelist only)

- Rice polished-grain metabolome + hyperspectral + mGWAS (Genome Biology 2025, doi:10.1186/s13059-025-03513-w).
- Plant metabolomics methods review (Plant Physiology 2025).
- Rice population metabolome and brown planthopper resistance (JIPB 2025).

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| Perfect OPLS-DA separation | Overfitting/small n | Permutation test and cross-validation |
| Empty pathway results | Unmapped metabolite IDs | Unify KEGG/ChEBI IDs |
| Biomarkers fail on new data | No independent validation | Hold-out or external cohort |
| Many strong correlations | Compositional artefacts | CLR/log-ratio transform |

## 6. Script map

`scripts/metabolomics/01-07`, `scripts/qtl/02,03`.
