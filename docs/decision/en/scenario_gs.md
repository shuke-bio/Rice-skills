# Scenario: Phenomics and genomic selection (GS)

## 1. Scope and prerequisites

- Input: phenotypes (multi-environment/years, preferably replicated) + genotypes.
- Requirements: strict sample ID matching; phenotypes corrected for block/environment (BLUP).
- Goals: heritability, prediction accuracy, candidate ranking.

## 2. Three checks first

1. Intersect phenotype and genotype samples; record n before/after.
2. Phenotype distribution and outliers (`phenomics/01_phenotype_blup.R`).
3. Power assessment (`phenomics/05_power_analysis.R`).

## 3. Three tiers

### Minimal
1. `phenomics/01_phenotype_blup.R`: BLUP and heritability.
2. `phenomics/05_power_analysis.R`: power.
3. Sample matching table.

### Standard
1. `phenomics/03_gs_gblup.R`: VanRaden kinship + ridge regression + k-fold CV.
2. Multi-environment: `phenomics/04_ammi_gge.R` for stability.
3. Report mean +/- SD of repeated CV and predicted-observed correlation.

### Advanced
1. `phenomics/06_menet_gs.py`: MeNet multi-environment prediction (PyTorch).
2. Joint multi-environment/year modelling.
3. Compare against the GBLUP baseline (mandatory).

## 4. Sources (whitelist only)

- MeNet: mixed-effects deep network validated on 12 rice traits x 3 environments
  (Plant Communications 2025, doi:10.1016/j.xplc.2025.101620).
- GBLUP implementation: `phenomics/03_gs_gblup.R` (standard tool, no non-whitelisted citation).
- AMMI/GGE: `phenomics/04_ammi_gge.R`.

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| CV accuracy > 0.9 | Data leakage/random split | Group splits by family/environment |
| Fails in new environments | Strong G x E | Joint multi-environment model |
| Heritability near zero | Noise/small n | More replicates/environments |
| Deep model worse than GBLUP | Small n/overfitting | Keep GBLUP as baseline; shrink the model |

## 6. Script map

`scripts/phenomics/01,03,04,05,06`, `scripts/gwas/01`.
