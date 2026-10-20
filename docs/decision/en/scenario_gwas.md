# Scenario: GWAS (natural population, families, multi-parent)

## 1. Scope and prerequisites

- Input: genotypes (SNP/InDel/SV/PAV) + phenotypes (replicates/multi-environment) + population structure.
- Requirements: population structure must be corrected; reliable phenotypes (blocked design).
- Goals: associated loci, candidate intervals, allele effects.

## 2. Design first, model second

| Design | Suggested n | Primary method | Note |
|---|---|---|---|
| Natural population | >=200-300 | GLM/MLM (PCA+K), FarmCPU | Low power for rare variants |
| RIL/F2 | 150-300 families | Linkage mapping | Complements association |
| MAGIC/NAM/ROAM | >=300 families | Multi-parent QTL (two-stage) | Requires family structure |
| Case-control | >=200/group | Logistic + structure correction | Case/control ratio affects power |

## 3. Three tiers

### Minimal
1. `01_plink_qc.sh`: MAF/missingness/HWE filters.
2. `02_plink_to_matrix.R`: export the genotype matrix.
3. `05_plot_gwas.R`: QQ/Manhattan templates (validate the plotting chain).

### Standard
1. `03_association.R` (GLM/MLM with PCA covariates).
2. `04_candidate_loci.R` (lead SNP +/- LD interval).
3. `07_haplotype.R` (haplotype-phenotype).
4. `08_2d_association.R` (pairwise/epistasis screening).

### Advanced
1. `09_multilocus_qtn.R` (3VmrMLM: multi-locus QTN / QEI / epistasis).
2. `10_case_control_qtn.R` (binary traits; Linux/WSL).
3. `11_multiparent_qtl.R` (MAGIC/NAM/ROAM).
4. `06_meta_gwas.R` (multi-environment/multi-cohort).

## 4. Sources (whitelist only)

- 3VmrMLM: multi-locus QTN, QTN x environment and epistasis (Molecular Plant 2022, doi:10.1016/j.molp.2022.02.012).
- Fast3VmrMLM: genome-wide multi-locus scan (Plant Communications 2025, doi:10.1016/j.xplc.2025.101385).
- Graph pangenome GWAS: SV/PAV-aware association (Nature Communications 2025).
- TWAS: expression-informed candidate prioritisation (Nature Communications 2023).
- mGWAS + colocalization in rice (Genome Biology 2025, doi:10.1186/s13059-025-03513-w).

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| Inflated QQ (lambda >> 1) | Uncorrected structure | Add PCA/kinship; check family structure |
| Signals cluster at centromeres | Low recombination/repeats | Flag or exclude; not candidate loci |
| Many candidates in a wide LD block | Large LD interval | Add haplotype/expression/colocalization evidence |
| Unstable binary results | Imbalanced groups | Adjust ratio; report OR with CI |
| Contradictory environments | Real G x E | Use QEI models and report environment-specific effects |

## 6. Script map

`scripts/gwas/01-11`, `scripts/qtl/01-03`, `scripts/popgen/`.
