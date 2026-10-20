# Scenario: QTL and BSA

## 1. Scope and prerequisites

- Family designs: RIL/F2/MAGIC/NAM/ROAM; genetic map or dense markers required.
- BSA-seq: extreme bulks from two parents; sufficient individuals and depth.
- Requirements: reliable phenotypes (blocks/replicates); accurate pedigree.

## 2. Design selection

| Design | Suggested scale | Method |
|---|---|---|
| F2/RIL | >=150 individuals/families | Linkage analysis (LOD, permutation) |
| MAGIC/NAM/ROAM | >=300 families | Multi-parent QTL (two-stage) |
| BSA-seq | >=30 individuals/bulk, >=30x | Delta SNP-index / G prime |
| Natural population | >=200 | Association (see GWAS card) |

## 3. Three tiers

### Minimal
1. `qtl/01_qtl_mapping.R`: genetic map/linkage analysis.
2. Marker QC via `gwas/01_plink_qc.sh` when needed.
3. Check phenotype distribution and blocks.

### Standard
1. `qtl/01_qtl_mapping.R` with permutation-based LOD thresholds.
2. `qtl/04_qtl_region_snps.py`: interval SNP extraction and annotation.
3. `qtl/05_bsa_scan.py`: Delta SNP-index / G prime with simulated confidence lines.

### Advanced
1. `qtl/02_colocalization.R` / `03_coloc_bayesian.R` with GWAS/expression.
2. `gwas/07_haplotype.R`: candidate haplotypes.
3. `design/01_sgrna_design.py`: editing targets for validation.

## 4. Sources (whitelist only)

- Multi-locus QTN models (Molecular Plant 2022, doi:10.1016/j.molp.2022.02.012).
- Multi-parent population QTL (Plant Communications 2025).
- Rice QTL fine mapping and haplotypes (JIPB 2024).
- Pangenome SV association (Nature Communications 2025).

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| High LOD everywhere | Phenotype/pedigree error | Check pedigree, blocks and distribution |
| Very wide intervals | Low marker density | Add markers or switch to BSA |
| No difference between bulks | Insufficient pooling/wrong phenotype | Increase individuals/depth; re-phenotype |
| Too many candidates | Large interval | Haplotype + expression + editing validation |

## 6. Script map

`scripts/qtl/01-05`, `scripts/gwas/01,07`, `scripts/design/01`.
