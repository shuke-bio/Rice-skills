# Scenario: Epigenomics (WGBS / ATAC / ChIP)

## 1. Scope and prerequisites

- WGBS: differential methylation (DMR); report coverage and bisulfite conversion.
- ATAC/ChIP: differential accessibility/binding; report FRiP, TSS enrichment, replicate concordance.
- Requirements: >=2 biological replicates; matched material and stage; batch records.

## 2. Three tiers

### Minimal
1. `epigenomics/01_methylation_analysis.R` (WGBS) or `02_atac_chip_analysis.R` (ATAC/ChIP).
2. Sample correlation and coverage checks.
3. Collect batch/technical covariates.

### Standard
1. DMR/DAR detection with explicit thresholds and smoothing.
2. Region annotation (promoter/enhancer/gene body).
3. Joint analysis with RNA-seq when available.

### Advanced
1. 3D genome/chromatin interactions (Hi-C when available).
2. Colocalization with GWAS/QTL.
3. Multi-omics integration (MOFA/DIABLO).

## 3. Sources (whitelist only)

- Single-cell chromatin accessibility atlas in rice (Nature Plants 2025).
- 3D epigenome architecture and heading date (Genome Biology 2026).
- Integrated 3D genome/epigenome/transcriptome (NAR 2023).
- H3K9ac-marked 3D interactions under drought (Genome Biology 2024).

## 4. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| Insufficient coverage | Low sequencing depth | Sequence more or restrict to high-coverage regions |
| Low replicate correlation | Batch/material differences | Check collection and library batches |
| DMRs only in repeats | Ambiguous mapping | Filter repeats and report the rule |
| No TSS enrichment in ATAC | Library failure | Rebuild the library; do not force analysis |

## 5. Script map

`scripts/epigenomics/01_methylation_analysis.R`, `02_atac_chip_analysis.R`.
