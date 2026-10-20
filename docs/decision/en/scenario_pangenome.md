# Scenario: Pangenome and structural variation

## 1. Scope and prerequisites

- Input: multiple chromosome-level assemblies (>=5, ideally >=10) with unified annotation.
- Requirements: consistent annotation and gene-ID mapping; otherwise PAV matrices are not comparable.
- Goals: core/variable/private genes, SV/PAV, Pan-GWAS, genus-level comparison.

## 2. Can a pangenome be built?

| Check | Alternative if it fails |
|---|---|
| Assembly continuity (N50/BUSCO) | Improve assemblies or use SV calling only |
| Consistent annotation source | Re-annotate or restrict to SV |
| Defined reference genome | Fix IRGSP-1.0 or the relevant assembly |
| Representative sampling | Add diverse accessions |

## 3. Three tiers

### Minimal
1. `pangenome/01_pangenome_graph.sh` (minigraph/cactus/PGGB).
2. Unify annotation and gene-ID mapping.
3. Assembly quality assessment (BUSCO/continuity).

### Standard
1. `pangenome/02_sv_detect.sh` (minimap2 + SYRI).
2. `pangenome/03_pav_analysis.R` (core/variable/private).
3. `pangenome/04_pan_gwas.R` (Pan-GWAS).

### Advanced
1. `pangenome/05_oryza_pangenome_pipeline.sh` (genus-level Snakemake stage runner).
2. `pangenome/06_annosine.sh` (SINE annotation).
3. Joint analysis with expression/phenotypes.

## 4. Sources (whitelist only)

- Wild-cultivated rice pangenome with 145 chromosome-level assemblies (Nature 2025).
- Graph pangenome GWAS (Nature Communications 2025).
- TE insertion polymorphisms (Nature Communications 2025).
- PanOryza pangene catalog (Genome Research 2026).
- Genus-scale pangenome workflow (Genome Biology 2023; MIT upstream).

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| PAV matrix mostly missing | Inconsistent annotation | Re-annotate and remap IDs |
| Excess SVs | Assembly errors/chimeras | Filter low-support SVs; cross-validate callers |
| Too few core genes | Fragmented assemblies | Improve assemblies or justify relaxed thresholds |
| No Pan-GWAS signal | Small population/sparse PAV | Expand population; use SV genotypes |

## 6. Script map

`scripts/pangenome/01-06`, `scripts/gwas/03,09`.
