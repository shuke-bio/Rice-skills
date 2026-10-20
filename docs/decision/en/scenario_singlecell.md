# Scenario: Single-cell and spatial omics

## 1. Scope and prerequisites

- Single-cell: scRNA-seq/snRNA-seq for cell types, composition, trajectories, regulation.
- Spatial: Stereo-seq/Visium for spatial domains, region-specific expression, neighbourhoods.
- Requirements: biological samples are the unit of replication; sample-level QC before merging; batch/section records.

## 2. Can between-condition comparisons be made?

| Check | Alternative if it fails |
|---|---|
| >=2 biological samples per group | Descriptive atlas only; no between-group statistics |
| Interpretable sample-level QC (mito%, genes) | Re-derive thresholds per tissue and report |
| Complete batch/section metadata | Conservative integration; label as exploratory |
| Registered coordinates and images | Fix registration first; otherwise spatial claims are invalid |

## 3. Three tiers

### Minimal
1. `single_cell/01_scanpy_pipeline.py`: QC, normalization, HVG, PCA.
2. Report cells and median genes per sample.
3. Integration (Harmony/scVI via script options).

### Standard
1. Clustering + marker annotation (optionally label transfer from a reference atlas).
2. Cell-type composition differences at the sample level (not cell level).
3. Pseudobulk differential expression to avoid pseudo-replication.

### Advanced
1. Trajectory inference with an explicit root.
2. SCENIC-like regulatory networks (resource intensive).
3. Spatial mapping or multi-section integration.

## 4. Sources (whitelist only)

- Rice single-cell multi-omics atlas (Nature 2025): reference atlas and label transfer.
- Rice spatiotemporal lifecycle atlas (Cell 2026): spatial and developmental design.
- Single-cell chromatin accessibility atlas (Nature Plants 2025): scATAC and cross-species comparison.
- Rice root scRNA-seq + spatial under soil stress (Nature 2025).

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| Differences only in proportions | Dissociation/capture bias | Report capture efficiency; confirm by pseudobulk |
| Signal disappears after integration | Over-integration | Reduce integration strength |
| Opposite trajectory direction | Subjective root | Support with RNA velocity/markers |
| Spatial domains disagree with histology | Registration error | Re-align images and coordinates |

## 6. Script map

`scripts/single_cell/01_scanpy_pipeline.py`, `scripts/spatial/01_spatial_pipeline.py`.
