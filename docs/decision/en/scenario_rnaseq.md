# Scenario: Bulk RNA-seq

## 1. Scope and prerequisites

- Input: integer count matrix + sample sheet (group, batch, covariates).
- Requirements: biological replicates (n>=2, n>=3 preferred); matched genome/annotation version.
- Goals: differential expression, markers, pathways/regulation, integration with metabolome.

## 2. Design-to-model mapping

| Design | Example formula | Note |
|---|---|---|
| Treatment vs control | `~ batch + group` | Batch must be modelled, or technical effects become "biology" |
| Time course | `~ group + time + group:time` | Interpret the interaction, not separate per-time tests |
| Multi-environment | `~ environment + group + environment:group` | Environment as fixed effect when few levels |

## 3. Three tiers

### Minimal
1. `01_qc_filter.py`: low-count filter, sample correlation, PCA.
2. `02_deseq2.R`: differential table only when replicates exist; descriptive summary otherwise.
3. `07_plot_transcriptome.R`: PCA/correlation/volcano drafts.

### Standard
1. `01_qc_filter.py` -> `02_deseq2.R` (batch in the design).
2. `03_enrichment.R` (GO/KEGG) -> `04_gsea.R` (ranked gene sets).
3. `06_tf_analysis.R` -> `07_plot_transcriptome.R`.
4. Sensitivity: intersect edgeR and DESeq2 results.

### Advanced
1. `05_wgcna.R` (n>=15) for module-trait associations.
2. `integration/01_mofa.R` or `02_diablo.R` when matched metabolome exists.
3. `integration/05_network_integration.R` for gene-metabolite networks.
4. `qtl/02_colocalization.R` for eQTL/mQTL colocalization when genotypes exist.

## 4. Sources (whitelist only)

- DESeq2: negative binomial with shrinkage (Genome Biology 2014, doi:10.1186/s13059-014-0550-8).
- edgeR: exact/GLM framework (Bioinformatics 2010, doi:10.1093/bioinformatics/btp616).
- limma-voom: large samples or non-count platforms (Genome Biology 2014, doi:10.1186/gb-2014-15-2-r29).
- GSEA: ranked gene-set enrichment (PNAS 2005, doi:10.1073/pnas.0506580102).
- WGCNA: standard tool; non-whitelisted method paper is not cited as a design source.

## 5. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| Thousands of DE genes, no enrichment | Batch/quality dominates | Check PCA and batch crosstab; refit with batch |
| Perfect separation, no DE genes | Small n or high variance | Report effect sizes; avoid strong claims |
| Candidates disagree with literature | Annotation mismatch | Unify MSU/RAP-DB and convert IDs |
| Enrichment changes with background | Wrong background set | Use detected genes as background |

## 6. Script map

`scripts/transcriptomics/01-07`, `scripts/integration/01,02,05`, `scripts/qtl/02`.
