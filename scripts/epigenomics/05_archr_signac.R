#!/usr/bin/env Rscript
###############################################################################
# 05_archr_signac.R -- ArchR / Signac plan generator
#
# 说明：ArchR/Signac 依赖较重，默认只生成可执行 R 计划。
#
# 用法：
#   Rscript scripts/epigenomics/05_archr_signac.R --fragments fragments.tsv.gz --outdir results/epigenomics/05_archr
#   Rscript scripts/epigenomics/05_archr_signac.R --demo --outdir results/epigenomics/05_archr_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/epigenomics/05_archr"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
frag <- opt("fragments", "fragments.tsv.gz")
archr <- c(
  "library(ArchR)",
  sprintf("addArchRGenome('hg19')  # 水稻请用自定义 BSgenome/基因注释"),
  sprintf("inputFiles <- '%s'", frag),
  "ArrowFiles <- createArrowFiles(inputFiles=inputFiles, sampleNames='S1', minTSS=4, minFrags=1000)",
  "proj <- ArchRProject(ArrowFiles=ArrowFiles, outputDirectory='ArchRProj')",
  "proj <- addDoubletScores(proj); proj <- filterDoublets(proj)",
  "proj <- addIterativeLSI(proj); proj <- addClusters(proj, resolution=0.8)",
  "proj <- addUMAP(proj); saveArchRProject(proj)"
)
writeLines(archr, file.path(outdir, "archr_plan.R"))
signac <- c(
  "library(Signac); library(Seurat)",
  sprintf("counts <- Read10X('%s')", "filtered_feature_bc_matrix"),
  "obj <- CreateSeuratObject(counts=counts)",
  "obj <- RunTFIDF(obj); obj <- FindTopFeatures(obj, min.cutoff='q5')",
  "obj <- RunSVD(obj); obj <- RunUMAP(obj, reduction='lsi', dims=2:30)",
  "obj <- FindNeighbors(obj, reduction='lsi', dims=2:30); obj <- FindClusters(obj, resolution=0.8)"
)
writeLines(signac, file.path(outdir, "signac_plan.R"))
writeLines(c(paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
             "mode=plan only (ArchR/Signac not executed)"), file.path(outdir, "run_info.txt"))
cat("[archr] plans ->", outdir, "\n")
