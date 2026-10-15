"""01_qc_filter.py -- Transcriptome QC and low-expression gene filtering.

Title:    Transcriptome quality control and low-expression gene filtering
Description:
    Loads the raw count matrix and sample metadata, aligns the samples,
    filters lowly expressed genes (a gene is retained when its count reaches
    --min-counts in at least --min-samples samples), computes CPM and
    log2(CPM+1), and writes per-sample QC metrics plus a standard QC figure
    set: library size, detected genes, sample PCA on the top highly variable
    genes, and the sample-wise Spearman correlation heatmap. The scientific
    interpretation of every figure is documented in the comment directly
    above the corresponding plotting code.

Usage:
    python scripts/transcriptomics/01_qc_filter.py \
        --counts counts.csv --metadata metadata.csv \
        --outdir results/transcriptome/01_qc \
        --min-counts 10 --min-samples 3 --top-n-hv 500

Inputs:
    counts.csv     gene x sample integer count matrix (first column gene_id)
    metadata.csv   sample metadata (first column sample_id, required column
                   `group`, optional column `batch`)

Outputs (written to --outdir):
    counts_filtered.csv                        filtered integer counts
    qc_metrics.csv                             per-sample QC metrics
    qc_report.txt                              filtering summary + outlier warnings
    library_size.{pdf,svg,tiff}                sequencing depth per sample
    detected_genes.{pdf,svg,tiff}              detected genes per sample
    sample_pca.{pdf,svg,tiff}                  PC1/PC2 coloured by group
    sample_correlation_heatmap.{pdf,svg,tiff}  Spearman correlation heatmap
    run_info.txt                               provenance (parameters, versions, seed)

Notes:
    - Counts must be non-negative integers; missing values are not allowed.
    - Samples whose median pairwise Spearman correlation is below 0.7 are
      reported as potential outliers in qc_report.txt.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

from data_io import (  # noqa: E402
    RiceDataError,
    align_samples,
    check_group_sizes,
    load_expression_matrix,
    load_metadata,
    save_table,
    set_seed,
    setup_logging,
    write_run_info,
)
from plot_theme import (  # noqa: E402
    diverging,
    group_colors,
    save_nature_figure,
    setup_nature_style,
)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402

# Marker shapes used to encode the batch column (only when present).
BATCH_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
# A sample is flagged as a potential outlier below this median correlation.
OUTLIER_COR_CUTOFF = 0.7


def filter_low_expression(
    counts: pd.DataFrame, min_counts: int, min_samples: int, logger
) -> pd.DataFrame:
    """Keep genes with count >= min_counts in at least min_samples samples."""
    n_expressed = (counts >= min_counts).sum(axis=1)
    keep = n_expressed >= min_samples
    n_before, n_after = counts.shape[0], int(keep.sum())
    logger.info(
        "Low-expression filter: %d -> %d genes (count >= %d in >= %d samples)",
        n_before, n_after, min_counts, min_samples,
    )
    if n_after == 0:
        raise RiceDataError(
            f"No genes passed the low-expression filter (count >= {min_counts} "
            f"in >= {min_samples} samples).\n"
            "Hint: lower --min-counts / --min-samples or check the count matrix."
        )
    return counts.loc[keep]


def compute_log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Return log2(CPM + 1); library sizes are computed on the input counts."""
    library = counts.sum(axis=0)
    if (library == 0).any():
        empty = library[library == 0].index.tolist()
        raise RiceDataError(
            f"Sample(s) with zero library size: {empty}\n"
            "Hint: remove empty samples from the count matrix and metadata."
        )
    cpm = counts.div(library, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def sample_qc_metrics(
    counts: pd.DataFrame, log_cpm: pd.DataFrame, meta: pd.DataFrame, group_col: str
) -> pd.DataFrame:
    """Per-sample QC metrics: library size, detected genes, median Spearman r."""
    # log_cpm columns are samples, so corr() yields the sample x sample matrix.
    corr = log_cpm.corr(method="spearman")
    n = corr.shape[0]
    values = corr.to_numpy()
    off_diagonal = values[~np.eye(n, dtype=bool)].reshape(n, -1)
    median_corr = np.median(off_diagonal, axis=1)
    metrics = pd.DataFrame(
        {
            "sample_id": list(counts.columns),
            group_col: meta.loc[counts.columns, group_col].to_numpy(),
            "library_size": counts.sum(axis=0).to_numpy(dtype=int),
            "detected_genes": (counts > 0).sum(axis=0).to_numpy(dtype=int),
            "median_spearman": median_corr,
        }
    )
    if "batch" in meta.columns:
        metrics["batch"] = meta.loc[counts.columns, "batch"].to_numpy()
    metrics["outlier_flag"] = metrics["median_spearman"] < OUTLIER_COR_CUTOFF
    return metrics


def pca_scores(log_cpm: pd.DataFrame, top_n_hv: int, n_components: int = 10):
    """PCA on the top-N highly variable genes (by variance of log2 CPM)."""
    n_hv = int(min(top_n_hv, log_cpm.shape[0]))
    hv_genes = log_cpm.var(axis=0).sort_values(ascending=False).index[:n_hv]
    X = log_cpm[hv_genes].to_numpy()
    n_comp = int(min(n_components, X.shape[0] - 1, X.shape[1]))
    if n_comp < 2:
        raise RiceDataError(
            "PCA needs at least 3 samples and 2 genes; got "
            f"{X.shape[0]} samples and {X.shape[1]} genes.\n"
            "Hint: provide at least three samples per experimental group."
        )
    pca = PCA(n_components=n_comp)
    scores = pca.fit_transform(X)
    return pca, scores, hv_genes


def plot_library_size(metrics: pd.DataFrame, group_col: str, outdir: Path) -> None:
    """Figure: library size (total counts) per sample.

    Scientific conclusion: comparable library sizes across groups indicate
    balanced sequencing depth, so counts can be compared directly after the
    low-expression filter; a sample far below the group median (e.g. <50%)
    signals a failed library that should be re-sequenced or excluded.
    """
    setup_nature_style()
    groups = list(dict.fromkeys(metrics[group_col]))
    colors = dict(zip(groups, group_colors(len(groups))))
    fig, ax = plt.subplots()
    x = np.arange(len(metrics))
    ax.bar(x, metrics["library_size"], color=[colors[g] for g in metrics[group_col]])
    ax.set_xticks(x)
    ax.set_xticklabels(metrics["sample_id"], rotation=45, ha="right")
    ax.set_ylabel("Library size (total counts)")
    ax.set_title("Library size per sample")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors.values()]
    ax.legend(handles, list(colors), title=group_col, frameon=False, loc="best")
    save_nature_figure(fig, outdir / "library_size", width_mm=89, height_mm=60)
    plt.close(fig)


def plot_detected_genes(metrics: pd.DataFrame, group_col: str, outdir: Path) -> None:
    """Figure: number of detected genes (count > 0) per sample.

    Scientific conclusion: samples with markedly fewer detected genes than
    their group indicate poor RNA quality or shallow sequencing; the metric
    also confirms that the low-expression filter removed mostly genes with
    near-zero counts rather than expressed genes.
    """
    setup_nature_style()
    groups = list(dict.fromkeys(metrics[group_col]))
    colors = dict(zip(groups, group_colors(len(groups))))
    fig, ax = plt.subplots()
    x = np.arange(len(metrics))
    ax.bar(x, metrics["detected_genes"], color=[colors[g] for g in metrics[group_col]])
    ax.set_xticks(x)
    ax.set_xticklabels(metrics["sample_id"], rotation=45, ha="right")
    ax.set_ylabel("Detected genes (count > 0)")
    ax.set_title("Detected genes per sample")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors.values()]
    ax.legend(handles, list(colors), title=group_col, frameon=False, loc="best")
    save_nature_figure(fig, outdir / "detected_genes", width_mm=89, height_mm=60)
    plt.close(fig)


def plot_sample_pca(
    scores: np.ndarray,
    var_ratio: np.ndarray,
    metrics: pd.DataFrame,
    group_col: str,
    outdir: Path,
) -> None:
    """Figure: PC1/PC2 scatter, colour = group, shape = batch (if available).

    Scientific conclusion: samples separating by group along PC1/PC2 show that
    the treatment effect dominates the expression variance (good signal);
    clustering driven by batch instead would motivate adding batch to the
    DESeq2 design formula (02_deseq2.R --design '~ batch + group').
    """
    setup_nature_style()
    groups = list(dict.fromkeys(metrics[group_col]))
    colors = dict(zip(groups, group_colors(len(groups))))
    has_batch = "batch" in metrics.columns
    batches = list(dict.fromkeys(metrics["batch"])) if has_batch else []
    fig, ax = plt.subplots()
    seen_labels = set()
    for i in range(len(metrics)):
        group = metrics[group_col].iloc[i]
        if has_batch:
            batch = batches.index(metrics["batch"].iloc[i])
            marker = BATCH_MARKERS[batch % len(BATCH_MARKERS)]
            label = f"{group} / batch {metrics['batch'].iloc[i]}"
        else:
            marker = "o"
            label = group
        ax.scatter(
            scores[i, 0], scores[i, 1], c=[colors[group]], marker=marker,
            s=22, alpha=0.9, label=label if label not in seen_labels else None,
        )
        seen_labels.add(label)
    ax.set_xlabel(f"PC1 ({var_ratio[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var_ratio[1]:.1f}% variance)")
    ax.set_title("Sample PCA (log2 CPM, top HV genes)")
    ax.legend(frameon=False, loc="best", title="group" if not has_batch else "group / batch")
    save_nature_figure(fig, outdir / "sample_pca", width_mm=89, height_mm=70)
    plt.close(fig)


def plot_correlation_heatmap(log_cpm: pd.DataFrame, outdir: Path) -> None:
    """Figure: sample-wise Spearman correlation heatmap on log2 CPM.

    Scientific conclusion: high within-group correlations (r > 0.8) and lower
    between-group correlations confirm that biological replicates are
    reproducible and that the group effect is visible before modelling;
    a sample with low correlation to all others is an outlier candidate
    (also reported in qc_report.txt).
    """
    setup_nature_style()
    # log_cpm columns are samples, so corr() yields the sample x sample matrix.
    corr = log_cpm.corr(method="spearman")
    labels = list(corr.columns)
    n = len(labels)
    cmap = LinearSegmentedColormap.from_list("rice_diverging", diverging())
    fig, ax = plt.subplots()
    im = ax.imshow(corr.to_numpy(), cmap=cmap, vmin=-1.0, vmax=1.0)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    if n <= 12:  # annotate only readable matrix sizes
        for i in range(n):
            for j in range(n):
                ax.text(
                    j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center",
                    fontsize=max(4.0, 7.0 - 0.4 * n),
                )
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman correlation")
    ax.set_title("Sample correlation (log2 CPM)")
    save_nature_figure(
        fig, outdir / "sample_correlation_heatmap", width_mm=89, height_mm=80
    )
    plt.close(fig)


def write_qc_report(
    path: Path,
    args: argparse.Namespace,
    counts_raw: pd.DataFrame,
    counts_filtered: pd.DataFrame,
    metrics: pd.DataFrame,
    group_col: str,
) -> None:
    """Write qc_report.txt: filtering summary, library stats, outlier warnings."""
    outliers = metrics.loc[metrics["outlier_flag"], "sample_id"].tolist()
    sizes = metrics["library_size"]
    detected = metrics["detected_genes"]
    lines = [
        "Rice-skills -- transcriptome QC report",
        "Generated by scripts/transcriptomics/01_qc_filter.py (seed 42)",
        "",
        "Inputs:",
        f"  counts:   {args.counts}",
        f"  metadata: {args.metadata}",
        "",
        f"Samples: {len(metrics)} (groups: "
        + ", ".join(
            f"{g}={int((metrics[group_col] == g).sum())}"
            for g in dict.fromkeys(metrics[group_col])
        )
        + ")",
        f"Genes before filtering: {counts_raw.shape[0]}",
        f"Genes after filtering:  {counts_filtered.shape[0]} "
        f"(kept genes with count >= {args.min_counts} in >= {args.min_samples} samples)",
        "",
        f"Library size:    min={sizes.min()}, median={int(sizes.median())}, max={sizes.max()}",
        f"Detected genes:  min={detected.min()}, median={int(detected.median())}, "
        f"max={detected.max()}",
        "",
    ]
    if outliers:
        lines += [
            "WARNING: potential outlier sample(s) with median Spearman r < "
            f"{OUTLIER_COR_CUTOFF}: {', '.join(outliers)}",
            "  -> Inspect these samples (RNA quality, batch, contamination) before",
            "     differential expression analysis; consider excluding them.",
            "",
        ]
    else:
        lines += [
            f"Outlier check: no sample with median Spearman r < {OUTLIER_COR_CUTOFF}.",
            "",
        ]
    lines += [
        "Outputs:",
        "  counts_filtered.csv, qc_metrics.csv",
        "  library_size.*, detected_genes.*, sample_pca.*, sample_correlation_heatmap.*",
        "  (each figure saved as PDF + SVG + TIFF at 300 dpi)",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Transcriptome QC and low-expression gene filtering.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--counts", required=True,
                        help="Raw count matrix CSV (first column gene_id, rest = samples)")
    parser.add_argument("--metadata", required=True,
                        help="Sample metadata CSV (first column sample_id, group column)")
    parser.add_argument("--outdir", default="results/transcriptome/01_qc",
                        help="Output directory")
    parser.add_argument("--min-counts", type=int, default=10,
                        help="Count threshold for the low-expression filter")
    parser.add_argument("--min-samples", type=int, default=3,
                        help="Minimum number of samples reaching --min-counts")
    parser.add_argument("--group-col", default="group",
                        help="Metadata column defining the experimental groups")
    parser.add_argument("--top-n-hv", type=int, default=500,
                        help="Number of highly variable genes used for the PCA")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(log_file=outdir / "01_qc_filter.log")
    set_seed(42)

    logger.info("Loading count matrix: %s", args.counts)
    counts = load_expression_matrix(args.counts)
    logger.info("Loading metadata: %s", args.metadata)
    meta = load_metadata(args.metadata, required_cols=(args.group_col,))
    counts, meta = align_samples(counts, meta, matrix_name="count matrix")
    logger.info("Aligned %d genes x %d samples", counts.shape[0], counts.shape[1])
    check_group_sizes(meta, group_col=args.group_col, min_size=2)

    counts_filtered = filter_low_expression(
        counts, args.min_counts, args.min_samples, logger
    )
    log_cpm = compute_log_cpm(counts_filtered)
    metrics = sample_qc_metrics(counts_filtered, log_cpm, meta, args.group_col)

    save_table(counts_filtered, outdir / "counts_filtered.csv", index=True)
    save_table(metrics, outdir / "qc_metrics.csv", index=False)
    logger.info("Saved counts_filtered.csv and qc_metrics.csv")

    pca, scores, _ = pca_scores(log_cpm, args.top_n_hv)
    var_ratio = pca.explained_variance_ratio_ * 100.0
    logger.info(
        "PCA on top %d highly variable genes: PC1 = %.1f%%, PC2 = %.1f%% variance",
        min(args.top_n_hv, log_cpm.shape[0]), var_ratio[0], var_ratio[1],
    )

    plot_library_size(metrics, args.group_col, outdir)
    plot_detected_genes(metrics, args.group_col, outdir)
    plot_sample_pca(scores, var_ratio, metrics, args.group_col, outdir)
    plot_correlation_heatmap(log_cpm, outdir)
    logger.info("QC figures saved (PDF + SVG + TIFF)")

    write_qc_report(
        outdir / "qc_report.txt", args, counts, counts_filtered, metrics, args.group_col
    )
    write_run_info(
        outdir,
        "01_qc_filter.py",
        params={
            "counts": args.counts,
            "metadata": args.metadata,
            "group_col": args.group_col,
            "min_counts": args.min_counts,
            "min_samples": args.min_samples,
            "top_n_hv": args.top_n_hv,
        },
        packages=["numpy", "pandas", "scikit-learn", "matplotlib"],
    )
    logger.info("QC complete. Results in %s", outdir)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except RiceDataError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
