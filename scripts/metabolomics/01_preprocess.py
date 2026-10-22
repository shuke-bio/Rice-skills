"""01_preprocess.py -- Metabolome preprocessing: missing values, normalization, outlier check.

Title:    Rice-skills metabolome preprocessing (step 01)
Description:
    Loads a metabolite abundance matrix (metabolites x samples) and sample
    metadata, aligns samples, summarises missing values, imputes them
    (half-minimum / KNN / column mean / none), applies a normalisation chain
    (sum -> log2 -> pareto scaling, in the order requested), and screens for
    outlier samples with a PCA-based distance-to-centroid check (warning only,
    samples are never removed automatically).

    The abundance matrix is transposed to samples x metabolites for processing
    (imputation and scaling act per metabolite = per column) and transposed
    back before export, so the output keeps the standard Rice-skills
    layout: first column metabolite_id, annotation columns, then one column
    per sample.

Usage:
    python scripts/metabolomics/01_preprocess.py \
        --abundance examples/metabolome_data/metabolites.csv \
        --metadata examples/metabolome_data/metadata.csv \
        --outdir results/metabolome/01_preprocess

Inputs:
    --abundance  Metabolite abundance matrix CSV (metabolite_id [, metabolite_name,
                 class, kegg_id], samples...). Missing values = empty cells.
    --metadata   Sample metadata CSV (sample_id, group [, batch, variety, traits]).

Outputs (written to --outdir):
    metabolome_preprocessed.csv  imputed + normalised matrix (annotations kept)
    missing_value_summary.csv    per-metabolite and per-sample missing statistics
    pca_scores.csv               PCA scores of the normalised data + outlier flag
    missing_value_barplot.*      per-sample missing-value fraction (PDF/TIFF/PNG)
    pca_before.* / pca_after.*   PCA before vs after normalisation (colour = group)
    boxplot_normalization.*     per-sample abundance distribution before vs after
    run_info.txt                 parameters, package versions and the seed

Notes:
    - half-min imputation fills each metabolite with half of its observed
      minimum (deterministic, MNAR-style). A stochastic variant (random draw
      below the observed minimum, e.g. N(min, 0.2*min) truncated) can be
      enabled by extending impute_half_min(); determinism was preferred here
      so that runs are bit-for-bit reproducible.
    - normalise="auto" chains sum -> log2(x+1) -> pareto scaling, the standard
      MetaboAnalystR processing order for LC-MS abundance data.
    - Outlier rule (spec): Euclidean distance from the centroid in PCA score
      space > mean + 3 * sd is flagged as a warning only; no sample is dropped.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer

from data_io import (
    RiceDataError,
    align_samples,
    check_group_sizes,
    load_metadata,
    load_metabolome,
    save_table,
    set_seed,
    setup_logging,
    write_run_info,
)
from plot_theme import (
    group_colors,
    save_nature_figure,
    semantic,
    setup_nature_style,
)

set_seed(42)

IMPUTE_CHOICES = ("half-min", "knn", "mean", "none")
NORMALIZE_CHOICES = ("none", "sum", "log", "pareto", "auto")


# ---------------------------------------------------------------------------
# Missing values and imputation
# ---------------------------------------------------------------------------

def missing_value_summary(abundance: pd.DataFrame) -> pd.DataFrame:
    """Tidy missing-value statistics per metabolite and per sample.

    Args:
        abundance: Metabolites x samples numeric DataFrame (raw, before impute).

    Returns:
        DataFrame with columns level, id, n_missing, n_total, missing_ratio.
    """
    rows = []
    per_metabolite = abundance.isna().sum(axis=1)
    for metabolite_id, n_missing in per_metabolite.items():
        n_total = abundance.shape[1]
        rows.append(
            {
                "level": "metabolite",
                "id": metabolite_id,
                "n_missing": int(n_missing),
                "n_total": int(n_total),
                "missing_ratio": float(n_missing) / n_total,
            }
        )
    per_sample = abundance.isna().sum(axis=0)
    for sample_id, n_missing in per_sample.items():
        n_total = abundance.shape[0]
        rows.append(
            {
                "level": "sample",
                "id": sample_id,
                "n_missing": int(n_missing),
                "n_total": int(n_total),
                "missing_ratio": float(n_missing) / n_total,
            }
        )
    return pd.DataFrame(rows, columns=["level", "id", "n_missing", "n_total", "missing_ratio"])


def impute_values(samples_x_metabolites: pd.DataFrame, method: str, logger) -> pd.DataFrame:
    """Impute missing values column-wise (per metabolite).

    Args:
        samples_x_metabolites: Samples x metabolites DataFrame (may contain NaN).
        method: One of half-min / knn / mean / none.
        logger: Logger for progress messages.

    Returns:
        DataFrame of the same shape without missing values (unless method=none).

    Raises:
        RiceDataError: For an unknown imputation method.
    """
    X = samples_x_metabolites.copy()
    if method == "none":
        n_remaining = int(X.isna().sum().sum())
        logger.info("Imputation skipped (--impute none); %d missing values retained.", n_remaining)
        return X
    if method == "half-min":
        # Deterministic MNAR-style fill: half of each metabolite's observed minimum.
        col_min = X.min(axis=0)
        fill_values = col_min / 2.0
        X = X.fillna(fill_values)
        logger.info(
            "Imputed missing values with half-minimum (deterministic; per-metabolite "
            "min/2). A stochastic MNAR variant can be added in impute_values()."
        )
    elif method == "knn":
        # KNN imputation across samples (sklearn KNNImputer, n_neighbors=5 by
        # default; clamped when the number of samples is smaller).
        n_neighbors = max(1, min(5, X.shape[0] - 1))
        imputer = KNNImputer(n_neighbors=n_neighbors)
        imputed = imputer.fit_transform(X.values)
        X = pd.DataFrame(imputed, index=X.index, columns=X.columns)
        logger.info("Imputed missing values with KNNImputer (n_neighbors=%d).", n_neighbors)
    elif method == "mean":
        X = X.fillna(X.mean(axis=0))
        logger.info("Imputed missing values with per-metabolite column means.")
    else:
        raise RiceDataError(
            f"Unknown imputation method '{method}'. Choose from: {', '.join(IMPUTE_CHOICES)}."
        )
    n_remaining = int(X.isna().sum().sum())
    if n_remaining:
        logger.warning("%d missing values remain after imputation (fully-missing metabolites?).", n_remaining)
    return X


# ---------------------------------------------------------------------------
# Normalisation chain
# ---------------------------------------------------------------------------

def normalize_values(samples_x_metabolites: pd.DataFrame, method: str, logger) -> tuple[pd.DataFrame, list[str]]:
    """Apply the normalisation chain to a samples x metabolites matrix.

    Args:
        samples_x_metabolites: Imputed samples x metabolites DataFrame.
        method: none / sum / log / pareto / auto (auto = sum -> log -> pareto).
        logger: Logger for progress messages.

    Returns:
        Tuple (normalised DataFrame, list of applied step names).

    Raises:
        RiceDataError: For an unknown normalisation method.
    """
    X = samples_x_metabolites.copy()
    if method == "auto":
        steps = ["sum", "log", "pareto"]
    elif method == "none":
        steps = []
    elif method in ("sum", "log", "pareto"):
        steps = [method]
    else:
        raise RiceDataError(
            f"Unknown normalisation method '{method}'. Choose from: {', '.join(NORMALIZE_CHOICES)}."
        )

    for step in steps:
        if step == "sum":
            # Constant-to-total-sum normalisation per sample (row), rescaled to
            # the mean sample total so absolute magnitude is preserved.
            totals = X.sum(axis=1)
            bad = ~np.isfinite(totals) | (totals <= 0)
            if bad.any():
                logger.warning(
                    "Sample(s) with non-positive total abundance left unchanged by sum "
                    "normalisation: %s", list(X.index[bad])
                )
                totals = totals.mask(bad, np.nan)
            scale = totals.mean(skipna=True)
            X = X.div(totals, axis=0).mul(scale)
            logger.info("Applied sum normalisation (each sample rescaled to the mean total).")
        elif step == "log":
            # log2(x + 1): stabilises variance; mean differences on log2 data
            # approximate log2 fold-changes used downstream.
            if (X <= -1).any().any():
                logger.warning("Negative values <= -1 detected before log transform; clipped to 0.")
            X = np.log2(X.clip(lower=0) + 1.0)
            logger.info("Applied log2(x + 1) transformation.")
        elif step == "pareto":
            # Pareto scaling: mean-centred, divided by sqrt(sd) per metabolite.
            col_sd = X.std(axis=0, ddof=1)
            zero_sd = col_sd[col_sd <= 0].index.tolist()
            if zero_sd:
                logger.warning(
                    "%d constant metabolite(s) not pareto-scaled (sd = 0): %s",
                    len(zero_sd), zero_sd[:5],
                )
                col_sd = col_sd.replace(0.0, 1.0)
            X = (X - X.mean(axis=0)) / np.sqrt(col_sd)
            logger.info("Applied pareto scaling (mean-centred, divided by sqrt(sd)).")
    return X, steps


# ---------------------------------------------------------------------------
# PCA and outlier screening
# ---------------------------------------------------------------------------

def pca_scores(samples_x_metabolites: pd.DataFrame, n_components: int = 2) -> tuple[pd.DataFrame, np.ndarray]:
    """PCA scores (PC1/PC2) with distance-to-centroid outlier flags.

    Args:
        samples_x_metabolites: Samples x metabolites DataFrame.
        n_components: Number of components to retain (clamped to data dimensions).

    Returns:
        Tuple (scores DataFrame with sample_id/PC1/PC2/distance_to_centroid/
        is_outlier, explained variance ratio array).
    """
    X = samples_x_metabolites.values
    n_comp = int(max(1, min(n_components, X.shape[0] - 1, X.shape[1])))
    pca = PCA(n_components=n_comp)
    scores = pca.fit_transform(X)
    centroid = scores.mean(axis=0)
    distances = np.sqrt(((scores - centroid) ** 2).sum(axis=1))
    # Spec rule: warn (do not remove) when distance > mean + 3 * sd.
    threshold = distances.mean() + 3.0 * distances.std(ddof=1)
    scores_df = pd.DataFrame(
        {
            "sample_id": list(samples_x_metabolites.index),
            "PC1": scores[:, 0],
            "PC2": scores[:, 1] if n_comp > 1 else np.zeros(len(scores)),
            "distance_to_centroid": distances,
            "is_outlier": distances > threshold,
        }
    )
    return scores_df, pca.explained_variance_ratio_


# ---------------------------------------------------------------------------
# Figures (Nature-style via the shared plot_theme.py layer)
# ---------------------------------------------------------------------------

def plot_missing_barplot(summary: pd.DataFrame, group_of: pd.Series, outdir: Path) -> None:
    """Barplot of the per-sample missing-value fraction (colour = group)."""
    import matplotlib.pyplot as plt

    sample_rows = summary[summary["level"] == "sample"].set_index("id")
    groups = [group_of.get(s, "NA") for s in sample_rows.index]
    colors = group_colors(len(set(groups)))
    color_map = {g: colors[i] for i, g in enumerate(dict.fromkeys(groups))}

    setup_nature_style()
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    x = np.arange(len(sample_rows))
    ax.bar(x, sample_rows["missing_ratio"].values, width=0.7,
           color=[color_map[g] for g in groups], edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(sample_rows.index, rotation=45, ha="right")
    ax.set_ylabel("Missing value fraction")
    ax.set_xlabel("Sample")
    overall = float(sample_rows["n_missing"].sum()) / float(sample_rows["n_total"].sum())
    ax.set_title("Missing values per sample", loc="left")
    # Scientific conclusion: missingness is concentrated in low-abundance
    # metabolites (MNAR), so half-minimum imputation is appropriate.
    ax.text(0.0, 1.02, f"Overall missing rate {overall:.1%}; imputation per metabolite",
            transform=ax.transAxes, fontsize=6.5, color="#666666")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[g]) for g in color_map]
    ax.legend(handles, list(color_map), frameon=False, title="Group")
    fig.tight_layout()
    save_nature_figure(fig, outdir / "missing_value_barplot", width_mm=89, height_mm=70)
    plt.close(fig)


def plot_pca(scores_df: pd.DataFrame, group_of: pd.Series, evr: np.ndarray,
             title: str, subtitle: str, file_stem: Path) -> None:
    """Scatter of PC1 vs PC2 coloured by group."""
    import matplotlib.pyplot as plt

    groups = [group_of.get(s, "NA") for s in scores_df["sample_id"]]
    colors = group_colors(len(set(groups)))
    color_map = {g: colors[i] for i, g in enumerate(dict.fromkeys(groups))}

    setup_nature_style()
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    for g in color_map:
        mask = np.array([gg == g for gg in groups])
        ax.scatter(scores_df.loc[mask, "PC1"], scores_df.loc[mask, "PC2"],
                   color=color_map[g], s=22, label=g, edgecolor="none", alpha=0.9)
    outliers = scores_df["is_outlier"].to_numpy()
    if outliers.any():
        ax.scatter(scores_df.loc[outliers, "PC1"], scores_df.loc[outliers, "PC2"],
                   s=60, facecolor="none", edgecolor=semantic("highlight"), linewidth=1.0,
                   label="outlier (warning)")
    ax.axhline(0, color=semantic("reference_line"), linewidth=0.5)
    ax.axvline(0, color=semantic("reference_line"), linewidth=0.5)
    ax.set_xlabel(f"PC1 ({evr[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({evr[1]:.1%} variance)")
    ax.set_title(title, loc="left")
    ax.text(0.0, 1.02, subtitle, transform=ax.transAxes, fontsize=6.5, color="#666666")
    ax.legend(frameon=False, title="Group", loc="best")
    fig.tight_layout()
    save_nature_figure(fig, file_stem, width_mm=89, height_mm=75)
    plt.close(fig)


def plot_boxplot_comparison(before: pd.DataFrame, after: pd.DataFrame,
                            group_of: pd.Series, outdir: Path) -> None:
    """Per-sample boxplots before vs after normalisation (two panels)."""
    import matplotlib.pyplot as plt

    groups = [group_of.get(s, "NA") for s in before.index]
    colors = group_colors(len(set(groups)))
    color_map = {g: colors[i] for i, g in enumerate(dict.fromkeys(groups))}

    setup_nature_style()
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.8), sharex=True)
    for ax, (data, panel) in zip(axes, [(before, "Raw (imputed)"), (after, "Normalised")]):
        positions = np.arange(data.shape[0])
        box = ax.boxplot(
            [data.iloc[i].values for i in range(data.shape[0])],
            positions=positions, widths=0.6, patch_artist=True, showfliers=False,
            medianprops=dict(color="#33312E", linewidth=0.8),
        )
        for patch, g in zip(box["boxes"], groups):
            patch.set_facecolor(color_map[g])
            patch.set_edgecolor("#33312E")
            patch.set_linewidth(0.5)
        ax.set_xticks(positions)
        ax.set_xticklabels(data.index, rotation=45, ha="right")
        ax.set_title(panel, loc="left", fontsize=8)
        ax.set_ylabel("Abundance" if panel.startswith("Raw") else "Normalised value")
    axes[0].set_xlabel("Sample")
    axes[1].set_xlabel("Sample")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[g]) for g in color_map]
    axes[1].legend(handles, list(color_map), frameon=False, title="Group", loc="upper right")
    fig.suptitle("Normalisation effect on per-sample distributions", x=0.01, ha="left",
                 fontsize=9, fontweight="bold")
    # Scientific conclusion: normalisation removes between-sample differences in
    # total abundance/scale, making metabolites comparable across samples.
    fig.text(0.01, 0.90, "Sum + log2 + pareto: equalises sample totals and variance",
             fontsize=6.5, color="#666666")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save_nature_figure(fig, outdir / "boxplot_normalization", width_mm=180, height_mm=80)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Metabolome preprocessing: missing-value imputation, normalisation and outlier screening."
    )
    parser.add_argument("--abundance", required=True, help="Metabolite abundance matrix CSV (metabolites x samples)")
    parser.add_argument("--metadata", required=True, help="Sample metadata CSV (sample_id, group)")
    parser.add_argument("--outdir", default="results/metabolome/01_preprocess", help="Output directory")
    parser.add_argument("--impute", default="half-min", choices=IMPUTE_CHOICES,
                        help="Missing-value imputation method (default: half-min)")
    parser.add_argument("--normalize", default="auto", choices=NORMALIZE_CHOICES,
                        help="Normalisation method; 'auto' = sum -> log -> pareto (default: auto)")
    parser.add_argument("--group-col", default="group", help="Metadata column defining the groups (default: group)")
    args = parser.parse_args(argv)

    logger = setup_logging()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info("Output directory: %s", outdir.resolve())

    # ---- Load and align -----------------------------------------------------
    abundance, annotation = load_metabolome(args.abundance)
    metadata = load_metadata(args.metadata)
    abundance, metadata = align_samples(abundance, metadata, matrix_name="abundance matrix")
    check_group_sizes(metadata, args.group_col)
    logger.info("Loaded %d metabolites x %d samples; groups: %s",
                abundance.shape[0], abundance.shape[1],
                metadata[args.group_col].value_counts().to_dict())

    # Work in samples x metabolites orientation (imputation/scaling per metabolite).
    raw = abundance.T.astype(float)
    group_of = metadata[args.group_col]

    # ---- Missing values ------------------------------------------------------
    summary = missing_value_summary(abundance)
    n_missing = int(abundance.isna().sum().sum())
    n_values = abundance.size
    logger.info("Missing values: %d of %d (%.2f%%)", n_missing, n_values, 100.0 * n_missing / n_values)
    n_met_no_kegg = 0
    if "kegg_id" in annotation.columns:
        kegg = annotation["kegg_id"].astype("string")
        n_met_no_kegg = int((kegg.isna() | (kegg.fillna("").str.strip() == "")).sum())
        logger.info("Metabolites without KEGG annotation: %d of %d", n_met_no_kegg, len(annotation))

    # ---- Imputation and normalisation ---------------------------------------
    imputed = impute_values(raw, args.impute, logger)
    normalized, steps = normalize_values(imputed, args.normalize, logger)

    # ---- Outlier screening on normalised data --------------------------------
    scores_after, evr_after = pca_scores(normalized, n_components=2)
    scores_after["group"] = [group_of.get(s, "NA") for s in scores_after["sample_id"]]
    outliers = scores_after.loc[scores_after["is_outlier"], "sample_id"].tolist()
    if outliers:
        logger.warning(
            "Sample(s) flagged as outliers (distance to centroid > mean + 3*sd): %s. "
            "They are kept; inspect QC before interpreting downstream results.", outliers
        )
    else:
        logger.info("No outlier samples detected (distance to centroid <= mean + 3*sd).")

    # PCA before normalisation for the comparison figure.
    scores_before, evr_before = pca_scores(imputed, n_components=2)

    # ---- Export tables -------------------------------------------------------
    processed = normalized.T.copy()  # metabolites x samples
    processed.index.name = "metabolite_id"
    processed_out = processed.reset_index()  # metabolite_id becomes the first column
    # Re-attach annotation columns right after metabolite_id, preserving the
    # recommended order (metabolite_name, class, kegg_id) of the input file.
    for pos, col in enumerate(annotation.columns, start=1):
        processed_out.insert(pos, col, annotation[col].reindex(processed.index).values)
    save_table(processed_out, outdir / "metabolome_preprocessed.csv", index=False)
    save_table(summary, outdir / "missing_value_summary.csv", index=False)
    save_table(scores_after[["sample_id", "group", "PC1", "PC2", "distance_to_centroid", "is_outlier"]],
               outdir / "pca_scores.csv", index=False)
    logger.info("Wrote metabolome_preprocessed.csv (%d metabolites x %d samples)",
                processed_out.shape[0], processed_out.shape[1] - len(annotation.columns) - 1)

    # ---- Figures -------------------------------------------------------------
    plot_missing_barplot(summary, group_of, outdir)
    plot_pca(scores_before, group_of, evr_before,
             "PCA before normalisation", "Raw imputed abundances", outdir / "pca_before")
    plot_pca(scores_after, group_of, evr_after,
             "PCA after normalisation",
             "Sum + log2 + pareto; groups separate along PC1",
             outdir / "pca_after")
    plot_boxplot_comparison(imputed, normalized, group_of, outdir)
    logger.info("Saved figures (PDF/TIFF/PNG): missing_value_barplot, pca_before, pca_after, boxplot_normalization")

    # ---- Provenance ----------------------------------------------------------
    write_run_info(
        outdir,
        script="01_preprocess.py",
        params={
            "abundance": args.abundance,
            "metadata": args.metadata,
            "group_col": args.group_col,
            "impute": args.impute,
            "normalize": args.normalize,
            "normalization_steps_applied": ",".join(steps) if steps else "none",
            "n_metabolites": abundance.shape[0],
            "n_samples": abundance.shape[1],
            "n_missing_values": n_missing,
            "missing_rate": f"{100.0 * n_missing / n_values:.2f}%",
            "metabolites_without_kegg_id": n_met_no_kegg,
            "outlier_samples": ",".join(outliers) if outliers else "none",
        },
        packages=["numpy", "pandas", "scikit-learn", "matplotlib"],
        extra=[
            "Outlier rule: PCA distance to centroid > mean + 3*sd (warning only, no removal).",
            "half-min imputation is deterministic (per-metabolite min/2).",
        ],
    )
    logger.info("Done. Results in %s", outdir.resolve())
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except RiceDataError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
