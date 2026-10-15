"""01_preprocess.py -- Proteome preprocessing: missing values, normalisation, QC.

Title:    Rice-skills proteome preprocessing (step 01)
Description:
    Loads a protein abundance matrix (proteins x samples, label-free LFQ /
    iBAQ / TMT ratios) and sample metadata, then:
      1. aligns samples and reports missing-value structure;
      2. filters proteins detected in >= --min-valid of samples
         (proteomics convention: 70% is a common default);
      3. applies log2 transformation (auto/always/never) and per-sample
         median normalisation, with optional per-batch median centring;
      4. imputes the remaining missing values (half-min / KNN / none);
      5. produces QC figures: missing-value barplot, PCA (colour = group),
         sample-sample Spearman correlation heatmap, within-group CV.
    All figures follow the shared Nature figure specification.

Usage:
    python scripts/proteomics/01_preprocess.py \
        --abundance demo_data/proteomics/proteins.csv \
        --metadata demo_data/proteomics/metadata.csv \
        --outdir results/proteomics/01_preprocess

Inputs:
    --abundance   Protein abundance matrix CSV (first column = protein_id,
                  remaining columns = samples; numeric).
    --metadata    Sample metadata CSV (sample_id, group [, batch, ...]).

Outputs (written to --outdir):
    proteins_filtered.csv      filtered + normalised (+ imputed) matrix
    missing_summary.csv        per-protein / per-sample missing statistics
    pca_scores.csv             PCA scores (normalised data) + group
    missing_barplot.*          per-sample missing fraction (PDF/SVG/TIFF)
    pca_samples.*              PCA scatter coloured by group
    sample_correlation.*       Spearman correlation heatmap
    cv_distribution.*          within-group CV per protein
    run_info.txt               parameters, versions, seed

Notes:
    - Imputation "halfmin" fills each protein with half its observed minimum
      (deterministic); "knn" uses sklearn KNNImputer (k=5); "none" keeps NA
      (downstream scripts must then handle missingness).
    - Normalisation chain: log2 (optional) -> per-sample median centring ->
      optional per-batch median centring. No variance scaling (that is a
      downstream/clustering choice).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer

from data_io import (
    align_samples,
    load_expression_matrix,
    load_metadata,
    save_table,
    set_seed,
    setup_logging,
    write_run_info,
)
from plot_theme import group_colors, save_nature_figure, setup_nature_style


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Proteome preprocessing (missing values, normalisation, QC).")
    p.add_argument("--abundance", required=True, help="Protein abundance matrix CSV (protein_id + samples).")
    p.add_argument("--metadata", required=True, help="Sample metadata CSV (sample_id, group, ...).")
    p.add_argument("--outdir", default="results/proteomics/01_preprocess")
    p.add_argument("--min-valid", type=float, default=0.7,
                   help="Keep proteins detected in >= this fraction of samples (default 0.7).")
    p.add_argument("--log2", choices=["auto", "always", "never"], default="auto",
                   help="log2 transform: auto (if max > 50), always, never.")
    p.add_argument("--impute", choices=["halfmin", "knn", "none"], default="halfmin")
    p.add_argument("--batch-col", default=None, help="Optional metadata column for per-batch median centring.")
    p.add_argument("--group-col", default="group", help="Metadata column with sample groups (default 'group').")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    set_seed(args.seed)
    setup_logging()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_nature_style()

    mat = load_expression_matrix(args.abundance, id_col="protein_id",
                                 what="Protein abundance matrix")  # proteins x samples
    meta = load_metadata(args.metadata)
    mat, meta = align_samples(mat, meta, matrix_name="protein matrix")
    if args.group_col not in meta.columns:
        raise SystemExit(f"metadata is missing the group column '{args.group_col}'")
    samples = list(mat.columns)
    print(f"[01_preprocess] proteins={mat.shape[0]} x samples={mat.shape[1]}")

    # ---------------- missing values ----------------
    per_protein_valid = mat.notna().sum(axis=1)
    per_sample_missing = mat.isna().mean(axis=0)
    summary = pd.DataFrame({
        "protein_id": mat.index,
        "n_valid": per_protein_valid.values,
        "valid_frac": (per_protein_valid / len(samples)).values,
    })
    keep = summary["valid_frac"] >= args.min_valid
    print(f"[01_preprocess] filter >= {args.min_valid:.0%} valid: keep {keep.sum()} / {len(keep)} proteins")
    mat = mat.loc[keep.values].copy()
    if mat.shape[0] < 10:
        raise SystemExit("Fewer than 10 proteins remain after filtering; relax --min-valid.")

    # ---------------- transform + normalise ----------------
    x = mat.to_numpy(dtype=float)
    do_log = args.log2 == "always" or (args.log2 == "auto" and np.nanmax(x) > 50)
    if do_log:
        print("[01_preprocess] log2 transform applied")
        x = np.log2(np.where(x > 0, x, np.nan))
    # per-sample median centring
    col_med = np.nanmedian(x, axis=0)
    x = x - col_med
    # optional per-batch median centring
    if args.batch_col:
        if args.batch_col not in meta.columns:
            raise SystemExit(f"metadata is missing batch column '{args.batch_col}'")
        batch = meta[args.batch_col].astype(str).values
        for b in np.unique(batch):
            cols = np.where(batch == b)[0]
            if len(cols) == 0:
                continue
            bm = np.nanmedian(x[:, cols], axis=0, keepdims=True)
            grand = np.nanmedian(x[:, cols])
            x[:, cols] = x[:, cols] - bm + grand
        print(f"[01_preprocess] per-batch median centring on '{args.batch_col}'")
    mat_norm = pd.DataFrame(x, index=mat.index, columns=mat.columns)

    # ---------------- impute ----------------
    if args.impute == "halfmin":
        for i in range(mat_norm.shape[0]):
            row = mat_norm.iloc[i]
            if row.isna().any():
                # half-min: in log2 space, halving = -1; raw scale = /2
                fill = (np.nanmin(row.values) - 1.0) if do_log else (np.nanmin(row.values) / 2.0)
                mat_norm.iloc[i] = row.fillna(fill)
    elif args.impute == "knn":
        imp = KNNImputer(n_neighbors=5)
        arr = mat_norm.to_numpy()
        if np.isnan(arr).any():
            mat_norm.iloc[:, :] = imp.fit_transform(arr)
    print(f"[01_preprocess] imputation: {args.impute}")

    save_table(mat_norm, outdir / "proteins_filtered.csv", index=True, index_label="protein_id")
    save_table(summary, outdir / "missing_summary.csv")

    # ---------------- figures ----------------
    groups = meta[args.group_col].astype(str)
    order = list(dict.fromkeys(groups.values))
    colors = {g: c for g, c in zip(order, group_colors(len(order)))}

    # (a) per-sample missing barplot
    fig, ax = plt.subplots()
    ax.bar(range(len(samples)), 100 * per_sample_missing.values,
           color=[colors[g] for g in groups.loc[samples].values])
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=4)
    ax.set_ylabel("Missing values (%)")
    ax.set_title("Missing value fraction per sample")
    save_nature_figure(fig, outdir / "missing_barplot", width_mm=120, height_mm=70)
    plt.close(fig)

    # (b) PCA on normalised matrix
    X = mat_norm.T.to_numpy()  # samples x proteins
    pcs = PCA(n_components=min(10, X.shape[0], X.shape[1]))
    scores = pcs.fit_transform(X)
    ve = pcs.explained_variance_ratio_ * 100
    pca_df = pd.DataFrame({"sample_id": samples, "PC1": scores[:, 0], "PC2": scores[:, 1] if scores.shape[1] > 1 else 0,
                           "group": groups.loc[samples].values})
    pca_df.to_csv(outdir / "pca_scores.csv", index=False)
    fig, ax = plt.subplots()
    for g in order:
        sub = pca_df[pca_df["group"] == g]
        ax.scatter(sub["PC1"], sub["PC2"], s=12, alpha=0.85, label=g, color=colors[g])
    ax.set_xlabel(f"PC1 ({ve[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({ve[1]:.1f}%)" if len(ve) > 1 else "PC2")
    ax.legend(frameon=False, fontsize=5)
    ax.set_title("PCA of normalised proteome")
    save_nature_figure(fig, outdir / "pca_samples", width_mm=89, height_mm=75)
    plt.close(fig)

    # (c) sample-sample Spearman correlation heatmap
    corr = mat_norm.corr(method="spearman")
    fig, ax = plt.subplots()
    im = ax.imshow(corr.values, cmap="RdBu_r", vmin=np.nanmin(corr.values), vmax=1.0)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=3.5)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=3.5)
    ax.set_title("Sample-sample Spearman correlation")
    fig.colorbar(im, ax=ax, shrink=0.7, label="rho")
    save_nature_figure(fig, outdir / "sample_correlation", width_mm=100, height_mm=95)
    plt.close(fig)

    # (d) within-group CV distribution
    fig, ax = plt.subplots()
    cvs = {}
    for g in order:
        cols = [c for c in samples if groups.loc[c] == g]
        if len(cols) >= 3:
            sub = mat_norm[cols]
            mu = sub.mean(axis=1).abs()
            sd = sub.std(axis=1)
            cv = (sd / mu.replace(0, np.nan) * 100).dropna()
            cvs[g] = cv.values
    if cvs:
        try:
            ax.boxplot(list(cvs.values()), tick_labels=list(cvs.keys()), showfliers=False)
        except TypeError:  # matplotlib < 3.9
            ax.boxplot(list(cvs.values()), labels=list(cvs.keys()), showfliers=False)
        ax.set_ylabel("CV (%)")
        ax.set_title("Within-group coefficient of variation")
    save_nature_figure(fig, outdir / "cv_distribution", width_mm=80, height_mm=65)
    plt.close(fig)

    # ---------------- run info ----------------
    write_run_info(
        outdir,
        script="01_preprocess.py",
        params={
            "abundance": str(args.abundance), "metadata": str(args.metadata),
            "min_valid": args.min_valid, "log2": args.log2, "do_log2": bool(do_log),
            "impute": args.impute, "batch_col": args.batch_col,
            "n_proteins_final": int(mat_norm.shape[0]), "n_samples": len(samples),
            "seed": args.seed,
        },
    )
    print(f"[01_preprocess] done -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
