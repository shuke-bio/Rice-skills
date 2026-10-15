"""01_scanpy_pipeline.py -- single-cell RNA-seq pipeline (Scanpy).

Title:    Rice-skills single-cell pipeline (step 01)
Description:
    End-to-end single-cell RNA-seq analysis for rice (and other plants):
    QC filtering -> normalisation -> HVG -> PCA -> neighbours -> UMAP/tSNE
    -> Leiden clustering -> marker genes (Wilcoxon) -> optional reference
    marker scoring (e.g. RSTA cell-type markers) -> optional PAGA trajectory
    -> Nature-spec figures.

    Input can be:
      - a Cell Ranger output directory (barcodes/features/matrix.mtx.gz), or
      - a CSV/h5ad expression matrix (cells x genes or genes x cells; auto).

    The script is a template runner: if Scanpy is not installed it prints
    install instructions and exits. Compute-heavy steps run comfortably on a
    laptop for <=50k cells.

Usage:
    python scripts/single_cell/01_scanpy_pipeline.py \
        --input demo_data/single_cell/counts.csv \
        --outdir results/single_cell/01_pipeline \
        --min-genes 200 --max-mito 10 --resolution 0.8

    # with cell-type marker scoring (CSV: cell_type, gene columns -> use GMT-like)
    python scripts/single_cell/01_scanpy_pipeline.py --input ... \
        --markers marker_sets.csv --outdir ...

Inputs:
    --input      Cell Ranger dir, .h5ad, or expression matrix CSV
                 (first column = cell/barcode id; remaining columns = genes).
    --metadata   optional CSV (cell_id, sample/group) for batch colouring
    --markers    optional CSV: first column cell_type, remaining columns =
                 marker genes (blank cells allowed) — used for score_genes
    --min-genes  minimum genes per cell (default 200)
    --max-mito   maximum mitochondrial fraction % (default 10)
    --min-cells  minimum cells per gene (default 3)
    --n-hvg      highly variable genes (default 2000)
    --n-pcs      PCA components (default 30)
    --resolution Leiden resolution (default 0.8)
    --batch-key  metadata column for harmony batch integration (optional)
    --outdir     output directory (default results/single_cell/01_pipeline)
Outputs:
    qc_violin.* / qc_scatter.*      QC distributions before filtering
    umap_clusters.* / umap_qc.*     UMAP coloured by cluster / QC metrics
    marker_genes.csv                top markers per cluster (Wilcoxon)
    cluster_counts.csv              cells per cluster
    cell_scores.csv                 (optional) marker-set scores per cell
    paga_trajectory.*               (optional) PAGA connectivity map
    processed.h5ad                  processed AnnData object
    run_info.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

import numpy as np
import pandas as pd

from data_io import set_seed, setup_logging, write_run_info
from plot_theme import save_nature_figure, setup_nature_style


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Scanpy single-cell pipeline (rice).")
    p.add_argument("--input", required=True, help="Cell Ranger dir / .h5ad / counts CSV.")
    p.add_argument("--metadata", default=None)
    p.add_argument("--markers", default=None)
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--max-mito", type=float, default=10.0)
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--resolution", type=float, default=0.8)
    p.add_argument("--batch-key", default=None)
    p.add_argument("--outdir", default="results/single_cell/01_pipeline")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    set_seed(args.seed)
    setup_logging()

    try:
        import scanpy as sc
    except ImportError:
        print("Scanpy is not installed.")
        print("")
        print("Install (recommended, in a dedicated conda env):")
        print("  conda create -n sc python=3.11")
        print("  conda activate sc")
        print("  pip install scanpy  leidenalg  # plus harmony-pytorch for batch")
        print("")
        print("Optional downstream tools:")
        print("  pip install scvelo  cellrank   # RNA velocity / fate")
        print("  R: CellChat / NicheNet         # cell-cell communication")
        print("")
        print("The script will run once Scanpy is available. See")
        print("docs/single_cell_usage.md for the full workflow and rice resources.")
        return 0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_nature_style()
    sc.settings.verbosity = 1
    sc.settings.figdir = str(outdir)

    # ------------------------------------------------------------------
    # 1) Load
    # ------------------------------------------------------------------
    inp = str(args.input)
    if Path(inp).is_dir():
        adata = sc.read_10x_mtx(inp, var_names="gene_symbols", cache=True)
        adata.var_names_make_unique()
    elif inp.endswith(".h5ad"):
        adata = sc.read_h5ad(inp)
    else:
        df = pd.read_csv(inp, index_col=0)
        # orientation: samples x genes is standard; detect by gene-like names
        if df.shape[0] < df.shape[1]:
            df = df.T
        adata = sc.AnnData(df.astype(np.float32))
        adata.var_names_make_unique()
    print(f"[single_cell] loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    if args.metadata:
        meta = pd.read_csv(args.metadata)
        meta = meta.set_index(meta.columns[0])
        shared = [c for c in adata.obs_names if c in meta.index]
        if len(shared) < 10:
            print("WARNING: fewer than 10 cells matched the metadata; skipping annotation.")
        else:
            adata = adata[shared].copy()
            for col in meta.columns:
                adata.obs[col] = meta.loc[adata.obs_names, col].astype(str).values

    # ------------------------------------------------------------------
    # 2) QC
    # ------------------------------------------------------------------
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT_", "OSMT"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    fig_ok = False
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(8, 2.6))
        for ax, col, lab in zip(
                axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                ["genes / cell", "UMI / cell", "mito %"]):
            ax.violinplot(adata.obs[col].values, showextrema=False)
            ax.set_ylabel(lab)
            ax.set_xticks([])
        save_nature_figure(fig, outdir / "qc_violin", width_mm=180, height_mm=60)
        plt.close(fig)
        fig_ok = True
    except Exception as exc:  # noqa: BLE001
        print(f"(QC figure skipped: {exc})")

    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < args.max_mito].copy()
    print(f"[single_cell] QC: {n0} -> {adata.n_obs} cells "
          f"({adata.n_vars} genes; min_genes={args.min_genes}, "
          f"max_mito={args.max_mito}%)")

    # ------------------------------------------------------------------
    # 3) Normalise + HVG + PCA + neighbours
    # ------------------------------------------------------------------
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat")
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1))

    if args.batch_key and args.batch_key in adata.obs.columns:
        try:
            import scanpy.external as sce
            sce.pp.harmony_integrate(adata, args.batch_key)
            rep = "X_pca_harmony"
            print(f"[single_cell] harmony integration on '{args.batch_key}'")
        except ImportError:
            rep = "X_pca"
            print("(harmony not installed - using uncorrected PCA; pip install harmony-pytorch)")
    else:
        rep = "X_pca"

    sc.pp.neighbors(adata, use_rep=rep)
    sc.tl.umap(adata)

    # ------------------------------------------------------------------
    # 4) Clustering + markers
    # ------------------------------------------------------------------
    sc.tl.leiden(adata, resolution=args.resolution, key_added="cluster",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon")
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)
    marker_df.to_csv(outdir / "marker_genes.csv", index=False)
    counts = adata.obs["cluster"].value_counts().sort_index()
    counts.to_csv(outdir / "cluster_counts.csv", header=["n_cells"])
    print(f"[single_cell] clusters: {len(counts)} @ resolution {args.resolution}")

    # UMAP figures
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color="cluster", ax=ax, show=False, legend_loc="on data",
                   size=8, frameon=False, title="Clusters (Leiden)")
        save_nature_figure(fig, outdir / "umap_clusters", width_mm=89, height_mm=80)
        plt.close(fig)
        qc_cols = ["n_genes_by_counts", "pct_counts_mt"]
        sc.pl.umap(adata, color=qc_cols, show=False, size=8, frameon=False,
                   save="_qc.png")
        for c in qc_cols:
            src = outdir / f"umap{c}.png"
            if src.exists():
                src.rename(outdir / f"umap_qc_{c}.png")
    except Exception as exc:  # noqa: BLE001
        print(f"(UMAP figures skipped: {exc})")

    # ------------------------------------------------------------------
    # 5) Optional: marker-set scoring (e.g. RSTA cell-type markers)
    # ------------------------------------------------------------------
    if args.markers:
        mk = pd.read_csv(args.markers)
        sets = {}
        for _, row in mk.iterrows():
            genes = [g for g in row.iloc[1:].dropna().astype(str).tolist() if g]
            if genes:
                sets[str(row.iloc[0])] = genes
        scores = pd.DataFrame(index=adata.obs_names)
        for name, genes in sets.items():
            present = [g for g in genes if g in adata.raw.var_names]
            if len(present) < 3:
                print(f"  (marker set '{name}': only {len(present)} genes found; skipped)")
                continue
            sc.tl.score_genes(adata, present, score_name=f"score_{name}")
            scores[f"score_{name}"] = adata.obs[f"score_{name}"].values
        if scores.shape[1]:
            scores.to_csv(outdir / "cell_scores.csv")
            best = scores.idxmax(axis=1).str.replace("score_", "", regex=False)
            pd.DataFrame({"cell_id": scores.index, "best_marker_set": best.values}) \
                .to_csv(outdir / "cell_type_assignment.csv", index=False)
            print(f"[single_cell] marker sets scored: {scores.shape[1]}")
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                sc.pl.umap(adata, color=list(scores.columns[:6]), ax=ax, show=False,
                           ncols=3, size=6, frameon=False)
                save_nature_figure(fig, outdir / "umap_marker_scores",
                                   width_mm=180, height_mm=120)
                plt.close(fig)
            except Exception as exc:  # noqa: BLE001
                print(f"(marker score figure skipped: {exc})")

    # ------------------------------------------------------------------
    # 6) Optional: PAGA trajectory
    # ------------------------------------------------------------------
    try:
        sc.tl.paga(adata, groups="cluster")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        sc.pl.paga(adata, plot=True, ax=ax, show=False, frameon=False)
        sc.tl.umap(adata, init_pos="paga")
        save_nature_figure(fig, outdir / "paga_trajectory", width_mm=89, height_mm=80)
        plt.close(fig)
    except Exception as exc:  # noqa: BLE001
        print(f"(PAGA skipped: {exc})")

    adata.write_h5ad(outdir / "processed.h5ad")

    write_run_info(
        outdir,
        script="01_scanpy_pipeline.py",
        params={
            "input": inp, "min_genes": args.min_genes, "max_mito": args.max_mito,
            "n_hvg": args.n_hvg, "n_pcs": args.n_pcs, "resolution": args.resolution,
            "batch_key": args.batch_key or "none", "seed": args.seed,
            "n_cells_final": int(adata.n_obs), "n_genes_final": int(adata.n_vars),
            "n_clusters": int(len(counts)),
        },
        packages=["scanpy", "anndata"],
    )
    print(f"[single_cell] done -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
