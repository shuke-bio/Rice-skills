"""01_spatial_pipeline.py -- spatial transcriptomics pipeline (Stereo-seq / Visium).

Title:    Rice-skills spatial transcriptomics pipeline (step 01)
Description:
    End-to-end spatial analysis for bin/spot-based platforms:
      load (h5ad or matrix + coords CSV) -> QC -> normalise -> HVG -> PCA ->
      expression clustering (spatial domains) -> spatial neighbour graph ->
      spatially variable genes (Moran's I, squidpy) -> spatial figures
      (domain map, gene expression maps) -> optional neighbourhood enrichment.
    Scanpy is required; Squidpy enables the spatial-graph steps (Moran's I,
    neighbourhood enrichment). Without Squidpy the script still produces the
    expression-domain map and gene expression plots from spatial coordinates.

Usage:
    # Stereo-seq bin matrix + coordinates
    python scripts/spatial/01_spatial_pipeline.py \
        --input demo_data/spatial/counts.csv --coords demo_data/spatial/coords.csv \
        --outdir results/spatial/01_pipeline --resolution 0.5

    # Visium / any h5ad carrying obsm['spatial']
    python scripts/spatial/01_spatial_pipeline.py --input sample.h5ad --outdir ...

    # highlight specific genes on the tissue map
    python scripts/spatial/01_spatial_pipeline.py ... --genes G1,G2,G3

Inputs:
    --input      count matrix CSV (bin/cell x gene) or .h5ad
    --coords     coordinates CSV: barcode, x, y (required for CSV input)
    --genes      comma-separated genes drawn on the spatial map
    --min-genes  min genes per bin/spot (default 100)
    --resolution Leiden resolution for domains (default 0.5)
    --n-hvg      highly variable genes (default 2000)
    --outdir     output directory (default results/spatial/01_pipeline)
Outputs:
    domain_map.*            spatial scatter coloured by expression domain
    gene_maps.*             gene expression spatial maps (with --genes)
    svg_results.csv         Moran's I per gene (squidpy; top-N)
    spatial_domains.csv     barcode -> domain assignment
    neighbourhood_enrichment.*  (squidpy) domain adjacency enrichment
    processed.h5ad / run_info.txt
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
    p = argparse.ArgumentParser(description="Spatial transcriptomics pipeline (Stereo-seq / Visium).")
    p.add_argument("--input", required=True, help="Count matrix CSV or .h5ad.")
    p.add_argument("--coords", default=None, help="Coordinates CSV (barcode, x, y) for CSV input.")
    p.add_argument("--genes", default=None, help="Comma-separated genes to map spatially.")
    p.add_argument("--min-genes", type=int, default=100)
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--resolution", type=float, default=0.5)
    p.add_argument("--n-svg", type=int, default=2000, help="Genes tested for Moran's I (top HVG).")
    p.add_argument("--outdir", default="results/spatial/01_pipeline")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    set_seed(args.seed)
    setup_logging()

    try:
        import scanpy as sc
        import anndata as ad
    except ImportError:
        print("Scanpy is not installed (required).")
        print("")
        print("Install:")
        print("  conda create -n spatial python=3.11 -y && conda activate spatial")
        print("  pip install scanpy squidpy   # squidpy enables Moran's I & neighbourhood tests")
        print("")
        print("For Stereo-seq data conversion (gef -> matrix) see the SAW pipeline:")
        print("  https://github.com/STOmics/SAW   (official; fetched by the user at runtime)")
        print("See docs/spatial_usage.md for the full workflow and RSTA integration.")
        return 0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_nature_style()
    sc.settings.verbosity = 1

    # ------------------------------------------------------------------
    # 1) Load + attach coordinates
    # ------------------------------------------------------------------
    inp = str(args.input)
    if inp.endswith(".h5ad"):
        adata = sc.read_h5ad(inp)
        if "spatial" not in adata.obsm:
            print("WARNING: h5ad has no obsm['spatial']; spatial plots will be skipped.")
    else:
        if not args.coords:
            print("ERROR: --coords is required for CSV input (barcode, x, y).")
            return 1
        df = pd.read_csv(inp, index_col=0)
        if df.shape[0] < df.shape[1]:
            df = df.T
        adata = ad.AnnData(df.astype(np.float32))
        adata.var_names_make_unique()
        coords = pd.read_csv(args.coords)
        coords.columns = [c.lower() for c in coords.columns]
        coords = coords.set_index(coords.columns[0])
        shared = [c for c in adata.obs_names if c in coords.index]
        if len(shared) < 10:
            print("ERROR: fewer than 10 barcodes shared between matrix and coords.")
            return 1
        adata = adata[shared].copy()
        xy = coords.loc[adata.obs_names, ["x", "y"]].to_numpy(dtype=float)
        adata.obsm["spatial"] = xy
    print(f"[spatial] loaded: {adata.n_obs} bins/spots x {adata.n_vars} genes")

    # ------------------------------------------------------------------
    # 2) QC + expression domains
    # ------------------------------------------------------------------
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat")
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="domain",
                 flavor="igraph", n_iterations=2, directed=False)
    pd.DataFrame({"barcode": adata.obs_names, "domain": adata.obs["domain"].values}) \
        .to_csv(outdir / "spatial_domains.csv", index=False)
    n_dom = adata.obs["domain"].nunique()
    print(f"[spatial] expression domains: {n_dom} @ resolution {args.resolution}")

    xy = adata.obsm.get("spatial", None)

    # domain map
    try:
        import matplotlib.pyplot as plt
        if xy is not None:
            fig, ax = plt.subplots()
            doms = sorted(adata.obs["domain"].unique(), key=int)
            palette = sc.pl.palettes.default_20
            for i, d in enumerate(doms):
                m = (adata.obs["domain"] == d).values
                ax.scatter(xy[m, 0], xy[m, 1], s=3, label=f"domain {d}",
                           color=palette[i % len(palette)], linewidths=0)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_title("Spatial domains (expression clusters)")
            ax.legend(markerscale=3, fontsize=4, frameon=False, ncol=2)
            save_nature_figure(fig, outdir / "domain_map", width_mm=120, height_mm=100)
            plt.close(fig)
    except Exception as exc:  # noqa: BLE001
        print(f"(domain map skipped: {exc})")

    # gene maps
    if args.genes and xy is not None:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        try:
            import matplotlib.pyplot as plt
            n = len(genes)
            fig, axes = plt.subplots(1, n, figsize=(2.4 * n, 2.6), squeeze=False)
            for j, g in enumerate(genes):
                if g not in adata.raw.var_names:
                    axes[0][j].set_title(f"{g}\n(not found)", fontsize=6)
                    continue
                v = adata.raw[:, g].X
                v = np.asarray(v.todense()).ravel() if hasattr(v, "todense") else np.asarray(v).ravel()
                sca = axes[0][j].scatter(xy[:, 0], xy[:, 1], c=v, s=3, cmap="viridis", linewidths=0)
                axes[0][j].set_title(g, fontsize=6)
                axes[0][j].set_aspect("equal")
                axes[0][j].invert_yaxis()
                axes[0][j].set_xticks([])
                axes[0][j].set_yticks([])
                fig.colorbar(sca, ax=axes[0][j], shrink=0.8)
            save_nature_figure(fig, outdir / "gene_maps", width_mm=min(180, 60 * n), height_mm=70)
            plt.close(fig)
        except Exception as exc:  # noqa: BLE001
            print(f"(gene maps skipped: {exc})")

    # ------------------------------------------------------------------
    # 3) Spatial statistics (squidpy)
    # ------------------------------------------------------------------
    try:
        import squidpy as sq
        sq.gr.spatial_neighbors(adata)
        print("[spatial] squidpy: spatial neighbour graph built")

        # Moran's I on top HVGs
        genes_test = adata.var_names[: min(args.n_svg, adata.n_vars)].tolist()
        sq.gr.spatial_autocorr(adata, mode="moran", genes=genes_test,
                               n_perms=100, n_jobs=1)
        svg = adata.uns["moranI"].copy()
        svg.sort_values("I", ascending=False, inplace=True)
        svg.to_csv(outdir / "svg_results.csv")
        print(f"[spatial] SVG (Moran's I) computed for {len(svg)} genes; "
              f"top: {svg.index[0]} (I={svg.iloc[0]['I']:.3f})")

        # neighbourhood enrichment between domains
        try:
            sq.gr.nhood_enrichment(adata, cluster_key="domain")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            sq.pl.nhood_enrichment(adata, cluster_key="domain", ax=ax, show=False)
            save_nature_figure(fig, outdir / "neighbourhood_enrichment",
                               width_mm=89, height_mm=80)
            plt.close(fig)
        except Exception as exc:  # noqa: BLE001
            print(f"(neighbourhood enrichment skipped: {exc})")
    except ImportError:
        print("(squidpy not installed - Moran's I / neighbourhood tests skipped;")
        print(" pip install squidpy to enable them)")
    except Exception as exc:  # noqa: BLE001
        print(f"(spatial statistics skipped: {exc})")

    adata.write_h5ad(outdir / "processed.h5ad")

    write_run_info(
        outdir,
        script="01_spatial_pipeline.py",
        params={
            "input": inp, "resolution": args.resolution,
            "min_genes": args.min_genes, "n_hvg": args.n_hvg,
            "n_bins_final": int(adata.n_obs), "n_domains": int(n_dom),
            "seed": args.seed,
        },
        packages=["scanpy", "anndata", "squidpy"],
    )
    print(f"[spatial] done -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
