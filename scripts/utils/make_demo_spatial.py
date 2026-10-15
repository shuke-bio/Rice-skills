"""make_demo_spatial.py -- synthetic spatial transcriptomics demo generator.

Title:    Rice-skills spatial demo data generator
Description:
    Creates a small deterministic synthetic Stereo-seq-like dataset:
      2500 bins on a 50x50 grid, 4 spatial domains (quadrant layout) with
      domain-specific marker genes and a smooth spatial gradient.

        counts.csv        bin x gene count matrix
        coords.csv        barcode, x, y
        true_domains.csv  barcode, true_domain (ground truth)

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_spatial.py --outdir demo_data/spatial
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42


def build(args: argparse.Namespace) -> dict:
    rng = np.random.default_rng(args.seed)
    g = args.grid
    xs, ys = np.meshgrid(np.arange(g), np.arange(g))
    xs = xs.ravel().astype(float)
    ys = ys.ravel().astype(float)
    n = len(xs)
    barcodes = [f"BIN{i + 1:05d}" for i in range(n)]

    # quadrant domains + a stripe through the middle for a 4th class
    dom = np.where(xs < g / 2,
                   np.where(ys < g / 2, "Q1", "Q2"),
                   np.where(ys < g / 2, "Q3", "Q4"))
    # spatial gradient (e.g. developmental axis)
    grad = (xs / g) * 0.6 + (ys / g) * 0.4

    n_dom = 4
    n_marker = args.n_markers_per_domain
    n_gene_noise = args.n_genes - n_dom * n_marker
    marker_ids = [f"GeneD{d + 1}_{j + 1}" for d in range(n_dom) for j in range(n_marker)]
    noise_ids = [f"GeneN{i + 1:04d}" for i in range(n_gene_noise)]
    genes = marker_ids + noise_ids

    counts = np.zeros((n, len(genes)), dtype=int)
    # noise genes: baseline + gradient contribution for some
    base = rng.gamma(2.0, 0.7, n_gene_noise)
    for k in range(n_gene_noise):
        grad_eff = 6.0 if k % 20 == 0 else 0.0  # some genes are spatially graded
        mu = base[k] + grad_eff * grad
        counts[:, n_dom * n_marker + k] = rng.poisson(np.clip(mu, 0.01, None))
    # domain markers
    for d in range(n_dom):
        m = (dom == f"Q{d + 1}")
        for j in range(n_marker):
            gi = d * n_marker + j
            counts[m, gi] = rng.poisson(20.0 * rng.uniform(0.8, 1.3), size=int(m.sum()))
            counts[~m, gi] = rng.poisson(0.2, size=int((~m).sum()))

    df = pd.DataFrame(counts, columns=genes)
    df.insert(0, "barcode", barcodes)
    coords = pd.DataFrame({"barcode": barcodes, "x": xs, "y": ys})
    truth = pd.DataFrame({"barcode": barcodes, "true_domain": dom})
    # suggest marker genes to map
    return {"counts": df, "coords": coords, "truth": truth,
            "example_genes": [marker_ids[0], marker_ids[n_marker], marker_ids[2 * n_marker]]}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a synthetic spatial demo dataset.")
    parser.add_argument("--outdir", default="demo_data/spatial")
    parser.add_argument("--grid", type=int, default=50)
    parser.add_argument("--n-genes", type=int, default=600)
    parser.add_argument("--n-markers-per-domain", type=int, default=12)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = build(args)
    d["counts"].to_csv(outdir / "counts.csv", index=False)
    d["coords"].to_csv(outdir / "coords.csv", index=False)
    d["truth"].to_csv(outdir / "true_domains.csv", index=False)

    print(f"[make_demo_spatial] wrote {outdir}/ :")
    print(f"  counts.csv       ({d['counts'].shape[0]} bins x {d['counts'].shape[1] - 1} genes)")
    print(f"  coords.csv       (50x50 grid)")
    print(f"  true_domains.csv (4 quadrant domains)")
    print(f"  example genes for --genes: {','.join(d['example_genes'])}")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
