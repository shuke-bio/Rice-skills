"""make_demo_single_cell.py -- synthetic single-cell demo dataset generator.

Title:    Rice-skills single-cell demo data generator
Description:
    Creates a small deterministic synthetic scRNA-seq dataset:
      1500 cells x 1000 genes, 3 cell types (500 cells each), with
      cell-type marker genes (10x enriched) and realistic count noise
      (negative-binomial-like).

        counts.csv        cell x gene count matrix (sparse-like integers)
        true_labels.csv   cell_id, true_type (ground truth for validation)
        marker_sets.csv   cell_type + marker gene columns (for --markers)

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_single_cell.py --outdir demo_data/single_cell
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42
CELL_TYPES = ["meristem", "epidermis", "cortex"]


def build(args: argparse.Namespace) -> dict:
    rng = np.random.default_rng(args.seed)
    n_per = args.n_cells_per_type
    n_ct = len(CELL_TYPES)
    n_cells = n_per * n_ct
    n_marker = args.n_markers_per_type
    n_genes = args.n_genes
    n_other = n_genes - n_marker * n_ct

    cells = [f"CELL{i + 1:05d}" for i in range(n_cells)]
    types = np.repeat(CELL_TYPES, n_per)
    marker_ids = [f"GeneM{ct + 1}_{j + 1}" for ct in range(n_ct) for j in range(n_marker)]
    other_ids = [f"GeneN{i + 1:04d}" for i in range(n_other)]
    genes = marker_ids + other_ids

    counts = np.zeros((n_cells, n_genes), dtype=int)
    for ct in range(n_ct):
        rows = np.where(types == CELL_TYPES[ct])[0]
        # baseline expression
        base = rng.gamma(shape=2.0, scale=0.6, size=n_other)
        for i in rows:
            counts[i, n_marker * n_ct:] = rng.poisson(base * rng.uniform(0.7, 1.3))
        # marker genes (strong in this type)
        for j in range(n_marker):
            g = ct * n_marker + j
            counts[rows, g] = rng.poisson(25.0 * rng.uniform(0.7, 1.4), size=len(rows))
        # marker genes from other types (near zero)
        for c2 in range(n_ct):
            if c2 == ct:
                continue
            for j in range(n_marker):
                g = c2 * n_marker + j
                counts[rows, g] = rng.poisson(0.15, size=len(rows))

    df = pd.DataFrame(counts, columns=genes)
    df.insert(0, "cell_id", cells)
    labels = pd.DataFrame({"cell_id": cells, "true_type": types})
    marker_rows = []
    for ct in range(n_ct):
        row = {"cell_type": CELL_TYPES[ct]}
        for j in range(n_marker):
            row[f"m{j + 1}"] = marker_ids[ct * n_marker + j]
        marker_rows.append(row)
    markers = pd.DataFrame(marker_rows)
    return {"counts": df, "labels": labels, "markers": markers}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a synthetic single-cell demo dataset.")
    parser.add_argument("--outdir", default="demo_data/single_cell")
    parser.add_argument("--n-cells-per-type", type=int, default=500)
    parser.add_argument("--n-genes", type=int, default=1000)
    parser.add_argument("--n-markers-per-type", type=int, default=8)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = build(args)
    d["counts"].to_csv(outdir / "counts.csv", index=False)
    d["labels"].to_csv(outdir / "true_labels.csv", index=False)
    d["markers"].to_csv(outdir / "marker_sets.csv", index=False)

    print(f"[make_demo_single_cell] wrote {outdir}/ :")
    print(f"  counts.csv       ({d['counts'].shape[0]} cells x {d['counts'].shape[1] - 1} genes)")
    print(f"  true_labels.csv  (3 cell types, ground truth)")
    print(f"  marker_sets.csv  (3 sets x 8 markers, for --markers)")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
