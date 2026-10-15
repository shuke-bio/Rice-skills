"""make_demo_proteomics.py -- synthetic rice proteomics demo dataset generator.

Title:    Rice-skills proteomics demo data generator
Description:
    Creates a small, fully deterministic synthetic label-free proteomics
    dataset in the Rice-skills proteomics input format:

        proteins.csv         protein x sample abundance matrix (LFQ-like)
        metadata.csv         sample_id, group (control / salt), batch
        protein_gene_map.csv protein_id -> gene_id (RAP-DB style)

    1200 proteins, 24 samples (12 vs 12), with:
      - 80 truly differential proteins (log2FC ~ +/-1.5, group effect)
      - 10-15% missing values (MNAR-style: low-abundance proteins missing
        more often), the realistic label-free pattern
      - a mild batch effect on half the samples

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_proteomics.py --outdir demo_data/proteomics
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
    n_a, n_b = args.n_ctrl, args.n_salt
    n = n_a + n_b
    samples = [f"C{i + 1:02d}" for i in range(n_a)] + [f"S{i + 1:02d}" for i in range(n_b)]
    group = np.array(["control"] * n_a + ["salt"] * n_b)
    # batches interleaved WITHIN each group (not confounded with treatment)
    batch = np.array([["B1", "B2"][i % 2] for i in range(n)])

    n_prot = args.n_proteins
    n_de = args.n_de
    base = rng.normal(22, 1.2, n_prot)          # log2-scale baseline abundance
    lfc = np.zeros(n_prot)
    lfc[:n_de] = rng.choice([-1.0, 1.0], n_de) * rng.uniform(1.5, 2.5, n_de)

    effect = np.where(group == "salt", 1.0, 0.0)
    batch_eff = np.where(batch == "B2", 0.35, 0.0)
    prot_ids = [f"OsProt{i + 1:05d}" for i in range(n_prot)]
    mat = np.zeros((n_prot, n))
    for i in range(n_prot):
        mat[i] = (base[i] + lfc[i] * effect + batch_eff
                  + rng.normal(0, 0.35, n))

    # MNAR-style missingness: lower-abundance proteins missing more
    abund_rank = (base - base.min()) / (base.max() - base.min())
    miss_rate = np.clip(0.35 - 0.30 * abund_rank, 0.02, 0.35)
    mask = rng.random((n_prot, n)) < miss_rate[:, None]
    mat[mask] = np.nan
    # make sure each protein keeps >= 2 values
    for i in range(n_prot):
        if np.isnan(mat[i]).sum() > n - 2:
            keep = rng.choice(n, 2, replace=False)
            drop = np.setdiff1d(np.arange(n), keep)
            mat[i, drop] = np.nan
            mat[i, keep] = base[i] + lfc[i] * effect[keep] + batch_eff[keep]

    prot_df = pd.DataFrame(np.round(mat, 3), columns=samples)
    prot_df.insert(0, "protein_id", prot_ids)
    meta_df = pd.DataFrame({"sample_id": samples, "group": group, "batch": batch})
    map_df = pd.DataFrame({
        "protein_id": prot_ids,
        "gene_id": [f"Os{1 + i % 12:02d}g{(i + 1) * 100:07d}" for i in range(n_prot)],
    })

    truth = pd.DataFrame({"protein_id": prot_ids[:n_de],
                          "true_log2FC": np.round(lfc[:n_de], 3)})
    return {"proteins": prot_df, "metadata": meta_df, "map": map_df,
            "truth": truth, "n_de": n_de}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a synthetic proteomics demo dataset.")
    parser.add_argument("--outdir", default="demo_data/proteomics")
    parser.add_argument("--n-proteins", type=int, default=1200)
    parser.add_argument("--n-de", type=int, default=80)
    parser.add_argument("--n-ctrl", type=int, default=12)
    parser.add_argument("--n-salt", type=int, default=12)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = build(args)
    d["proteins"].to_csv(outdir / "proteins.csv", index=False)
    d["metadata"].to_csv(outdir / "metadata.csv", index=False)
    d["map"].to_csv(outdir / "protein_gene_map.csv", index=False)
    d["truth"].to_csv(outdir / "ground_truth.csv", index=False)

    print(f"[make_demo_proteomics] wrote {outdir}/ :")
    print(f"  proteins.csv  ({d['proteins'].shape[0]} proteins x {d['proteins'].shape[1] - 1} samples)")
    print(f"  metadata.csv  (2 groups, 2 batches)")
    print(f"  protein_gene_map.csv / ground_truth.csv ({d['n_de']} true DE proteins)")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
