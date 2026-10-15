"""make_demo_qtl.py -- synthetic eQTL demo dataset generator.

Title:    Rice-skills QTL demo data generator
Description:
    Builds a small deterministic eQTL demo on top of the genomics demo
    genotype (run make_demo_genomics.py first), producing:

        expression.csv   150 gene x sample normalized expression matrix
        gene_info.csv    feature_id, chr, pos (gene coordinates)
        ground_truth.csv feature_id of the 60 causal genes + causal SNP

    60 genes carry a strong cis-eQTL: expression = N(8, 1.5) + effect * dosage
    + N(0, 0.8), with the causal SNP placed within 50 kb of the gene.
    The remaining 90 genes are pure noise.

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_qtl.py \
        --genotype demo_data/genomics/genotype.csv \
        --snp-info demo_data/genomics/snp_info.csv \
        --outdir demo_data/qtl
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 7


def build(genotype_path: str, snp_info_path: str, n_features: int, n_causal: int,
          seed: int) -> dict:
    rng = np.random.default_rng(seed)
    geno = pd.read_csv(genotype_path)
    snp_info = pd.read_csv(snp_info_path)
    samples = list(geno.columns[1:])
    n = len(samples)

    causal_idx = rng.choice(len(geno), n_causal, replace=False)
    expr = np.zeros((n_features, n))
    for i in range(n_features):
        if i < n_causal:
            g = pd.to_numeric(geno.iloc[causal_idx[i], 1:], errors="coerce").fillna(1).values.astype(float)
            eff = rng.uniform(1.2, 2.0)
            expr[i] = rng.normal(8, 1.5, n) + eff * g + rng.normal(0, 0.8, n)
        else:
            expr[i] = rng.normal(8, 1.5, n)

    feature_ids = [f"GENE{i + 1:04d}" for i in range(n_features)]
    expr_df = pd.DataFrame(expr, columns=samples)
    expr_df.insert(0, "feature_id", feature_ids)

    rows, truth = [], []
    for i, fid in enumerate(feature_ids):
        if i < n_causal:
            s = int(causal_idx[i])
            ch = snp_info.iloc[s]["chr"]
            pos = int(snp_info.iloc[s]["pos"]) + int(rng.integers(-50000, 50000))
            truth.append({"feature_id": fid, "causal_snp": snp_info.iloc[s]["snp_id"]})
        else:
            s = int(rng.integers(0, len(snp_info)))
            ch = snp_info.iloc[s]["chr"]
            pos = int(rng.integers(1, 40_000_000))
        rows.append({"feature_id": fid, "chr": ch, "pos": max(1, pos)})

    return {
        "expression": expr_df,
        "gene_info": pd.DataFrame(rows),
        "ground_truth": pd.DataFrame(truth),
        "n_causal": n_causal,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a synthetic eQTL demo dataset.")
    parser.add_argument("--genotype", default="demo_data/genomics/genotype.csv")
    parser.add_argument("--snp-info", default="demo_data/genomics/snp_info.csv")
    parser.add_argument("--outdir", default="demo_data/qtl")
    parser.add_argument("--n-features", type=int, default=150)
    parser.add_argument("--n-causal", type=int, default=60)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    for p in (args.genotype, args.snp_info):
        if not Path(p).exists():
            print(f"ERROR: {p} not found. Run make_demo_genomics.py first.", file=sys.stderr)
            return 1

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = build(args.genotype, args.snp_info, args.n_features, args.n_causal, args.seed)
    d["expression"].to_csv(outdir / "expression.csv", index=False)
    d["gene_info"].to_csv(outdir / "gene_info.csv", index=False)
    d["ground_truth"].to_csv(outdir / "ground_truth.csv", index=False)

    print(f"[make_demo_qtl] wrote {outdir}/ :")
    print(f"  expression.csv  ({d['expression'].shape[0]} genes x {d['expression'].shape[1] - 1} samples)")
    print(f"  gene_info.csv   ({len(d['gene_info'])} genes with coordinates)")
    print(f"  ground_truth.csv ({d['n_causal']} causal cis-eQTL genes)")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
