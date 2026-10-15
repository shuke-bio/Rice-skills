"""make_demo_pangenome.py -- synthetic pangenome PAV demo dataset generator.

Title:    Rice-skills pangenome demo data generator
Description:
    Creates a small, fully deterministic synthetic PAV (presence/absence
    variation) dataset in the Rice-skills pangenome input format:

        pav_matrix.csv     marker x sample presence matrix (0/1, ~3% NA)
        marker_info.csv    marker_id, chr, pos
        phenotype.csv      sample_id + two traits (continuous + categorical)

    The data mimic a two-subpopulation panel (A-like / B-like) with:
      - core markers        (always present)
      - variable markers    (frequency gradient, some subgroup-differentiated)
      - private markers     (single sample)
      - 5 causal markers    (presence associated with the continuous trait)
    so that PAV classification, structure PCA and association find
    plausible signals.

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_pangenome.py --outdir demo_data/pangenome
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
    n_a = args.n_subgroup_a
    n_b = args.n_subgroup_b
    n = n_a + n_b
    samples = [f"A{i + 1:03d}" for i in range(n_a)] + [f"B{i + 1:03d}" for i in range(n_b)]
    group = np.array(["A"] * n_a + ["B"] * n_b)

    n_markers = args.n_markers
    n_core = int(n_markers * 0.30)
    n_private = int(n_markers * 0.05)
    n_causal = 3
    n_var = n_markers - n_core - n_private - n_causal

    chr_ids = [f"chr{i}" for i in range(1, 13)]
    marker_ids, chrs, poss = [], [], []
    for i in range(n_markers):
        c = chr_ids[i % 12]
        marker_ids.append(f"MKR{i + 1:04d}")
        chrs.append(c)
        poss.append(int(rng.integers(1, 43_000_000)))
    poss = np.array(poss)

    pav = np.zeros((n_markers, n), dtype=int)
    # core: present everywhere
    pav[:n_core, :] = 1
    # subgroup-differentiated variable markers
    for i in range(n_core, n_core + n_var):
        base = rng.uniform(0.05, 0.95)
        if i % 3 == 0:  # A-enriched
            p = np.where(group == "A", min(0.95, base + 0.35), max(0.05, base - 0.35))
        elif i % 3 == 1:  # B-enriched
            p = np.where(group == "B", min(0.95, base + 0.35), max(0.05, base - 0.35))
        else:  # neutral
            p = np.full(n, base)
        pav[i, :] = (rng.random(n) < p).astype(int)
    # private: exactly one sample
    for k in range(n_private):
        i = n_core + n_var + k
        pav[i, :] = 0
        pav[i, int(rng.integers(0, n))] = 1
    # causal markers: presence raises the trait
    causal_idx = np.arange(n_core + n_var + n_private, n_markers)
    for i in causal_idx:
        pav[i, :] = (rng.random(n) < 0.5).astype(int)

    # 2% missing cells
    mask = rng.random((n_markers, n)) < 0.02
    pav_float = pav.astype(float)
    pav_float[mask] = np.nan

    # phenotype: continuous trait driven by causal markers (dosage-like sum)
    causal_sum = pav[causal_idx, :].sum(axis=0)
    trait = 100 + 15.0 * causal_sum + rng.normal(0, 5.0, n)
    # a categorical trait (e.g. subgroup) for demonstrating Fisher tests
    subtype = np.where(rng.random(n) < 0.5, "early", "late")

    pav_df = pd.DataFrame(pav_float, columns=samples)
    pav_df.insert(0, "marker_id", marker_ids)
    info_df = pd.DataFrame({"marker_id": marker_ids, "chr": chrs, "pos": poss})
    pheno_df = pd.DataFrame({"sample_id": samples,
                             "grain_yield": np.round(trait, 2),
                             "maturity_group": subtype})

    return {"pav": pav_df, "info": info_df, "pheno": pheno_df,
            "n_causal": n_causal, "causal": [marker_ids[i] for i in causal_idx]}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate a synthetic pangenome PAV demo dataset."
    )
    parser.add_argument("--outdir", default="demo_data/pangenome")
    parser.add_argument("--n-subgroup-a", type=int, default=30)
    parser.add_argument("--n-subgroup-b", type=int, default=30)
    parser.add_argument("--n-markers", type=int, default=200)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    d = build(args)
    d["pav"].to_csv(outdir / "pav_matrix.csv", index=False)
    d["info"].to_csv(outdir / "marker_info.csv", index=False)
    d["pheno"].to_csv(outdir / "phenotype.csv", index=False)

    print(f"[make_demo_pangenome] wrote {outdir}/ :")
    print(f"  pav_matrix.csv  ({d['pav'].shape[0]} markers x {d['pav'].shape[1] - 1} samples)")
    print(f"  marker_info.csv ({len(d['info'])} markers, 12 chromosomes)")
    print(f"  phenotype.csv   ({len(d['pheno'])} samples, 2 traits)")
    print(f"  {d['n_causal']} causal markers: {', '.join(d['causal'])}")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
