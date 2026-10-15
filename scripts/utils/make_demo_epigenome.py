"""make_demo_epigenome.py -- synthetic epigenome demo dataset generator.

Title:    Rice-skills epigenome demo data generator
Description:
    Creates deterministic synthetic datasets for the epigenomics module:

        methylation_regions.csv  region x sample methylation fractions
                                 (3000 regions, 200 true DMRs)
        peak_counts.csv          peak x sample fragment counts
                                 (4000 peaks, 250 true differential peaks)
        metadata.csv             sample_id, group (control / salt)

    Values mimic a two-condition experiment (10 vs 10) with realistic noise:
    methylation ~ Beta-like variation (sd 0.05 on the fraction scale);
    peak counts ~ negative-binomial-like overdispersion.

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_epigenome.py --outdir demo_data/epigenome
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
    chrs = [f"Chr{i:02d}" for i in range(1, 13)]

    # ---- methylation ----
    n_reg = args.n_meth_regions
    n_dmr = args.n_dmr
    base = rng.uniform(0.08, 0.92, n_reg)
    delta = np.zeros(n_reg)
    delta[:n_dmr] = rng.choice([-1, 1], n_dmr) * rng.uniform(0.12, 0.30, n_dmr)
    meth = np.zeros((n_reg, n))
    for i in range(n_reg):
        eff = np.where(group == "salt", delta[i], 0.0)
        meth[i] = np.clip(base[i] + eff + rng.normal(0, 0.045, n), 0.001, 0.999)
    reg_ids = [f"R{i + 1:05d}" for i in range(n_reg)]
    starts = rng.integers(1, 40_000_000, n_reg)
    reg_df = pd.DataFrame({
        "region_id": reg_ids,
        "chr": [chrs[i % 12] for i in range(n_reg)],
        "start": starts,
        "end": starts + 200,
    })
    for j, s in enumerate(samples):
        reg_df[s] = np.round(meth[:, j], 4)

    # ---- peaks ----
    n_peak = args.n_peaks
    n_dar = args.n_dar
    pk_base = rng.integers(40, 600, n_peak).astype(float)
    lfc = np.zeros(n_peak)
    lfc[:n_dar] = rng.choice([-1, 1], n_dar) * rng.uniform(1.4, 2.6, n_dar)
    disp = 0.15  # overdispersion
    counts = np.zeros((n_peak, n), dtype=int)
    for i in range(n_peak):
        mu = pk_base[i] * np.where(group == "salt", 2.0 ** lfc[i], 1.0)
        lam = rng.gamma(shape=1 / disp, scale=mu * disp)  # NB via gamma-poisson
        counts[i] = rng.poisson(lam)
    pk_ids = [f"P{i + 1:05d}" for i in range(n_peak)]
    pstarts = rng.integers(1, 40_000_000, n_peak)
    peak_df = pd.DataFrame({
        "peak_id": pk_ids,
        "chr": [chrs[i % 12] for i in range(n_peak)],
        "start": pstarts,
        "end": pstarts + rng.integers(200, 800, n_peak),
    })
    for j, s in enumerate(samples):
        peak_df[s] = counts[:, j]

    meta_df = pd.DataFrame({"sample_id": samples, "group": group})
    truth_meth = pd.DataFrame({"region_id": reg_ids[:n_dmr], "true_delta": np.round(delta[:n_dmr], 3)})
    truth_peak = pd.DataFrame({"peak_id": pk_ids[:n_dar], "true_log2FC": np.round(lfc[:n_dar], 3)})
    return {"meth": reg_df, "peaks": peak_df, "meta": meta_df,
            "truth_meth": truth_meth, "truth_peak": truth_peak,
            "n_dmr": n_dmr, "n_dar": n_dar}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate synthetic epigenome demo data.")
    parser.add_argument("--outdir", default="demo_data/epigenome")
    parser.add_argument("--n-meth-regions", type=int, default=3000)
    parser.add_argument("--n-dmr", type=int, default=200)
    parser.add_argument("--n-peaks", type=int, default=4000)
    parser.add_argument("--n-dar", type=int, default=250)
    parser.add_argument("--n-ctrl", type=int, default=10)
    parser.add_argument("--n-salt", type=int, default=10)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = build(args)
    d["meth"].to_csv(outdir / "methylation_regions.csv", index=False)
    d["peaks"].to_csv(outdir / "peak_counts.csv", index=False)
    d["meta"].to_csv(outdir / "metadata.csv", index=False)
    d["truth_meth"].to_csv(outdir / "ground_truth_dmr.csv", index=False)
    d["truth_peak"].to_csv(outdir / "ground_truth_dar.csv", index=False)

    print(f"[make_demo_epigenome] wrote {outdir}/ :")
    print(f"  methylation_regions.csv  ({d['meth'].shape[0]} regions, {d['n_dmr']} true DMRs)")
    print(f"  peak_counts.csv          ({d['peaks'].shape[0]} peaks, {d['n_dar']} true differential)")
    print(f"  metadata.csv / ground_truth_*.csv")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
