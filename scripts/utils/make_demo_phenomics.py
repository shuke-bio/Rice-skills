"""make_demo_phenomics.py -- synthetic phenomics demo dataset generator.

Title:    Rice-skills phenomics demo data generator
Description:
    Creates deterministic synthetic datasets for the phenomics module:

        phenotype_long.csv  line x environment trial table (60 lines x 3 envs,
                            2 traits: grain_yield, plant_height; includes
                            a handful of injected outliers)
        spectra.csv         canopy reflectance at 6 bands (470/550/670/730/
                            800/850 nm) for 40 plots, salt vs control

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_phenomics.py --outdir demo_data/phenomics
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
    n_line = args.n_lines
    envs = ["E1_lowland", "E2_saline", "E3_drought"]
    lines = [f"L{i + 1:03d}" for i in range(n_line)]

    g_yield = rng.normal(0, 9.0, n_line)
    g_height = rng.normal(0, 4.0, n_line)
    env_yield = {"E1_lowland": 0.0, "E2_saline": -6.0, "E3_drought": -4.5}
    env_height = {"E1_lowland": 0.0, "E2_saline": -3.0, "E3_drought": -1.5}

    rows = []
    for li, line in enumerate(lines):
        for env in envs:
            yield_v = 120 + g_yield[li] + env_yield[env] + rng.normal(0, 6.5)
            height_v = 95 + g_height[li] + env_height[env] + rng.normal(0, 3.0)
            rows.append({"line_id": line, "env": env,
                         "grain_yield": round(yield_v, 2),
                         "plant_height": round(height_v, 2)})
    pheno = pd.DataFrame(rows)

    # inject realistic outliers (data-entry style errors)
    idx = rng.choice(len(pheno), 5, replace=False)
    for i in idx:
        pheno.loc[i, "grain_yield"] = round(pheno.loc[i, "grain_yield"] * rng.uniform(1.8, 2.3), 2)

    # canopy spectra: salt reduces NIR, raises red band
    n_plot = args.n_plots
    groups = ["control", "salt"] * (n_plot // 2)
    bands = [470, 550, 670, 730, 800, 850]
    spec = {"sample_id": [f"P{i + 1:03d}" for i in range(n_plot)],
            "group": groups}
    for b in bands:
        base = {470: 0.045, 550: 0.09, 670: 0.06, 730: 0.24, 800: 0.42, 850: 0.46}[b]
        vals = base + rng.normal(0, base * 0.08, n_plot)
        salt_eff = {470: +0.010, 550: +0.005, 670: +0.012, 730: -0.02, 800: -0.05, 850: -0.05}[b]
        vals = np.where(np.array(groups) == "salt", vals + salt_eff, vals)
        spec[f"X{b}"] = np.round(np.clip(vals, 0.001, None), 4)
    spectra = pd.DataFrame(spec)

    return {"pheno": pheno, "spectra": spectra}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate synthetic phenomics demo data.")
    parser.add_argument("--outdir", default="demo_data/phenomics")
    parser.add_argument("--n-lines", type=int, default=60)
    parser.add_argument("--n-plots", type=int, default=40)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = build(args)
    d["pheno"].to_csv(outdir / "phenotype_long.csv", index=False)
    d["spectra"].to_csv(outdir / "spectra.csv", index=False)

    print(f"[make_demo_phenomics] wrote {outdir}/ :")
    print(f"  phenotype_long.csv  ({d['pheno'].shape[0]} records = 60 lines x 3 envs, 2 traits, 5 outliers)")
    print(f"  spectra.csv         ({d['spectra'].shape[0]} plots x 6 bands)")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
