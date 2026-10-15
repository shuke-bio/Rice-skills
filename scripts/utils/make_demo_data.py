"""make_demo_data.py -- generate a small synthetic rice multi-omics demo dataset.

Title:    Rice-skills synthetic demo data generator
Description:
    Creates a small, fully deterministic synthetic dataset in the Rice-skills
    input format (count matrix, metabolite abundance matrix, sample metadata) so
    that the pipeline can be smoke-tested without any real data. The data mimic
    a salt-treatment experiment on rice with two batches, batch effects,
    missing values in the metabolome and a subset of differential
    genes/metabolites.

    IMPORTANT: the generated data are random and biologically meaningless.
    Use them only to verify that the pipeline runs end-to-end.

Usage:
    python scripts/utils/make_demo_data.py --outdir demo_data

Outputs (written to --outdir):
    counts.csv        gene x sample integer count matrix (rows = LOC_Os... genes)
    metabolites.csv   metabolite x sample abundance matrix (with annotations)
    metadata.csv      sample metadata (sample_id, group, batch, variety, traits)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42

GENE_CHROMOSOMES = [f"{i:02d}" for i in range(1, 13)]
METABOLITE_CLASSES = [
    "Flavonoids", "Alkaloids", "Lipids", "Amino acids",
    "Organic acids", "Carbohydrates", "Nucleotides", "Terpenoids",
]
# A small pool of real KEGG compound IDs (rice central metabolism / secondary
# metabolism); metabolites not in the pool get an empty kegg_id.
KEGG_COMPOUND_POOL = [
    "C00031", "C00022", "C00024", "C00026", "C00036", "C00037", "C00041",
    "C00042", "C00049", "C00058", "C00062", "C00064", "C00074", "C00078",
    "C00079", "C00082", "C00091", "C00122", "C00157", "C00158", "C00197",
    "C00227", "C00279", "C00311", "C00315", "C00418", "C00534", "C00631",
    "C00719", "C00803", "C00877", "C00900", "C01185", "C01236", "C01401",
    "C01571", "C02059", "C02291", "C02504", "C02916", "C03287", "C03617",
]


def build_metadata(groups: list[str], n_per_group: int, n_batches: int) -> pd.DataFrame:
    """Build sample metadata with batch, variety and two phenotypic traits."""
    rows = []
    for g_idx, group in enumerate(groups):
        for rep in range(1, n_per_group + 1):
            sample = f"{group[:3].upper()}_{rep}"
            batch = (g_idx + rep) % n_batches + 1
            # Traits differ by group with biological noise (for WGCNA demo).
            tiller = (12.0 if group == "control" else 8.0) + np.random.normal(0, 1.0)
            na = (0.31 if group == "control" else 1.87) + np.random.normal(0, 0.15)
            rows.append(
                {
                    "sample_id": sample,
                    "group": group,
                    "batch": batch,
                    "variety": "Nipponbare",
                    "tiller_number": round(float(tiller), 2),
                    "Na_content": round(float(na), 3),
                }
            )
    return pd.DataFrame(rows)


def build_counts(
    n_genes: int,
    samples: list[str],
    groups: list[str],
    n_per_group: int,
    de_fraction: float = 0.15,
    batch_fraction: float = 0.10,
) -> pd.DataFrame:
    """Build a gene x sample integer count matrix with group and batch effects."""
    gene_ids = [
        f"LOC_Os{np.random.choice(GENE_CHROMOSOMES)}g{10000 + i:05d}"
        for i in range(n_genes)
    ]
    gene_ids = list(dict.fromkeys(gene_ids))  # ensure uniqueness
    n_genes = len(gene_ids)

    lib_sizes = np.random.lognormal(mean=np.log(1.5e7), sigma=0.15, size=len(samples))
    base = np.random.lognormal(mean=6.0, sigma=1.2, size=n_genes)  # per-gene abundance

    is_de = np.random.rand(n_genes) < de_fraction
    log2fc = np.zeros(n_genes)
    log2fc[is_de] = np.random.normal(0.0, 1.6, size=is_de.sum())  # signed effect

    batch_col = []
    group_of = []
    for g_idx, group in enumerate(groups):
        group_of.extend([group] * n_per_group)
        batch_col.extend([(g_idx + r) % 2 + 1 for r in range(1, n_per_group + 1)])

    counts = np.zeros((n_genes, len(samples)), dtype=int)
    for j, (group, batch) in enumerate(zip(group_of, batch_col)):
        group_effect = np.where(group == groups[1], 2.0 ** log2fc, 1.0)
        batch_effect = np.where(
            np.random.rand(n_genes) < batch_fraction,
            1.0 + 0.3 * (batch - 1.5),
            1.0,
        )
        lam = base * group_effect * batch_effect * lib_sizes[j] / 1.5e7
        counts[:, j] = np.random.poisson(np.maximum(lam, 0.1))

    df = pd.DataFrame(counts, index=gene_ids, columns=samples)
    df.index.name = "gene_id"
    return df


def build_metabolites(
    n_metabolites: int,
    samples: list[str],
    groups: list[str],
    n_per_group: int,
    de_fraction: float = 0.20,
    missing_fraction: float = 0.08,
) -> pd.DataFrame:
    """Build a metabolite x sample abundance matrix with annotations and NAs."""
    met_ids = [f"M{i + 1:04d}" for i in range(n_metabolites)]
    names = [f"Compound-{i + 1:04d}" for i in range(n_metabolites)]
    classes = list(np.random.choice(METABOLITE_CLASSES, size=n_metabolites))
    kegg = [
        KEGG_COMPOUND_POOL[i % len(KEGG_COMPOUND_POOL)]
        if i % 3 == 0
        else ""
        for i in range(n_metabolites)
    ]

    base = np.random.lognormal(mean=12.0, sigma=1.0, size=n_metabolites)
    is_de = np.random.rand(n_metabolites) < de_fraction
    effect = np.zeros(n_metabolites)
    effect[is_de] = np.random.normal(0.0, 0.8, size=is_de.sum())

    group_of = []
    for g_idx, group in enumerate(groups):
        group_of.extend([group] * n_per_group)

    data = np.zeros((n_metabolites, len(samples)))
    for j, group in enumerate(group_of):
        fold = np.where(group == groups[1], np.exp(effect), 1.0)
        data[:, j] = base * fold * np.exp(np.random.normal(0, 0.15, n_metabolites))

    # Missing values: lower-abundance metabolites are more likely to be missing.
    miss_prob = missing_fraction * np.exp(-(base - base.min()) / (base.std() + 1e-9))
    miss_prob = np.clip(miss_prob, 0, 0.5)
    mask = np.random.rand(n_metabolites, len(samples)) < miss_prob[:, None]
    data[mask] = np.nan

    df = pd.DataFrame(data, index=met_ids, columns=samples)
    df.insert(0, "kegg_id", kegg)
    df.insert(0, "class", classes)
    df.insert(0, "metabolite_name", names)
    df.index.name = "metabolite_id"
    return df


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate a small synthetic rice multi-omics demo dataset."
    )
    parser.add_argument("--outdir", default="demo_data", help="Output directory")
    parser.add_argument("--n-genes", type=int, default=2000)
    parser.add_argument("--n-metabolites", type=int, default=300)
    parser.add_argument("--n-per-group", type=int, default=3)
    parser.add_argument("--groups", default="control,salt",
                        help="Comma-separated group names (first = reference)")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    np.random.seed(args.seed)
    groups = [g.strip() for g in args.groups.split(",") if g.strip()]
    if len(groups) < 2:
        parser.error("--groups must contain at least two groups")

    meta = build_metadata(groups, args.n_per_group, n_batches=2)
    samples = meta["sample_id"].tolist()

    counts = build_counts(args.n_genes, samples, groups, args.n_per_group)
    metabolites = build_metabolites(
        args.n_metabolites, samples, groups, args.n_per_group
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(outdir / "counts.csv")
    metabolites.to_csv(outdir / "metabolites.csv")
    meta.to_csv(outdir / "metadata.csv", index=False)

    n_missing = int(metabolites[samples].isna().sum().sum())
    print(f"Demo data written to: {outdir.resolve()}")
    print(f"  counts.csv       : {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"  metabolites.csv  : {metabolites.shape[0]} metabolites x {len(samples)} samples"
          f" ({n_missing} missing values)")
    print(f"  metadata.csv     : {len(meta)} samples, groups = {groups}")
    print("\nNext steps:")
    print(f"  python scripts/transcriptomics/01_qc_filter.py \\")
    print(f"      --counts {outdir}/counts.csv --metadata {outdir}/metadata.csv")
    print(f"  python scripts/metabolomics/01_preprocess.py \\")
    print(f"      --abundance {outdir}/metabolites.csv --metadata {outdir}/metadata.csv")
    print("\nNOTE: synthetic data are for pipeline testing only, not for biology.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
