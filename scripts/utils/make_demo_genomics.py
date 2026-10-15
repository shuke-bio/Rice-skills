"""make_demo_genomics.py -- synthetic rice genomics demo dataset generator.

Title:    Rice-skills genomics demo data generator
Description:
    Creates a small, fully deterministic synthetic population-genomics / GWAS
    dataset in the Rice-skills genomics input format:

        genotype.csv    SNP x sample dosage matrix (0/1/2; rows = SNPs)
        snp_info.csv    snp_id, chr, pos
        phenotype.csv   sample_id + three traits
        annotation.gff3 minimal gene models covering the causal SNPs

    The data mimic a two-subpopulation panel (indica-like A / japonica-like B,
    plus admixed individuals) with ~10% differentiated SNPs and four causal
    loci, so that every downstream step (PCA, LD decay, Fst/pi, GWAS,
    candidate-gene overlap) finds a plausible signal.

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_genomics.py --outdir demo_data/genomics
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42
TRAITS = ["salt_tolerance_index", "grain_weight", "plant_height"]


def build(args: argparse.Namespace) -> dict:
    rng = np.random.default_rng(args.seed)
    n_a, n_b, n_mix = args.n_indica, args.n_japonica, args.n_admixed
    samples = ([f"A{i + 1:03d}" for i in range(n_a)]
               + [f"B{i + 1:03d}" for i in range(n_b)]
               + [f"M{i + 1:03d}" for i in range(n_mix)])
    n = len(samples)
    m = args.n_snps
    chr_ids = [f"Chr{i:02d}" for i in range(1, 13)]
    per_chr = int(np.ceil(m / len(chr_ids)))

    # --- SNP positions and allele frequencies ----------------------------------
    snp_ids, snp_chr, snp_pos = [], [], []
    freqs_a = []
    for chrom in chr_ids:
        pos = np.sort(rng.integers(1, 40_000_000, size=per_chr))
        for k, p in enumerate(pos):
            snp_ids.append(f"s{len(snp_ids) + 1:05d}")
            snp_chr.append(chrom)
            snp_pos.append(int(p))
        freqs_a.append(rng.uniform(0.05, 0.95, size=per_chr))
    # trim to exactly m SNPs (per-chromosome counts use ceil)
    snp_ids = snp_ids[:m]
    snp_chr = snp_chr[:m]
    snp_pos = snp_pos[:m]
    freqs_a = np.concatenate(freqs_a)[:m]
    freqs_b = freqs_a.copy()

    # ~12% differentiated SNPs (structure signal)
    diff_idx = rng.choice(m, int(m * 0.12), replace=False)
    shift = rng.uniform(0.25, 0.50, size=len(diff_idx)) * rng.choice([-1, 1], len(diff_idx))
    freqs_b[diff_idx] = np.clip(freqs_a[diff_idx] + shift, 0.05, 0.95)

    # --- genotypes (dosages 0/1/2) ---------------------------------------------
    geno = np.zeros((m, n), dtype=float)
    for j in range(m):
        geno[j, :n_a] = rng.binomial(2, freqs_a[j], n_a)
        geno[j, n_a:n_a + n_b] = rng.binomial(2, freqs_b[j], n_b)
        p_mix = (freqs_a[j] + freqs_b[j]) / 2
        geno[j, n_a + n_b:] = rng.binomial(2, p_mix, n_mix)

    # --- local LD structure ------------------------------------------------------
    # Real panels carry linkage disequilibrium; without it the LD-decay figure
    # would be a flat line. For ~40% of SNPs (per chromosome, in order), copy the
    # previous SNP's dosage in 85% of individuals to create short-range LD blocks.
    snp_chr_arr = np.array(snp_chr)
    for j in range(1, m):
        if snp_chr_arr[j] != snp_chr_arr[j - 1]:
            continue
        if rng.random() < 0.40:
            copy_mask = rng.random(n) < 0.85
            geno[j, copy_mask] = geno[j - 1, copy_mask]

    # missing genotypes (0.5%)
    miss = rng.random((m, n)) < 0.005
    geno[miss] = np.nan

    # --- causal loci ------------------------------------------------------------
    # pick well-behaved common SNPs (MAF in [0.15, 0.45]) as causal
    af = np.nanmean(geno, axis=1) / 2
    maf = np.minimum(af, 1 - af)
    candidates = np.where((maf > 0.15) & (maf < 0.45))[0]
    causal = rng.choice(candidates, 4, replace=False)
    causal_chr = [snp_chr[i] for i in causal]

    # --- phenotypes -------------------------------------------------------------
    def dosage(j: int) -> np.ndarray:
        d = geno[j].copy()
        d[np.isnan(d)] = np.nanmean(d)
        return d

    subpop_effect = np.array([1.0] * n_a + [-1.0] * n_b + [0.0] * n_mix)
    salt = (0.8 * dosage(causal[0]) + 0.55 * dosage(causal[1])
            + rng.normal(0, 0.6, n))
    grain = 0.9 * dosage(causal[2]) + rng.normal(0, 0.7, n)
    height = 5.0 * subpop_effect + 0.6 * dosage(causal[3]) + rng.normal(0, 0.8, n)

    pheno = pd.DataFrame({
        "sample_id": samples,
        "group": (["indica"] * n_a + ["japonica"] * n_b + ["admixed"] * n_mix),
        TRAITS[0]: np.round(salt, 3),
        TRAITS[1]: np.round(grain, 3),
        TRAITS[2]: np.round(height, 3),
    })

    # --- annotation (GFF3) covering each causal SNP ----------------------------
    rng_ann = np.random.default_rng(args.seed + 1)
    gene_rows = []
    gene_no = 0
    for chrom in chr_ids:
        n_genes = 25
        starts = np.sort(rng_ann.integers(1, 39_000_000, size=n_genes))
        for s in starts:
            gene_no += 1
            gene_rows.append((chrom, int(s), int(s) + 2000, f"DEMO_g{gene_no:05d}"))
    # guarantee one gene around every causal SNP
    for gi, j in enumerate(causal, start=1):
        gene_no += 1
        s = max(1, int(snp_pos[j]) - 500)
        gene_rows.append((snp_chr[j], s, s + 3000, f"CAUSAL_g{gi:02d}_{snp_ids[j]}"))
    ann = pd.DataFrame(gene_rows, columns=["chr", "start", "end", "gene_id"])
    ann["score"] = "."
    ann["strand"] = "+"
    ann["phase"] = "."
    ann["attrs"] = ("ID=" + ann["gene_id"] + ";Name=" + ann["gene_id"])

    return {
        "geno": pd.DataFrame(geno, index=snp_ids, columns=samples).rename_axis("snp_id"),
        "snp_info": pd.DataFrame({"snp_id": snp_ids, "chr": snp_chr, "pos": snp_pos}),
        "pheno": pheno,
        "ann": ann[["chr", "score", "start", "end", "score", "strand", "phase", "attrs"]]
                   .assign(col3="gene")[["chr", "col3", "start", "end", "score", "strand",
                                         "phase", "attrs"]],
        "causal": [snp_ids[j] for j in causal],
        "causal_chr": causal_chr,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate a synthetic rice genomics demo dataset (GWAS + popgen)."
    )
    parser.add_argument("--outdir", default="demo_data/genomics")
    parser.add_argument("--n-indica", type=int, default=140)
    parser.add_argument("--n-japonica", type=int, default=140)
    parser.add_argument("--n-admixed", type=int, default=20)
    parser.add_argument("--n-snps", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    out = build(args)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out["geno"].to_csv(outdir / "genotype.csv", float_format="%.0f", na_rep="")
    out["snp_info"].to_csv(outdir / "snp_info.csv", index=False)
    out["pheno"].to_csv(outdir / "phenotype.csv", index=False)
    # GFF3 body: chr \t gene \t start \t end \t . \t + \t . \t attrs
    with open(outdir / "annotation.gff3", "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        for row in out["ann"].itertuples(index=False):
            fh.write(f"{row.chr}\tDEMO\tgene\t{row.start}\t{row.end}\t.\t+\t.\t{row.attrs}\n")

    # --- demo gene for the design module (sgRNA / primer design) ----------------
    # A single 2,400 bp gene with three exons and two introns, so that qPCR
    # exon-junction primers and CRISPR targets can be demonstrated end-to-end.
    rng_gene = np.random.default_rng(args.seed + 7)
    seq = "".join(rng_gene.choice(list("ACGT"), size=2400))
    with open(outdir / "demo_gene.fa", "w", encoding="utf-8") as fh:
        fh.write(">DEMO_Os01g99999 demo gene (3 exons)\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    with open(outdir / "demo_gene.gff3", "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        fh.write("Chr01\tDEMO\tgene\t1\t2400\t.\t+\t.\tID=DEMO_Os01g99999;Name=DEMO_Os01g99999\n")
        exon_spans = [(1, 600), (751, 1350), (1531, 2400)]
        cds_spans = [(1, 600), (751, 1350), (1531, 2400)]
        for idx, (s, e) in enumerate(exon_spans, start=1):
            fh.write(f"Chr01\tDEMO\texon\t{s}\t{e}\t.\t+\t.\tID=DEMO_Os01g99999.exon{idx}\n")
        for idx, (s, e) in enumerate(cds_spans, start=1):
            fh.write(f"Chr01\tDEMO\tCDS\t{s}\t{e}\t.\t+\t.\tID=DEMO_Os01g99999.cds{idx}\n")

    n_total = args.n_indica + args.n_japonica + args.n_admixed
    print(f"Genomics demo data written to: {outdir.resolve()}")
    print(f"  genotype.csv   : {args.n_snps} SNPs x {n_total} samples "
          "(dosages 0/1/2, ~0.5% missing)")
    print(f"  snp_info.csv   : 12 chromosomes")
    print(f"  phenotype.csv  : {', '.join(TRAITS)}")
    print(f"  annotation.gff3: gene models incl. causal loci")
    print(f"  causal SNPs    : {', '.join(out['causal'])}")
    print("\nNext steps:")
    print("  Rscript scripts/gwas/03_association.R \\")
    print(f"      --genotype {outdir}/genotype.csv --snp-info {outdir}/snp_info.csv \\")
    print(f"      --phenotype {outdir}/phenotype.csv --trait salt_tolerance_index")
    print("  Rscript scripts/popgen/01_pca_structure.R \\")
    print(f"      --genotype {outdir}/genotype.csv --snp-info {outdir}/snp_info.csv")
    print("\nNOTE: synthetic data are for pipeline testing only, not for biology.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
