"""01_promoter_kmer.py -- promoter k-mer enrichment for candidate gene sets.

Title:    Rice-skills cis-element screen (step 01)
Description:
    Detects over-represented k-mer motifs in the promoters of a candidate
    gene set (e.g. DEGs, GWAS candidates, TF targets):
      1. extracts upstream promoter sequences from the genome using GFF3
         coordinates (strand-aware; plus strand takes [TSS-up, TSS),
         minus strand takes the reverse complement);
      2. counts k-mers (both strands collapsed to canonical orientation);
      3. tests each k-mer with Fisher's exact test on "genes whose promoter
         contains the k-mer" vs background promoters, BH-corrected;
      4. reports fold enrichment and the top hits. Top k-mers can then be
         matched against known TF motifs via TOMTOM/PlantTFDB (see
         docs/tf_motif_usage notes in the module documentation).
    Pure Python (scipy only) — no MEME/FIMO required for the screen.

Usage:
    python scripts/tf/01_promoter_kmer.py \
        --genome genome.fa --gff annotation.gff3 \
        --genes candidate_genes.csv --outdir results/tf/01_promoter \
        --up 2000 --k 6

Inputs:
    --genome   genome FASTA (chromosome names must match GFF seqid)
    --gff      GFF3 with gene/mRNA features (attribute ID= or Name= used)
    --genes    CSV with a gene id column (first column used; header optional
               single column)
    --background-id-col  optional: same file may hold all genes; background =
               all annotated genes minus candidates (default behaviour)
    --up       promoter length upstream of TSS/ATG (default 2000)
    --k        k-mer length (default 6)
    --max-bg   maximum background genes sampled (default 3000)
    --outdir   output directory (default results/tf/01_promoter)
Outputs:
    promoters.fasta          extracted promoter sequences (all genes used)
    kmer_enrichment.csv      kmer, n_cand, freq_cand, n_bg, freq_bg,
                             fold_enrichment, p, fdr
    top_kmers_bar.*          top enriched k-mers bar plot
    run_info.txt
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

import numpy as np
import pandas as pd
from scipy import stats

from data_io import set_seed, write_run_info
from plot_theme import save_nature_figure, setup_nature_style, group_colors

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta_indexed(path: str) -> tuple[dict, list]:
    """Load a genome FASTA into memory (chromosomes keyed by header token)."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line.upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs, list(seqs.keys())


def parse_gff_genes(path: str) -> dict:
    """Return gene_id -> (seqid, start, end, strand) from GFF3 gene/mRNA lines."""
    genes = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] not in ("gene", "mRNA"):
                continue
            attrs = f[8]
            gid = None
            for part in attrs.split(";"):
                part = part.strip()
                if part.startswith("ID="):
                    gid = part[3:].split(",")[0]
                    break
            if gid is None:
                for part in attrs.split(";"):
                    if part.strip().startswith("Name="):
                        gid = part.strip()[5:].split(",")[0]
                        break
            if gid is None:
                continue
            if f[2] == "mRNA" and gid in genes:
                continue
            genes[gid] = (f[0], int(f[3]), int(f[4]), f[6])
    return genes


def extract_promoters(genome: dict, genes: dict, up: int) -> tuple[dict, list]:
    seqs, skipped = {}, []
    for gid, (seqid, start, end, strand) in genes.items():
        chrom = genome.get(seqid)
        if chrom is None:
            skipped.append((gid, "chromosome not found"))
            continue
        if strand == "+":
            tss = start
            seg_start, seg_end = max(1, tss - up), tss - 1
        else:
            tss = end
            seg_start, seg_end = tss + 1, min(len(chrom), tss + up)
        if seg_end <= seg_start:
            skipped.append((gid, "promoter out of range"))
            continue
        s = chrom[seg_start - 1:seg_end]
        if strand == "-":
            s = s.translate(COMP)[::-1]
        if len(s) < up // 2:
            skipped.append((gid, "promoter too short"))
            continue
        seqs[gid] = s
    return seqs, skipped


def revcomp_canon(kmer: str) -> str:
    rc = kmer.translate(COMP)[::-1]
    return min(kmer, rc)


def kmer_profile(seq_dict: dict, k: int) -> Counter:
    c = Counter()
    for s in seq_dict.values():
        seen = set()
        for i in range(len(s) - k + 1):
            km = revcomp_canon(s[i:i + k])
            seen.add(km)
        for km in seen:
            c[km] += 1
    return c


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Promoter k-mer enrichment screen.")
    p.add_argument("--genome", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--genes", required=True, help="CSV: candidate gene ids (first column)")
    p.add_argument("--up", type=int, default=2000)
    p.add_argument("--k", type=int, default=6)
    p.add_argument("--max-bg", type=int, default=3000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--outdir", default="results/tf/01_promoter")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    set_seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_nature_style()

    genome, _ = read_fasta_indexed(args.genome)
    print(f"[promoter_kmer] genome chromosomes: {len(genome)}")
    genes = parse_gff_genes(args.gff)
    print(f"[promoter_kmer] genes in GFF: {len(genes)}")

    gdf = pd.read_csv(args.genes, header=None, names=["gene_id"])
    candidates = [str(g) for g in gdf["gene_id"].dropna().unique() if str(g) in genes]
    if not candidates:
        raise SystemExit("No candidate gene ids matched the GFF. Check id format (ID=/Name=).")
    bg_pool = sorted(set(genes) - set(candidates))
    if len(bg_pool) > args.max_bg:
        import random
        rng = random.Random(args.seed)
        bg_pool = rng.sample(bg_pool, args.max_bg)
    print(f"[promoter_kmer] candidates: {len(candidates)} | background: {len(bg_pool)}")

    promoters, skipped = extract_promoters(genome, genes, args.up)
    with open(outdir / "promoters.fasta", "w") as fh:
        for gid, s in promoters.items():
            fh.write(f">{gid}\n{s}\n")
    if skipped:
        print(f"  (skipped {len(skipped)} genes, e.g. {skipped[0]})")

    cand_use = [g for g in candidates if g in promoters]
    bg_use = [g for g in bg_pool if g in promoters]
    cand_prof = kmer_profile({g: promoters[g] for g in cand_use}, args.k)
    bg_prof = kmer_profile({g: promoters[g] for g in bg_use}, args.k)

    n_cand, n_bg = len(cand_use), len(bg_use)
    rows = []
    for km in set(cand_prof) | set(bg_prof):
        a = cand_prof[km]              # candidate promoters containing km
        b = max(n_cand - a, 0)
        c = bg_prof[km]                # background promoters containing km
        dd = max(n_bg - c, 0)
        table = [[a, b], [c, dd]]
        try:
            _, p = stats.fisher_exact(table, alternative="greater")
        except ValueError:
            p = 1.0
        fold = (a / n_cand) / max(c / n_bg, 1e-9) if c else float("inf") if a else 0.0
        rows.append({"kmer": km, "n_cand": a, "freq_cand": round(a / n_cand, 4),
                     "n_bg": c, "freq_bg": round(c / n_bg, 4),
                     "fold_enrichment": round(fold, 3), "p": p})
    if not rows:
        raise SystemExit("No k-mers counted; check sequence extraction.")
    df = pd.DataFrame(rows)
    df["fdr"] = np.minimum(1.0, stats.false_discovery_control(df["p"].values, method="bh"))
    df = df.sort_values("p").reset_index(drop=True)
    df.to_csv(outdir / "kmer_enrichment.csv", index=False)
    print(f"[promoter_kmer] k-mers tested: {len(df)}")

    import matplotlib.pyplot as plt
    top = df.head(15)
    fig, ax = plt.subplots()
    ax.barh(range(len(top))[::-1], top["fold_enrichment"].clip(upper=20),
            color=group_colors(3)[0], edgecolor="gray", linewidth=0.3)
    ax.set_yticks(range(len(top))[::-1])
    ax.set_yticklabels(top["kmer"], fontsize=5, family="monospace")
    ax.set_xlabel("Fold enrichment (capped at 20)")
    ax.set_title(f"Top {len(top)} enriched {args.k}-mers in candidate promoters")
    save_nature_figure(fig, outdir / "top_kmers_bar", width_mm=100, height_mm=90)
    plt.close(fig)

    write_run_info(outdir, script="01_promoter_kmer.py",
                   params={"genome": str(args.genome), "gff": str(args.gff),
                           "n_candidates": len(cand_use), "n_background": len(bg_use),
                           "up": args.up, "k": args.k, "seed": args.seed})
    print(f"[promoter_kmer] done -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
