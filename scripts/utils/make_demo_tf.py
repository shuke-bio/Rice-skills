"""make_demo_tf.py -- synthetic promoter-motif demo dataset generator.

Title:    Rice-skills TF motif demo data generator
Description:
    Creates a small synthetic genome + GFF3 + candidate gene list for the
    promoter k-mer enrichment screen:
      300 genes along a 500 kb chromosome, 40 candidates whose promoters
      carry an injected 6-mer motif (ACGCGT, an E-box-like element).

        genome.fa            2 chromosomes x 500 kb
        annotation.gff3      gene features
        candidate_genes.csv  40 candidate ids

    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_tf.py --outdir demo_data/tf
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a synthetic promoter-motif demo dataset.")
    parser.add_argument("--outdir", default="demo_data/tf")
    parser.add_argument("--n-genes", type=int, default=300)
    parser.add_argument("--n-candidates", type=int, default=40)
    parser.add_argument("--motif", default="ACGCGT")
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)

    rng = random.Random(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seqs = {"Chr1": "".join(rng.choice("ACGT") for _ in range(500_000)),
            "Chr2": "".join(rng.choice("ACGT") for _ in range(500_000))}
    genes, gids = [], []
    for i in range(args.n_genes):
        start = 5000 + i * 1500
        gid = f"OsG{i + 1:04d}"
        strand = "+" if i % 2 == 0 else "-"
        genes.append((gid, "Chr1", start, start + 1000, strand))
        gids.append(gid)

    cand = gids[: args.n_candidates]
    for gid in cand:
        gi = gids.index(gid)
        gid_, seqid, start, end, strand = genes[gi]
        pos = start - rng.randint(100, 1900)      # upstream position (plus-strand coords)
        idx = pos - 1
        if strand == "-":
            # minus-strand promoter sits downstream of the gene end
            pos = end + rng.randint(100, 1900)
            idx = pos - 1
        if 0 <= idx < len(seqs["Chr1"]) - len(args.motif):
            inject = args.motif if strand == "+" else args.motif.translate(
                str.maketrans("ACGT", "TGCA"))[::-1]
            seqs["Chr1"] = seqs["Chr1"][:idx] + inject + seqs["Chr1"][idx + len(args.motif):]

    with open(outdir / "genome.fa", "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
    with open(outdir / "annotation.gff3", "w") as fh:
        fh.write("##gff-version 3\n")
        for gid, seqid, start, end, strand in genes:
            fh.write(f"{seqid}\tDEMO\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                     f"ID={gid};Name={gid}\n")
    pd_cand = pd.DataFrame({"gene_id": cand})
    pd_cand.to_csv(outdir / "candidate_genes.csv", index=False)

    print(f"[make_demo_tf] wrote {outdir}/ :")
    print(f"  genome.fa (2 x 500 kb), annotation.gff3 ({args.n_genes} genes),")
    print(f"  candidate_genes.csv ({args.n_candidates} candidates, motif '{args.motif}' injected)")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0


if __name__ == "__main__":
    import pandas as pd  # deferred import
    sys.exit(main())
