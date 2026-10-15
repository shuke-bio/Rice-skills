"""make_demo_evolution.py -- synthetic KA/KS demo dataset generator.

Title:    Rice-skills evolution demo data generator
Description:
    Simulates 30 orthologous CDS pairs from random ancestors: each branch
    accumulates 300 substitutions with ~25% non-synonymous (true dN/dS ~
    0.11), producing:
        cds.fasta            60 CDS sequences (30 pairs)
        ortholog_pairs.csv   gene1, gene2
        ground_truth.csv     realised substitution counts per pair
    IMPORTANT: random synthetic data — for pipeline testing only, not biology.

Usage:
    python scripts/utils/make_demo_evolution.py --outdir demo_data/evolution
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import pandas as pd

CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"
          if a + b + c not in ("TAA", "TAG", "TGA")]

# standard genetic code (condensed; identical to scripts/evolution/01_ka_ks.py)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TGT": "C", "TGC": "C",
    "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
    "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I",
    "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T",
    "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
    "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V",
    "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAT": "D",
    "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a synthetic KA/KS demo dataset.")
    parser.add_argument("--outdir", default="demo_data/evolution")
    parser.add_argument("--n-pairs", type=int, default=30)
    parser.add_argument("--n-codons", type=int, default=400)
    parser.add_argument("--n-events", type=int, default=300)
    parser.add_argument("--ns-frac", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)

    rng = random.Random(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def rand_cds(n_codons: int) -> str:
        seq = ""
        while len(seq) < n_codons * 3:
            c = rng.choice(CODONS)
            if CODON_TABLE.get(c, "*") != "*":
                seq += c
        return seq

    # codon table: module-level CODON_TABLE (mirrors scripts/evolution/01_ka_ks.py)

    def evolve(cds_seq: str, n_events: int, ns_frac: float):
        codons = [cds_seq[i:i + 3] for i in range(0, len(cds_seq), 3)]
        n_syn = n_ns = 0
        for _ in range(n_events):
            ci = rng.randrange(len(codons))
            codon = codons[ci]
            pos = rng.randrange(3)
            alt = [b for b in "TCAG" if b != codon[pos]]
            muts = [codon[:pos] + b + codon[pos + 1:] for b in alt]
            syn = [m for m in muts if CODON_TABLE.get(m, "X") == CODON_TABLE[codon]]
            ns = [m for m in muts if m in CODON_TABLE and CODON_TABLE[m] not in ("*", CODON_TABLE[codon])]
            if ns and rng.random() < ns_frac:
                codons[ci] = rng.choice(ns)
                n_ns += 1
            elif syn:
                codons[ci] = rng.choice(syn)
                n_syn += 1
        return "".join(codons), n_syn, n_ns

    pairs, cds_seqs = [], {}
    for i in range(1, args.n_pairs + 1):
        anc = rand_cds(args.n_codons)
        b1, s1, n1 = evolve(anc, args.n_events, args.ns_frac)
        b2, s2, n2 = evolve(anc, args.n_events, args.ns_frac)
        g1, g2 = f"GENE_A{i:02d}", f"GENE_B{i:02d}"
        cds_seqs[g1], cds_seqs[g2] = b1, b2
        pairs.append({"gene1": g1, "gene2": g2,
                      "true_syn": s1 + s2, "true_ns": n1 + n2})

    with open(outdir / "cds.fasta", "w") as fh:
        for gid, seq in cds_seqs.items():
            fh.write(f">{gid}\n{seq}\n")
    pd.DataFrame([{"gene1": p["gene1"], "gene2": p["gene2"]} for p in pairs]).to_csv(
        outdir / "ortholog_pairs.csv", index=False)
    pd.DataFrame(pairs).to_csv(outdir / "ground_truth.csv", index=False)
    print(f"[make_demo_evolution] wrote {outdir}/ :")
    print(f"  cds.fasta ({len(cds_seqs)} CDS), ortholog_pairs.csv ({len(pairs)} pairs),")
    print(f"  ground_truth.csv (non-synonymous event fraction {args.ns_frac};")
    print("  note: realised dN/dS also depends on the synonymous-site fraction)")
    print("  NOTE: synthetic data for pipeline testing only.")
    return 0




if __name__ == "__main__":
    sys.exit(main())
