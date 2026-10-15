"""01_ka_ks.py -- pairwise KA/KS (dN/dS) for orthologous CDS pairs.

Title:    Rice-skills gene-family evolution step
Description:
    Computes non-synonymous (Ka/dN) and synonymous (Ks/dS) substitution rates
    for a set of orthologous CDS pairs:
      1. protein-level pairwise alignment (Needleman-Wunsch, BLOSUM62) per
         pair — or an external protein alignment (aligned FASTA) can be
         supplied with --aligned to skip this step;
      2. codon-aware back-translation of the alignment (pal2nal logic);
      3. Nei-Gojobori (1986) counting with Jukes-Cantor correction:
           Ka = JC(dN_sites), Ks = JC(dS_sites), reported with Ka/Ks.
    Interpretation: Ka/Ks < 1 purifying selection; ~1 neutral; > 1 positive
    selection. Standard gene-family evolution summary.

Usage:
    python scripts/evolution/01_ka_ks.py \
        --pairs demo_data/evolution/ortholog_pairs.csv \
        --cds demo_data/evolution/cds.fasta \
        --outdir results/evolution/01_ka_ks

    # with a pre-aligned protein FASTA (pal2nal-style, aligned ids = pair ids):
    python scripts/evolution/01_ka_ks.py --pairs pairs.csv --cds cds.fasta \
        --protein aligned_proteins.fasta --outdir results/evolution/01_ka_ks

Inputs:
    --pairs    CSV with columns gene1, gene2 (ortholog/paralog pair ids).
    --cds      FASTA of CDS sequences (ids = gene ids; length % 3 == 0).
    --protein  optional FASTA of ALIGNED protein sequences (same ids);
               when given, the CDS are mapped through it (recommended for
               divergent pairs; otherwise NW alignment is computed here).
Outputs:
    ka_ks.csv        gene1, gene2, S, N, Sd, Nd, Ka, Ks, Ka_Ks
    ka_ks_scatter.*  Ka vs Ks scatter with Ka=Ks and Ka=0.1*Ks guides
    ks_histogram.*   Ks distribution (synonymous divergence ~ age proxy)
    run_info.txt

Notes:
    - NG counting with Jukes-Cantor is the fast classic estimator (biased
      upward for distant pairs, Ks > ~1); for publication-grade estimates on
      divergent families use PAML (codeml, runmode=-2) — the pairs file this
      script consumes is exactly what codeml needs.
    - Upstream: get ortholog pairs by reciprocal-best-BLAST (RBH) between
      the two proteomes, then align proteins with MAFFT before this script.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

from data_io import set_seed, write_run_info
from plot_theme import save_nature_figure, setup_nature_style, group_colors

B62 = None  # lazy-loaded BLOSUM62 (built below to avoid biopython dependency)


def _blosum62() -> dict:
    """Minimal BLOSUM62 matrix (common pairs; default -4 for the rest)."""
    global B62
    if B62 is not None:
        return B62
    rows = "ARNDCQEGHILKMFPSTWYV"
    # compact standard BLOSUM62 (symmetric)
    vals = {
        ("A", "A"): 4, ("R", "R"): 5, ("N", "N"): 6, ("D", "D"): 6, ("C", "C"): 9,
        ("Q", "Q"): 5, ("E", "E"): 5, ("G", "G"): 6, ("H", "H"): 8, ("I", "I"): 4,
        ("L", "L"): 4, ("K", "K"): 5, ("M", "M"): 5, ("F", "F"): 6, ("P", "P"): 7,
        ("S", "S"): 4, ("T", "T"): 5, ("W", "W"): 11, ("Y", "Y"): 7, ("V", "V"): 4,
        ("A", "R"): -1, ("A", "N"): -2, ("A", "D"): -2, ("A", "C"): 0, ("A", "Q"): -1,
        ("A", "E"): -1, ("A", "G"): 0, ("A", "H"): -2, ("A", "I"): -1, ("A", "L"): -1,
        ("A", "K"): -1, ("A", "M"): -1, ("A", "F"): -2, ("A", "P"): -1, ("A", "S"): 1,
        ("A", "T"): 0, ("A", "W"): -3, ("A", "Y"): -2, ("A", "V"): 0,
        ("R", "N"): 0, ("R", "D"): -2, ("R", "C"): -3, ("R", "Q"): 1, ("R", "E"): 0,
        ("R", "G"): -2, ("R", "H"): 0, ("R", "I"): -3, ("R", "L"): -2, ("R", "K"): 2,
        ("R", "M"): -1, ("R", "F"): -3, ("R", "P"): -2, ("R", "S"): -1, ("R", "T"): -1,
        ("R", "W"): -3, ("R", "Y"): -2, ("R", "V"): -2,
        ("N", "D"): 1, ("N", "C"): -3, ("N", "Q"): 0, ("N", "E"): 0, ("N", "G"): 0,
        ("N", "H"): 1, ("N", "I"): -3, ("N", "L"): -3, ("N", "K"): 0, ("N", "M"): -2,
        ("N", "F"): -3, ("N", "P"): -2, ("N", "S"): 1, ("N", "T"): 0, ("N", "W"): -4,
        ("N", "Y"): -2, ("N", "V"): -3,
        ("D", "C"): -3, ("D", "Q"): 0, ("D", "E"): 2, ("D", "G"): -1, ("D", "H"): -1,
        ("D", "I"): -3, ("D", "L"): -4, ("D", "K"): -1, ("D", "M"): -3, ("D", "F"): -3,
        ("D", "P"): -1, ("D", "S"): 0, ("D", "T"): -1, ("D", "W"): -4, ("D", "Y"): -3,
        ("D", "V"): -3,
        ("C", "Q"): -3, ("C", "E"): -4, ("C", "G"): -3, ("C", "H"): -3, ("C", "I"): -1,
        ("C", "L"): -1, ("C", "K"): -3, ("C", "M"): -1, ("C", "F"): -2, ("C", "P"): -3,
        ("C", "S"): -1, ("C", "T"): -1, ("C", "W"): -2, ("C", "Y"): -2, ("C", "V"): -1,
        ("Q", "E"): 2, ("Q", "G"): -2, ("Q", "H"): 0, ("Q", "I"): -3, ("Q", "L"): -2,
        ("Q", "K"): 1, ("Q", "M"): 0, ("Q", "F"): -3, ("Q", "P"): -1, ("Q", "S"): 0,
        ("Q", "T"): -1, ("Q", "W"): -2, ("Q", "Y"): -1, ("Q", "V"): -2,
        ("E", "G"): -2, ("E", "H"): 0, ("E", "I"): -3, ("E", "L"): -3, ("E", "K"): 1,
        ("E", "M"): -2, ("E", "F"): -3, ("E", "P"): 0, ("E", "S"): 0, ("E", "T"): -1,
        ("E", "W"): -3, ("E", "Y"): -2, ("E", "V"): -2,
        ("G", "H"): -2, ("G", "I"): -4, ("G", "L"): -4, ("G", "K"): -2, ("G", "M"): -3,
        ("G", "F"): -3, ("G", "P"): -2, ("G", "S"): 0, ("G", "T"): -2, ("G", "W"): -2,
        ("G", "Y"): -3, ("G", "V"): -3,
        ("H", "I"): -3, ("H", "L"): -3, ("H", "K"): -1, ("H", "M"): -2, ("H", "F"): -1,
        ("H", "P"): -2, ("H", "S"): -1, ("H", "T"): -2, ("H", "W"): -2, ("H", "Y"): 2,
        ("H", "V"): -3,
        ("I", "L"): 2, ("I", "K"): -3, ("I", "M"): 1, ("I", "F"): 0, ("I", "P"): -3,
        ("I", "S"): -2, ("I", "T"): -1, ("I", "W"): -3, ("I", "Y"): -1, ("I", "V"): 3,
        ("L", "K"): -2, ("L", "M"): 2, ("L", "F"): 0, ("L", "P"): -3, ("L", "S"): -2,
        ("L", "T"): -1, ("L", "W"): -2, ("L", "Y"): -1, ("L", "V"): 1,
        ("K", "M"): -1, ("K", "F"): -3, ("K", "P"): -1, ("K", "S"): 0, ("K", "T"): -1,
        ("K", "W"): -3, ("K", "Y"): -2, ("K", "V"): -2,
        ("M", "F"): 0, ("M", "P"): -2, ("M", "S"): -1, ("M", "T"): -1, ("M", "W"): -1,
        ("M", "Y"): -1, ("M", "V"): 1,
        ("F", "P"): -4, ("F", "S"): -2, ("F", "T"): -2, ("F", "W"): 1, ("F", "Y"): 3,
        ("F", "V"): -1,
        ("P", "S"): -1, ("P", "T"): -1, ("P", "W"): -4, ("P", "Y"): -3, ("P", "V"): -2,
        ("S", "T"): 1, ("S", "W"): -4, ("S", "Y"): -2, ("S", "V"): -2,
        ("T", "W"): -5, ("T", "Y"): -2, ("T", "V"): 0,
        ("W", "Y"): 2, ("W", "V"): -3,
        ("Y", "V"): -1,
    }
    B62 = {}
    for (a, b), v in vals.items():
        B62[(a, b)] = v
        B62[(b, a)] = v
    for a in rows:
        for b in rows:
            if (a, b) not in B62:
                B62[(a, b)] = -4
    return B62


def read_fasta(path: str) -> dict:
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
    return seqs


def needleman_wunsch(a: str, b: str) -> tuple[str, str]:
    """Global protein alignment with BLOSUM62, affine-free (linear gap -8)."""
    m, n = len(a), len(b)
    gap = -8
    score = np.zeros((m + 1, n + 1), dtype=np.int32)
    ptr = np.zeros((m + 1, n + 1), dtype=np.int8)  # 1 diag, 2 up, 3 left
    score[1:, 0] = np.arange(1, m + 1) * gap
    ptr[1:, 0] = 2
    score[0, 1:] = np.arange(1, n + 1) * gap
    ptr[0, 1:] = 3
    bl = _blosum62()
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            d = score[i - 1, j - 1] + bl.get((a[i - 1], b[j - 1]), -4)
            u = score[i - 1, j] + gap
            l = score[i, j - 1] + gap
            best = max(d, u, l)
            score[i, j] = best
            ptr[i, j] = 1 if best == d else (2 if best == u else 3)
    i, j = m, n
    aa, bb = [], []
    while i > 0 or j > 0:
        p = ptr[i, j]
        if p == 1:
            aa.append(a[i - 1]); bb.append(b[j - 1]); i -= 1; j -= 1
        elif p == 2:
            aa.append(a[i - 1]); bb.append("-"); i -= 1
        else:
            aa.append("-"); bb.append(b[j - 1]); j -= 1
    return "".join(reversed(aa)), "".join(reversed(bb))


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

BASES = "TCAG"


def syn_potential(codon: str, pos: int) -> int:
    """Number of the 3 alternative bases at `pos` keeping the amino acid."""
    aa0 = CODON_TABLE.get(codon, "X")
    n = 0
    for b in BASES:
        if b == codon[pos]:
            continue
        alt = codon[:pos] + b + codon[pos + 1:]
        if CODON_TABLE.get(alt, "X") == aa0:
            n += 1
    return n


def jukes_cantor(p: float) -> float:
    """Jukes-Cantor corrected distance per site (p = observed divergence)."""
    p = min(max(p, 0.0), 0.74)  # guard the JC domain
    return -0.75 * np.log(max(1e-12, 1 - 4.0 * p / 3.0))


def ka_ks_ng(cds1: str, cds2: str) -> dict | None:
    """Nei-Gojobori counting on a codon-aligned CDS pair."""
    if len(cds1) != len(cds2) or len(cds1) % 3 != 0 or len(cds1) == 0:
        return None
    S = N = Sd = Nd = 0.0
    for i in range(0, len(cds1), 3):
        c1, c2 = cds1[i:i + 3], cds2[i:i + 3]
        if "-" in c1 or "-" in c2:
            continue
        if c1 not in CODON_TABLE or c2 not in CODON_TABLE:
            continue
        # degenerate sites per codon (average of the two)
        s1 = sum(syn_potential(c1, p) for p in range(3)) / 3.0
        s2 = sum(syn_potential(c2, p) for p in range(3)) / 3.0
        S += (s1 + s2) / 2.0
        N += 3.0 - (s1 + s2) / 2.0
        if c1 == c2:
            continue
        diff = [p for p in range(3) if c1[p] != c2[p]]
        w = 1.0 / len(diff)
        for p in diff:
            if CODON_TABLE[c1] == CODON_TABLE[c2]:
                Sd += w  # full-codon difference is synonymous (rare)
            else:
                # classify by the single substitution at this position when
                # the rest matches; multi-step diffs counted as non-syn
                alt = c2[:p] + c1[p] + c2[p + 1:]
                if CODON_TABLE.get(alt, "X") == CODON_TABLE[c2]:
                    Sd += w
                else:
                    Nd += w
    if S <= 0 or N <= 0:
        return None
    ka = jukes_cantor(Nd / N)
    ks = jukes_cantor(Sd / S)
    ka_ks = ka / ks if ks > 0 else float("inf")
    return {"S": round(S, 2), "N": round(N, 2), "Sd": round(Sd, 2),
            "Nd": round(Nd, 2), "Ka": round(ka, 5), "Ks": round(ks, 5),
            "Ka_Ks": round(ka_ks, 4) if ks > 0 else np.nan}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pairwise KA/KS for orthologous CDS pairs.")
    p.add_argument("--pairs", required=True, help="CSV: gene1, gene2")
    p.add_argument("--cds", required=True, help="CDS FASTA")
    p.add_argument("--protein", default=None, help="optional ALIGNED protein FASTA")
    p.add_argument("--outdir", default="results/evolution/01_ka_ks")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    set_seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_nature_style()

    pairs = pd.read_csv(args.pairs)
    if not {"gene1", "gene2"}.issubset(pairs.columns):
        raise SystemExit("--pairs must contain columns gene1, gene2")
    cds = read_fasta(args.cds)
    prot_aligned = read_fasta(args.protein) if args.protein else None

    rows = []
    for _, r in pairs.iterrows():
        g1, g2 = str(r["gene1"]), str(r["gene2"])
        if g1 not in cds or g2 not in cds:
            print(f"  (skip {g1}/{g2}: CDS missing)")
            continue
        c1, c2 = cds[g1], cds[g2]
        if prot_aligned and g1 in prot_aligned and g2 in prot_aligned:
            pa, pb = prot_aligned[g1], prot_aligned[g2]
            c1a, c2a, k = [], [], 0
            for aa1, aa2 in zip(pa, pb):
                if aa1 == "-":
                    c2a.append("---"); continue
                if aa2 == "-":
                    c1a.append("---"); continue
                c1a.append(c1[k * 3:(k + 1) * 3])
                c2a.append(c2[k * 3:(k + 1) * 3])
                k += 1
            aln1, aln2 = "".join(c1a), "".join(c2a)
        else:
            pa, pb = needleman_wunsch(
                "".join(CODON_TABLE.get(c1[i:i + 3], "X") for i in range(0, len(c1), 3)),
                "".join(CODON_TABLE.get(c2[i:i + 3], "X") for i in range(0, len(c2), 3)))
            # map back: rebuild codon alignment through the protein alignment
            c1a, c2a, k1, k2 = [], [], 0, 0
            for aa1, aa2 in zip(pa, pb):
                if aa1 != "-":
                    c1a.append(c1[k1 * 3:(k1 + 1) * 3] if aa1 != "-" else "---")
                    k1 += 1 if aa1 != "-" else 0
                else:
                    c1a.append("---")
                if aa2 != "-":
                    c2a.append(c2[k2 * 3:(k2 + 1) * 3] if aa2 != "-" else "---")
                    k2 += 1 if aa2 != "-" else 0
                else:
                    c2a.append("---")
            aln1, aln2 = "".join(c1a), "".join(c2a)
        res = ka_ks_ng(aln1, aln2)
        if res is None:
            print(f"  (skip {g1}/{g2}: alignment/counting failed)")
            continue
        rows.append({"gene1": g1, "gene2": g2, **res})

    if not rows:
        print("No pairs produced results.")
        return 1
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "ka_ks.csv", index=False)
    print(f"[01_ka_ks] pairs computed: {len(df)}")
    print(df[["gene1", "gene2", "Ka", "Ks", "Ka_Ks"]].head(10).to_string(index=False))

    # figures
    import matplotlib.pyplot as plt
    ok = df[(df["Ks"] > 0) & (df["Ks"] < 2)]
    fig, ax = plt.subplots()
    ax.scatter(ok["Ks"], ok["Ka"], s=14, alpha=0.8, color=group_colors(3)[0])
    xs = np.linspace(0, max(ok["Ks"].max(), 0.5), 50)
    ax.plot(xs, xs, ls="--", lw=0.7, color="gray", label="Ka = Ks (neutral)")
    ax.plot(xs, 0.1 * xs, ls=":", lw=0.7, color="gray", label="Ka = 0.1 Ks")
    ax.set_xlabel("Ks (synonymous substitutions)")
    ax.set_ylabel("Ka (non-synonymous substitutions)")
    ax.legend(frameon=False, fontsize=5)
    ax.set_title("Ka / Ks of orthologous pairs")
    save_nature_figure(fig, outdir / "ka_ks_scatter", width_mm=89, height_mm=80)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(df["Ks"].replace([np.inf, -np.inf], np.nan).dropna(), bins=30,
            color=group_colors(3)[1], edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Ks")
    ax.set_ylabel("Number of pairs")
    ax.set_title("Ks distribution (divergence-time proxy)")
    save_nature_figure(fig, outdir / "ks_histogram", width_mm=89, height_mm=70)
    plt.close(fig)

    write_run_info(outdir, script="01_ka_ks.py",
                   params={"pairs": str(args.pairs), "n_pairs": len(df),
                           "protein_alignment": bool(args.protein), "seed": args.seed})
    print(f"[01_ka_ks] done -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
