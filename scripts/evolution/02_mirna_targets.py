"""02_mirna_targets.py -- Local miRNA target prediction (psRNATarget-style).

Title:    Rice-skills miRNA target prediction (v3.2, local implementation)
Description:
    Dependency-free local implementation of the psRNATarget-style
    seed + pairing scoring scheme for plant miRNA targets:

      - scan every 21-nt window of each transcript (strand-aware: the
        target site is the reverse complement of the miRNA);
      - pair miRNA position i with window position (21 - i + 1) of the
        alignment (anti-parallel); G:U pairs count as half matches;
      - SEED (miRNA nt 2-8) mismatches are penalised heavily (1.0),
        3' region mismatches lightly (0.25), bulges are not modelled
        (this is a heuristic, NOT the full RNAplex energy model);
      - classification (psRNATarget logic, simplified):
          cleavage   : seed perfect (0 mismatch) AND miRNA nt 10-11 pair
                       perfectly to the target  -> slice the transcript
          inhibition : score <= --max-score but not cleavage-class
                       (translational repression candidate)
      - a single combined score (weighted mismatch sum, lower = better)
        ranks candidate sites; report sites with score <= --max-score.

    This is a screening tool: for publication-grade predictions follow up
    with psRNATarget (Web) or TargetFinder, and degradome/PARE validation
    (see docs/small_rna_splicing_usage.md).

Usage:
    python scripts/evolution/02_mirna_targets.py \
        --mirna mirna.fasta --transcripts demo_data/evolution/cds.fasta \
        --max-score 4.0 --outdir results/evolution/02_mirna_targets
    # self-test (no inputs)
    python scripts/evolution/02_mirna_targets.py --demo --outdir results/evolution/02_mirna_targets_demo

Inputs:
    --mirna       FASTA of mature miRNA sequences (one or more)
    --transcripts FASTA of target transcripts/CDS (spliced sequences)

Outputs:
    mirna_targets.csv   mirna, target, position(1-based on transcript),
                        score, mode (cleavage/inhibition), seed_mm, 3p_mm,
                        gu_pairs, site (21nt target alignment window)
    mirna_summary.txt   counts per miRNA, demo self-test result

Notes:
    - Only A/C/G/U handled; N in target windows disqualifies the site.
    - Score scale is heuristic (not an energy); thresholds are relative.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from pathlib import Path

COMP = str.maketrans("ACGUNacgun", "UGCAAugcuu".replace("u", "U"))


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return [(id, seq)] keeping order; ids stripped at first whitespace."""
    out: list[tuple[str, str]] = []
    name = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        out.append((name, "".join(chunks).upper()))
    return out


def pair_ok(a: str, b: str) -> bool:
    """True if a (miRNA nt, in RNA) pairs with b (target RNA nt)."""
    a, b = a.upper(), b.upper()
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    if b == comp.get(a):
        return True
    gu = {("G", "U"), ("U", "G")}
    return (a, b) in gu


def pair_cost(a: str, b: str) -> float:
    """0 perfect, 0.5 G:U, 1 mismatch."""
    if pair_ok(a, b):
        return 0.0 if (a.upper(), b.upper()) in {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")} else 0.5
    return 1.0


def score_site(mirna_rna: str, window: str, seed_mm_max: int = 2) -> dict | None:
    """Score one 21nt window against the miRNA (anti-parallel).
    mirna_rna position i pairs with window position (21 - i)."""
    if "N" in window:
        return None
    L = 21
    seed_mm = 0
    three_mm = 0.0
    gu = 0
    for i in range(1, L + 1):
        m = mirna_rna[i - 1]
        t = window[L - i]  # anti-parallel pairing
        cost = pair_cost(m, t)
        if cost == 0.5:
            gu += 1
        if 2 <= i <= 8:  # seed
            seed_mm += int(cost == 1.0)
            three_mm += 0.5 * cost if cost > 0 else 0.0
        else:  # 3'/5' supplementary region
            three_mm += 0.25 * cost
    if seed_mm > seed_mm_max:
        return None
    # cleavage requires perfect seed and perfect pairing at miRNA nt 10-11
    cleav_pair = pair_cost(mirna_rna[9], window[L - 10]) == 0 and \
                 pair_cost(mirna_rna[10], window[L - 11]) == 0
    cleavage = (seed_mm == 0) and cleav_pair
    score = seed_mm * 1.0 + three_mm + 0.1 * gu
    return {"score": score, "seed_mm": seed_mm, "three_mm": round(three_mm, 2),
            "gu_pairs": gu, "mode": "cleavage" if cleavage else "inhibition"}


def predict(mirnas: list[tuple[str, str]], targets: list[tuple[str, str]],
            max_score: float) -> list[dict]:
    rows: list[dict] = []
    for mid, mseq in mirnas:
        m_rna = mseq.upper().replace("T", "U")
        if len(m_rna) < 21:
            continue
        core = m_rna[:21]
        for tid, tseq in targets:
            t_u = tseq.upper().replace("T", "U")
            for pos in range(0, len(t_u) - 20):
                window = t_u[pos:pos + 21]
                s = score_site(core, window)
                if s is None or s["score"] > max_score:
                    continue
                rows.append({"mirna": mid, "target": tid,
                             "position": pos + 1, "score": round(s["score"], 2),
                             "mode": s["mode"], "seed_mm": s["seed_mm"],
                             "three_mm": s["three_mm"], "gu_pairs": s["gu_pairs"],
                             "site": window})
    rows.sort(key=lambda r: (r["mirna"], r["score"], r["target"]))
    return rows


# ------------------------------------------------------------------- demo
DEMO_MIRNA = ">osa_demo_mir\nUGAGGAAAGAUGCAAGCUCAG"
# target RNA (mRNA strand): revcomp of miRNA = CUGAGCUUGCAUCUUCCUCA (21nt)
# embedded at position 31 of a random-ish transcript
def demo_transcripts() -> list[tuple[str, str]]:
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    mir = DEMO_MIRNA.split("\n")[1]
    rc = "".join(comp[c] for c in reversed(mir))  # perfect target site
    filler = "ACGU" * 15
    t_perfect = filler[:30] + rc + filler[30 + 21:]          # cleavage class
    # seed-shifted site: one seed mismatch (position 5 of miRNA)
    bad = list(rc)
    mi_pos5 = mir[5 - 1]
    for cand_b in "ACGU":
        if cand_b != bad[21 - 5] and not pair_ok(mi_pos5, cand_b):
            bad[21 - 5] = cand_b  # force a real seed mismatch
            break
    t_seedmm = filler[:30] + "".join(bad) + filler[30 + 21:]
    t_none = "ACGU" * 10                                     # no site
    return [("TGT_perfect", t_perfect), ("TGT_seedmm", t_seedmm), ("TGT_none", t_none)]


# ----------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mirna", default=None)
    ap.add_argument("--transcripts", default=None)
    ap.add_argument("--max-score", type=float, default=4.0)
    ap.add_argument("--outdir", default="results/evolution/02_mirna_targets")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        mirnas = read_fasta_lines(DEMO_MIRNA.split("\n"))
        targets = demo_transcripts()
        source = "demo (synthetic)"
    else:
        if not (args.mirna and args.transcripts):
            ap.error("--mirna and --transcripts required (or --demo)")
        mirnas = read_fasta(Path(args.mirna))
        targets = read_fasta(Path(args.transcripts))
        source = f"{args.mirna} vs {args.transcripts}"

    rows = predict(mirnas, targets, args.max_score)
    out_csv = outdir / "mirna_targets.csv"
    if rows:
        with out_csv.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    else:
        out_csv.write_text("mirna,target,position,score,mode,seed_mm,three_mm,gu_pairs,site\n",
                           encoding="utf-8")

    n_cleav = sum(1 for r in rows if r["mode"] == "cleavage")
    stats = [
        f"source: {source}",
        f"timestamp: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"miRNAs: {len(mirnas)}; transcripts: {len(targets)}; max_score: {args.max_score}",
        f"sites reported: {len(rows)} (cleavage-class: {n_cleav})",
    ]
    (outdir / "mirna_summary.txt").write_text("\n".join(stats) + "\n", encoding="utf-8")
    print("\n".join("[02m] " + s for s in stats))
    print(f"[02m] wrote {out_csv}")

    if args.demo:
        per = {(r["mirna"], r["target"]): r for r in rows}
        p = per.get(("osa_demo_mir", "TGT_perfect"))
        s = per.get(("osa_demo_mir", "TGT_seedmm"))
        assert p is not None and p["mode"] == "cleavage" and p["score"] == 0.0, p
        assert s is not None and s["mode"] == "inhibition" and s["seed_mm"] == 1, s
        assert ("osa_demo_mir", "TGT_none") not in per
        print("[02m] demo self-test PASSED (perfect=cleavage score 0; "
              "seed-mismatch=inhibition; no-site target absent)")
    return 0


def read_fasta_lines(lines: list[str]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    name = None
    chunks: list[str] = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        out.append((name, "".join(chunks).upper()))
    return out


if __name__ == "__main__":
    sys.exit(main())
