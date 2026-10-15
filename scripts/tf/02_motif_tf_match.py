"""02_motif_tf_match.py -- Match enriched promoter k-mers to known TF motifs.

Title:    Rice-skills motif -> TF matching (v3.2 step 02)
Description:
    Takes the enriched k-mers from 01_promoter_kmer.py (or any motif/k-mer
    list) and matches them against a library of known TF position weight
    matrices (PWMs) to propose candidate trans-acting TFs.

    Matching = best log-odds score of the k-mer (both strands, sliding
    window) against each PWM under a 0-order background:
        log-odds(i, base) = log2( (PWM[i, base] + pc) / (bg[base] + pc) )
    summed over motif positions; A PWM hit is reported per k-mer with its
    normalised score (fraction of the PWM's max attainable log-odds), and
    TFs are ranked by how many top k-mers they explain.

    PWM library format: MEME motif-file (letter-probability matrix) — the
    export format of PlantPAN 4.0, Cis-BP, and TOMTOM-compatible databases.
    Download the rice TF motif set from PlantPAN (data used with its own
    terms) and pass --pwm; see docs/skill_routing.md (motif-analysis row).

    NOTE: this is a screening match, not a TOMTOM E-value test. Follow up
    significant k-mers with TOMTOM/PlantPAN for publication-grade claims.

Usage:
    python scripts/tf/02_motif_tf_match.py \
        --kmers results/tf/01_promoter_kmer/kmer_enrichment.csv \
        --pwm data/motifs/oryza_plantpan.meme \
        --top 20 --outdir results/tf/02_motif_match
    # self-test (no inputs)
    python scripts/tf/02_motif_tf_match.py --demo --outdir results/tf/02_motif_match_demo

Inputs:
    --kmers  CSV with a 'kmer' column (01_promoter_kmer.py output), or TXT
             with one sequence per line. Sequences must be DNA (A/C/G/T).
    --pwm    MEME-format PWM library; motif name line should look like
             'MOTIF <accession> <TF name>' (TF name optional).

Outputs:
    motif_matches.csv   kmer, motif_id, tf_name, strand, logodds,
                        norm_score (0-1)
    tf_ranking.csv      tf_name, n_kmers_hit, best_norm_score (ranked)
    motif_match_summary.txt

Notes:
    - Requires numpy. Background = uniform 25% by default (--bg to override).
    - Demo PWMs are synthetic; real rice libraries: PlantPAN 4.0 / Cis-BP.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from pathlib import Path

import numpy as np

BASES = "ACGT"


# ------------------------------------------------------------- MEME parser
def parse_meme(path: Path) -> list[dict]:
    """Parse a MEME motif file into [{id, name, pwm: (w x 4) float array}]."""
    motifs: list[dict] = []
    cur: dict | None = None
    in_matrix = False
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("MOTIF"):
                parts = s.split()
                cur = {"id": parts[1] if len(parts) > 1 else f"motif{len(motifs)+1}",
                       "name": parts[2] if len(parts) > 2 else parts[1] if len(parts) > 1 else "",
                       "rows": []}
                in_matrix = False
            elif s.startswith("letter-probability matrix"):
                in_matrix = True
            elif in_matrix and cur is not None:
                if s.startswith(("URL", " ", "\t")) or (s and (s[0].isdigit() or s[0] == "0")):
                    vals = [float(x) for x in s.split()[:4]]
                    if len(vals) == 4:
                        cur["rows"].append(vals)
                elif not s:
                    if cur["rows"]:
                        motifs.append({"id": cur["id"], "name": cur["name"],
                                       "pwm": np.array(cur["rows"], dtype=float)})
                    cur = None
                    in_matrix = False
    if cur is not None and cur["rows"]:
        motifs.append({"id": cur["id"], "name": cur["name"],
                       "pwm": np.array(cur["rows"], dtype=float)})
    return motifs


def read_kmers(path: Path) -> list[str]:
    p = Path(path)
    if p.suffix.lower() == ".csv":
        df = list(csv.DictReader(p.open("r", encoding="utf-8")))
        if df and "kmer" in df[0]:
            return [str(r["kmer"]).upper() for r in df]
        key = next(iter(df[0])) if df else ""
        return [str(r[key]).upper() for r in df]
    out = []
    for line in p.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            out.append(s.split(",")[0].split("\t")[0].upper())
    return out


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def pwm_logodds(pwm: np.ndarray, bg: np.ndarray, pc: float = 0.1) -> np.ndarray:
    return np.log2((pwm + pc) / (bg + pc)).sum(axis=0) if False else \
        np.log2((pwm + pc) / (bg + pc))  # (w x 4)


def match_kmer(kmer: str, motifs: list[dict], bg: np.ndarray) -> list[dict]:
    """Best match of one k-mer against every motif (both strands, sliding)."""
    out: list[dict] = []
    seqs = {"+": kmer, "-": revcomp(kmer)}
    for m in motifs:
        pwm = m["pwm"]
        if not pwm.size or pwm.shape[0] < 2:
            continue
        lo = pwm_logodds(pwm, bg)          # (w x 4)
        max_possible = lo.max(axis=1).sum()
        best = None
        for strand, seq in seqs.items():
            for off in range(0, max(1, len(seq) - pwm.shape[0] + 1)):
                win = seq[off:off + pwm.shape[0]]
                if len(win) < pwm.shape[0]:
                    continue
                idx = [BASES.find(c) for c in win]
                if -1 in idx:
                    continue
                score = float(lo[np.arange(pwm.shape[0]), idx].sum())
                norm = score / max_possible if max_possible > 0 else 0.0
                if best is None or norm > best["norm_score"]:
                    best = {"motif_id": m["id"], "tf_name": m["name"] or m["id"],
                            "strand": strand, "logodds": round(score, 3),
                            "norm_score": round(norm, 4)}
        if best:
            out.append(best)
    out.sort(key=lambda r: -r["norm_score"])
    return out


# ------------------------------------------------------------------- demo
DEMO_MEME = """MEME version 4

MOTIF DEMO1 OsDREB_like
letter-probability matrix: alength= 4 w= 6 nsites= 100
0.85 0.05 0.05 0.05
0.05 0.05 0.85 0.05
0.05 0.05 0.85 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.85 0.05
0.85 0.05 0.05 0.05

MOTIF DEMO2 OsTBP_like
letter-probability matrix: alength= 4 w= 5 nsites= 100
0.9 0.0333 0.0333 0.0333
0.9 0.0333 0.0333 0.0333
0.0333 0.9 0.0333 0.0333
0.0333 0.0333 0.9 0.0333
0.0333 0.0333 0.0333 0.9

MOTIF DEMO3 OsGT_like
letter-probability matrix: alength= 4 w= 4 nsites= 100
0.0333 0.0333 0.9 0.0333
0.0333 0.9 0.0333 0.0333
0.9 0.0333 0.0333 0.0333
0.0333 0.9 0.0333 0.0333
"""

# DEMO1 consensus AGGAGA, DEMO2 AACGT, DEMO3 GCAC (from the matrices above);
# each consensus is embedded in a 7-mer (sliding window finds it)
DEMO_KMERS = ["AGGAGAA", "TAACGTA", "GGCACGC", "TTACCGC", "GGATCCG", "CATTAGC"]


# ----------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--kmers", default=None)
    ap.add_argument("--pwm", default=None)
    ap.add_argument("--bg", default=None, help="comma A,C,G,T background freqs")
    ap.add_argument("--norm-cut", type=float, default=0.6,
                    help="report matches with norm_score >= this")
    ap.add_argument("--outdir", default="results/tf/02_motif_match")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bg = np.array([0.25, 0.25, 0.25, 0.25]) if not args.bg else \
        np.array([float(x) for x in args.bg.split(",")])

    if args.demo:
        pwm_path = outdir / "demo_motifs.meme"
        pwm_path.write_text(DEMO_MEME, encoding="utf-8")
        kmers = DEMO_KMERS
        source = "demo (synthetic motifs/kmers)"
    else:
        if not (args.kmers and args.pwm):
            ap.error("--kmers and --pwm required (or --demo)")
        pwm_path = Path(args.pwm)
        kmers = read_kmers(Path(args.kmers))
        source = f"{args.kmers} vs {args.pwm}"

    motifs = parse_meme(pwm_path)
    rows: list[dict] = []
    for kmer in kmers:
        best = match_kmer(kmer, motifs, bg)
        if best and best[0]["norm_score"] >= args.norm_cut:
            r = dict(best[0])
            r.update({"kmer": kmer})
            rows.append(r)
    rows.sort(key=lambda r: -r["norm_score"])

    out_csv = outdir / "motif_matches.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=["kmer"] + list(rows[0].keys() - {"kmer"}))
            w.writeheader()
            w.writerows(rows)
        else:
            fh.write("kmer,motif_id,tf_name,strand,logodds,norm_score\n")

    tf_rank: dict[str, dict] = {}
    for r in rows:
        d = tf_rank.setdefault(r["tf_name"], {"tf_name": r["tf_name"],
                                              "n_kmers_hit": 0, "best_norm_score": 0.0})
        d["n_kmers_hit"] += 1
        d["best_norm_score"] = max(d["best_norm_score"], r["norm_score"])
    ranking = sorted(tf_rank.values(), key=lambda d: (-d["n_kmers_hit"], -d["best_norm_score"]))
    with (outdir / "tf_ranking.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["tf_name", "n_kmers_hit", "best_norm_score"])
        w.writeheader()
        w.writerows(ranking)

    stats = [
        f"source: {source}",
        f"timestamp: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"motifs: {len(motifs)}; kmers: {len(kmers)}; norm-cut: {args.norm_cut}",
        f"kmers with a PWM match: {len(rows)}",
        f"top TF: {ranking[0]['tf_name'] if ranking else '-'}",
        "follow up with TOMTOM/PlantPAN for E-values (this is a screening match)",
    ]
    (outdir / "motif_match_summary.txt").write_text("\n".join(stats) + "\n", encoding="utf-8")
    print("\n".join("[02t] " + s for s in stats))
    print(f"[02t] wrote {out_csv}")

    if args.demo:
        by_kmer = {r["kmer"]: r for r in rows}
        # DEMO1 consensus ~ AGCGTA(C) ; DEMO2 TATAA ; DEMO3 GTACG
        assert by_kmer["AGGAGAA"]["motif_id"] == "DEMO1", by_kmer.get("AGGAGAA")
        assert by_kmer["TAACGTA"]["motif_id"] == "DEMO2", by_kmer.get("TAACGTA")
        assert by_kmer["GGCACGC"]["motif_id"] == "DEMO3", by_kmer.get("GGCACGC")
        for noise in ("GGATCCG", "CATTAGC"):
            assert noise not in by_kmer, f"noise k-mer matched: {noise}"
        print("[02t] demo self-test PASSED (3 consensus k-mers -> correct motifs; "
              "2 noise k-mers unmatched)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
