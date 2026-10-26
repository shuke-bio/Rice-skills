"""01_sgrna_design.py -- rice CRISPR sgRNA target design (SpCas9, pYLCRISPR/Cas9-ready).

Title:    Rice-skills sgRNA target-design step
Description:
    Designs SpCas9 sgRNA targets (20 nt protospacer + NGG PAM) for a rice gene
    sequence, scores them with the selection principles of the standard
    pYLCRISPR/Cas9 protocol (Ma et al. 2015; 刘耀光课题组):

      1. editing window prefers the CDS (when a GFF is supplied);
      2. GC content 40-75 % preferred;
      3. reject protospacers with >= 4 consecutive T (Pol III terminator);
      4. reject protospacers pairing >= 7 consecutive bases with the sgRNA
         scaffold tail (G-U wobble allowed);
      5. prefer upstream/mid-gene positions (first two-thirds of the CDS);
      6. flag BsaI sites (GGTCTC / GAGACC) inside protospacer+PAM (they would
         interfere with the Golden-Gate construction);
      7. local specificity check: count near-matches (<= 3 mismatches + NGG)
         elsewhere in the supplied sequence; genome-wide off-target analysis
         should be done with CRISPR-P 2.0 or Cas-OFFinder (see docs).

    It then pairs targets for dual-target knockout (>= --min-gap bp apart,
    by default both in the CDS, similar GC) and generates the pYLCRISPR/Cas9
    construction primers (gRT#+/U6aT1-/U6bT2- and the SP-L1/SP-R verification
    pair) plus the shared primer table.

Usage:
    python scripts/design/01_sgrna_design.py \
        --fasta demo_data/genomics/demo_gene.fa \
        --gff demo_data/genomics/demo_gene.gff3 \
        --outdir results/design/01_sgrna --n-pairs 3

Outputs:
    sgrna_candidates.csv         all candidates with every metric + verdict
    sgrna_recommended_pairs.csv  recommended dual-target pairs
    construction_primers.csv     gRT#+/U6#T#- + SP-L1/SP-R + shared primers
    sgrna_map.{pdf,tiff,png}     gene model with target positions
    run_info.txt
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

from data_io import RiceDataError, SEED, save_table, setup_logging, write_run_info  # noqa: E402
from plot_theme import (  # noqa: E402
    load_palette, nature_palette, save_nature_figure, semantic, setup_nature_style,
)

LOGGER = setup_logging("Rice-skills")

# ---------------------------------------------------------------------------
# Protocol constants (pYLCRISPR/Cas9, Ma et al. Mol Plant 2015; 刘耀光课题组文档)
# ---------------------------------------------------------------------------
SCAFFOLD_TAIL = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
GRT_PLUS_SUFFIX = "gttttagagctagaaat"          # appended after the target in gRT#+
U6A_MINUS_SUFFIX = "Cggcagccaagccagca"          # appended after revcomp(target) in U6aT1-
U6B_MINUS_SUFFIX = "Caacacaagcggcagc"           # appended after revcomp(target) in U6bT2-
SHARED_PRIMERS = [
    ("U-F", "CTCCGTTTTACCTGTGGAATCG", "OsU6 启动子扩增正向"),
    ("gR-R", "CGGAGGAAAATTCCATCCAC", "sgRNA 骨架扩增反向"),
    ("Pps-GGL", "TTCAGAggtctcTctcgACTAGTATGGAATCGGCAGCAAAGG", "第二轮 PCR 位点 B-L"),
    ("Pgs-GG2", "AGCGTGggtctcGtcagggTCCATCCACTCCAAGCTC", "第二轮 PCR 位点 2"),
    ("Pps-GG2", "TTCAGAggtctcTctgacacTGGAATCGGCAGCAAAGG", "第二轮 PCR 位点 2"),
    ("Pgs-GGR", "AGCGTGggtctcGaccgACGCGTATCCATCCACTCCAAGCTC", "第二轮 PCR 位点 B-R"),
    ("PB-L", "GCGCGCgGTctcGCTCGACTAGTATGG", "连接产物扩增正向"),
    ("PB-R", "GCGCGCggtctcTACCGACGCGTATCC", "连接产物扩增反向"),
    ("SP-L1", "GCGGTGTCATCTATGTTACTAG", "菌液 PCR / 验证正向"),
    ("SP-R", "CCGACATAGATGCAATAACTTC", "菌液 PCR / 验证反向"),
    ("F-Casjc", "CTGACGCTAACCTCGACAAG", "载体骨架引物 F"),
    ("R-Nos13", "CCGATCTAGTAACATAGATGACACC", "载体骨架引物 R"),
]
COMPLEMENT = {"A": "U", "T": "U", "U": "A", "G": ["C", "U"], "C": "G"}


def revcomp(seq: str) -> str:
    """Reverse complement of a DNA sequence (keeps upper case)."""
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def read_fasta(path: Path) -> tuple[str, str]:
    """Read a (single-record) FASTA file; returns (name, sequence)."""
    if not path.is_file():
        raise RiceDataError(f"FASTA not found: {path}\nHint: pass a gene/genomic sequence file.")
    name, chunks = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                break  # first record only
            name = line[1:].strip().split()[0]
        elif line.strip():
            chunks.append(line.strip())
    if not chunks:
        raise RiceDataError(f"No sequence found in FASTA: {path}")
    seq = "".join(chunks).upper()
    bad = set(seq) - set("ACGTN")
    if bad:
        raise RiceDataError(
            f"Unexpected characters in FASTA: {sorted(bad)}\n"
            "Hint: expected a plain nucleotide sequence (ACGTN)."
        )
    return (name or "sequence"), seq


def read_gff_cds(path: Path, feature: str = "CDS") -> tuple[list[tuple[int, int]], int]:
    """Return (cds_spans, offset): spans are 1-based inclusive relative to the FASTA.

    The offset is the coordinate of FASTA position 1 in the GFF coordinate
    system minus 1, i.e. spans_relative = (start - offset, end - offset).
    """
    if not path.is_file():
        raise RiceDataError(
            f"GFF not found: {path}\n"
            "Hint: pass the GFF that matches the FASTA coordinate system, or omit --gff."
        )
    spans, min_start = [], None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != feature:
            continue
        s, e = int(f[3]), int(f[4])
        spans.append((s, e))
        min_start = s if min_start is None else min(min_start, s)
    if not spans:
        raise RiceDataError(
            f"No '{feature}' features found in {path}.\n"
            "Hint: use a GFF3 with gene/CDS/exon records (e.g. RAP-DB or MSU annotation)."
        )
    offset = (min_start or 1) - 1
    rel = [(s - offset, e - offset) for s, e in sorted(spans)]
    return rel, offset


def in_spans(pos: int, spans: list[tuple[int, int]]) -> bool:
    """True when 1-based position lies in any span."""
    return any(s <= pos <= e for s, e in spans)


def gc_content(seq: str) -> float:
    """GC fraction in percent."""
    gc = sum(1 for b in seq if b in "GC")
    return 100.0 * gc / max(len(seq), 1)


def max_t_run(seq: str) -> int:
    """Length of the longest consecutive T run."""
    return max((len(m) for m in re.findall("T+", seq)), default=0)


def max_scaffold_pairing(protospacer: str, scaffold: str = SCAFFOLD_TAIL) -> int:
    """Longest consecutive complementary stretch between protospacer and scaffold.

    Base-pairing rules include the G-U wobble (RNA-RNA duplex).
    """
    target = protospacer.replace("T", "U")
    sca = scaffold.replace("T", "U")
    best = 0
    for offset in range(-len(sca) + 1, len(target)):
        run = 0
        for i, b in enumerate(target):
            j = i - offset
            if 0 <= j < len(sca):
                sb = sca[j]
                partners = COMPLEMENT.get(b, [])
                if isinstance(partners, str):
                    partners = [partners]
                if sb in partners:
                    run += 1
                    best = max(best, run)
                else:
                    run = 0
        # small optimisation is unnecessary for 20-nt windows
    return best


def local_offtargets(protospacer: str, seq: str, pam: str = "GG", max_mm: int = 3) -> tuple[int, int]:
    """Count near-match sites (<= max_mm mismatches, with PAM) elsewhere in seq.

    Returns (n_sites, best_mm). Checks both strands; the on-target site itself
    (0 mismatches) is excluded.
    """
    pat = re.compile(f"(?=([ACGT]{{{len(protospacer)}}}{pam}))")
    n_sites, best = 0, 99
    for strand_seq in (seq, revcomp(seq)):
        for m in pat.finditer(strand_seq):
            candidate = m.group(1)[: len(protospacer)]
            mm = sum(1 for a, b in zip(candidate, protospacer) if a != b)
            if mm == 0:
                continue  # on-target (or perfect repeat: counted below as mm=0 elsewhere)
            if mm <= max_mm:
                n_sites += 1
                best = min(best, mm)
    # exact repeats of the protospacer with PAM (other copies) are worth flagging
    exact_extra = sum(
        len(re.findall(f"(?={re.escape(protospacer)}{pam})", s)) for s in (seq, revcomp(seq))
    ) - 1
    if exact_extra > 0:
        n_sites += exact_extra
        best = min(best, 1)
    return n_sites, best


def bsaI_warning(protospacer_pam: str) -> bool:
    """True when a BsaI recognition site (GGTCTC / GAGACC) occurs in the window."""
    return ("GGTCTC" in protospacer_pam) or ("GAGACC" in protospacer_pam)


def find_targets(seq: str, name: str) -> list[dict]:
    """Scan both strands for SpCas9 targets (20 nt + NGG)."""
    n = len(seq)
    targets = []
    pat = re.compile(r"(?=([ACGT]{20}GG))")
    for m in pat.finditer(seq):
        targets.append({
            "protospacer": m.group(1)[:20], "pam": m.group(1)[20:22],
            "strand": "+", "start": m.start() + 1, "end": m.start() + 20,
        })
    rc = revcomp(seq)
    for m in pat.finditer(rc):
        start_rc = m.start()
        # The 20-mer found on the reverse-complement scan IS the guide sequence
        # (as it appears in the sgRNA, 5'->3' on the minus strand) - no further
        # reverse complement must be applied here.
        targets.append({
            "protospacer": m.group(1)[:20], "pam": m.group(1)[20:22],
            "strand": "-", "start": n - start_rc - 20 + 1, "end": n - start_rc,
        })
    for t in targets:
        t["seq_name"] = name
    return targets


def score_target(t: dict, cds_spans: list[tuple[int, int]] | None,
                 seq: str, min_gc: float, max_gc: float, max_poly_t: int) -> dict:
    """Apply the protocol selection principles; returns the target row."""
    ps = t["protospacer"]
    gc = gc_content(ps)
    t_run = max_t_run(ps)
    pair_len = max_scaffold_pairing(ps)
    center = (t["start"] + t["end"]) / 2
    in_cds = bool(cds_spans) and in_spans(int(center), cds_spans) if cds_spans else None
    # relative position inside the CDS (0 = 5' end) when available
    if cds_spans:
        cds_start, cds_end = cds_spans[0][0], cds_spans[-1][1]
        rel = min(max((center - cds_start) / max(cds_end - cds_start, 1), 0.0), 1.0)
    else:
        rel = min(max(center / max(len(seq), 1), 0.0), 1.0)
    off_n, off_best = local_offtargets(ps, seq)
    bsa = bsaI_warning(ps + t["pam"])

    # --- transparent composite score (0-100) --------------------------------
    score = 0.0
    score += max(0.0, 40.0 - 2.0 * abs(gc - 60.0))          # GC near 60 % ideal
    seed_gc = gc_content(ps[14:])                            # PAM-proximal seed
    score += max(0.0, 20.0 - 2.0 * abs(seed_gc - 55.0))      # seed region GC
    score += 30.0 * (1.0 - rel)                              # upstream/mid-gene
    if ps[-1] in "GC":
        score += 5.0                                         # 3' GC clamp
    if t_run >= 4:
        score -= 25.0
    if pair_len >= 7:
        score -= 25.0
    if off_n > 0:
        score -= min(15.0, 5.0 * off_n)
    if bsa:
        score -= 10.0
    score = round(float(np.clip(score, 0.0, 100.0)), 1)

    fails = []
    if not (min_gc <= gc <= max_gc):
        fails.append(f"GC {gc:.0f}% outside {min_gc:.0f}-{max_gc:.0f}%")
    if t_run >= 4:
        fails.append(f"{t_run} consecutive T (Pol III terminator risk)")
    if pair_len >= 7:
        fails.append(f"pairs {pair_len} bp with sgRNA scaffold (>=7)")
    if cds_spans is not None and in_cds is False:
        fails.append("outside CDS (protocol prefers CDS)")
    if off_n > 0 and off_best <= 2:
        fails.append(f"local near-match with {off_best} mismatch(es)")

    return {
        **t, "gc_pct": round(gc, 1), "max_t_run": t_run,
        "scaffold_pair_bp": pair_len, "in_cds": in_cds,
        "cds_rel_pos": round(rel, 3), "local_offtarget_hits": off_n,
        "best_local_mm": None if off_best == 99 else off_best,
        "bsaI_site": bsa, "on_score": score,
        "pass": len(fails) == 0, "reject_reasons": "; ".join(fails),
    }


def pair_targets(cand: pd.DataFrame, min_gap: int, max_gap: int | None,
                 n_pairs: int) -> list[tuple[int, int]]:
    """Greedy dual-target pairing on the passing candidates (best-scored first).

    Targets within 20 bp of an already-used target are treated as consumed, so
    recommended pairs are genuinely distinct (no 1-bp-shifted duplicates).
    """
    ok = cand[cand["pass"]].sort_values("on_score", ascending=False).reset_index(drop=True)
    used: set[int] = set()
    used_pos: list[int] = []

    def near(pos: int, window: int = 20) -> bool:
        return any(abs(pos - p) <= window for p in used_pos)

    pairs: list[tuple[int, int]] = []
    for i in range(len(ok)):
        if i in used:
            continue
        a = ok.iloc[i]
        if near(int(a["start"])):
            continue
        for j in range(i + 1, len(ok)):
            if j in used:
                continue
            b = ok.iloc[j]
            if near(int(b["start"])):
                continue
            gap = abs(int(a["start"]) - int(b["start"]))
            if gap < min_gap:
                continue
            if max_gap is not None and gap > max_gap:
                continue
            if abs(float(a["gc_pct"]) - float(b["gc_pct"])) > 15:
                continue
            pairs.append((i, j))
            used.update({i, j})
            used_pos.extend([int(a["start"]), int(b["start"])])
            if len(pairs) >= n_pairs:
                return pairs
            break
    return pairs


def construction_primers(a: dict, b: dict) -> pd.DataFrame:
    """gRT#+/U6#T#- primers for the two targets (G-start rule applied)."""
    def core(ps: str) -> str:
        # protocol: targets starting with G use the N19 without this G;
        # targets not starting with G use the full N20.
        return ps[1:] if ps.startswith("G") else ps

    c1, c2 = core(a["protospacer"]), core(b["protospacer"])
    rows = [
        ("gRT1+", f"{c1}{GRT_PLUS_SUFFIX}", "靶点1 + sgRNA骨架起始（PCR1a 正向）"),
        ("U6aT1-", f"{revcomp(c1)}{U6A_MINUS_SUFFIX}", "靶点1 互补 + U6a 端（PCR1b 反向）"),
        ("gRT2+", f"{c2}{GRT_PLUS_SUFFIX}", "靶点2 + sgRNA骨架起始（PCR2a 正向）"),
        ("U6bT2-", f"{revcomp(c2)}{U6B_MINUS_SUFFIX}", "靶点2 互补 + U6b 端（PCR2b 反向）"),
    ]
    return pd.DataFrame(rows, columns=["primer", "sequence_5to3", "purpose"])


def plot_map(seq_len: int, exons: list[tuple[int, int]], cand: pd.DataFrame,
             pairs: list[tuple[int, int]], out_stem: str) -> None:
    """Gene model with target arrows (publication figure, PDF+TIFF+PNG)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    palette = load_palette()
    setup_nature_style()
    fig, ax = plt.subplots(figsize=(7.08, 2.0))
    # gene spine
    ax.plot([1, seq_len], [0.5, 0.5], color="#B4B2A9", lw=0.8, zorder=1)
    for (s, e) in exons:
        ax.add_patch(Rectangle((s, 0.36), e - s, 0.28, facecolor="#B5D4F4",
                               edgecolor="#185FA5", lw=0.5, zorder=2))
    # targets
    ok = cand[cand["pass"]]
    others = cand[~cand["pass"]]
    ax.scatter(others["start"], np.full(len(others), 0.08), s=6, color="#C9C4BE",
               zorder=3, label="rejected target")
    for _, row in ok.iterrows():
        col = semantic("up") if row["strand"] == "+" else semantic("down")
        if row["strand"] == "+":
            ax.annotate("", xy=(row["end"], 0.82), xytext=(row["start"], 0.82),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0), zorder=4)
        else:  # 5'->3' on the minus strand points right-to-left in gene coordinates
            ax.annotate("", xy=(row["start"], 0.18), xytext=(row["end"], 0.18),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0), zorder=4)
    if len(ok) > 0:
        ax.scatter([], [], marker=">", color=semantic("up"), label="+ strand target")
        ax.scatter([], [], marker="<", color=semantic("down"), label="- strand target")
    # recommended pairs: shade the region
    for k, (i, j) in enumerate(pairs):
        a, b = ok.iloc[i], ok.iloc[j]
        x0, x1 = sorted([float(a["start"]), float(b["end"])])
        ax.axvspan(x0, x1, color="#E8A33D", alpha=0.12, zorder=0)
        ax.text((x0 + x1) / 2, 1.06, f"pair {k + 1}", ha="center", fontsize=7,
                color="#854F0B")
    ax.set_xlim(0, seq_len * 1.02)
    ax.set_ylim(-0.05, 1.15)
    ax.set_yticks([])
    ax.set_xlabel("Position in sequence (bp)")
    ax.set_title("sgRNA target map")
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.legend(loc="lower right", fontsize=7, ncol=3, frameon=False)
    fig.tight_layout()
    save_nature_figure(fig, out_stem, width_mm=180, height_mm=52)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Design SpCas9 sgRNA targets (pYLCRISPR/Cas9-ready) for a rice gene."
    )
    parser.add_argument("--fasta", required=True, help="Gene/genomic sequence (FASTA)")
    parser.add_argument("--gff", default=None, help="GFF3 with CDS features (optional)")
    parser.add_argument("--outdir", default="results/design/01_sgrna")
    parser.add_argument("--n-pairs", type=int, default=3, help="Recommended dual pairs")
    parser.add_argument("--min-gap", type=int, default=100,
                        help="Minimum distance between the two targets (bp)")
    parser.add_argument("--max-gap", type=int, default=None,
                        help="Maximum distance between the two targets (bp)")
    parser.add_argument("--min-gc", type=float, default=40.0)
    parser.add_argument("--max-gc", type=float, default=80.0)
    parser.add_argument("--genome", default=None,
                        help="Optional genome FASTA for genome-wide off-target search (delegates to 04_offtarget.py)")
    parser.add_argument("--max-mm", type=int, default=3,
                        help="Max mismatches for --genome off-target search")
    args = parser.parse_args(argv)

    name, seq = read_fasta(Path(args.fasta))
    cds_spans = None
    exons = []
    if args.gff:
        cds_spans, offset = read_gff_cds(Path(args.gff))
        try:
            exon_spans, _ = read_gff_cds(Path(args.gff), feature="exon")
            exons = exon_spans
        except RiceDataError:
            exons = cds_spans
        LOGGER.info("CDS features: %d (relative spans %s)", len(cds_spans), cds_spans[:3])

    targets = find_targets(seq, name)
    LOGGER.info("Scanned %d bp: %d candidate protospacers (both strands).", len(seq), len(targets))
    if not targets:
        raise RiceDataError(
            "No NGG-PAM targets found in the supplied sequence.\n"
            "Hint: check that the FASTA is a real gene sequence (and that the PAM is NGG)."
        )

    rows = [score_target(t, cds_spans, seq, args.min_gc, args.max_gc, max_poly_t=3)
            for t in targets]
    cand = pd.DataFrame(rows).sort_values(
        ["pass", "on_score"], ascending=[False, False]).reset_index(drop=True)
    n_pass = int(cand["pass"].sum())
    LOGGER.info("Candidates passing all filters: %d / %d.", n_pass, len(cand))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(cand, outdir / "sgrna_candidates.csv")

    # Optional genome-wide off-target search (delegates to 04_offtarget.py)
    if args.genome:
        top = cand[cand["pass"]].head(50).copy()
        if len(top):
            top["sgrna_id"] = [f"sg{i+1}" for i in range(len(top))]
            top[["sgrna_id", "protospacer"]].rename(
                columns={"protospacer": "sequence"}).to_csv(
                outdir / "offtarget_input.csv", index=False)
            cmd = [sys.executable, str(Path(__file__).with_name("04_offtarget.py")),
                   "--sgrna", str(outdir / "offtarget_input.csv"),
                   "--genome", args.genome, "--outdir", str(outdir / "offtarget"),
                   "--max-mm", str(args.max_mm)]
            LOGGER.info("Genome-wide off-target: %s", " ".join(cmd))
            subprocess.run(cmd, check=True)
            summ_path = outdir / "offtarget" / "offtarget_summary.csv"
            if summ_path.exists():
                summ = pd.read_csv(summ_path)
                summ = summ.rename(columns={"sgrna": "protospacer"})
                cand = cand.merge(
                    summ[["protospacer", "total_hits", "mm0", "mm1", "mm2", "mm3plus"]],
                    on="protospacer", how="left")
                save_table(cand, outdir / "sgrna_candidates.csv")
                LOGGER.info("Merged genome-wide off-target counts into sgrna_candidates.csv")
    else:
        LOGGER.warning("No --genome: only local (in-sequence) off-target counts were computed.")

    pairs = pair_targets(cand, args.min_gap, args.max_gap, args.n_pairs)
    ok = cand[cand["pass"]].reset_index(drop=True)
    pair_rows = []
    for k, (i, j) in enumerate(pairs, start=1):
        a, b = ok.iloc[i], ok.iloc[j]
        pair_rows.append({
            "pair": f"pair{k}",
            "target1": a["protospacer"], "t1_start": int(a["start"]), "t1_strand": a["strand"],
            "t1_score": a["on_score"], "t1_gc": a["gc_pct"],
            "target2": b["protospacer"], "t2_start": int(b["start"]), "t2_strand": b["strand"],
            "t2_score": b["on_score"], "t2_gc": b["gc_pct"],
            "gap_bp": abs(int(a["start"]) - int(b["start"])),
            "combined_score": round(float(a["on_score"]) + float(b["on_score"]), 1),
        })
    pairs_df = pd.DataFrame(pair_rows)
    save_table(pairs_df, outdir / "sgrna_recommended_pairs.csv")

    primer_frames = []
    for k, (i, j) in enumerate(pairs, start=1):
        a, b = ok.iloc[i], ok.iloc[j]
        pf = construction_primers(a.to_dict(), b.to_dict())
        pf.insert(0, "pair", f"pair{k}")
        primer_frames.append(pf)
    if primer_frames:
        prim_df = pd.concat(primer_frames, ignore_index=True)
    else:
        prim_df = pd.DataFrame(columns=["pair", "primer", "sequence_5to3", "purpose"])
    shared = pd.DataFrame(SHARED_PRIMERS, columns=["primer", "sequence_5to3", "purpose"])
    shared.insert(0, "pair", "shared")
    prim_all = pd.concat([prim_df, shared], ignore_index=True)
    save_table(prim_all, outdir / "construction_primers.csv")

    try:
        plot_map(len(seq), exons if exons else [(1, len(seq))], cand, pairs,
                 str(outdir / "sgrna_map"))
    except ImportError:
        LOGGER.warning("matplotlib not available - the target map figure was skipped "
                       "(install matplotlib to enable it). Tables are unaffected.")

    write_run_info(outdir, "01_sgrna_design.py",
                   params={"fasta": args.fasta, "gff": args.gff or "none",
                           "n_candidates": len(cand), "n_pass": n_pass,
                           "n_pairs": len(pairs), "min_gap": args.min_gap,
                           "gc_range": f"{args.min_gc}-{args.max_gc}"},
                   packages=("pandas", "numpy", "matplotlib"),
                   extra=["Selection principles follow the pYLCRISPR/Cas9 protocol "
                          "(Ma et al. 2015); construction-primer adaptors transcribed "
                          "from the official protocol document - verify against the "
                          "latest official version before ordering."])

    print(f"Candidates: {len(cand)} scanned, {n_pass} passed all filters.")
    if len(pairs_df) > 0:
        print(pairs_df.to_string(index=False))
    else:
        print("No dual-target pair satisfied the constraints; try --min-gap smaller "
              "or relax --max-gc/--min-gc.")
    print(f"\nResults in {outdir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
