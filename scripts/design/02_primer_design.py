"""02_primer_design.py -- PCR / qPCR / validation primer design (primer3 engine).

Title:    Rice-skills primer-design step
Description:
    Designs primer pairs with Primer3 (via primer3-py) for three common rice
    lab workflows:

      qpcr        exon-junction-spanning qPCR primers: the genomic sequence is
                  spliced according to the GFF exon records and pairs whose
                  amplicon crosses at least one exon-exon junction are kept
                  (so genomic DNA contamination does not amplify);
      validation  primers amplifying a CRISPR target site +/- --flank bp for
                  Sanger validation of edited plants;
      generic     primers amplifying any region of interest (or anywhere).

    Every pair is reported with Tm, GC, hairpin/homodimer/heterodimer Tm,
    end stability and the Primer3 penalty score, plus a coverage figure.

Usage:
    python scripts/design/02_primer_design.py --mode qpcr \
        --fasta demo_data/genomics/demo_gene.fa --gff demo_data/genomics/demo_gene.gff3 \
        --outdir results/design/02_primers

    python scripts/design/02_primer_design.py --mode validation \
        --fasta demo_data/genomics/demo_gene.fa \
        --target-start 300 --target-end 330 --flank 100 \
        --outdir results/design/02_primers_validation

Notes:
    - primer3-py is licensed GPL-2.0 and used here strictly as an external
      dependency:  pip install primer3-py
    - For NCBI Primer-BLAST-style genome-wide specificity, BLAST the designed
      sequences against the rice genome (or use the NCBI web service); this
      script reports within-sequence specificity only.
Outputs:
    primers.csv            ranked primer pairs with full quality metrics
    primer_report.txt      human-readable summary
    primer_map.{pdf,svg,tiff}  positions of the designed pairs
    run_info.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

from data_io import RiceDataError, save_table, setup_logging, write_run_info  # noqa: E402

LOGGER = setup_logging("Rice-skills")


def require_primer3():
    """Import primer3-py with a friendly error message."""
    try:
        import primer3  # noqa: PLC0415

        return primer3
    except ImportError as exc:  # noqa: BLE001
        raise RiceDataError(
            "primer3-py is required for primer design but is not installed.\n"
            "Hint: pip install primer3-py\n"
            "      (primer3-py is GPL-2.0; it is used as an external engine only.)"
        ) from exc


def read_fasta(path: Path) -> tuple[str, str]:
    """Read a single-record FASTA; returns (name, sequence)."""
    if not path.is_file():
        raise RiceDataError(f"FASTA not found: {path}")
    name, chunks = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].strip().split()[0]
        elif line.strip():
            chunks.append(line.strip())
    if not chunks:
        raise RiceDataError(f"No sequence found in FASTA: {path}")
    return (name or "sequence"), "".join(chunks).upper()


def read_exons(path: Path) -> list[tuple[int, int]]:
    """Read exon spans (1-based inclusive) from a GFF3."""
    if not path.is_file():
        raise RiceDataError(f"GFF not found: {path}")
    spans = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 9 and f[2] == "exon":
            spans.append((int(f[3]), int(f[4])))
    if not spans:
        raise RiceDataError(
            f"No 'exon' features in {path}.\nHint: use a GFF3 with exon/CDS records."
        )
    return sorted(spans)


def splice(seq: str, exons: list[tuple[int, int]]) -> tuple[str, list[int], list[int]]:
    """Splice exons from the genomic sequence.

    Returns (spliced_seq, junctions, spliced_to_genomic):
      junctions: 1-based positions in the spliced sequence where an exon ENDS
                 (i.e. the last base before the next exon starts);
      spliced_to_genomic: for each spliced position (0-based), the 1-based
                 genomic position (needed to map primer coordinates back).
    """
    parts, mapping, junctions = [], [], []
    pos = 0
    for i, (s, e) in enumerate(exons):
        chunk = seq[s - 1 : e]
        parts.append(chunk)
        mapping.extend(range(s, e + 1))
        pos += len(chunk)
        if i < len(exons) - 1:
            junctions.append(pos)  # 1-based position of this exon's last base
    return "".join(parts), junctions, mapping


def design(primer3, template: str, target: tuple[int, int] | None,
           product_range: tuple[int, int], n_return: int, opt_tm: float) -> pd.DataFrame:
    """Run Primer3 and return the pairs as a DataFrame.

    ``target`` is a 1-based inclusive (start, end) region that the amplicon
    must cover (converted to Primer3's 0-based SEQUENCE_TARGET internally).
    """
    seq_args = {
        "SEQUENCE_ID": "template",
        "SEQUENCE_TEMPLATE": template,
        "PRIMER_OPT_SIZE": 21, "PRIMER_MIN_SIZE": 19, "PRIMER_MAX_SIZE": 24,
        "PRIMER_OPT_TM": opt_tm, "PRIMER_MIN_TM": opt_tm - 3, "PRIMER_MAX_TM": opt_tm + 3,
        "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 62.0,
        "PRIMER_MAX_POLY_X": 4, "PRIMER_GC_CLAMP": 1,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_range)],
        "PRIMER_NUM_RETURN": n_return,
    }
    if target is not None:
        start0 = max(int(target[0]) - 1, 0)
        length = max(int(target[1]) - int(target[0]) + 1, 1)
        seq_args["SEQUENCE_TARGET"] = [start0, length]
    global_args = {
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_SALT_DIVALENT": 1.5,
        "PRIMER_DNTP_CONC": 0.2,
    }
    res = primer3.design_primers(seq_args, global_args)
    if "PRIMER_PAIR_NUM_RETURNED" not in res or res["PRIMER_PAIR_NUM_RETURNED"] == 0:
        detail = res.get("PRIMER_ERROR", "no pair satisfied the constraints")
        LOGGER.warning("Primer3 returned no pair (%s).", detail)
        return pd.DataFrame()
    rows = []
    for i in range(res["PRIMER_PAIR_NUM_RETURNED"]):
        f_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        r_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        f_pos, f_len = res[f"PRIMER_LEFT_{i}"]
        r_pos, r_len = res[f"PRIMER_RIGHT_{i}"]
        rows.append({
            "pair_rank": i + 1,
            "forward_seq": f_seq, "reverse_seq": r_seq,
            "forward_pos": int(f_pos) + 1, "reverse_pos": int(r_pos) + 1,
            "forward_tm": round(float(res[f"PRIMER_LEFT_{i}_TM"]), 1),
            "reverse_tm": round(float(res[f"PRIMER_RIGHT_{i}_TM"]), 1),
            "forward_gc": round(float(res[f"PRIMER_LEFT_{i}_GC_PERCENT"]), 1),
            "reverse_gc": round(float(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"]), 1),
            "product_size": int(res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]),
            "pair_penalty": round(float(res[f"PRIMER_PAIR_{i}_PENALTY"]), 3),
            "f_hairpin_tm": round(float(primer3.calc_hairpin_tm(f_seq)), 1),
            "r_hairpin_tm": round(float(primer3.calc_hairpin_tm(r_seq)), 1),
            "f_homodimer_tm": round(float(primer3.calc_homodimer_tm(f_seq)), 1),
            "r_homodimer_tm": round(float(primer3.calc_homodimer_tm(r_seq)), 1),
            "heterodimer_tm": round(float(primer3.calc_heterodimer_tm(f_seq, r_seq)), 1),
            "f_end_stability": round(float(res[f"PRIMER_LEFT_{i}_END_STABILITY"]), 3),
            "r_end_stability": round(float(res[f"PRIMER_RIGHT_{i}_END_STABILITY"]), 3),
        })
    return pd.DataFrame(rows)


def plot_primer_map(seq_name: str, exons: list[tuple[int, int]], pairs: pd.DataFrame,
                    target: tuple[int, int] | None, out_stem: str,
                    mapping: list[int] | None = None, spliced: bool = False) -> None:
    """Coverage figure: gene model + primer positions + amplicon spans."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    from plot_theme import load_palette, setup_nature_style  # noqa: PLC0415
    from plot_theme import save_nature_figure  # noqa: PLC0415

    palette = load_palette()
    setup_nature_style()
    total = max(e for _, e in exons) if not spliced else max(mapping)
    fig, ax = plt.subplots(figsize=(7.08, 1.9))
    # gene model (genomic view)
    ax.plot([exons[0][0], exons[-1][1]], [0.5, 0.5], color="#B4B2A9", lw=0.8)
    for s, e in exons:
        ax.add_patch(Rectangle((s, 0.42), e - s, 0.16, facecolor="#B5D4F4",
                               edgecolor="#185FA5", lw=0.5))
    if target is not None and not spliced:
        ax.axvspan(target[0], target[1], color="#E8A33D", alpha=0.18)
    # pairs
    for _, row in pairs.iterrows():
        f, r = row["forward_pos"], row["reverse_pos"]
        if spliced and mapping is not None:
            f = mapping[min(max(f - 1, 0), len(mapping) - 1)]
            r = mapping[min(max(r - 1, 0), len(mapping) - 1)]
        y = 0.72 + 0.09 * int(row["pair_rank"])
        ax.annotate("", xy=(f + 18, y), xytext=(f, y),
                    arrowprops=dict(arrowstyle="-|>", color="#185FA5", lw=0.9))
        ax.annotate("", xy=(r - 18, y), xytext=(r, y),
                    arrowprops=dict(arrowstyle="-|>", color="#A32D2D", lw=0.9))
        ax.plot([f, r], [y, y], color="#8A8578", lw=0.5, ls=":")
        ax.text(r + 10, y, f"P{int(row['pair_rank'])} ({row['product_size']} bp)",
                fontsize=6.5, va="center", color="#5F5E5A")
    ax.set_xlim(exons[0][0] - 50, exons[-1][1] + 320)
    ax.set_ylim(0.3, 0.72 + 0.09 * (len(pairs) + 1.2))
    ax.set_yticks([])
    view = "spliced transcript" if spliced else "genomic"
    ax.set_xlabel(f"Position ({view}, bp)")
    ax.set_title(f"Primer coverage - {seq_name}")
    ax.spines[["left", "right", "top"]].set_visible(False)
    fig.tight_layout()
    save_nature_figure(fig, out_stem, width_mm=180, height_mm=48)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="PCR / qPCR / validation primer design (primer3).")
    parser.add_argument("--mode", required=True, choices=["qpcr", "validation", "generic"])
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gff", default=None, help="GFF3 with exons (required for --mode qpcr)")
    parser.add_argument("--target-start", type=int, default=None, help="1-based region start")
    parser.add_argument("--target-end", type=int, default=None, help="1-based region end (inclusive)")
    parser.add_argument("--flank", type=int, default=100, help="validation mode: bp around target")
    parser.add_argument("--outdir", default="results/design/02_primers")
    parser.add_argument("--n-return", type=int, default=5)
    args = parser.parse_args(argv)

    primer3 = require_primer3()
    seq_name, seq = read_fasta(Path(args.fasta))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.mode == "qpcr":
        if not args.gff:
            raise RiceDataError("--mode qpcr requires --gff (exon records for splicing).")
        exons = read_exons(Path(args.gff))
        spliced_seq, junctions, mapping = splice(seq, exons)
        LOGGER.info("Spliced transcript: %d bp, %d junction(s) at %s.",
                    len(spliced_seq), len(junctions), junctions)
        # Design per junction: forcing the amplicon to cover each junction's
        # position guarantees qPCR pairs that span an exon-exon boundary
        # (genomic DNA contamination then cannot amplify).
        frames = []
        for j in junctions:
            dfj = design(primer3, spliced_seq, (max(1, j - 2), j + 3),
                         (80, 200), args.n_return, 60.0)
            if len(dfj) > 0:
                dfj["junction_target"] = j
                frames.append(dfj)
        if not frames:
            LOGGER.warning("No junction-covering pair found; falling back to a "
                           "global design (spanning flag may be False).")
            pairs = design(primer3, spliced_seq, None, (80, 200), args.n_return, 60.0)
        else:
            pairs = pd.concat(frames, ignore_index=True)
            pairs = pairs.drop_duplicates(subset=["forward_seq", "reverse_seq"]).reset_index(drop=True)
        if len(pairs) == 0:
            raise RiceDataError(
                "No qPCR primer pair found.\n"
                "Hint: relax by extending the sequence, or check that the GFF exons match the FASTA."
            )
        # keep pairs whose amplicon (forward 5' end .. reverse 5' end, 1-based)
        # contains at least one exon-exon junction position
        keep = []
        for _, row in pairs.iterrows():
            lo, hi = sorted([int(row["forward_pos"]), int(row["reverse_pos"])])
            if any(lo <= j <= hi for j in junctions):
                keep.append(row["pair_rank"])
        spanning = pairs[pairs["pair_rank"].isin(keep)].copy()
        if len(spanning) == 0:
            LOGGER.warning("No designed pair crosses an exon junction; reporting all pairs "
                           "with spanning_junction = False.")
            pairs["spanning_junction"] = False
            out_pairs = pairs
        else:
            spanning["spanning_junction"] = True
            dropped = pairs[~pairs["pair_rank"].isin(keep)].copy()
            dropped["spanning_junction"] = False
            out_pairs = pd.concat([spanning, dropped], ignore_index=True)
        # map spliced positions back to genomic coordinates
        out_pairs["forward_genomic_pos"] = out_pairs["forward_pos"].apply(
            lambda p: mapping[min(max(p - 1, 0), len(mapping) - 1)])
        out_pairs["reverse_genomic_pos"] = out_pairs["reverse_pos"].apply(
            lambda p: mapping[min(max(p - 1, 0), len(mapping) - 1)])
        n_span = int(out_pairs["spanning_junction"].sum())
        LOGGER.info("qPCR pairs: %d total, %d span an exon junction.", len(out_pairs), n_span)
        try:
            plot_primer_map(seq_name, exons, spanning if len(spanning) else out_pairs,
                            None, str(outdir / "primer_map"), mapping=mapping, spliced=True)
        except ImportError:
            LOGGER.warning("matplotlib not available - the coverage figure was skipped.")
        template_desc = f"spliced transcript ({len(spliced_seq)} bp, {len(exons)} exons)"
    else:
        target = None
        if args.mode == "validation":
            if args.target_start is None or args.target_end is None:
                raise RiceDataError(
                    "--mode validation requires --target-start/--target-end (1-based)."
                )
            lo = max(1, args.target_start - args.flank)
            hi = min(len(seq), args.target_end + args.flank)
            target = (lo, hi)
            span = hi - lo + 1
            product_range = (max(span + 20, 250), span + 300)
            LOGGER.info("Validation amplicon region: %d-%d (%d bp); product range %s.",
                        lo, hi, span, product_range)
        else:
            if args.target_start is not None and args.target_end is not None:
                target = (max(1, args.target_start), min(len(seq), args.target_end))
            product_range = (150, 1000)
        pairs = design(primer3, seq, target, product_range, args.n_return, 58.0)
        if len(pairs) == 0:
            raise RiceDataError(
                "No primer pair found.\n"
                "Hint: enlarge the sequence window or adjust --flank / region."
            )
        exons = [(1, len(seq))]
        try:
            if args.gff:
                exons = read_exons(Path(args.gff))
        except RiceDataError:
            exons = [(1, len(seq))]
        try:
            plot_primer_map(seq_name, exons, pairs, target, str(outdir / "primer_map"))
        except ImportError:
            LOGGER.warning("matplotlib not available - the coverage figure was skipped.")
        template_desc = f"genomic sequence ({len(seq)} bp)"

    save_table(out_pairs if args.mode == "qpcr" else pairs, outdir / "primers.csv")

    report = [
        "Rice-skills primer-design report",
        f"mode     : {args.mode}",
        f"template : {template_desc}",
        f"pairs    : {len(out_pairs) if args.mode == 'qpcr' else len(pairs)}",
    ]
    if args.mode == "qpcr":
        report.append(f"junction-spanning pairs: {int(out_pairs['spanning_junction'].sum())}")
    report.append("")
    show = out_pairs if args.mode == "qpcr" else pairs
    for _, row in show.head(10).iterrows():
        report.append(
            f"P{int(row['pair_rank'])}  F:{row['forward_seq']}  R:{row['reverse_seq']}  "
            f"Tm {row['forward_tm']}/{row['reverse_tm']}  GC {row['forward_gc']}/{row['reverse_gc']}  "
            f"product {row['product_size']} bp  penalty {row['pair_penalty']}"
        )
    report += ["", "Quality guide: hairpin/homodimer Tm below the annealing temperature is "
                   "acceptable; |Tm(F)-Tm(R)| <= 2 C is preferred; heterodimer Tm < product "
                   "annealing - 10 C is comfortable."]
    (outdir / "primer_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    write_run_info(outdir, "02_primer_design.py",
                   params={"mode": args.mode, "fasta": args.fasta, "gff": args.gff or "none",
                           "target_start": args.target_start, "target_end": args.target_end,
                           "flank": args.flank, "n_return": args.n_return},
                   packages=("primer3-py", "pandas", "matplotlib"))
    print("\n".join(report))
    print(f"\nResults in {outdir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
