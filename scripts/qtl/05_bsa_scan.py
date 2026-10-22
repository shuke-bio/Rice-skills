"""05_bsa_scan.py -- Genome-wide BSA / QTL-seq scan (delta SNP-index + G').

Title:    Rice-skills genome-wide BSA scan (QTL-seq upstream of 04)
Description:
    Genome-wide bulk segregant analysis scan from a two-bulk sequencing VCF
    (each bulk = pooled individuals; per-sample AD fields carry the pooled
    allele counts). Completes the BSA chain: 05 scans the whole genome for
    candidate regions, 04_qtl_region_snps.py then zooms into a peak interval
    for per-SNP annotation.

    The pipeline implements the two classic statistics (no external R
    package needed):
      1. per-SNP SNP-index for each bulk:
             index = ALT_reads / (REF_reads + ALT_reads)
         delta(index) = index_high - index_low        (QTL-seq, Takagi 2013)
      2. sliding-window aggregation (--window / --step):
             mean delta-index and a G' surrogate
             G_window = sum_i 2 * (index_high_i - 0.5)**2 scaled by window
         significance of mean delta-index from an analytical simulation:
         under H0 the per-SNP index difference has SD
             sqrt(0.25 * (1/n_high + 1/n_low)),
         and window means are compared against the --simulate-fold
         simulation quantiles at --ci (default 0.99). Windows whose |mean
         delta| exceeds the threshold are merged into candidate regions.
      3. outputs: per-window table, candidate-region table, and a two-panel
         genome-wide figure (delta SNP-index and G'), plus run_info.

Usage:
    python scripts/qtl/05_bsa_scan.py \
        --vcf bulk.vcf --bulk-high HighPool --bulk-low LowPool \
        --window 2000000 --step 100000 --ci 0.99 \
        --outdir results/qtl/05_bsa_scan
    # self-test on synthetic bulks with a planted QTL (no inputs needed)
    python scripts/qtl/05_bsa_scan.py --demo --outdir results/qtl/05_bsa_scan_demo

Inputs:
    --vcf         VCF with two pooled samples (plain or .gz); per-sample AD
                  required (--use-gt falls back to genotype dosage indices)
    --bulk-high / --bulk-low   sample names of the two extreme bulks in the VCF

Outputs (--outdir):
    bsa_windows.csv      chrom, start, end, n_snp, mean_delta, abs_mean_delta,
                         mean_index_high, mean_index_low, Gp, threshold
    bsa_peaks.csv        merged candidate regions above the CI threshold
    bsa_scan.png/pdf     genome-wide delta-index and G' tracks (matplotlib)
    run_info.txt

Notes:
    - This is a screening tool aligned with the QTLseqr / QTL-seq methods;
      for publication-grade analysis cross-check with the QTLseqr R package
      (see docs/tools_navigation.md).
    - Chromosomes are processed independently; window size/step should be
      chosen from marker density (e.g. 2 Mb / 100 kb for resequencing data).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:  # pragma: no cover
    HAS_MPL = False


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


# ----------------------------------------------------------------- core
def read_bulk_snps(vcf_path: Path, high: str, low: str, use_gt: bool):
    """Return per-chromosome sorted arrays of (pos, n_high, n_low, i_high, i_low)."""
    per_chrom: dict[str, list[tuple]] = defaultdict(list)
    with open_maybe_gz(vcf_path) as fh:
        ih = il = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                if high not in samples or low not in samples:
                    raise SystemExit(f"bulk samples not in VCF. available: {samples}")
                ih, il = 9 + samples.index(high), 9 + samples.index(low)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            fmt = p[8].split(":")
            if "AD" not in fmt:
                continue
            ia, ib = fmt.index("AD"), None
            gta = None
            d_high = dict(zip(fmt, p[ih].split(":")))
            d_low = dict(zip(fmt, p[il].split(":")))
            if use_gt or d_high.get("AD", ".") in (".", ""):
                try:
                    gh = [int(x) for x in d_high.get("GT", "./.").replace("|", "/").split("/") if x != "."]
                    gl = [int(x) for x in d_low.get("GT", "./.").replace("|", "/").split("/") if x != "."]
                except ValueError:
                    continue
                if not gh or not gl:
                    continue
                i_h = sum(1 for x in gh if x > 0) / len(gh)
                i_l = sum(1 for x in gl if x > 0) / len(gl)
                n_h = n_l = 2
            else:
                try:
                    ah = [int(x) for x in d_high["AD"].split(",") if x != "."]
                    al = [int(x) for x in d_low["AD"].split(",") if x != "."]
                except ValueError:
                    continue
                if len(ah) < 2 or len(al) < 2:
                    continue
                nh = sum(ah); nl = sum(al)
                if nh < 1 or nl < 1:
                    continue
                i_h = sum(ah[1:]) / nh
                i_l = sum(al[1:]) / nl
                n_h, n_l = nh, nl
            per_chrom[p[0]].append((int(p[1]), n_h, n_l, i_h, i_l))
    for c in per_chrom:
        per_chrom[c].sort(key=lambda r: r[0])
    return per_chrom


def window_scan(records, window: int, step: int):
    """Sliding-window delta-index and G' over one chromosome's SNP records."""
    if not records:
        return []
    pos = np.array([r[0] for r in records])
    nh = np.array([r[1] for r in records], dtype=float)
    nl = np.array([r[2] for r in records], dtype=float)
    ih = np.array([r[3] for r in records], dtype=float)
    il = np.array([r[4] for r in records], dtype=float)
    d = ih - il
    # per-SNP SD of the index difference under H0
    sd = np.sqrt(0.25 * (1.0 / np.maximum(nh, 1) + 1.0 / np.maximum(nl, 1)))
    out = []
    start = int(pos.min())
    end_all = int(pos.max())
    s = start
    while s <= end_all:
        e = s + window
        m = (pos >= s) & (pos < e)
        n = int(m.sum())
        if n >= 1:
            md = float(d[m].mean())
            thr = float(1.96 * sd[m].mean())  # placeholder; CI set below
            out.append({
                "start": s, "end": min(e, end_all + 1) - 1, "n_snp": n,
                "mean_delta": round(md, 4),
                "abs_mean_delta": round(abs(md), 4),
                "mean_index_high": round(float(ih[m].mean()), 4),
                "mean_index_low": round(float(il[m].mean()), 4),
                "Gp": round(float((2.0 * (ih[m] - 0.5) ** 2).sum() / n), 4),
                "_sd_mean": float(sd[m].mean()),
            })
        s += step
    return out


def simulate_threshold(sd_means, ci: float, reps: int, seed: int = 42):
    """CI threshold for the window mean delta from N(0, sd_mean) simulation."""
    rng = np.random.default_rng(seed)
    sds = np.array([w["_sd_mean"] for w in sd_means])
    if sds.size == 0:
        return 0.0
    med = float(np.median(sds))
    draws = np.abs(rng.normal(0.0, med, size=(reps,)))
    return float(np.quantile(draws, ci))


def merge_peaks(windows, thr: float):
    peaks = []
    cur = None
    for w in sorted(windows, key=lambda x: x["start"]):
        if w["abs_mean_delta"] >= thr:
            if cur is None:
                cur = dict(start=w["start"], end=w["end"], n_snp=w["n_snp"],
                           peak_abs=w["abs_mean_delta"], peak_Gp=w["Gp"],
                           max_abs=w["abs_mean_delta"])
            else:
                cur["end"] = w["end"]
                cur["n_snp"] += w["n_snp"]
                cur["max_abs"] = max(cur["max_abs"], w["abs_mean_delta"])
                cur["peak_Gp"] = max(cur["peak_Gp"], w["Gp"])
        else:
            if cur is not None:
                peaks.append(cur)
                cur = None
    if cur is not None:
        peaks.append(cur)
    return peaks


# ----------------------------------------------------------------- demo
def make_demo(outdir: Path, chrom_len=1_200_000, n_per_chrom=240, seed=7):
    """Two bulks x 3 chromosomes; a QTL on Chr2 [400k, 700k] drives
    delta-index = 0.6; elsewhere 0.0. AD depth ~40 per bulk."""
    rng = np.random.default_rng(seed)
    vcf = outdir / "demo_bulk.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHighPool\tLowPool",
    ]
    qtl = {"Chr2": (400_000, 700_000)}
    for chrom in ("Chr1", "Chr2", "Chr3"):
        positions = sorted(rng.choice(np.arange(50_000, chrom_len - 50_000),
                                      size=n_per_chrom, replace=False))
        for pos in positions:
            lo, hi = qtl.get(chrom, (10**12, 10**12))
            in_q = lo <= pos <= hi
            idx_h = 0.5 + (0.6 if in_q and rng.random() < 0.9 else rng.normal(0, 0.06))
            idx_l = 0.5 + (0.0 if in_q else rng.normal(0, 0.06))
            idx_h = min(max(idx_h, 0.02), 0.98)
            idx_l = min(max(idx_l, 0.02), 0.98)
            depth = 40
            ah, al = int(round(depth * (1 - idx_h))), int(round(depth * (1 - idx_l)))
            bh, bl = depth - ah, depth - al
            lines.append(f"{chrom}\t{pos}\t.\tA\tG\t.\t.\t.\tGT:AD\t0/1:{ah},{bh}\t0/1:{al},{bl}")
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return vcf, qtl


# ----------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vcf", default=None)
    ap.add_argument("--bulk-high", dest="high", default=None)
    ap.add_argument("--bulk-low", dest="low", default=None)
    ap.add_argument("--window", type=int, default=200_000)
    ap.add_argument("--step", type=int, default=50_000)
    ap.add_argument("--ci", type=float, default=0.99)
    ap.add_argument("--use-gt", action="store_true", help="use GT dosage when AD missing")
    ap.add_argument("--outdir", default="results/qtl/05_bsa_scan")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        vcf, planted = make_demo(outdir)
        high, low = "HighPool", "LowPool"
        window, step = 100_000, 25_000
    else:
        if not (args.vcf and args.high and args.low):
            ap.error("--vcf --bulk-high --bulk-low required (or --demo)")
        vcf, planted = Path(args.vcf), {}
        high, low = args.high, args.low
        window, step = args.window, args.step

    per_chrom = read_bulk_snps(vcf, high, low, args.use_gt)
    if not per_chrom:
        raise SystemExit("no usable SNP records (need AD fields or --use-gt)")

    all_windows: dict[str, list[dict]] = {}
    flat_sd = []
    for chrom in sorted(per_chrom):
        ws = window_scan(per_chrom[chrom], window, step)
        all_windows[chrom] = ws
        flat_sd.extend(w["_sd_mean"] for w in ws)
    thr = simulate_threshold([{"_sd_mean": s} for s in flat_sd], args.ci, reps=4000)
    print(f"[05b] window {window} bp / step {step}; CI {args.ci} threshold on |delta| = {thr:.4f}")

    rows, peak_rows = [], []
    for chrom in sorted(all_windows):
        for w in all_windows[chrom]:
            rows.append([chrom, w["start"], w["end"], w["n_snp"], w["mean_delta"],
                         w["abs_mean_delta"], w["mean_index_high"], w["mean_index_low"],
                         w["Gp"], round(thr, 4)])
        for pk in merge_peaks(all_windows[chrom], thr):
            peak_rows.append([chrom, pk["start"], pk["end"], pk["n_snp"],
                              round(pk["max_abs"], 4), round(pk["peak_Gp"], 4)])

    with (outdir / "bsa_windows.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["chrom", "start", "end", "n_snp", "mean_delta", "abs_mean_delta",
                    "mean_index_high", "mean_index_low", "Gp", "threshold"])
        w.writerows(rows)
    with (outdir / "bsa_peaks.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["chrom", "start", "end", "n_snp", "max_abs_delta", "peak_Gp"])
        w.writerows(peak_rows)

    if HAS_MPL:
        chroms = sorted(all_windows)
        fig, axes = plt.subplots(2, 1, figsize=(11, 6.4), sharex=False)
        xoff = 0
        ticks, labels = [], []
        for chrom in chroms:
            ws = all_windows[chrom]
            xs = [xoff + (w["start"] + w["end"]) / 2 for w in ws]
            dh = [w["mean_delta"] for w in ws]
            gp = [w["Gp"] for w in ws]
            axes[0].plot(xs, dh, lw=1.1, color="#4C72B0")
            axes[1].plot(xs, gp, lw=1.1, color="#55A868")
            axes[0].axhline(thr, color="#C44E52", ls="--", lw=0.8)
            axes[0].axhline(-thr, color="#C44E52", ls="--", lw=0.8)
            ticks.append(xoff + (max(w["end"] for w in ws)) / 2)
            labels.append(chrom)
            xoff += max(w["end"] for w in ws) + 200_000
        axes[0].set_ylabel("mean ΔSNP-index")
        axes[1].set_ylabel("G' (window)")
        axes[1].set_xlabel("chromosome")
        for ax in axes:
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
        fig.tight_layout()
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(outdir / "bsa_scan.pdf") as pdf:
            pdf.savefig(fig)
        fig.savefig(outdir / "bsa_scan.tiff", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / "bsa_scan.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    info = [
        f"timestamp: {datetime.now().isoformat(timespec='seconds')}",
        f"samples: {high} vs {low}; window {window}; step {step}; CI {args.ci}",
        f"SNPs used: {sum(len(v) for v in per_chrom.values())}; windows: {len(rows)}",
        f"threshold: {thr:.4f}; candidate regions: {len(peak_rows)}",
    ]
    (outdir / "run_info.txt").write_text("\n".join(info) + "\n", encoding="utf-8")
    print("\n".join("[05b] " + s for s in info))
    for pk in peak_rows[:8]:
        print(f"[05b]   peak {pk[0]}:{pk[1]}-{pk[2]} n={pk[3]} |delta|={pk[4]} G'={pk[5]}")

    if args.demo:
        hits = [p for p in peak_rows if p[0] == "Chr2" and p[1] < 700_000 and p[2] > 400_000]
        assert hits, f"planted Chr2 QTL not detected: {peak_rows}"
        off = [p for p in peak_rows if p[0] in ("Chr1", "Chr3")]
        assert not off, f"false-positive peaks on no-QTL chromosomes: {off}"
        print("[05b] demo self-test PASSED (Chr2:400-700k recovered; Chr1/Chr3 clean)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
