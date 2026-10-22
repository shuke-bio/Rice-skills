"""03_qpcr_analysis.py -- qPCR relative quantification (2^-ddCt) + RNA-seq cross-check.

Title:    Rice-skills qPCR quantification (design/validation module 9)
Description:
    Standard wet-lab validation layer for the design/transcriptome workflow:
    turns raw qPCR Ct/Cp values into publication-ready relative expression:

      1. read a long-format Ct table (gene, group, bio_rep, ct); technical
         replicates are averaged within (gene, group, bio_rep);
      2. dCt = Ct(target) - Ct(reference gene, same group/bio_rep)
         (--ref-gene; every sample must provide it);
      3. ddCt = dCt - mean(dCt of the calibrator group) per gene
         (--calibrator, usually the control); relative expression = 2^-ddCt;
      4. per-gene summary (mean / SD / SEM / n) and figures: one bar plot per
         gene (mean +/- SD with individual points, calibrator = 1) plus an
         optional all-genes grid; colour scheme from config/nature_palette.json
         when present;
      5. optional cross-check against transcriptome data (--rna-table
         gene,group,value in FPKM/TPM/counts): values are relativised with the
         same calibrator, then Pearson/Spearman correlation + log2 scatter
         with r/p annotation (the standard "qPCR validates RNA-seq" figure).

    Missing values: "Undetermined"/"NA"/empty are dropped with a warning;
    a sample missing its reference Ct is skipped.

Usage:
    python scripts/design/03_qpcr_analysis.py \
        --ct-table ct_values.csv --ref-gene OsUBQ5 --calibrator Control \
        --outdir results/design/03_qpcr
    # with RNA-seq cross-check
    python scripts/design/03_qpcr_analysis.py --ct-table ct.csv \
        --ref-gene OsUBQ5 --calibrator Control \
        --rna-table rnaseq.csv --rna-format fpkm \
        --outdir results/design/03_qpcr
    # self-test (synthetic wells with known fold changes)
    python scripts/design/03_qpcr_analysis.py --demo --outdir results/design/03_qpcr_demo

Inputs:
    --ct-table   CSV: gene, group, bio_rep, ct   (one row per well; technical
                 replicates share gene/group/bio_rep; the reference gene
                 appears as a normal `gene` row, named by --ref-gene)
    --ref-gene   reference gene name (e.g. OsUBQ5 / OsActin)
    --calibrator group name used as calibrator (relative expression = 1)
    --rna-table  optional CSV: gene, group, value (mean transcriptome value)
    --rna-format label only (fpkm/tpm/counts) — count data are log2(x+1)
                 transformed before correlation scatter

Outputs (--outdir):
    ddct_table.csv             per-bio-replicate ct/dct/ddct/rel_expression
    rel_expression_summary.csv gene, group, n, mean, sd, sem
    qpcr_<gene>.pdf/tiff/png   per-gene bar plot (mean +/- SD, points)
    qpcr_all_genes_grid.pdf    optional grid (auto when <= 24 genes)
    correlation_qpcr_rnaseq.csv + scatter (when --rna-table given)
    run_info.txt

Notes:
    - This is the quantification/visualisation step only; primer design lives
      in 02_primer_design.py, and the Ct data itself comes from the qPCR
      machine export (any vendor; convert to the long format above).
    - Significance testing is intentionally NOT implemented here: for a
      two-group fold-change test use the Welch t / Mann-Whitney on
      2^-ddCt values (see scripts/phenomics/05_power_analysis.R for sizing).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
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

BAD = {"", "na", "nan", "undetermined", "none", "-"}


def load_palette() -> list[str]:
    """Shared Nature-style palette (fallback to a sensible default)."""
    default = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]
    for cand in (Path("config/nature_palette.json"),):
        if cand.exists():
            try:
                data = json.loads(cand.read_text(encoding="utf-8"))
                cols = data.get("colors") or data.get("palette") or []
                if cols and isinstance(cols[0], str):
                    return cols
            except Exception:
                pass
    return default


def read_ct(path: Path):
    """Return {(gene, group, bio_rep): [ct, ...]} from a long-format table."""
    wells: dict[tuple, list[float]] = defaultdict(list)
    with path.open("r", encoding="utf-8-sig", errors="replace") as fh:
        reader = csv.DictReader(fh)
        need = {"gene", "group", "bio_rep", "ct"}
        if not need <= {c.lower() for c in (reader.fieldnames or [])}:
            raise SystemExit(f"ct table needs columns {sorted(need)} (got {reader.fieldnames})")
        cols = {c.lower(): c for c in reader.fieldnames}
        for row in reader:
            gene = (row[cols["gene"]] or "").strip()
            group = (row[cols["group"]] or "").strip()
            rep = (row[cols["bio_rep"]] or "").strip()
            raw = (row[cols["ct"]] or "").strip()
            if not gene or not group or raw.lower() in BAD:
                continue
            try:
                ct = float(raw)
            except ValueError:
                continue
            wells[(gene, group, rep)].append(ct)
    return wells


def compute_ddct(wells, ref_gene: str, calibrator: str):
    """Return (detail rows, summary rows, groups, genes)."""
    genes = sorted({g for (g, _, _) in wells if g != ref_gene})
    groups = sorted({grp for (_, grp, _) in wells if grp})
    if calibrator not in groups:
        raise SystemExit(f"calibrator '{calibrator}' not among groups {groups}")
    if not any(g == ref_gene for (g, _, _) in wells):
        raise SystemExit(f"reference gene '{ref_gene}' not found in ct table")

    # dCt per (gene, group, bio_rep)
    dct = {}
    warned = set()
    for (gene, group, rep), cts in wells.items():
        if gene == ref_gene:
            continue
        ref_cts = wells.get((ref_gene, group, rep))
        if not ref_cts:
            warned.add((group, rep))
            continue
        dct[(gene, group, rep)] = float(np.mean(cts)) - float(np.mean(ref_cts))

    detail, summary = [], []
    for gene in genes:
        calib = [dct[(gene, calibrator, r)] for r in
                 sorted({rep for (g, grp, rep) in dct if g == gene and grp == calibrator})
                 if (gene, calibrator, r) in dct]
        if not calib:
            print(f"[03q] WARNING: no calibrator dCt for {gene}, skipped")
            continue
        cal_mean = float(np.mean(calib))
        for group in groups:
            reps = sorted({rep for (g, grp, rep) in dct if g == gene and grp == group})
            vals = []
            for rep in reps:
                key = (gene, group, rep)
                if key not in dct:
                    continue
                ddct = dct[key] - cal_mean
                rel = 2.0 ** (-ddct)
                vals.append(rel)
                detail.append([gene, group, rep, round(float(np.mean(wells[(gene, group, rep)])), 3),
                               round(float(np.mean(wells[(ref_gene, group, rep)])), 3),
                               round(dct[key], 4), round(ddct, 4), round(rel, 6)])
            if vals:
                m = float(np.mean(vals)); s = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
                sem = s / math.sqrt(len(vals)) if vals else 0.0
                summary.append([gene, group, len(vals), round(m, 6), round(s, 6), round(sem, 6)])
    return detail, summary, groups, [g for g in genes]


def plot_gene(gene, summary, groups, palette, outdir: Path, ymax=None):
    rows = {r[1]: r for r in summary if r[0] == gene}
    xs = [i for i, g in enumerate(groups) if g in rows]
    means = [rows[g][3] for g in groups if g in rows]
    sds = [rows[g][4] for g in groups if g in rows]
    fig, ax = plt.subplots(figsize=(1.2 * len(xs) + 1.6, 3.2))
    ax.bar(xs, means, yerr=sds, capsize=3, width=0.62,
           color=[palette[i % len(palette)] for i in xs],
           edgecolor="black", linewidth=0.6,
           error_kw=dict(elinewidth=0.8, ecolor="black"))
    for i, g in zip(xs, [g for g in groups if g in rows]):
        ax.text(i, -0.04 * (ymax or max(means) or 1), g, ha="center", va="top",
                fontsize=9, rotation=0)
    ax.set_xticks([])
    ax.set_ylabel(f"{gene}\nrelative expression (2^-ddCt)")
    ax.set_ylim(bottom=-0.18 * (ymax or max(means) or 1))
    if ymax:
        ax.set_ylim(top=ymax)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / f"qpcr_{gene}.pdf")
    fig.savefig(outdir / f"qpcr_{gene}.tiff", dpi=300)
    fig.savefig(outdir / f"qpcr_{gene}.png", dpi=300)
    plt.close(fig)


def plot_grid(genes, summary, groups, palette, outdir: Path):
    n = len(genes)
    ncol = min(4, n)
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.8 * nrow), squeeze=False)
    for k, gene in enumerate(genes):
        ax = axes[k // ncol][k % ncol]
        rows = {r[1]: r for r in summary if r[0] == gene}
        xs = [i for i, g in enumerate(groups) if g in rows]
        means = [rows[g][3] for g in groups if g in rows]
        sds = [rows[g][4] for g in groups if g in rows]
        ax.bar(xs, means, yerr=sds, capsize=2, width=0.62,
               color=[palette[i % len(palette)] for i in xs],
               edgecolor="black", linewidth=0.5,
               error_kw=dict(elinewidth=0.7, ecolor="black"))
        ax.set_xticks(xs)
        ax.set_xticklabels([g for g in groups if g in rows], fontsize=7)
        ax.set_title(gene, fontsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.tight_layout()
    fig.savefig(outdir / "qpcr_all_genes_grid.pdf")
    fig.savefig(outdir / "qpcr_all_genes_grid.tiff", dpi=300)
    fig.savefig(outdir / "qpcr_all_genes_grid.png", dpi=300)
    plt.close(fig)


def rnaseq_crosscheck(summary, groups, calibrator, rna_path: Path, rna_format: str, outdir: Path):
    """Relativise transcriptome values (same calibrator) and correlate."""
    rna = {}
    with rna_path.open("r", encoding="utf-8-sig", errors="replace") as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        key = next((cols[k] for k in ("value", "expression", "fpkm", "tpm", "counts") if k in cols), None)
        if not ("gene" in cols and "group" in cols and key):
            raise SystemExit("rna table needs columns gene, group, value")
        for row in reader:
            try:
                rna[(row[cols["gene"]].strip(), row[cols["group"]].strip())] = float(row[key])
            except (ValueError, TypeError):
                continue
    # relativise per gene: calibrator group = 1
    genes = sorted({g for (g, _) in rna})
    rows = []
    for gene in genes:
        cal = rna.get((gene, calibrator))
        if not cal or cal <= 0:
            continue
        for (g, grp), v in rna.items():
            if g != gene:
                continue
            rel = v / cal
            if rna_format == "counts":
                rel = math.log2(v + 1)
            rows.append((gene, grp, rel))
    rna_rel = {(g, grp): v for (g, grp, v) in rows}
    paired = []
    for r in summary:
        gene, grp, n, m = r[0], r[1], r[2], r[3]
        rv = rna_rel.get((gene, grp))
        if rv is None:
            continue
        if rna_format == "counts":
            qv = math.log2(m + 1)
        else:
            qv = m
        paired.append((gene, grp, qv, rv))
    if len(paired) >= 3:
        q = np.array([p[2] for p in paired]); t = np.array([p[3] for p in paired])
        from scipy.stats import pearsonr, spearmanr
        rp, pp = pearsonr(q, t)
        rs, ps = spearmanr(q, t)
        with (outdir / "correlation_qpcr_rnaseq.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["gene", "group", "qpcr_value", "rnaseq_value"])
            for p in paired:
                w.writerow([p[0], p[1], round(p[2], 6), round(p[3], 6)])
            w.writerow([])
            w.writerow(["n", len(paired), "pearson_r", round(rp, 4), "pearson_p", round(pp, 6)])
            w.writerow(["", "", "spearman_rho", round(rs, 4), "spearman_p", round(ps, 6)])
        print(f"[03q] qPCR vs RNA-seq: n={len(paired)} Pearson r={rp:.4f} (p={pp:.2g}); "
              f"Spearman rho={rs:.4f}")
        if HAS_MPL:
            fig, ax = plt.subplots(figsize=(3.6, 3.4))
            ax.scatter(q, t, s=18, color="#4C72B0", edgecolor="black", linewidth=0.4)
            lim = [0, max(q.max(), t.max()) * 1.08]
            ax.plot(lim, lim, ls="--", lw=0.7, color="grey")
            ax.set_xlabel("qPCR value"); ax.set_ylabel("RNA-seq value")
            ax.set_title(f"r = {rp:.3f}, p = {pp:.2g}  (n = {len(paired)})", fontsize=9)
            fig.tight_layout()
            fig.savefig(outdir / "qpcr_vs_rnaseq_scatter.pdf")
            fig.savefig(outdir / "qpcr_vs_rnaseq_scatter.tiff", dpi=300)
            fig.savefig(outdir / "qpcr_vs_rnaseq_scatter.png", dpi=300)
            plt.close(fig)
        return {"n": len(paired), "pearson_r": float(rp), "pearson_p": float(pp),
                "spearman_rho": float(rs), "spearman_p": float(ps)}
    print("[03q] not enough paired points for correlation (need >= 3)")
    return None


# ------------------------------------------------------------------- demo
def make_demo(outdir: Path):
    """Synthetic Ct table with known fold changes:
    GeneA: Treat = 4x Ctrl (ddCt = -2); GeneB: Treat = 2x (ddCt = -1);
    GeneC: Treat = 0.5x (ddCt = +1); ref gene flat at 22."""
    rows = [["gene", "group", "bio_rep", "ct"]]
    effects = {"GeneA": -2.0, "GeneB": -1.0, "GeneC": 1.0}
    for gene, dd in effects.items():
        for rep in (1, 2, 3):
            rows.append([gene, "Ctrl", rep, f"{18.0 + (rep - 2) * 0.15:.2f}"])
            rows.append([gene, "Treat", rep, f"{18.0 + dd + (rep - 2) * 0.15:.2f}"])
    for rep in (1, 2, 3):
        for _ in range(2):  # two technical replicates
            rows.append(["RefGene", "Ctrl", rep, "22.0"])
            rows.append(["RefGene", "Treat", rep, "22.0"])
    for gene, dd in effects.items():
        for rep in (1, 2, 3):
            for _ in range(2):
                rows.append([gene, "Ctrl", rep, f"{18.0 + (rep - 2) * 0.15:.2f}"])
                rows.append([gene, "Treat", rep, f"{18.0 + dd + (rep - 2) * 0.15:.2f}"])
    ct = outdir / "demo_ct.csv"
    with ct.open("w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)
    # RNA-seq values aligned with the same fold changes
    rna = outdir / "demo_rnaseq.csv"
    with rna.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "group", "value"])
        for gene, dd in effects.items():
            w.writerow([gene, "Ctrl", 100.0])
            w.writerow([gene, "Treat", 100.0 * (2 ** (-dd))])
    return ct, rna


# ------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ct-table", default=None)
    ap.add_argument("--ref-gene", dest="ref_gene", default=None)
    ap.add_argument("--calibrator", default=None)
    ap.add_argument("--rna-table", default=None)
    ap.add_argument("--rna-format", dest="rna_format", choices=["fpkm", "tpm", "counts"], default="fpkm")
    ap.add_argument("--outdir", default="results/design/03_qpcr")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        ct_path, rna_path = make_demo(outdir)
        ref, calib = "RefGene", "Ctrl"
    else:
        if not (args.ct_table and args.ref_gene and args.calibrator):
            ap.error("--ct-table --ref-gene --calibrator required (or --demo)")
        ct_path, rna_path = Path(args.ct_table), (Path(args.rna_table) if args.rna_table else None)
        ref, calib = args.ref_gene, args.calibrator

    wells = read_ct(ct_path)
    detail, summary, groups, genes = compute_ddct(wells, ref, calib)
    if not detail:
        raise SystemExit("no computable dCt values (check ref gene / bio_rep pairing)")

    with (outdir / "ddct_table.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "group", "bio_rep", "ct_target_mean", "ct_ref_mean",
                    "dct", "ddct", "rel_expression"])
        w.writerows(detail)
    with (outdir / "rel_expression_summary.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "group", "n", "mean", "sd", "sem"])
        w.writerows(summary)

    if HAS_MPL:
        palette = load_palette()
        ymax = max(r[3] + r[4] for r in summary) * 1.18
        for gene in genes:
            plot_gene(gene, summary, groups, palette, outdir, ymax=ymax)
        if 1 < len(genes) <= 24:
            plot_grid(genes, summary, groups, palette, outdir)

    corr = None
    if rna_path:
        corr = rnaseq_crosscheck(summary, groups, calib, rna_path, args.rna_format, outdir)

    info = [
        f"timestamp: {datetime.now().isoformat(timespec='seconds')}",
        f"reference gene: {ref}; calibrator: {calib}",
        f"genes: {len(genes)}; groups: {groups}; bio-rep rows: {len(detail)}",
        f"figures: {'yes' if HAS_MPL else 'no (matplotlib missing)'}",
    ]
    if corr:
        info.append(f"RNA-seq cross-check: n={corr['n']} r={corr['pearson_r']:.3f}")
    (outdir / "run_info.txt").write_text("\n".join(info) + "\n", encoding="utf-8")
    print("\n".join("[03q] " + s for s in info))

    if args.demo:
        # known fold changes: Treat = 4x / 2x / 0.5x of Ctrl
        got = {(r[0], r[1]): r[3] for r in summary}
        for gene, expect in (("GeneA", 4.0), ("GeneB", 2.0), ("GeneC", 0.5)):
            obs = got.get((gene, "Treat"))
            assert obs is not None and abs(obs - expect) < 0.05, \
                f"{gene} Treat expected {expect}, got {obs}"
            # calibrator group mean is ~1 (individual bio-reps scatter around it)
            assert abs(got.get((gene, "Ctrl"), 0) - 1.0) < 0.01, f"{gene} calibrator != 1"
        assert corr and corr["pearson_r"] > 0.99, f"RNA-seq cross-check too low: {corr}"
        print("[03q] demo self-test PASSED (2^-ddCt = 4x/2x/0.5x exact; calibrator=1; "
              "RNA-seq r > 0.99)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
