"""05_interpro_domains.py -- InterProScan TSV parsing + differential protein
domain enrichment.

Title:    Rice-skills proteome domain enrichment (step 05)
Description:
    Parses an InterProScan output TSV and tests protein structural-domain
    sets for over-representation among differential proteins (module 14
    output), the protein-level analogue of GO/KEGG over-representation:

      1. parse InterProScan5 TSV columns:
         1 Protein accession, 2 MD5, 5 Analysis, 6 Signature accession,
         7 Signature description, 12 InterPro accession, 13 InterPro
         description (GO terms in col 14 optional);
      2. build domain -> protein sets (InterPro accession level by default;
         --level signature for e.g. PFAM-specific views);
      3. take the differential protein list (02_diff_protein.R output, or
         any one-column/first-column list);
      4. hypergeometric over-representation per domain (universe = proteins
         with >=1 domain annotation), Benjamini-Hochberg FDR;
      5. write ranked results + optional barplot of top domains.

    All-ALL: if a background (measured) protein list is available pass
    --background to restrict the universe to detected proteins — strongly
    recommended for MS data (missingness is not random).

Usage:
    python scripts/proteomics/05_interpro_domains.py \
        --interpro results/proteomics/interproscan.tsv \
        --de results/proteomics/02_diffpro/de_proteins.csv \
        --outdir results/proteomics/05_domains

    # self-test with synthetic data (no input files needed)
    python scripts/proteomics/05_interpro_domains.py --demo --outdir results/proteomics/05_domains_demo

Inputs:
    --interpro   InterProScan5 TSV (as produced by interproscan.sh -T or the
                 EBI online service; comment lines starting with '#' ignored)
    --de         Differential proteins: CSV/TXT, protein_id in the first
                 column (a header line is auto-detected)
    --background Optional CSV/TXT of all measured protein_ids (first column)

Outputs (written to --outdir):
    domain_enrichment.csv  domain_id, name, n_de, n_universe, fold,
                           p_hyper, p_adj (BH), top terms first
    domain_barplot.*       top-N enriched domains (-log10 FDR)
    run_info.txt           parameters, versions, seed

Notes:
    - Hypergeometric via scipy.stats.hypergeom; domains with n_de < 2 are
      skipped (unstable); BH FDR across tested domains only.
    - InterProScan is GPL; this script only READS its TSV output — no
      InterProScan code is used or bundled.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:  # pragma: no cover
    HAS_MPL = False


# --------------------------------------------------------------- parsing
def parse_interpro_tsv(path: Path, level: str) -> dict[str, set[str]]:
    """domain_id -> set of protein ids.

    level='interpro' uses column 12 (InterPro accession, cross-analysis);
    level='signature' uses column 6 (database-specific signature, e.g. PFAM).
    """
    domains: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 13:
                continue
            prot = c[0].strip()
            if level == "signature":
                dom_id, desc = c[5].strip(), (c[6].strip() if len(c) > 6 else "")
            else:
                dom_id, desc = c[11].strip(), (c[12].strip() if len(c) > 12 else "")
                if dom_id == "-" or not dom_id:
                    # fall back to the analysis signature at InterPro level '-'
                    dom_id, desc = c[5].strip(), (c[6].strip() if len(c) > 6 else "")
            if not dom_id or dom_id == "-":
                continue
            domains.setdefault(dom_id, set()).add(prot)
            # remember description once
            dom_desc.setdefault(dom_id, desc or dom_id)
    return domains


dom_desc: dict[str, str] = {}


def read_id_list(path: Path) -> list[str]:
    """Read one protein id per line; auto-detect a single header line."""
    ids: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        first = True
        for line in fh:
            line = line.strip()
            if not line:
                continue
            head = line.split(",")[0].split("\t")[0].strip()
            if first:
                first = False
                if head.lower().replace("_", "").replace("-", "") in {
                    "proteinid", "protein", "proteinid", "id", "protein_id", "accession",
                }:
                    continue
            ids.append(head)
    return ids


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment."""
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n, dtype=float)
    adj[order] = np.clip(ranked, 0, 1)
    return adj


# ------------------------------------------------------------------- demo
def write_demo_inputs(tmp: Path) -> tuple[Path, Path, Path]:
    """Synthetic InterProScan TSV + DE list with a planted signal:
    domain IPR0DEMO1 is over-represented among DE proteins."""
    rng = np.random.default_rng(42)
    n_all, n_de = 400, 60
    all_ids = [f"OsProt{i:05d}" for i in range(1, n_all + 1)]
    de_ids = [f"OsProt{i:05d}" for i in range(1, n_de + 1)]
    # 4 domains; IPR0DEMO1 planted in 80% of DE and 10% of background
    dom_rows: list[str] = []
    for i, pid in enumerate(all_ids):
        if pid in de_ids:
            hit = "IPR0DEMO1" if rng.random() < 0.8 else rng.choice(["IPR0DEMO2", "IPR0DEMO3"])
        else:
            hit = "IPR0DEMO1" if rng.random() < 0.1 else rng.choice(["IPR0DEMO2", "IPR0DEMO3", "IPR0DEMO4"])
        dom_rows.append("\t".join([pid, f"md5{i}", "Pfam", f"PF{i:05d}", "demo analysis",
                                   hit, f"demo domain {hit[-1]}", "-", "-", "-", "-", hit, f"demo domain {hit[-1]}", "-"]))
    # give a second domain to half the proteins (realistic multi-domain)
    for i, pid in enumerate(all_ids):
        if i % 2 == 0:
            dom_rows.append("\t".join([pid, f"md5{i}", "PANTHER", f"PTHR{i:05d}", "demo analysis",
                                       "IPR0DEMO4", "demo domain 4", "-", "-", "-", "-", "IPR0DEMO4", "demo domain 4", "-"]))
    ipr = tmp / "interproscan_demo.tsv"
    ipr.write_text("\n".join(dom_rows) + "\n", encoding="utf-8")
    de = tmp / "de_proteins_demo.csv"
    de.write_text("protein_id\n" + "\n".join(de_ids) + "\n", encoding="utf-8")
    bg = tmp / "background_demo.csv"
    bg.write_text("protein_id\n" + "\n".join(all_ids) + "\n", encoding="utf-8")
    return ipr, de, bg


# ------------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--interpro", default=None)
    ap.add_argument("--de", default=None)
    ap.add_argument("--background", default=None)
    ap.add_argument("--outdir", default="results/proteomics/05_domains")
    ap.add_argument("--level", choices=["interpro", "signature"], default="interpro")
    ap.add_argument("--min-de", type=int, default=2)
    ap.add_argument("--top-plot", type=int, default=15)
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        ipr_path, de_path, bg_path = write_demo_inputs(outdir)
        args.interpro, args.de, args.background = str(ipr_path), str(de_path), str(bg_path)

    if not (args.interpro and args.de):
        ap.error("--interpro and --de are required (or use --demo)")

    dom_desc.clear()
    domains = parse_interpro_tsv(Path(args.interpro), args.level)
    de_ids = set(read_id_list(Path(args.de)))
    universe_ids = set(read_id_list(Path(args.background))) if args.background else set(domains)
    if not universe_ids:
        universe_ids = set().union(*domains.values()) if domains else set()

    n_universe = len(universe_ids)
    results: list[dict] = []
    for dom_id, prots in domains.items():
        overlap = prots & de_ids
        n_dom = len(prots & universe_ids)
        n_de_hit = len(overlap)
        if n_de_hit < args.min_de or n_dom == 0:
            continue
        # P(X >= n_de_hit) with X ~ Hypergeom(N=n_universe, K=n_dom, n=n_de)
        p = hypergeom.sf(n_de_hit - 1, n_universe, n_dom, len(de_ids))
        results.append({
            "domain_id": dom_id,
            "description": dom_desc.get(dom_id, dom_id),
            "n_de": n_de_hit,
            "n_de_total": len(de_ids),
            "n_universe_domain": n_dom,
            "n_universe": n_universe,
            "fold": (n_de_hit / len(de_ids)) / (n_dom / n_universe) if n_dom else np.nan,
            "p_hyper": float(p),
        })

    if not results:
        print("[05] no domain passed the min-de threshold — nothing to report")
        return 1

    pvals = np.array([r["p_hyper"] for r in results])
    padj = bh_adjust(pvals)
    for r, pa in zip(results, padj):
        r["p_adj"] = float(pa)
    results.sort(key=lambda r: (r["p_adj"], r["p_hyper"]))

    out_csv = outdir / "domain_enrichment.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)

    n_sig = sum(1 for r in results if r["p_adj"] < 0.05)
    stats = [
        f"timestamp: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"level: {args.level}; universe: {n_universe} proteins; DE: {len(de_ids)}",
        f"domains tested (n_de >= {args.min_de}): {len(results)}",
        f"significant at BH FDR < 0.05: {n_sig}",
        f"top: {results[0]['domain_id']} ({results[0]['description']}) "
        f"fold={results[0]['fold']:.2f} p_adj={results[0]['p_adj']:.3g}",
    ]
    (outdir / "run_info.txt").write_text("\n".join(stats) + "\n", encoding="utf-8")
    print("\n".join("[05] " + s for s in stats))
    print(f"[05] wrote {out_csv}")

    # demo self-test: planted domain must rank first and be significant
    if args.demo:
        top = results[0]
        assert top["domain_id"] == "IPR0DEMO1", f"planted domain not top: {top}"
        assert top["fold"] > 2 and top["p_adj"] < 0.01, f"signal too weak: {top}"
        print("[05] demo self-test PASSED (planted domain top-ranked, significant)")

    if HAS_MPL:
        import sys as _s
        _s.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
        top = results[: args.top_plot]
        labels = [f"{r['domain_id']}\n{r['description'][:28]}" for r in top][::-1]
        vals = [-np.log10(max(r["p_adj"], 1e-300)) for r in top][::-1]
        fig, ax = plt.subplots(figsize=(6.0, 0.4 * len(top) + 1.2))
        ax.barh(labels, vals, color="#4C72B0")
        ax.axvline(-np.log10(0.05), color="#C44E52", ls="--", lw=0.8)
        ax.set_xlabel("-log10(BH FDR)")
        ax.set_ylabel("")
        ax.tick_params(labelsize=7)
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"domain_barplot.{ext}", bbox_inches="tight")
        plt.close(fig)
        print(f"[05] wrote {outdir / 'domain_barplot.pdf'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
