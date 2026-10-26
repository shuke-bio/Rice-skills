#!/usr/bin/env python3
###############################################################################
# 06_ptm_ppi.py -- PTM enrichment + STRING PPI network summary
#
# 输入：
#   --proteins CSV: protein_id, gene_id, modification, site, log2FC, p
#   --string   optional TSV: protein1, protein2, score（STRING links）
# 输出：
#   ptm_enrichment.csv（每种修饰的 Fisher 富集）
#   ppi_edges.csv / ppi_hubs.csv（STRING 过滤后的边与 hub）
#
# 用法：
#   python scripts/proteomics/06_ptm_ppi.py --proteins ptm.csv --string string.tsv \
#       --outdir results/proteomics/06_ptm_ppi
#   python scripts/proteomics/06_ptm_ppi.py --demo --outdir results/proteomics/06_ptm_ppi_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

def bh(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(ranked, 1)
    return out

def make_demo(outdir, n=300, edges=400):
    rng = np.random.default_rng(42)
    mods = ["phospho", "ubiquitin", "acetyl", "methyl", "glyco"]
    rows = []
    for i in range(n):
        mod = mods[i % len(mods)]
        diff = (mod == "phospho" and rng.random() < 0.35) or rng.random() < 0.05
        rows.append([f"P{i:04d}", f"G{i:04d}", mod, f"S{rng.integers(1,900)}",
                     rng.normal(1.8 if diff else 0, 0.5),
                     rng.uniform(1e-6, 0.01) if diff else rng.uniform(0.05, 0.9)])
    df = pd.DataFrame(rows, columns=["protein_id","gene_id","modification","site","log2FC","p"])
    df.to_csv(os.path.join(outdir, "demo_proteins.csv"), index=False)
    e = pd.DataFrame({
        "protein1": [f"P{rng.integers(0,n):04d}" for _ in range(edges)],
        "protein2": [f"P{rng.integers(0,n):04d}" for _ in range(edges)],
        "score": rng.integers(200, 999, edges)})
    e = e[e.protein1 != e.protein2]
    e.to_csv(os.path.join(outdir, "demo_string.tsv"), sep="\t", index=False)
    return (os.path.join(outdir,"demo_proteins.csv"), os.path.join(outdir,"demo_string.tsv"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteins"); ap.add_argument("--string")
    ap.add_argument("--min-score", type=int, default=400)
    ap.add_argument("--outdir", default="results/proteomics/06_ptm_ppi")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.proteins, a.string = make_demo(a.outdir)
    if not a.proteins: sys.exit("Need --proteins or --demo")
    df = pd.read_csv(a.proteins)
    need = {"protein_id","modification","log2FC","p"}
    if not need <= set(df.columns): sys.exit("proteins needs protein_id, modification, log2FC, p")
    df["differential"] = (df["p"] < 0.05) & (df["log2FC"].abs() > 1)
    rows = []
    total_diff = int(df["differential"].sum())
    total_nondiff = int((~df["differential"]).sum())
    for mod, sub in df.groupby("modification"):
        a_ = int(sub["differential"].sum())
        b_ = int(len(sub) - a_)
        c_ = total_diff - a_
        d_ = total_nondiff - b_
        try:
            or_, p_ = fisher_exact([[a_, b_], [c_, d_]])
        except Exception:
            or_, p_ = np.nan, np.nan
        rows.append([mod, a_, b_, c_, d_, or_, p_])
    enr = pd.DataFrame(rows, columns=["modification","mod_diff","mod_nondiff","other_diff","other_nondiff","odds_ratio","p"])
    enr["padj"] = bh(enr["p"].fillna(1).to_numpy())
    enr.to_csv(os.path.join(a.outdir, "ptm_enrichment.csv"), index=False)
    if a.string and os.path.exists(a.string):
        ppi = pd.read_csv(a.string, sep="\t")
        ppi.columns = [c.lower() for c in ppi.columns]
        if {"protein1","protein2","score"} <= set(ppi.columns):
            ppi = ppi[ppi["score"] >= a.min_score]
            ppi.to_csv(os.path.join(a.outdir, "ppi_edges.csv"), index=False)
            deg = pd.concat([ppi["protein1"], ppi["protein2"]]).value_counts().reset_index()
            deg.columns = ["protein_id","degree"]
            deg.head(50).to_csv(os.path.join(a.outdir, "ppi_hubs.csv"), index=False)
        else:
            print("(STRING file missing protein1/protein2/score; skipped PPI)")
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"proteins={len(df)}\nmodifications={df['modification'].nunique()}\n")
        f.write(f"min_string_score={a.min_score}\n")
        f.write("ptm_test=Fisher exact; BH\n")
    print(f"[ptm] proteins={len(df)} modifications={df['modification'].nunique()}")
    print(enr.to_string(index=False))
    print(f"[ptm] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
