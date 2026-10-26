#!/usr/bin/env python3
###############################################################################
# 08_causal_network.py -- gene-metabolite association network (+ optional MR edges)
#
# 输入：
#   --x CSV: sample_id, gene1, gene2, ...（基因表达）
#   --y CSV: sample_id, met1, met2, ...（代谢物）
#   --mr  可选 CSV: gene_id, metabolite_id, beta, p（外部 MR 结果，作为因果边）
# 方法：
#   Pearson 相关 + BH；|r| >= 阈值保留；可选合并 MR 边并标注 evidence。
#
# 说明：相关边不等于因果；MR 边才作为因果证据，需自行提供。
#
# 用法：
#   python scripts/integration/08_causal_network.py --x genes.csv --y mets.csv \
#       --mr mr_edges.csv --outdir results/integration/08_network
#   python scripts/integration/08_causal_network.py --demo --outdir results/integration/08_network_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

def bh(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(ranked, 1)
    return out

def make_demo(outdir, n=120):
    rng = np.random.default_rng(42)
    genes = [f"G{i:02d}" for i in range(20)]
    mets = [f"M{i:02d}" for i in range(15)]
    X = rng.normal(6, 1, (n, len(genes)))
    Y = rng.normal(5, 1, (n, len(mets)))
    # module 1: G0-G3 -> M0-M2; module 2: G10-G12 -> M7-M8
    for i in range(3):
        Y[:, i] += 0.9 * X[:, 0] + rng.normal(0, 0.3, n)
    for i in range(2):
        Y[:, 7+i] += 0.8 * X[:, 10] + rng.normal(0, 0.3, n)
    xdf = pd.DataFrame(X, columns=genes); xdf.insert(0, "sample_id", [f"S{i:03d}" for i in range(n)])
    ydf = pd.DataFrame(Y, columns=mets); ydf.insert(0, "sample_id", xdf["sample_id"])
    xdf.to_csv(os.path.join(outdir, "demo_genes.csv"), index=False)
    ydf.to_csv(os.path.join(outdir, "demo_metabolites.csv"), index=False)
    return os.path.join(outdir, "demo_genes.csv"), os.path.join(outdir, "demo_metabolites.csv")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x"); ap.add_argument("--y"); ap.add_argument("--mr")
    ap.add_argument("--min-r", type=float, default=0.5)
    ap.add_argument("--outdir", default="results/integration/08_network")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.x, a.y = make_demo(a.outdir)
    if not a.x or not a.y: sys.exit("Need --x and --y, or --demo")
    x = pd.read_csv(a.x); y = pd.read_csv(a.y)
    if "sample_id" not in x.columns or "sample_id" not in y.columns: sys.exit("need sample_id")
    d = x.merge(y, on="sample_id", suffixes=("_g", "_m"))
    gcols = [c for c in x.columns if c != "sample_id"]
    mcols = [c for c in y.columns if c != "sample_id"]
    rows = []
    for g in gcols:
        for m in mcols:
            r, p = pearsonr(d[g], d[m])
            if np.isfinite(r):
                rows.append([g, m, float(r), float(p)])
    res = pd.DataFrame(rows, columns=["gene_id", "metabolite_id", "r", "p"])
    res["padj"] = bh(res["p"].to_numpy())
    res["evidence"] = "correlation"
    if a.mr and os.path.exists(a.mr):
        mr = pd.read_csv(a.mr)
        mr = mr.rename(columns={"gene_id": "gene_id", "metabolite_id": "metabolite_id"})
        res = res.merge(mr[["gene_id", "metabolite_id", "beta", "p"]].rename(
            columns={"beta": "mr_beta", "p": "mr_p"}), on=["gene_id", "metabolite_id"], how="outer")
        res["evidence"] = np.where(res["mr_beta"].notna(), "MR", "correlation")
    res = res[(res["r"].abs() >= a.min_r) & (res["padj"] < 0.05) | (res["evidence"] == "MR")]
    res = res.sort_values("r", key=lambda s: s.abs(), ascending=False)
    res.to_csv(os.path.join(a.outdir, "network_edges.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"samples={len(d)}\ngenes={len(gcols)}\nmetabolites={len(mcols)}\n")
        f.write(f"min_r={a.min_r}\nedges={len(res)}\n")
        f.write("note=correlation edges are not causal; MR edges require external estimates\n")
    print(f"[net] genes={len(gcols)} metabolites={len(mcols)} edges={len(res)}")
    print(res.head(6).to_string(index=False))
    print(f"[net] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
