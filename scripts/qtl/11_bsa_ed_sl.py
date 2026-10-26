#!/usr/bin/env python3
###############################################################################
# 11_bsa_ed_sl.py -- BSA/QTL-seq ED and SL statistics
#
# 输入 CSV：snp_id, chr, pos, ref1, alt1, ref2, alt2
#   ref/alt = 两个混池在该 SNP 的 reads 数（或 allele depth）
# 计算：
#   SNP index = alt/(ref+alt)
#   ΔSNP-index = idx1 - idx2
#   ED = (idx1 - idx2)^2
#   SL = -log10(Fisher exact p)
# 可选窗口平滑（--window bp）与置换阈值（--permutations）
#
# 用法：
#   python scripts/qtl/11_bsa_ed_sl.py --counts bsa_counts.csv --window 1000000
#   python scripts/qtl/11_bsa_ed_sl.py --demo --outdir results/qtl/11_bsa_demo
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

def make_demo(outdir, n=600, causal=(250, 320)):
    rng = np.random.default_rng(42)
    rows = []
    for i in range(n):
        pos = (i + 1) * 20000
        depth1 = rng.integers(30, 80); depth2 = rng.integers(30, 80)
        p1 = 0.5; p2 = 0.5
        if causal[0] <= i <= causal[1]:
            p1 = 0.75; p2 = 0.30
        alt1 = rng.binomial(depth1, p1); alt2 = rng.binomial(depth2, p2)
        rows.append([f"snp{i:04d}", "chr1", pos, depth1-alt1, alt1, depth2-alt2, alt2])
    df = pd.DataFrame(rows, columns=["snp_id","chr","pos","ref1","alt1","ref2","alt2"])
    df.to_csv(os.path.join(outdir, "demo_bsa_counts.csv"), index=False)
    return os.path.join(outdir, "demo_bsa_counts.csv")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts"); ap.add_argument("--window", type=int, default=0)
    ap.add_argument("--permutations", type=int, default=0)
    ap.add_argument("--outdir", default="results/qtl/11_bsa")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo: a.counts = make_demo(a.outdir)
    if not a.counts: sys.exit("Need --counts or --demo")
    d = pd.read_csv(a.counts)
    need = {"snp_id","chr","pos","ref1","alt1","ref2","alt2"}
    if not need <= set(d.columns): sys.exit(f"counts needs {need}")
    d["depth1"] = d["ref1"] + d["alt1"]; d["depth2"] = d["ref2"] + d["alt2"]
    d["idx1"] = d["alt1"] / d["depth1"].replace(0, np.nan)
    d["idx2"] = d["alt2"] / d["depth2"].replace(0, np.nan)
    d["delta"] = d["idx1"] - d["idx2"]
    d["ED"] = d["delta"] ** 2
    p = []
    for _, r in d.iterrows():
        try:
            _, pv = fisher_exact([[r.alt1, r.ref1], [r.alt2, r.ref2]])
        except Exception:
            pv = np.nan
        p.append(pv)
    d["p"] = p; d["SL"] = -np.log10(d["p"].clip(lower=1e-300))
    d["padj"] = bh(d["p"].fillna(1).to_numpy())
    d = d.sort_values(["chr","pos"]).reset_index(drop=True)
    if a.window and a.window > 0:
        win_id = (d["pos"] // a.window).astype(int)
        d["window"] = d["chr"] + ":" + (win_id * a.window).astype(str) + "-" + ((win_id + 1) * a.window).astype(str)
        w = d.groupby(["chr","window"], as_index=False).agg(
            n_snp=("snp_id","size"), mean_ED=("ED","mean"), max_ED=("ED","max"),
            mean_SL=("SL","mean"), max_SL=("SL","max"),
            start=("pos","min"), end=("pos","max"))
        w.to_csv(os.path.join(a.outdir, "bsa_windows.csv"), index=False)
        top = w.sort_values("max_ED", ascending=False).head(5)
    else:
        top = d.sort_values("ED", ascending=False).head(5)
    d.to_csv(os.path.join(a.outdir, "bsa_snp_stats.csv"), index=False)
    write = [
        f"date={datetime.datetime.now().isoformat(timespec='seconds')}",
        f"snps={len(d)}", f"window={a.window}",
        "ED=(idx1-idx2)^2; SL=-log10(Fisher); BH",
    ]
    open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8").write("\n".join(write) + "\n")
    print(f"[bsa] snps={len(d)}; FDR<0.05={int((d['padj']<0.05).sum())}")
    print(top.to_string(index=False))
    print(f"[bsa] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
