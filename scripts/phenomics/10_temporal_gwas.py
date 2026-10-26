#!/usr/bin/env python3
###############################################################################
# 10_temporal_gwas.py -- time-resolved GWAS (per time point + growth rate)
#
# 输入：
#   --genotype CSV: sample_id, snp1, ...（0/1/2）
#   --phenotype CSV 长表: sample_id, time, trait
#   --covariates 可选
# 方法：
#   每个时间点做 PCA 校正 OLS；再对每个样本拟合 trait~time 斜率，做斜率 GWAS。
#
# 用法：
#   python scripts/phenomics/10_temporal_gwas.py --genotype geno.csv \
#       --phenotype pheno_long.csv --outdir results/phenomics/10_temporal
#   python scripts/phenomics/10_temporal_gwas.py --demo --outdir results/phenomics/10_temporal_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

def bh(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(ranked, 1)
    return out

def scan(y, G, Z):
    m = G.shape[1]
    beta = np.full(m, np.nan); se = np.full(m, np.nan); p = np.full(m, np.nan)
    for j in range(m):
        D = np.column_stack([Z, G[:, j]])
        coef, *_ = np.linalg.lstsq(D, y, rcond=None)
        df = D.shape[0] - D.shape[1]
        rss = float(((y - D @ coef) ** 2).sum())
        try:
            inv = np.linalg.inv(D.T @ D)
        except np.linalg.LinAlgError:
            continue
        se[j] = np.sqrt(inv[-1, -1] * rss / df)
        beta[j] = coef[-1]
        if se[j] > 0: p[j] = 2 * stats.t.sf(abs(beta[j] / se[j]), df)
    return beta, se, p

def make_demo(outdir, n=80, m=200, times=(0, 7, 14, 21, 28)):
    rng = np.random.default_rng(42)
    G = rng.binomial(2, 0.3, size=(n, m)).astype(float)
    maf = np.minimum(G.mean(0)/2, 1-G.mean(0)/2)
    G = G[:, maf >= 0.05]
    samples = [f"S{i:03d}" for i in range(n)]
    causal = 20
    rows = []
    for s, sample in enumerate(samples):
        base = rng.normal(0, 1)
        for t in times:
            eff = 0.05 * t  # effect grows with time
            y = base + eff * G[s, causal] + rng.normal(0, 0.5)
            rows.append([sample, t, y])
    geno = pd.DataFrame(G, columns=[f"snp{j:04d}" for j in range(G.shape[1])])
    geno.insert(0, "sample_id", samples)
    pheno = pd.DataFrame(rows, columns=["sample_id", "time", "trait"])
    geno.to_csv(os.path.join(outdir, "demo_genotype.csv"), index=False)
    pheno.to_csv(os.path.join(outdir, "demo_phenotype_long.csv"), index=False)
    return (os.path.join(outdir, "demo_genotype.csv"), os.path.join(outdir, "demo_phenotype_long.csv"),
            f"snp{causal:04d}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genotype"); ap.add_argument("--phenotype"); ap.add_argument("--covariates")
    ap.add_argument("--npc", type=int, default=3)
    ap.add_argument("--outdir", default="results/phenomics/10_temporal")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    truth = None
    if a.demo:
        a.genotype, a.phenotype, truth = make_demo(a.outdir)
    if not a.genotype or not a.phenotype: sys.exit("Need --genotype and --phenotype, or --demo")
    g = pd.read_csv(a.genotype); ph = pd.read_csv(a.phenotype)
    if "sample_id" not in g.columns or not {"sample_id","time","trait"} <= set(ph.columns):
        sys.exit("genotype needs sample_id; phenotype needs sample_id,time,trait")
    snps = [c for c in g.columns if c != "sample_id"]
    d = ph.merge(g, on="sample_id", how="inner")
    if d.empty: sys.exit("No shared samples")
    G = d[snps].to_numpy(float)
    colmean = np.nanmean(G, 0); inds = np.where(np.isnan(G)); G[inds] = np.take(colmean, inds[1])
    maf = np.minimum(G.mean(0)/2, 1-G.mean(0)/2); keep = maf >= 0.05
    G = G[:, keep]; snps = [s for s, k in zip(snps, keep) if k]
    Z = np.ones((len(d), 1))
    if a.npc > 0:
        pcs = PCA(n_components=min(a.npc, G.shape[0]-1, G.shape[1]), random_state=0).fit_transform(G)
        Z = np.column_stack([Z, pcs])
    rows = []
    for t, sub in d.groupby("time"):
        idx = sub.index.to_numpy()
        # d is merged with default index; groupby retains index into d
        y = d.loc[idx, "trait"].to_numpy(float)
        b, s, p = scan(y, G[idx], Z[idx])
        for j, snp in enumerate(snps):
            if np.isfinite(p[j]): rows.append([snp, t, b[j], s[j], p[j]])
    res = pd.DataFrame(rows, columns=["snp_id", "time", "beta", "se", "p"])
    res["padj"] = bh(res["p"].to_numpy())
    res.to_csv(os.path.join(a.outdir, "temporal_gwas.csv"), index=False)
    # growth rate per sample
    slopes = []
    for sid, sub in ph.groupby("sample_id"):
        if sub["time"].nunique() >= 2:
            sl = np.polyfit(sub["time"], sub["trait"], 1)[0]
            slopes.append([sid, sl])
    sdf = pd.DataFrame(slopes, columns=["sample_id", "slope"]).merge(g, on="sample_id", how="inner")
    Gs = sdf[snps].to_numpy(float)
    colmean = np.nanmean(Gs, 0); inds = np.where(np.isnan(Gs)); Gs[inds] = np.take(colmean, inds[1])
    Zs = np.ones((len(sdf), 1))
    if a.npc > 0:
        pcs = PCA(n_components=min(a.npc, Gs.shape[0]-1, Gs.shape[1]), random_state=0).fit_transform(Gs)
        Zs = np.column_stack([Zs, pcs])
    b, s, p = scan(sdf["slope"].to_numpy(float), Gs, Zs)
    gr = pd.DataFrame({"snp_id": snps, "beta": b, "se": s, "p": p})
    gr["padj"] = bh(gr["p"].to_numpy())
    gr.to_csv(os.path.join(a.outdir, "growth_rate_gwas.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"samples={g.shape[0]}\nsnps={len(snps)}\ntimes={ph['time'].nunique()}\n")
        f.write("model=per-time OLS with PCs; growth-rate slope GWAS\n")
    print(f"[temporal] snps={len(snps)} times={ph['time'].nunique()} tests={len(res)}")
    print(res.sort_values("p").head(5).to_string(index=False))
    if truth:
        top = res.sort_values("p").iloc[0]
        print(f"[temporal] demo causal {truth} top at time {top['time']}: {top['snp_id']==truth}")
    print(f"[temporal] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
