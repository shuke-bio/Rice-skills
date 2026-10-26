#!/usr/bin/env python3
###############################################################################
# 13_burden_test.py -- rare-variant gene-based burden test
#
# 输入：
#   --genotype CSV: sample_id, snp1, ...（0/1/2）
#   --snp-info CSV: snp_id, gene_id（可含 maf）
#   --phenotype CSV: sample_id, trait
#   --covariates 可选
# 方法：
#   MAF <= --max-maf 的 SNP 按 gene 求和为 burden -> PCA 校正 OLS -> BH
#
# 用法：
#   python scripts/gwas/13_burden_test.py --genotype geno.csv --snp-info info.csv \
#       --phenotype phe.csv --outdir results/gwas/13_burden
#   python scripts/gwas/13_burden_test.py --demo --outdir results/gwas/13_burden_demo
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

def make_demo(outdir, n=300, genes=40, per_gene=10):
    rng = np.random.default_rng(42)
    rows = []
    for g in range(genes):
        for j in range(per_gene):
            rows.append([f"snp{g:02d}_{j:02d}", f"G{g:02d}"])
    info = pd.DataFrame(rows, columns=["snp_id", "gene_id"])
    X = rng.binomial(2, 0.005, size=(n, len(info))).astype(float)  # rare
    # make 5 causal genes with 3 variants each slightly more frequent
    causal = [1, 5, 9, 15, 22]
    for g in causal:
        idx = info.index[info.gene_id == f"G{g:02d}"][:3]
        X[:, idx] = rng.binomial(2, 0.015, size=(n, 3))
    y = rng.normal(0, 1, n)
    for g in causal:
        idx = info.index[info.gene_id == f"G{g:02d}"]
        y += 0.8 * X[:, idx].sum(1)
    geno = pd.DataFrame(X, columns=info.snp_id); geno.insert(0, "sample_id", [f"S{i:03d}" for i in range(n)])
    phe = pd.DataFrame({"sample_id": geno["sample_id"], "trait": y})
    geno.to_csv(os.path.join(outdir, "demo_genotype.csv"), index=False)
    info.to_csv(os.path.join(outdir, "demo_snp_info.csv"), index=False)
    phe.to_csv(os.path.join(outdir, "demo_phenotype.csv"), index=False)
    return (os.path.join(outdir,"demo_genotype.csv"), os.path.join(outdir,"demo_snp_info.csv"),
            os.path.join(outdir,"demo_phenotype.csv"), [f"G{g:02d}" for g in causal])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genotype"); ap.add_argument("--snp-info"); ap.add_argument("--phenotype")
    ap.add_argument("--max-maf", type=float, default=0.02)
    ap.add_argument("--min-variants", type=int, default=2)
    ap.add_argument("--npc", type=int, default=3)
    ap.add_argument("--outdir", default="results/gwas/13_burden")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    truth = None
    if a.demo:
        a.genotype, a.snp_info, a.phenotype, truth = make_demo(a.outdir)
    if not all([a.genotype, a.snp_info, a.phenotype]): sys.exit("Need --genotype --snp-info --phenotype, or --demo")
    g = pd.read_csv(a.genotype); info = pd.read_csv(a.snp_info); ph = pd.read_csv(a.phenotype)
    if not {"snp_id","gene_id"} <= set(info.columns): sys.exit("snp-info needs snp_id,gene_id")
    d = g.merge(ph, on="sample_id", how="inner")
    snps = [c for c in g.columns if c != "sample_id"]
    X = d[snps].to_numpy(float)
    colmean = np.nanmean(X, 0); inds = np.where(np.isnan(X)); X[inds] = np.take(colmean, inds[1])
    maf = np.minimum(X.mean(0)/2, 1-X.mean(0)/2)
    keep = maf <= a.max_maf
    X = X[:, keep]; snps = [s for s, k in zip(snps, keep) if k]
    info = info[info.snp_id.isin(snps)]
    gene_map = dict(zip(info.snp_id, info.gene_id))
    genes = {}
    for j, s in enumerate(snps):
        gid = gene_map.get(s)
        if gid: genes.setdefault(gid, []).append(j)
    genes = {k: v for k, v in genes.items() if len(v) >= a.min_variants}
    if not genes: sys.exit("No gene has enough rare variants")
    burden = np.column_stack([X[:, idx].sum(1) for idx in genes.values()])
    Z = np.ones((len(d), 1))
    if a.npc > 0:
        pcs = PCA(n_components=min(a.npc, X.shape[0]-1, X.shape[1]), random_state=0).fit_transform(X)
        Z = np.column_stack([Z, pcs])
    y = d["trait"].to_numpy(float)
    rows = []
    for j, gid in enumerate(genes):
        D = np.column_stack([Z, burden[:, j]])
        coef, *_ = np.linalg.lstsq(D, y, rcond=None)
        df = D.shape[0] - D.shape[1]
        rss = float(((y - D @ coef) ** 2).sum())
        try: inv = np.linalg.inv(D.T @ D)
        except np.linalg.LinAlgError: continue
        se = np.sqrt(inv[-1,-1] * rss / df)
        p = 2 * stats.t.sf(abs(coef[-1]/se), df) if se > 0 else np.nan
        rows.append([gid, len(genes[gid]), float(coef[-1]), float(se), float(p)])
    res = pd.DataFrame(rows, columns=["gene_id","n_rare","beta","se","p"])
    res["padj"] = bh(res["p"].to_numpy())
    res = res.sort_values("p")
    res.to_csv(os.path.join(a.outdir, "burden_results.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"samples={len(d)}\nsnps={len(snps)}\ngenes={len(genes)}\n")
        f.write(f"max_maf={a.max_maf}\nmin_variants={a.min_variants}\n")
    print(f"[burden] genes={len(genes)}; FDR<0.05={int((res['padj']<0.05).sum())}")
    print(res.head(6).to_string(index=False))
    if truth:
        hit = sum(1 for x in truth if x in set(res[res.padj<0.05].gene_id))
        print(f"[burden] demo causal recovered: {hit}/{len(truth)}")
    print(f"[burden] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
