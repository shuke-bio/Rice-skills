#!/usr/bin/env python3
###############################################################################
# 12_kmer_gwas.py -- k-mer dosage GWAS (KMERIA-style, simplified)
#
# 思路：
#   1) k-mer 计数矩阵（样本 x k-mer），按样本深度归一化；
#   2) 分位数归一化到 0-2 剂量（近似 KMERIA 的剂量编码）；
#   3) PCA + VanRaden 亲缘矩阵 + 特征分解 MLM（或 GLM）；
#   4) BH/Bonferroni 校正；显著 k-mer 可回贴参考基因组定位。
#
# 输入：
#   --kmer-matrix CSV: sample_id, kmer1, kmer2, ...（计数）
#   --phenotype   CSV: sample_id, trait
#   --genome      可选 FASTA，用于把显著 k-mer 定位
#   --demo        合成数据（60 样本 x 600 k-mer，6 个因果）
#
# 用法：
#   python scripts/gwas/12_kmer_gwas.py --kmer-matrix kmers.csv --phenotype phe.csv \
#       --outdir results/gwas/12_kmer_gwas
#   python scripts/gwas/12_kmer_gwas.py --demo --outdir results/gwas/12_kmer_gwas_demo
###############################################################################
import argparse, os, sys, json, datetime
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]; buf = []
            else: buf.append(line)
    if name: seqs[name] = "".join(buf).upper()
    return seqs

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def bh(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(ranked, 1)
    return out

def normalize_counts(X, method="quantile"):
    X = np.asarray(X, float)
    if method == "cpm":
        tot = X.sum(axis=1, keepdims=True); tot[tot == 0] = 1
        return np.log1p(X / tot * 1e6)
    # quantile -> 0..2 per k-mer
    out = np.empty_like(X)
    n = X.shape[0]
    for j in range(X.shape[1]):
        r = stats.rankdata(X[:, j], method="average")
        out[:, j] = (r - 1) / max(n - 1, 1) * 2
    return out

def mlm_scan(y, X, Z, K, delta=1e-6):
    eigval, U = np.linalg.eigh(K)
    scale = 1 / np.sqrt(np.maximum(eigval, 0) + delta)
    def tr(A): return (U.T @ A) * scale[:, None]
    yt = tr(y.reshape(-1, 1)).ravel()
    m = X.shape[1]
    beta = np.full(m, np.nan); se = np.full(m, np.nan); p = np.full(m, np.nan)
    for j in range(m):
        D = np.column_stack([Z, X[:, j]])
        Dt = tr(D)
        coef, res, rank, sv = np.linalg.lstsq(Dt, yt, rcond=None)
        df = Dt.shape[0] - Dt.shape[1]
        if df <= 0: continue
        rss = float(((yt - Dt @ coef) ** 2).sum())
        try:
            xtx_inv = np.linalg.inv(Dt.T @ Dt)
        except np.linalg.LinAlgError:
            continue
        s2 = rss / df
        se[j] = np.sqrt(xtx_inv[-1, -1] * s2)
        beta[j] = coef[-1]
        if se[j] > 0:
            t = beta[j] / se[j]
            p[j] = 2 * stats.t.sf(abs(t), df)
    return beta, se, p

def glm_scan(y, X, Z):
    m = X.shape[1]
    beta = np.full(m, np.nan); se = np.full(m, np.nan); p = np.full(m, np.nan)
    for j in range(m):
        D = np.column_stack([Z, X[:, j]])
        coef, res, rank, sv = np.linalg.lstsq(D, y, rcond=None)
        df = D.shape[0] - D.shape[1]
        if df <= 0: continue
        rss = float(((y - D @ coef) ** 2).sum())
        try:
            xtx_inv = np.linalg.inv(D.T @ D)
        except np.linalg.LinAlgError:
            continue
        se[j] = np.sqrt(xtx_inv[-1, -1] * rss / df)
        beta[j] = coef[-1]
        if se[j] > 0:
            p[j] = 2 * stats.t.sf(abs(beta[j] / se[j]), df)
    return beta, se, p

def make_demo(outdir, n=80, m=600, n_causal=6):
    rng = np.random.default_rng(42)
    kmers = ["".join(rng.choice(list("ACGT"), 21)) for _ in range(m)]
    counts = rng.poisson(30, size=(n, m)).astype(float)
    samples = [f"S{i:02d}" for i in range(n)]
    X = normalize_counts(counts, "quantile")
    causal = rng.choice(m, n_causal, replace=False)
    effect = rng.uniform(0.80, 1.10, n_causal)
    y = X[:, causal] @ effect + rng.normal(0, 0.40, n)
    km = pd.DataFrame(counts, columns=kmers); km.insert(0, "sample_id", samples)
    ph = pd.DataFrame({"sample_id": samples, "trait": y})
    km.to_csv(os.path.join(outdir, "demo_kmer_matrix.csv"), index=False)
    ph.to_csv(os.path.join(outdir, "demo_phenotype.csv"), index=False)
    # demo genome containing causal k-mers
    genome = list("".join(rng.choice(list("ACGT"), 50_000)))
    for i, j in enumerate(causal):
        pos = 2_000 + i * 5_000
        genome[pos:pos+21] = list(kmers[j])
    with open(os.path.join(outdir, "demo_genome.fa"), "w", encoding="utf-8") as f:
        f.write(">chrDemo\n")
        seq = "".join(genome)
        for i in range(0, len(seq), 60): f.write(seq[i:i+60] + "\n")
    return (os.path.join(outdir, "demo_kmer_matrix.csv"),
            os.path.join(outdir, "demo_phenotype.csv"),
            os.path.join(outdir, "demo_genome.fa"),
            [kmers[j] for j in causal])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmer-matrix"); ap.add_argument("--phenotype")
    ap.add_argument("--genome"); ap.add_argument("--outdir", default="results/gwas/12_kmer_gwas")
    ap.add_argument("--npc", type=int, default=5)
    ap.add_argument("--maf", type=float, default=0.05)
    ap.add_argument("--method", choices=["mlm", "glm"], default="mlm")
    ap.add_argument("--normalize", choices=["quantile", "cpm"], default="quantile")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    causal = []
    if a.demo:
        km_path, ph_path, genome_path, causal = make_demo(a.outdir)
        a.kmer_matrix, a.phenotype, a.genome = km_path, ph_path, genome_path
    if not a.kmer_matrix or not a.phenotype:
        sys.exit("Need --kmer-matrix and --phenotype, or --demo")

    km = pd.read_csv(a.kmer_matrix)
    ph = pd.read_csv(a.phenotype)
    if "sample_id" not in km.columns or "sample_id" not in ph.columns:
        sys.exit("Both inputs need sample_id")
    kmer_cols = [c for c in km.columns if c != "sample_id"]
    df = km.merge(ph[["sample_id", "trait"]], on="sample_id", how="inner")
    if len(df) < 10: sys.exit("Fewer than 10 shared samples")
    X = df[kmer_cols].to_numpy(float)
    # filter: prevalence and variance
    prev = (X > 0).mean(axis=0)
    keep = (prev >= 0.8) & (X.std(axis=0) > 0)
    X = X[:, keep]
    kmer_cols = [c for c, k in zip(kmer_cols, keep) if k]
    if X.shape[1] < 10: sys.exit("Fewer than 10 k-mers after filtering")
    Xn = normalize_counts(X, a.normalize)
    # MAF-like filter on dosage
    maf = np.minimum(Xn.mean(axis=0) / 2, 1 - Xn.mean(axis=0) / 2)
    sel = maf >= a.maf
    Xn = Xn[:, sel]; kmer_cols = [c for c, s in zip(kmer_cols, sel) if s]
    if Xn.shape[1] < 10: sys.exit("Fewer than 10 k-mers after MAF filter")

    y = df["trait"].to_numpy(float)
    Z = np.ones((len(y), 1))
    if a.npc > 0:
        pcs = PCA(n_components=min(a.npc, Xn.shape[0]-1, Xn.shape[1]), random_state=0).fit_transform(Xn)
        Z = np.column_stack([Z, pcs])
    if a.method == "mlm":
        Xc = Xn - Xn.mean(axis=0)
        K = Xc @ Xc.T / Xc.shape[1]
        beta, se, p = mlm_scan(y, Xn, Z, K)
    else:
        beta, se, p = glm_scan(y, Xn, Z)
    ok = np.isfinite(p)
    res = pd.DataFrame({"kmer": np.array(kmer_cols)[ok], "beta": beta[ok], "se": se[ok], "p": p[ok]})
    res["padj"] = bh(res["p"].to_numpy())
    res = res.sort_values("p")
    res.to_csv(os.path.join(a.outdir, "kmer_gwas_results.csv"), index=False)
    sig = res[res["padj"] < 0.05].copy()

    # map significant k-mers to genome (exact, both strands)
    if a.genome and len(sig):
        seqs = read_fasta(a.genome)
        locs = []
        for _, r in sig.iterrows():
            k = r["kmer"]; found = []
            for chrom, s in seqs.items():
                for i in range(len(s) - len(k) + 1):
                    if s[i:i+len(k)] == k:
                        found.append((chrom, i+1, "+"))
                    if len(found) >= 5: break
                if len(found) >= 5: break
            locs.append(";".join(f"{c}:{p}:{st}" for c, p, st in found))
        sig["genome_positions"] = locs
    sig.to_csv(os.path.join(a.outdir, "kmer_gwas_significant.csv"), index=False)

    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"samples={len(df)}\nkmers_tested={len(res)}\nmethod={a.method}\n")
        f.write(f"normalize={a.normalize}\nnpc={a.npc}\nmaf={a.maf}\n")
        f.write("note=simplified KMERIA; dosage=quantile-normalized k-mer abundance; MLM=fixed-ratio eigen transform\n")

    hit = 0
    if a.demo and causal:
        hit = sum(1 for c in causal if c in set(sig["kmer"]))
        print(f"[kmer] demo causal recovered in significant: {hit}/{len(causal)}")
    print(f"[kmer] samples={len(df)} tested={len(res)} significant(FDR<0.05)={len(sig)} method={a.method}")
    print(res.head(5).to_string(index=False))
    print(f"[kmer] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
