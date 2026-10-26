#!/usr/bin/env python3
###############################################################################
# 08_tf_prioritize.py -- interpretable TF prioritization (PURE-like)
#
# 输入：
#   --expression CSV: sample_id, gene1, gene2, ...（log/TPM）
#   --deg        CSV: gene_id, log2FC, padj（差异基因表）
#   --tfs        CSV: gene_id（候选 TF 列表）
#   --motif      可选 CSV: tf_id, target_gene, score（启动子 motif 匹配）
#   --binding    可选 CSV: tf_id, target_gene（ChIP/DAP-seq 证据）
#
# 评分：
#   共表达（top 靶基因平均 |r|）+ motif 命中比例 + 结合证据比例，
#   各分量 z-score 后加权（默认 0.5/0.3/0.2），输出贡献与排名。
#
# 用法：
#   python scripts/transcriptomics/08_tf_prioritize.py --expression expr.csv \
#       --deg deg.csv --tfs tfs.csv --motif motif.csv --binding binding.csv \
#       --outdir results/transcriptomics/08_tf_prioritize
#   python scripts/transcriptomics/08_tf_prioritize.py --demo --outdir results/transcriptomics/08_tf_prioritize_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, zscore

def make_demo(outdir, n=60, n_genes=120, n_tf=8):
    rng = np.random.default_rng(42)
    tfs = [f"TF{i}" for i in range(1, n_tf + 1)]
    genes = tfs + [f"G{i:03d}" for i in range(n_genes - n_tf)]
    # baseline expression
    X = rng.normal(6, 1.2, size=(n, n_genes))
    truth = {}
    # 3 true regulators, each drives 15 targets
    for ti in range(3):
        targets = rng.choice([g for g in genes if g not in tfs], 15, replace=False)
        eff = rng.uniform(0.8, 1.2, 15)
        tf_expr = rng.normal(6, 1.0, n)
        X[:, genes.index(tfs[ti])] = tf_expr
        for g, e in zip(targets, eff):
            X[:, genes.index(g)] += e * (tf_expr - tf_expr.mean())
        truth[tfs[ti]] = set(targets)
    expr = pd.DataFrame(X, columns=genes); expr.insert(0, "sample_id", [f"S{i:02d}" for i in range(n)])
    # DEGs: targets + noise
    deg_genes = set()
    for s in truth.values(): deg_genes |= s
    deg_genes |= set(rng.choice(genes, 25, replace=False))
    deg = pd.DataFrame({"gene_id": sorted(deg_genes),
                        "log2FC": rng.normal(1.2, 0.4, len(deg_genes)),
                        "padj": rng.uniform(1e-6, 0.04, len(deg_genes))})
    tfs_df = pd.DataFrame({"gene_id": tfs})
    # motif/binding evidence only for true TFs (with some noise)
    motif_rows, bind_rows = [], []
    for ti in range(3):
        for g in truth[tfs[ti]]:
            if rng.random() < 0.8: motif_rows.append([tfs[ti], g, rng.uniform(5, 15)])
            if rng.random() < 0.7: bind_rows.append([tfs[ti], g])
    motif = pd.DataFrame(motif_rows, columns=["tf_id", "target_gene", "score"])
    bind = pd.DataFrame(bind_rows, columns=["tf_id", "target_gene"])
    expr.to_csv(os.path.join(outdir, "demo_expression.csv"), index=False)
    deg.to_csv(os.path.join(outdir, "demo_deg.csv"), index=False)
    tfs_df.to_csv(os.path.join(outdir, "demo_tfs.csv"), index=False)
    motif.to_csv(os.path.join(outdir, "demo_motif.csv"), index=False)
    bind.to_csv(os.path.join(outdir, "demo_binding.csv"), index=False)
    return ([os.path.join(outdir, f) for f in
             ["demo_expression.csv", "demo_deg.csv", "demo_tfs.csv", "demo_motif.csv", "demo_binding.csv"]], truth)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expression"); ap.add_argument("--deg"); ap.add_argument("--tfs")
    ap.add_argument("--motif"); ap.add_argument("--binding")
    ap.add_argument("--outdir", default="results/transcriptomics/08_tf_prioritize")
    ap.add_argument("--top-targets", type=int, default=50)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    truth = {}
    if a.demo:
        paths, truth = make_demo(a.outdir)
        a.expression, a.deg, a.tfs, a.motif, a.binding = paths
    if not a.expression or not a.deg or not a.tfs:
        sys.exit("Need --expression, --deg, --tfs, or --demo")

    expr = pd.read_csv(a.expression)
    deg = pd.read_csv(a.deg)
    tfs = pd.read_csv(a.tfs)
    if "sample_id" not in expr.columns: sys.exit("expression needs sample_id")
    genes = [c for c in expr.columns if c != "sample_id"]
    if "gene_id" not in deg.columns or "gene_id" not in tfs.columns: sys.exit("deg/tfs need gene_id")
    deg_genes = [g for g in deg["gene_id"] if g in genes]
    if len(deg_genes) < 5: sys.exit("Fewer than 5 DEGs present in expression matrix")
    tf_list = [t for t in tfs["gene_id"] if t in genes]
    if not tf_list: sys.exit("No candidate TFs present in expression matrix")

    X = expr[genes].to_numpy(float)
    gidx = {g: i for i, g in enumerate(genes)}
    # co-expression: top targets by |r|
    co = {}
    for tf in tf_list:
        x = X[:, gidx[tf]]
        rs = []
        for g in deg_genes:
            if g == tf: continue
            r = pearsonr(x, X[:, gidx[g]])[0]
            if np.isfinite(r): rs.append((g, abs(r), r))
        rs.sort(key=lambda t: t[1], reverse=True)
        top = rs[:a.top_targets]
        co[tf] = {"targets": {t[0] for t in top},
                  "coexp": float(np.mean([t[1] for t in top])) if top else 0.0,
                  "n": len(top)}

    motif_map, bind_map = {}, {}
    if a.motif and os.path.exists(a.motif):
        m = pd.read_csv(a.motif)
        for tf, sub in m.groupby("tf_id"):
            motif_map[tf] = set(sub["target_gene"])
    if a.binding and os.path.exists(a.binding):
        b = pd.read_csv(a.binding)
        for tf, sub in b.groupby("tf_id"):
            bind_map[tf] = set(sub["target_gene"])

    rows = []
    for tf in tf_list:
        targets = co[tf]["targets"]
        n = max(len(targets), 1)
        motif_frac = len(targets & motif_map.get(tf, set())) / n
        bind_frac = len(targets & bind_map.get(tf, set())) / n
        rows.append({"tf_id": tf, "n_targets": len(targets),
                     "coexp": co[tf]["coexp"], "motif_frac": motif_frac,
                     "binding_frac": bind_frac})
    res = pd.DataFrame(rows)
    # z-score components; constant -> 0
    for c in ["coexp", "motif_frac", "binding_frac"]:
        sd = res[c].std(ddof=0)
        res[c + "_z"] = (res[c] - res[c].mean()) / sd if sd > 0 else 0.0
    res["combined_score"] = (0.5 * res["coexp_z"] + 0.3 * res["motif_frac_z"] + 0.2 * res["binding_frac_z"])
    res = res.sort_values("combined_score", ascending=False).reset_index(drop=True)
    res["rank"] = np.arange(1, len(res) + 1)
    res.to_csv(os.path.join(a.outdir, "tf_ranking.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"expression={a.expression}\ndeg={a.deg}\ntfs={a.tfs}\n")
        f.write(f"motif={a.motif or 'none'}\nbinding={a.binding or 'none'}\n")
        f.write("score=0.5*coexp_z + 0.3*motif_frac_z + 0.2*binding_frac_z\n")
    print("[tf] top candidates:")
    print(res[["rank", "tf_id", "n_targets", "coexp", "motif_frac", "binding_frac", "combined_score"]].head(6).to_string(index=False))
    if truth:
        hits = sum(1 for t in truth if t in set(res.head(3)["tf_id"]))
        print(f"[tf] demo truth recovered in top3: {hits}/{len(truth)}")
    print(f"[tf] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
