#!/usr/bin/env python3
###############################################################################
# 02_integration_annotation.py -- single-cell integration, clustering, annotation
#
# 输入（CSV，避免强制 h5ad 依赖）：
#   --matrix   cells x genes（原始计数或 log 表达）
#   --metadata cell_id, batch(可选), cell_type(可选，用于评估)
#   --markers  可选 CSV: cell_type, gene（标记基因）
# 方法：
#   CPM/log1p -> HVG -> scale -> PCA -> 批次中心化（harmonypy 可用则 Harmony）
#   -> KMeans/Leiden -> marker 打分注释；scrublet 可选
#
# 用法：
#   python scripts/single_cell/02_integration_annotation.py --matrix expr.csv \
#       --metadata meta.csv --markers markers.csv --outdir results/single_cell/02_integration
#   python scripts/single_cell/02_integration_annotation.py --demo --outdir results/single_cell/02_integration_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score

def make_demo(outdir, n=300, g=500, batches=3):
    rng = np.random.default_rng(42)
    types = ["Meristem", "Epidermis", "Vascular"]
    markers = {
        "Meristem": [f"G{i:03d}" for i in range(0, 20)],
        "Epidermis": [f"G{i:03d}" for i in range(20, 40)],
        "Vascular": [f"G{i:03d}" for i in range(40, 60)],
    }
    genes = [f"G{i:03d}" for i in range(g)]
    labels = rng.choice(types, n)
    batch = rng.choice([f"B{i}" for i in range(batches)], n)
    X = rng.negative_binomial(2, 0.3, size=(n, g)).astype(float)
    for i, t in enumerate(labels):
        for mg in markers[t]:
            X[i, genes.index(mg)] += rng.integers(20, 60)
    # batch effect on first 50 genes
    for b in range(batches):
        idx = np.where(batch == f"B{b}")[0]
        X[np.ix_(idx, np.arange(0, 50))] += rng.integers(0, 15, size=(len(idx), 50))
    mat = pd.DataFrame(X, columns=genes); mat.insert(0, "cell_id", [f"C{i:03d}" for i in range(n)])
    meta = pd.DataFrame({"cell_id": mat["cell_id"], "batch": batch, "cell_type": labels})
    mk = pd.DataFrame([[t, g] for t, gs in markers.items() for g in gs], columns=["cell_type", "gene"])
    mat.to_csv(os.path.join(outdir, "demo_matrix.csv"), index=False)
    meta.to_csv(os.path.join(outdir, "demo_metadata.csv"), index=False)
    mk.to_csv(os.path.join(outdir, "demo_markers.csv"), index=False)
    return (os.path.join(outdir, "demo_matrix.csv"), os.path.join(outdir, "demo_metadata.csv"),
            os.path.join(outdir, "demo_markers.csv"), meta)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix"); ap.add_argument("--metadata"); ap.add_argument("--markers")
    ap.add_argument("--npc", type=int, default=30)
    ap.add_argument("--n-clusters", type=int, default=3)
    ap.add_argument("--batch-key", default="batch")
    ap.add_argument("--outdir", default="results/single_cell/02_integration")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    truth = None
    if a.demo:
        m, meta, mk, truth = make_demo(a.outdir)
        a.matrix, a.metadata, a.markers = m, meta, mk
    if not a.matrix or not a.metadata:
        sys.exit("Need --matrix and --metadata, or --demo")

    mat = pd.read_csv(a.matrix)
    meta = pd.read_csv(a.metadata)
    if "cell_id" not in mat.columns or "cell_id" not in meta.columns:
        sys.exit("matrix/metadata need cell_id")
    meta = meta.set_index("cell_id").loc[mat["cell_id"]].reset_index()
    genes = [c for c in mat.columns if c != "cell_id"]
    X = mat[genes].to_numpy(float)
    # CPM/log1p if likely counts
    if X.max() > 50:
        tot = X.sum(1, keepdims=True); tot[tot == 0] = 1
        X = np.log1p(X / tot * 1e6)
    # HVG
    var = X.var(0)
    hvg = np.argsort(var)[-min(2000, len(var)):]
    X = X[:, hvg]; genes = [genes[i] for i in hvg]
    X = StandardScaler().fit_transform(X)
    npc = min(a.npc, X.shape[0] - 1, X.shape[1])
    pcs = PCA(n_components=npc, random_state=0).fit_transform(X)
    # batch centering or harmony
    batch = meta[a.batch_key].astype(str).to_numpy() if a.batch_key in meta.columns else None
    method = "none"
    if batch is not None:
        try:
            import harmonypy
            ho = harmonypy.run_harmony(pcs, meta, [a.batch_key])
            pcs = np.asarray(ho.Z_corr).T
            method = "harmony"
        except Exception:
            for b in np.unique(batch):
                idx = batch == b
                pcs[idx] -= pcs[idx].mean(0)
            method = "batch-centering"
    km = KMeans(n_clusters=a.n_clusters, n_init=20, random_state=0).fit(pcs)
    clusters = km.labels_
    # marker annotation
    cell_type = [f"cluster_{c}" for c in clusters]
    if a.markers and os.path.exists(a.markers):
        mk = pd.read_csv(a.markers)
        gene_idx = {g: i for i, g in enumerate(genes)}
        scores = []
        for c in range(a.n_clusters):
            cells = clusters == c
            row = {}
            for t, sub in mk.groupby("cell_type"):
                gs = [gene_idx[g] for g in sub["gene"] if g in gene_idx]
                row[t] = float(np.mean(X[np.ix_(cells, gs)])) if gs else -np.inf
            scores.append(row)
        sdf = pd.DataFrame(scores)
        best = sdf.idxmax(1).to_dict()
        cell_type = [best[c] for c in clusters]
        sdf.to_csv(os.path.join(a.outdir, "cluster_marker_scores.csv"), index=False)
    out = pd.DataFrame({"cell_id": mat["cell_id"], "cluster": clusters, "cell_type": cell_type})
    if batch is not None: out["batch"] = batch
    out.to_csv(os.path.join(a.outdir, "cell_metadata.csv"), index=False)
    # top markers per cluster (simple mean difference)
    rows = []
    for c in range(a.n_clusters):
        cells = clusters == c
        diff = X[cells].mean(0) - X[~cells].mean(0)
        for j in np.argsort(diff)[-10:][::-1]:
            rows.append([c, genes[j], float(diff[j])])
    pd.DataFrame(rows, columns=["cluster", "gene", "score"]).to_csv(
        os.path.join(a.outdir, "cluster_markers.csv"), index=False)
    ari = adjusted_rand_score(truth["cell_type"], cell_type) if truth is not None else None
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"cells={X.shape[0]}\ngenes_used={X.shape[1]}\npcs={npc}\n")
        f.write(f"integration={method}\nclusters={a.n_clusters}\n")
        if ari is not None: f.write(f"demo_ARI={ari:.3f}\n")
    print(f"[sc] cells={X.shape[0]} genes={X.shape[1]} integration={method} clusters={a.n_clusters}")
    if ari is not None: print(f"[sc] demo ARI vs truth: {ari:.3f}")
    print(f"[sc] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
