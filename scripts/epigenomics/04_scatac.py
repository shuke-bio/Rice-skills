#!/usr/bin/env python3
###############################################################################
# 04_scatac.py -- single-cell ATAC lightweight pipeline (TF-IDF + LSI + clustering)
#
# 输入：
#   --matrix CSV: cell_id, peak1, peak2, ...（peak 计数）
#   --metadata 可选: cell_id, cell_type（用于评估）
# 方法：
#   TF-IDF -> TruncatedSVD (LSI) -> KMeans -> cluster marker peaks
#
# 用法：
#   python scripts/epigenomics/04_scatac.py --matrix scatac.csv --metadata meta.csv \
#       --n-clusters 8 --outdir results/epigenomics/04_scatac
#   python scripts/epigenomics/04_scatac.py --demo --outdir results/epigenomics/04_scatac_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

def make_demo(outdir, n=200, p=500):
    rng = np.random.default_rng(42)
    types = ["A", "B", "C"]
    labels = rng.choice(types, n)
    X = rng.poisson(2, size=(n, p)).astype(float)
    for i, t in enumerate(labels):
        idx = types.index(t)
        X[i, idx*40:(idx+1)*40] += rng.integers(10, 30, 40)
    mat = pd.DataFrame(X, columns=[f"peak{i:03d}" for i in range(p)])
    mat.insert(0, "cell_id", [f"C{i:03d}" for i in range(n)])
    meta = pd.DataFrame({"cell_id": mat["cell_id"], "cell_type": labels})
    mat.to_csv(os.path.join(outdir, "demo_matrix.csv"), index=False)
    meta.to_csv(os.path.join(outdir, "demo_metadata.csv"), index=False)
    return os.path.join(outdir, "demo_matrix.csv"), os.path.join(outdir, "demo_metadata.csv"), meta

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix"); ap.add_argument("--metadata")
    ap.add_argument("--n-clusters", type=int, default=3)
    ap.add_argument("--outdir", default="results/epigenomics/04_scatac")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    truth = None
    if a.demo:
        a.matrix, a.metadata, truth = make_demo(a.outdir)
    if not a.matrix: sys.exit("Need --matrix or --demo")
    mat = pd.read_csv(a.matrix)
    if "cell_id" not in mat.columns: sys.exit("matrix needs cell_id")
    peaks = [c for c in mat.columns if c != "cell_id"]
    X = mat[peaks].to_numpy(float)
    # TF-IDF
    tf = TfidfTransformer(norm="l2", use_idf=True, smooth_idf=True)
    Xt = tf.fit_transform(X)
    k = min(30, Xt.shape[0] - 1, Xt.shape[1] - 1)
    lsi = TruncatedSVD(n_components=k, random_state=0).fit_transform(Xt)
    km = KMeans(n_clusters=a.n_clusters, n_init=20, random_state=0).fit(lsi[:, 1:])
    clusters = km.labels_
    out = pd.DataFrame({"cell_id": mat["cell_id"], "cluster": clusters})
    out.to_csv(os.path.join(a.outdir, "cell_clusters.csv"), index=False)
    # markers: cluster mean - rest mean
    rows = []
    for c in range(a.n_clusters):
        idx = clusters == c
        diff = X[idx].mean(0) - X[~idx].mean(0)
        for j in np.argsort(diff)[-10:][::-1]:
            rows.append([c, peaks[j], float(diff[j])])
    pd.DataFrame(rows, columns=["cluster", "peak", "score"]).to_csv(
        os.path.join(a.outdir, "cluster_marker_peaks.csv"), index=False)
    pd.DataFrame(lsi[:, :2], columns=["LSI1", "LSI2"]).assign(cell_id=mat["cell_id"]).to_csv(
        os.path.join(a.outdir, "lsi_coordinates.csv"), index=False)
    ari = adjusted_rand_score(truth["cell_type"], clusters) if truth is not None else None
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"cells={X.shape[0]}\npeaks={X.shape[1]}\nlsi={k}\nclusters={a.n_clusters}\n")
        if ari is not None: f.write(f"demo_ARI={ari:.3f}\n")
    print(f"[scatac] cells={X.shape[0]} peaks={X.shape[1]} clusters={a.n_clusters}")
    if ari is not None: print(f"[scatac] demo ARI vs truth: {ari:.3f}")
    print(f"[scatac] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
