#!/usr/bin/env python3
###############################################################################
# 02_deconvolution.py -- spatial spot deconvolution with NNLS signatures
#
# 输入：
#   --sc       CSV: cell_id, gene1, ...（scRNA 参考）
#   --labels   CSV: cell_id, cell_type
#   --spatial  CSV: spot_id, gene1, ...（空间表达）
#   --coords   可选 CSV: spot_id, x, y
#
# 方法：
#   选标记基因（类型间差异大）-> 类型平均表达作签名 -> 每个 spot 做 NNLS
#   -> 归一化为比例。是 cell2location/SPADE 的轻量替代，不是概率模型。
#
# 用法：
#   python scripts/spatial/02_deconvolution.py --sc sc.csv --labels labels.csv \
#       --spatial spots.csv --coords coords.csv --outdir results/spatial/02_deconv
#   python scripts/spatial/02_deconvolution.py --demo --outdir results/spatial/02_deconv_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from scipy.optimize import nnls

def make_demo(outdir, n_cells=300, n_spots=120, g=200):
    rng = np.random.default_rng(42)
    types = ["Mesophyll", "Epidermis", "Vascular"]
    genes = [f"G{i:03d}" for i in range(g)]
    labels = rng.choice(types, n_cells)
    X = rng.negative_binomial(2, 0.3, size=(n_cells, g)).astype(float)
    for i, t in enumerate(labels):
        idx = types.index(t)
        X[i, idx*20:(idx+1)*20] += rng.integers(20, 60, 20)
    sc = pd.DataFrame(X, columns=genes); sc.insert(0, "cell_id", [f"C{i:03d}" for i in range(n_cells)])
    lab = pd.DataFrame({"cell_id": sc["cell_id"], "cell_type": labels})
    # spots: mixture of 2-3 types
    prop = rng.dirichlet([2, 2, 2], n_spots)
    S = np.zeros((n_spots, g))
    for k in range(n_spots):
        for ti, t in enumerate(types):
            cells = X[labels == t]
            S[k] += prop[k, ti] * cells.mean(0)
    S = rng.poisson(np.maximum(S, 0)).astype(float)
    sp = pd.DataFrame(S, columns=genes); sp.insert(0, "spot_id", [f"P{i:03d}" for i in range(n_spots)])
    coords = pd.DataFrame({"spot_id": sp["spot_id"], "x": rng.uniform(0, 100, n_spots), "y": rng.uniform(0, 100, n_spots)})
    sc.to_csv(os.path.join(outdir, "demo_sc.csv"), index=False)
    lab.to_csv(os.path.join(outdir, "demo_labels.csv"), index=False)
    sp.to_csv(os.path.join(outdir, "demo_spots.csv"), index=False)
    coords.to_csv(os.path.join(outdir, "demo_coords.csv"), index=False)
    pd.DataFrame(prop, columns=types).assign(spot_id=sp["spot_id"]).to_csv(
        os.path.join(outdir, "demo_truth.csv"), index=False)
    return ("demo_sc.csv", "demo_labels.csv", "demo_spots.csv", "demo_coords.csv", types)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sc"); ap.add_argument("--labels"); ap.add_argument("--spatial"); ap.add_argument("--coords")
    ap.add_argument("--outdir", default="results/spatial/02_deconv")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    truth = None; types = None
    if a.demo:
        *files, types = make_demo(a.outdir)
        a.sc, a.labels, a.spatial, a.coords = [os.path.join(a.outdir, f) for f in files[:4]]
        truth = pd.read_csv(os.path.join(a.outdir, "demo_truth.csv")).set_index("spot_id")
    if not a.sc or not a.labels or not a.spatial:
        sys.exit("Need --sc, --labels, --spatial, or --demo")
    sc = pd.read_csv(a.sc); lab = pd.read_csv(a.labels); sp = pd.read_csv(a.spatial)
    if "cell_id" not in sc.columns or "cell_id" not in lab.columns or "spot_id" not in sp.columns:
        sys.exit("sc/labels need cell_id; spatial needs spot_id")
    lab = lab.set_index("cell_id").loc[sc["cell_id"]].reset_index()
    genes = [g for g in sc.columns if g != "cell_id" and g in set(sp.columns)]
    X = sc[genes].to_numpy(float)
    if X.max() > 50:
        tot = X.sum(1, keepdims=True); tot[tot == 0] = 1
        X = np.log1p(X / tot * 1e6)
    types = sorted(lab["cell_type"].unique())
    means = np.vstack([X[(lab["cell_type"] == t).to_numpy()].mean(0) for t in types]).T  # genes x types
    # marker selection: top 300 genes by max type mean / overall
    ratio = means.max(1) / (means.mean(1) + 1e-9)
    keep = np.argsort(ratio)[-min(300, len(ratio)):]
    sig = means[keep, :]
    S = sp[genes].to_numpy(float)[:, keep]
    if S.max() > 50:
        tot = S.sum(1, keepdims=True); tot[tot == 0] = 1
        S = np.log1p(S / tot * 1e6)
    props = np.zeros((S.shape[0], len(types)))
    for i in range(S.shape[0]):
        w, _ = nnls(sig, S[i])
        props[i] = w / w.sum() if w.sum() > 0 else w
    res = pd.DataFrame(props, columns=types); res.insert(0, "spot_id", sp["spot_id"])
    if a.coords and os.path.exists(a.coords):
        co = pd.read_csv(a.coords)
        res = res.merge(co, on="spot_id", how="left")
    res.to_csv(os.path.join(a.outdir, "spot_proportions.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"spots={len(res)}\ntypes={len(types)}\ngenes_used={len(keep)}\n")
        f.write("method=NNLS on type-mean signatures; not a probabilistic model\n")
    if truth is not None:
        cors = [np.corrcoef(res[t], truth.loc[res['spot_id'], t])[0, 1] for t in types]
        print("[deconv] demo correlation with truth:", ", ".join(f"{t}={c:.2f}" for t, c in zip(types, cors)))
    print(f"[deconv] spots={len(res)} types={len(types)} -> {a.outdir}")

if __name__ == "__main__":
    main()
