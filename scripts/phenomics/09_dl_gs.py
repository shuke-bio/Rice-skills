#!/usr/bin/env python3
###############################################################################
# 09_dl_gs.py -- deep-learning genomic selection baseline
#
# 方法：
#   ridge      : 岭回归（GBLUP 等价基线）
#   elasticnet : 弹性网
#   mlp        : 多层感知机（sklearn）
#   cnn        : 1D CNN（PyTorch，可选；未安装则跳过）
#
# 说明：
#   - 这是"可跑的 DL 基线"，不是 CLCNet/GP-WAITER 的官方实现；
#     外部模型权重可后续通过 --external 接口接入；
#   - 交叉验证按样本划分；报告 Pearson r / MSE / MAE。
#
# 用法：
#   python scripts/phenomics/09_dl_gs.py --genotype geno.csv --phenotype phe.csv \
#       --method all --folds 5 --outdir results/phenomics/09_dl_gs
#   python scripts/phenomics/09_dl_gs.py --demo --outdir results/phenomics/09_dl_gs_demo
###############################################################################
import argparse, os, sys, json, datetime, warnings
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, ElasticNet
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error

warnings.filterwarnings("ignore")

def make_demo(outdir, n=500, m=500, n_causal=20):
    rng = np.random.default_rng(42)
    G = rng.binomial(2, 0.3, size=(n, m)).astype(float)
    maf = np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2)
    G = G[:, maf >= 0.05]
    causal = rng.choice(G.shape[1], n_causal, replace=False)
    eff = rng.uniform(0.30, 0.50, n_causal)
    y = G[:, causal] @ eff + rng.normal(0, 0.6, n)
    # add 5 epistatic pairs so nonlinear models have signal to learn
    for i in range(5):
        a, b = causal[i], causal[i + 5]
        y = y + 0.35 * G[:, a] * G[:, b]
    samples = [f"S{i:03d}" for i in range(n)]
    geno = pd.DataFrame(G, columns=[f"snp{j:04d}" for j in range(G.shape[1])])
    geno.insert(0, "sample_id", samples)
    phe = pd.DataFrame({"sample_id": samples, "trait": y})
    geno.to_csv(os.path.join(outdir, "demo_genotype.csv"), index=False)
    phe.to_csv(os.path.join(outdir, "demo_phenotype.csv"), index=False)
    return os.path.join(outdir, "demo_genotype.csv"), os.path.join(outdir, "demo_phenotype.csv")

def load(genotype, phenotype, maf_cut=0.05):
    g = pd.read_csv(genotype); p = pd.read_csv(phenotype)
    if "sample_id" not in g.columns or "sample_id" not in p.columns:
        sys.exit("Both files need sample_id")
    trait = [c for c in p.columns if c != "sample_id"][0]
    d = g.merge(p[["sample_id", trait]], on="sample_id", how="inner")
    if len(d) < 30: sys.exit("Fewer than 30 shared samples")
    X = d.drop(columns=["sample_id", trait]).to_numpy(float)
    # mean impute
    colmean = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(colmean, inds[1])
    maf = np.minimum(X.mean(0) / 2, 1 - X.mean(0) / 2)
    keep = maf >= maf_cut
    X = X[:, keep]
    y = d[trait].to_numpy(float)
    ok = np.isfinite(y)
    return X[ok], y[ok]

def metrics(y, pred):
    r = pearsonr(y, pred)[0] if len(y) > 2 and np.std(pred) > 0 else np.nan
    return r, mean_squared_error(y, pred), mean_absolute_error(y, pred)

def fit_predict(method, Xtr, ytr, Xte):
    if method == "ridge":
        model = Ridge(alpha=1.0)
    elif method == "elasticnet":
        model = ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=5000)
    elif method == "mlp":
        model = MLPRegressor(hidden_layer_sizes=(256, 128), activation="relu",
                             alpha=1e-3, learning_rate_init=5e-4, max_iter=1200,
                             early_stopping=True, n_iter_no_change=30, random_state=0)
    elif method == "cnn":
        try:
            import torch
            import torch.nn as nn
        except Exception:
            return None, "PyTorch not installed"
        class CNN(nn.Module):
            def __init__(self, m):
                super().__init__()
                self.net = nn.Sequential(
                    nn.Conv1d(1, 16, 7, padding=3), nn.ReLU(), nn.MaxPool1d(2),
                    nn.Conv1d(16, 32, 5, padding=2), nn.ReLU(), nn.AdaptiveAvgPool1d(1),
                    nn.Flatten(), nn.Linear(32, 1))
            def forward(self, x): return self.net(x).squeeze(-1)
        torch.manual_seed(0)
        model = CNN(Xtr.shape[1])
        opt = torch.optim.Adam(model.parameters(), lr=1e-3)
        lossf = torch.nn.MSELoss()
        Xt = torch.tensor(Xtr[:, None, :], dtype=torch.float32)
        yt = torch.tensor(ytr, dtype=torch.float32)
        model.train()
        for _ in range(60):
            opt.zero_grad(); loss = lossf(model(Xt), yt); loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            pred = model(torch.tensor(Xte[:, None, :], dtype=torch.float32)).numpy()
        return pred, None
    else:
        return None, f"unknown method {method}"
    model.fit(Xtr, ytr)
    return model.predict(Xte), None

def run_method(method, X, y, folds, seed=0):
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    rows, preds = [], []
    for i, (tr, te) in enumerate(kf.split(X), 1):
        sc = StandardScaler().fit(X[tr])
        Xtr, Xte = sc.transform(X[tr]), sc.transform(X[te])
        if method == "mlp":
            from sklearn.decomposition import PCA
            pca = PCA(n_components=min(100, Xtr.shape[0] - 1, Xtr.shape[1]), random_state=0)
            Xtr = pca.fit_transform(Xtr); Xte = pca.transform(Xte)
        ys = StandardScaler().fit(y[tr].reshape(-1, 1))
        pred, err = fit_predict(method, Xtr, ys.transform(y[tr].reshape(-1, 1)).ravel(), Xte)
        if pred is None:
            return None, err
        pred = ys.inverse_transform(np.asarray(pred).reshape(-1, 1)).ravel()
        r, mse, mae = metrics(y[te], pred)
        rows.append([method, i, len(tr), len(te), r, mse, mae])
        preds.extend([[method, i, j, float(y[te][j]), float(pred[j])] for j in range(len(te))])
    return rows, preds

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genotype"); ap.add_argument("--phenotype")
    ap.add_argument("--method", default="all", choices=["ridge", "elasticnet", "mlp", "cnn", "all"])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--maf", type=float, default=0.05)
    ap.add_argument("--outdir", default="results/phenomics/09_dl_gs")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        g, p = make_demo(a.outdir)
        a.genotype, a.phenotype = g, p
    if not a.genotype or not a.phenotype:
        sys.exit("Need --genotype and --phenotype, or --demo")
    X, y = load(a.genotype, a.phenotype, a.maf)
    print(f"[dl_gs] samples={len(y)} snps={X.shape[1]} folds={a.folds}")

    methods = ["ridge", "elasticnet", "mlp", "cnn"] if a.method == "all" else [a.method]
    all_rows, all_preds, skipped = [], [], []
    for m in methods:
        rows, preds = run_method(m, X, y, a.folds)
        if rows is None:
            skipped.append(f"{m}: {preds}")
            print(f"[dl_gs] skip {m}: {preds}")
            continue
        all_rows.extend(rows); all_preds.extend(preds)
        arr = np.array([[r[4], r[5], r[6]] for r in rows], float)
        print(f"[dl_gs] {m:10s} r={np.nanmean(arr[:,0]):.3f} MSE={np.nanmean(arr[:,1]):.3f} MAE={np.nanmean(arr[:,2]):.3f}")
    if not all_rows: sys.exit("No method succeeded")

    res = pd.DataFrame(all_rows, columns=["method", "fold", "n_train", "n_test", "pearson_r", "mse", "mae"])
    res.to_csv(os.path.join(a.outdir, "dl_gs_results.csv"), index=False)
    pd.DataFrame(all_preds, columns=["method", "fold", "sample_index", "observed", "predicted"]).to_csv(
        os.path.join(a.outdir, "cv_predictions.csv"), index=False)
    summary = res.groupby("method")[["pearson_r", "mse", "mae"]].mean().reset_index()
    summary.to_csv(os.path.join(a.outdir, "dl_gs_summary.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"genotype={a.genotype}\nphenotype={a.phenotype}\n")
        f.write(f"samples={len(y)}\nsnps={X.shape[1]}\nfolds={a.folds}\n")
        f.write(f"methods={','.join(methods)}\nskipped={'; '.join(skipped) or 'none'}\n")
        f.write("note=DL baseline; not official CLCNet/GP-WAITER weights\n")
    print(summary.to_string(index=False))
    print(f"[dl_gs] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
