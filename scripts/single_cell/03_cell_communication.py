#!/usr/bin/env python3
###############################################################################
# 03_cell_communication.py -- ligand-receptor communication scoring
#
# 输入：
#   --matrix CSV: cell_id, gene1, gene2, ...（counts/log）
#   --labels CSV: cell_id, cell_type
#   --lr     可选 CSV: ligand, receptor（默认内置一组植物 LR 候选）
# 方法：
#   细胞类型平均表达 -> score = mean(ligand, source) * mean(receptor, target)
#   -> 标签置换检验 -> BH
#
# 说明：
#   植物 LR 数据库远不如哺乳动物完整；内置列表仅作演示，正式分析请用
#   CellPhoneDB 风格的自定义 LR 表或文献整理的作物 LR 对。
#
# 用法：
#   python scripts/single_cell/03_cell_communication.py --matrix expr.csv \
#       --labels labels.csv --lr lr.csv --outdir results/single_cell/03_ccc
#   python scripts/single_cell/03_cell_communication.py --demo --outdir results/single_cell/03_ccc_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from scipy.stats import rankdata

BUILTIN_LR = [
    ("CLV3", "CLV1"), ("CLV3", "CLV2"), ("CLE40", "CLV1"),
    ("EPF1", "ERECTA"), ("EPFL4", "ERECTA"), ("IDA", "HAESA"),
    ("TDIF", "TDR"), ("PSK", "PSKR1"), ("RALF1", "FER"),
    ("CEP1", "CEPR1"), ("TGFB1", "TGFBR1"), ("IL6", "IL6R"),
]

def bh(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(ranked, 1)
    return out

def make_demo(outdir, n=240, g=300):
    rng = np.random.default_rng(42)
    types = ["Meristem", "Vascular", "Epidermis"]
    labels = rng.choice(types, n)
    genes = [f"G{i:03d}" for i in range(g)] + ["CLV3", "CLV1", "EPF1", "ERECTA"]
    X = rng.negative_binomial(2, 0.3, size=(n, len(genes))).astype(float)
    for i, t in enumerate(labels):
        if t == "Meristem": X[i, genes.index("CLV3")] += 60
        if t == "Vascular": X[i, genes.index("CLV1")] += 60
    mat = pd.DataFrame(X, columns=genes); mat.insert(0, "cell_id", [f"C{i:03d}" for i in range(n)])
    lab = pd.DataFrame({"cell_id": mat["cell_id"], "cell_type": labels})
    mat.to_csv(os.path.join(outdir, "demo_matrix.csv"), index=False)
    lab.to_csv(os.path.join(outdir, "demo_labels.csv"), index=False)
    return os.path.join(outdir, "demo_matrix.csv"), os.path.join(outdir, "demo_labels.csv")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix"); ap.add_argument("--labels"); ap.add_argument("--lr")
    ap.add_argument("--permutations", type=int, default=200)
    ap.add_argument("--outdir", default="results/single_cell/03_ccc")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.matrix, a.labels = make_demo(a.outdir)
    if not a.matrix or not a.labels:
        sys.exit("Need --matrix and --labels, or --demo")
    mat = pd.read_csv(a.matrix); lab = pd.read_csv(a.labels)
    if "cell_id" not in mat.columns or "cell_id" not in lab.columns:
        sys.exit("matrix/labels need cell_id")
    lab = lab.set_index("cell_id").loc[mat["cell_id"]].reset_index()
    genes = [c for c in mat.columns if c != "cell_id"]
    X = mat[genes].to_numpy(float)
    if X.max() > 50:
        tot = X.sum(1, keepdims=True); tot[tot == 0] = 1
        X = np.log1p(X / tot * 1e6)
    types = sorted(lab["cell_type"].unique())
    means = pd.DataFrame({t: X[(lab["cell_type"] == t).to_numpy()].mean(0) for t in types}, index=genes)
    lr = pd.read_csv(a.lr) if a.lr and os.path.exists(a.lr) else pd.DataFrame(BUILTIN_LR, columns=["ligand", "receptor"])
    rows = []
    rng = np.random.default_rng(0)
    labels_arr = lab["cell_type"].to_numpy()
    for _, r in lr.iterrows():
        lig, rec = r["ligand"], r["receptor"]
        if lig not in means.index or rec not in means.index: continue
        for s in types:
            for t in types:
                if s == t: continue
                obs = means.loc[lig, s] * means.loc[rec, t]
                null = []
                for _ in range(a.permutations):
                    perm = rng.permutation(labels_arr)
                    ms = X[perm == s, genes.index(lig)].mean()
                    mt = X[perm == t, genes.index(rec)].mean()
                    null.append(ms * mt)
                p = (np.sum(np.array(null) >= obs) + 1) / (a.permutations + 1)
                rows.append([s, t, lig, rec, obs, p])
    res = pd.DataFrame(rows, columns=["source", "target", "ligand", "receptor", "score", "p"])
    res["padj"] = bh(res["p"].to_numpy())
    res = res.sort_values("score", ascending=False)
    res.to_csv(os.path.join(a.outdir, "cell_communication.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"cells={len(lab)}\ntypes={len(types)}\nlr_pairs={len(lr)}\npermutations={a.permutations}\n")
        f.write("score=mean(ligand,source)*mean(receptor,target); plant LR list is illustrative\n")
    print(f"[ccc] types={types}; interactions tested={len(res)}; FDR<0.05={int((res['padj']<0.05).sum())}")
    print(res.head(6).to_string(index=False))
    print(f"[ccc] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
