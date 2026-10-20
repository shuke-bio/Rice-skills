#!/usr/bin/env python3
###############################################################################
# 02_check_data.py -- 分析前的数据体检（不依赖 pandas）
#
# 检查：维度、重复 ID、缺失率、零方差特征、文库大小、PCA 离群、
#       分组样本量、批次×分组交叉表（混杂检测）、多组学样本重叠。
#
# 用法：
#   python scripts/plan/02_check_data.py --matrix counts.csv --metadata meta.csv \
#       --group-col group --batch-col batch --outdir results/qc
#   python scripts/plan/02_check_data.py --demo
#
# 输入格式（默认 samples_rows）：
#   counts.csv : 第一列 sample_id，其余列 feature（基因/代谢物）
#   meta.csv   : 第一列 sample_id，至少含 group；可选 batch
#   若矩阵是 features x samples，用 --orientation features_rows
#
# 输出：qc_report.md、issues.json、sample_flags.csv
###############################################################################
import argparse, csv, json, os, sys

try:
    import numpy as np
except Exception:
    sys.exit("需要 numpy：pip install numpy")

def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    if not rows:
        sys.exit(f"空文件: {path}")
    header, body = rows[0], [r for r in rows[1:] if any(x.strip() for x in r)]
    return header, body

def to_matrix(path, orientation="samples_rows"):
    header, body = read_csv(path)
    sample_ids = [r[0].strip() for r in body]
    feats = [h.strip() for h in header[1:]]
    data = np.full((len(body), len(feats)), np.nan, dtype=float)
    for i, r in enumerate(body):
        for j, v in enumerate(r[1:]):
            v = (v or "").strip()
            if v == "" or v.upper() in ("NA", "NAN", "NULL"):
                continue
            try:
                data[i, j] = float(v)
            except ValueError:
                pass
    if orientation == "features_rows":
        data = data.T
        sample_ids, feats = feats, sample_ids
    return sample_ids, feats, data

def read_meta(path, group_col="group", batch_col="batch"):
    header, body = read_csv(path)
    idx = {h.strip(): k for k, h in enumerate(header)}
    if group_col not in idx:
        sys.exit(f"metadata 缺少列: {group_col}")
    out = {}
    for r in body:
        sid = r[0].strip()
        out[sid] = {"group": r[idx[group_col]].strip(),
                    "batch": r[idx[batch_col]].strip() if batch_col in idx else ""}
    return out

def pca_scores(X, k=2):
    X = np.nan_to_num(X, nan=0.0)
    X = X - X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True); sd[sd == 0] = 1
    X = X / sd
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    return U[:, :k] * S[:k]

def check(matrix, meta, second=None):
    ids, feats, X = matrix
    issues = []
    def add(level, code, msg, detail=None):
        issues.append({"level": level, "code": code, "message": msg, "detail": detail or {}})

    if len(set(ids)) != len(ids):
        add("error", "duplicate_samples", "样本 ID 有重复", {"n": len(ids) - len(set(ids))})
    if len(set(feats)) != len(feats):
        add("warning", "duplicate_features", "特征名有重复", {"n": len(feats) - len(set(feats))})

    miss_s = np.isnan(X).mean(axis=1)
    miss_f = np.isnan(X).mean(axis=0)
    for i in np.where(miss_s > 0.2)[0]:
        add("error", "sample_missing", f"样本缺失率 >20%：{ids[i]}", {"sample": ids[i], "rate": round(float(miss_s[i]), 3)})
    if (miss_f > 0.2).any():
        add("warning", "feature_missing", f"{int((miss_f>0.2).sum())} 个特征缺失率 >20%", {"n": int((miss_f>0.2).sum())})

    var = np.nanvar(X, axis=0)
    const = int((var == 0).sum())
    if const:
        add("warning", "constant_features", f"{const} 个零方差特征", {"n": const})

    lib = np.nansum(X, axis=1)
    med = float(np.median(lib)) if len(lib) else 0
    if med > 0:
        for i in np.where(lib < 0.25 * med)[0]:
            add("error", "low_library", f"文库大小 <25% 中位数：{ids[i]}", {"sample": ids[i], "lib": round(float(lib[i]), 1)})

    if len(ids) >= 5 and X.shape[1] >= 2:
        top = np.argsort(var)[-min(500, X.shape[1]):]
        pc = pca_scores(X[:, top])
        for j, name in enumerate(["PC1", "PC2"]):
            if pc.shape[1] <= j: continue
            z = (pc[:, j] - pc[:, j].mean()) / (pc[:, j].std() or 1)
            for i in np.where(np.abs(z) > 3)[0]:
                add("warning", "pca_outlier", f"{name} 离群：{ids[i]} (z={z[i]:.1f})", {"sample": ids[i], "pc": name})

    groups, batches = {}, {}
    for sid in ids:
        m = meta.get(sid, {})
        groups.setdefault(m.get("group", "NA"), []).append(sid)
        if m.get("batch"): batches.setdefault(m.get("batch"), []).append(sid)
    if groups:
        for g, s in groups.items():
            if len(s) < 2:
                add("error", "group_too_small", f"分组 {g} 只有 {len(s)} 个样本", {"group": g, "n": len(s)})
        if batches:
            conf = []
            for g, s in groups.items():
                bs = {meta.get(x, {}).get("batch", "") for x in s}
                if len(bs) == 1 and all(x for x in bs):
                    conf.append((g, list(bs)[0]))
            if conf:
                add("error", "batch_confounded", "批次与分组完全混杂", {"pairs": conf})

    if second:
        ids2 = set(second[0])
        overlap = len(set(ids) & ids2)
        add("info", "multiomics_overlap", f"与第二层数据重叠样本 {overlap} 个", {"overlap": overlap, "n1": len(ids), "n2": len(ids2)})

    return {"n_samples": len(ids), "n_features": len(feats),
            "missing_rate": round(float(np.isnan(X).mean()), 4),
            "issues": issues,
            "sample_flags": {ids[i]: {"missing": round(float(miss_s[i]), 3),
                                      "library": round(float(lib[i]), 1)}
                             for i in range(len(ids))}}

def render_report(r):
    L = ["# 数据体检报告", "",
         f"- 样本数：{r['n_samples']}｜特征数：{r['n_features']}｜总体缺失率：{r['missing_rate']:.2%}", ""]
    errs = [x for x in r["issues"] if x["level"] == "error"]
    warns = [x for x in r["issues"] if x["level"] == "warning"]
    infos = [x for x in r["issues"] if x["level"] == "info"]
    if errs:
        L += ["## ⛔ 必须先处理", ""] + [f"- **{x['code']}**：{x['message']}" for x in errs] + [""]
    if warns:
        L += ["## ⚠️ 需要说明", ""] + [f"- **{x['code']}**：{x['message']}" for x in warns] + [""]
    if infos:
        L += ["## ℹ️ 信息", ""] + [f"- {x['message']}" for x in infos] + [""]
    if not errs and not warns:
        L += ["## ✅ 未发现阻断性问题", ""]
    L += ["---", "下一步：把 `n_per_group`、`replicates`、`batch` 等字段填入数据契约，"
          "再运行 `scripts/plan/01_make_plan.py` 生成三档方案。"]
    return "\n".join(L)

def demo(outdir):
    rng = np.random.default_rng(42)
    n, m = 40, 300
    ids = [f"S{i:02d}" for i in range(n)]
    feats = [f"G{j}" for j in range(m)]
    X = rng.normal(100, 15, size=(n, m))
    X[:, :15] = 5.0                       # 零方差
    X[3, :] = np.nan                       # 缺失样本
    X[3, :180] = np.nan                    # 60% 缺失
    X[7, 15:] *= 0.05                      # 低文库（保留常量特征块）
    X[11, 15:] *= 6                        # 离群（保留常量特征块）
    groups = ["control"] * 20 + ["salt"] * 20
    batches = ["B1"] * 20 + ["B2"] * 20     # 与分组完全混杂
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "matrix.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sample_id"] + feats)
        for i in range(n):
            w.writerow([ids[i]] + ["" if np.isnan(v) else f"{v:.3f}" for v in X[i]])
    with open(os.path.join(outdir, "meta.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sample_id", "group", "batch"])
        for i in range(n):
            w.writerow([ids[i], groups[i], batches[i]])
    r = check(to_matrix(os.path.join(outdir, "matrix.csv")),
              read_meta(os.path.join(outdir, "meta.csv")))
    codes = {x["code"] for x in r["issues"]}
    expect = {"sample_missing", "low_library", "constant_features",
              "pca_outlier", "batch_confounded"}
    missing = expect - codes
    assert not missing, f"demo 未检出: {missing}"
    open(os.path.join(outdir, "qc_report.md"), "w", encoding="utf-8").write(render_report(r))
    open(os.path.join(outdir, "issues.json"), "w", encoding="utf-8").write(
        json.dumps(r, ensure_ascii=False, indent=2))
    print(f"DEMO OK: {len(r['issues'])} issues -> {outdir}/qc_report.md")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix"); ap.add_argument("--metadata")
    ap.add_argument("--second-matrix"); ap.add_argument("--orientation", default="samples_rows")
    ap.add_argument("--group-col", default="group"); ap.add_argument("--batch-col", default="batch")
    ap.add_argument("--outdir", default="results/qc"); ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        demo(os.path.join(a.outdir, "demo")); return
    if not a.matrix:
        sys.exit("请提供 --matrix（或 --demo）")
    mat = to_matrix(a.matrix, a.orientation)
    meta = read_meta(a.metadata, a.group_col, a.batch_col) if a.metadata else {}
    second = to_matrix(a.second_matrix)[:2] if a.second_matrix else None
    r = check(mat, meta, second)
    os.makedirs(a.outdir, exist_ok=True)
    open(os.path.join(a.outdir, "qc_report.md"), "w", encoding="utf-8").write(render_report(r))
    open(os.path.join(a.outdir, "issues.json"), "w", encoding="utf-8").write(
        json.dumps(r, ensure_ascii=False, indent=2))
    with open(os.path.join(a.outdir, "sample_flags.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sample_id", "missing_rate", "library"])
        for sid, v in r["sample_flags"].items():
            w.writerow([sid, v["missing"], v["library"]])
    print(f"issues={len(r['issues'])} -> {a.outdir}/qc_report.md")

if __name__ == "__main__":
    main()
