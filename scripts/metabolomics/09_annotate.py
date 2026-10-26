#!/usr/bin/env python3
###############################################################################
# 09_annotate.py -- merge SIRIUS and GNPS metabolite annotations
#
# 输入：
#   --features CSV: feature_id, mz, rt
#   --sirius   CSV: feature_id, formula, adduct, score（SIRIUS 输出整理后）
#   --gnps     CSV: feature_id, library, score（GNPS 库匹配）
# 输出：
#   annotated_features.csv：每个 feature 的最佳 SIRIUS / GNPS 注释
#
# 说明：本脚本只做"结果合并与筛选"，不运行 SIRIUS/GNPS；两者需自行安装/上传。
#
# 用法：
#   python scripts/metabolomics/09_annotate.py --features features.csv \
#       --sirius sirius.csv --gnps gnps.csv --outdir results/metabolomics/09_annotate
#   python scripts/metabolomics/09_annotate.py --demo --outdir results/metabolomics/09_annotate_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd

def make_demo(outdir, n=100):
    rng = np.random.default_rng(42)
    feat = pd.DataFrame({"feature_id": [f"F{i:03d}" for i in range(n)],
                         "mz": rng.uniform(100, 1000, n).round(4),
                         "rt": rng.uniform(0.5, 20, n).round(2)})
    sir = feat.sample(30, random_state=1).copy()
    sir["formula"] = [f"C{int(rng.integers(5,30))}H{int(rng.integers(6,50))}O{int(rng.integers(1,10))}" for _ in range(len(sir))]
    sir["adduct"] = rng.choice(["[M+H]+", "[M-H]-", "[M+Na]+"], len(sir))
    sir["score"] = rng.uniform(0.4, 0.99, len(sir)).round(3)
    gnps = feat.sample(20, random_state=2).copy()
    gnps["library"] = [f"GNPS_{i}" for i in range(len(gnps))]
    gnps["score"] = rng.uniform(0.3, 0.98, len(gnps)).round(3)
    feat.to_csv(os.path.join(outdir, "demo_features.csv"), index=False)
    sir[["feature_id", "formula", "adduct", "score"]].to_csv(os.path.join(outdir, "demo_sirius.csv"), index=False)
    gnps[["feature_id", "library", "score"]].to_csv(os.path.join(outdir, "demo_gnps.csv"), index=False)
    return [os.path.join(outdir, f) for f in ["demo_features.csv", "demo_sirius.csv", "demo_gnps.csv"]]

def best(df, prefix):
    if df is None or not len(df): return None
    df = df.sort_values("score", ascending=False).drop_duplicates("feature_id")
    return df.rename(columns={c: f"{prefix}_{c}" for c in df.columns if c != "feature_id"})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features"); ap.add_argument("--sirius"); ap.add_argument("--gnps")
    ap.add_argument("--min-sirius", type=float, default=0.5)
    ap.add_argument("--min-gnps", type=float, default=0.3)
    ap.add_argument("--outdir", default="results/metabolomics/09_annotate")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.features, a.sirius, a.gnps = make_demo(a.outdir)
    if not a.features:
        sys.exit("Need --features, or --demo")
    feat = pd.read_csv(a.features)
    if "feature_id" not in feat.columns: sys.exit("features need feature_id")
    sir = pd.read_csv(a.sirius) if a.sirius and os.path.exists(a.sirius) else None
    gnp = pd.read_csv(a.gnps) if a.gnps and os.path.exists(a.gnps) else None
    if sir is not None:
        sir = sir[sir["score"] >= a.min_sirius]
        sir = best(sir, "sirius")
    if gnp is not None:
        gnp = gnp[gnp["score"] >= a.min_gnps]
        gnp = best(gnp, "gnps")
    out = feat.copy()
    if sir is not None: out = out.merge(sir, on="feature_id", how="left")
    if gnp is not None: out = out.merge(gnp, on="feature_id", how="left")
    out["annotated"] = False
    if "sirius_formula" in out.columns:
        out.loc[out["sirius_formula"].notna(), "annotated"] = True
    if "gnps_library" in out.columns:
        out.loc[out["gnps_library"].notna(), "annotated"] = True
    out.to_csv(os.path.join(a.outdir, "annotated_features.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"features={len(feat)}\nsirius_hits={0 if sir is None else len(sir)}\n")
        f.write(f"gnps_hits={0 if gnp is None else len(gnp)}\n")
        f.write(f"min_sirius={a.min_sirius}\nmin_gnps={a.min_gnps}\n")
    print(f"[annot] features={len(feat)} annotated={int(out['annotated'].sum())} -> {a.outdir}")

if __name__ == "__main__":
    main()
