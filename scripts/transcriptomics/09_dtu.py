#!/usr/bin/env python3
###############################################################################
# 09_dtu.py -- differential transcript usage (DTU) from isoform counts
#
# 输入：
#   --counts CSV: sample_id, isoform1, isoform2, ...（counts）
#   --info   CSV: isoform_id, gene_id
#   --groups CSV: sample_id, group（两组）
# 方法：
#   每个基因聚合同组 isoform counts -> 卡方检验（isoform x group）
#   -> BH；报告主 isoform 比例变化。适用于 bulk/空间 isoform 数据。
#
# 用法：
#   python scripts/transcriptomics/09_dtu.py --counts iso_counts.csv \
#       --info iso_info.csv --groups groups.csv --outdir results/transcriptomics/09_dtu
#   python scripts/transcriptomics/09_dtu.py --demo --outdir results/transcriptomics/09_dtu_demo
###############################################################################
import argparse, os, sys, datetime
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

def bh(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(ranked, 1)
    return out

def make_demo(outdir, n=6):
    rng = np.random.default_rng(42)
    genes = ["G1", "G2", "G3", "G4"]
    isos = {g: [f"{g}.{i}" for i in range(1, 4)] for g in genes}
    samples = [f"S{i}" for i in range(n)]
    groups = ["A"] * (n // 2) + ["B"] * (n - n // 2)
    rows = []
    for s, grp in zip(samples, groups):
        row = {"sample_id": s}
        for g in genes:
            base = rng.integers(50, 100, 3)
            if g == "G2" and grp == "B":
                base = np.array([10, 80, 10])  # DTU
            row.update({iso: int(v) for iso, v in zip(isos[g], base)})
        rows.append(row)
    counts = pd.DataFrame(rows)
    info = pd.DataFrame([[iso, g] for g in genes for iso in isos[g]], columns=["isoform_id", "gene_id"])
    grp = pd.DataFrame({"sample_id": samples, "group": groups})
    counts.to_csv(os.path.join(outdir, "demo_counts.csv"), index=False)
    info.to_csv(os.path.join(outdir, "demo_info.csv"), index=False)
    grp.to_csv(os.path.join(outdir, "demo_groups.csv"), index=False)
    return (os.path.join(outdir, "demo_counts.csv"), os.path.join(outdir, "demo_info.csv"),
            os.path.join(outdir, "demo_groups.csv"), "G2")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts"); ap.add_argument("--info"); ap.add_argument("--groups")
    ap.add_argument("--min-count", type=int, default=50)
    ap.add_argument("--outdir", default="results/transcriptomics/09_dtu")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    truth = None
    if a.demo:
        a.counts, a.info, a.groups, truth = make_demo(a.outdir)
    if not a.counts or not a.info or not a.groups:
        sys.exit("Need --counts, --info, --groups, or --demo")
    cnt = pd.read_csv(a.counts); info = pd.read_csv(a.info); grp = pd.read_csv(a.groups)
    if "sample_id" not in cnt.columns or "sample_id" not in grp.columns: sys.exit("need sample_id")
    info = info.drop_duplicates("isoform_id")
    grp = grp.set_index("sample_id").loc[cnt["sample_id"]].reset_index()
    groups = sorted(grp["group"].unique())
    if len(groups) != 2: sys.exit("DTU currently supports exactly two groups")
    iso_cols = [c for c in cnt.columns if c != "sample_id" and c in set(info["isoform_id"])]
    if len(iso_cols) < 2: sys.exit("Fewer than 2 isoforms present")
    iso2gene = dict(zip(info["isoform_id"], info["gene_id"]))
    rows, prop_rows = [], []
    for g, sub in info[info["isoform_id"].isin(iso_cols)].groupby("gene_id"):
        isos = [i for i in sub["isoform_id"] if i in iso_cols]
        if len(isos) < 2: continue
        tab = []
        for gr in groups:
            idx = grp["group"] == gr
            tab.append(cnt.loc[idx, isos].sum(0).to_numpy(float))
        tab = np.vstack(tab)

        if tab.sum() < a.min_count: continue
        try:
            chi2, p, dof, exp = chi2_contingency(tab + 0.5)
        except Exception:
            continue
        props = tab / tab.sum(1, keepdims=True)
        # major isoform shift
        top = int(np.argmax(tab.sum(0)))
        delta = props[1, top] - props[0, top]
        rows.append([g, len(isos), int(tab.sum()), float(chi2), float(p), isos[top], float(props[0, top]), float(props[1, top]), float(delta)])
        for i, iso in enumerate(isos):
            prop_rows.append([g, iso] + list(props[:, i]))
    if not rows: sys.exit("No gene passed filters")
    res = pd.DataFrame(rows, columns=["gene_id", "n_isoforms", "total_count", "chi2", "p",
                                      "major_isoform", "prop_group1", "prop_group2", "delta_prop"])
    res["padj"] = bh(res["p"].to_numpy())
    res = res.sort_values("p")
    res.to_csv(os.path.join(a.outdir, "dtu_results.csv"), index=False)
    prop = pd.DataFrame(prop_rows, columns=["gene_id", "isoform"] + [f"prop_{g}" for g in groups])
    prop.to_csv(os.path.join(a.outdir, "isoform_proportions.csv"), index=False)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"samples={len(cnt)}\ngenes_tested={len(res)}\ngroups={groups}\n")
        f.write("test=chi-square on isoform x group counts; BH\n")
    print(f"[dtu] genes tested={len(res)}; FDR<0.05={int((res['padj']<0.05).sum())}")
    print(res.head(5).to_string(index=False))
    if truth: print(f"[dtu] demo causal gene {truth} significant: {bool(res.loc[res.gene_id==truth,'padj'].iloc[0] < 0.05)}")
    print(f"[dtu] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
