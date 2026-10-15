"""07_grn_boost.py -- TF regulatory network inference (GRNBoost2-style).

Title:    Rice-skills gene regulatory network inference (v3.2 step 07)
Description:
    Infers candidate TF -> target regulatory edges from an expression matrix
    (samples x genes) using stochastic-gradient-boosting regression with
    early stopping -- the same algorithmic core as GRNBoost2 (Arboreto),
    implemented on scikit-learn so no Java/arboreto install is needed.
    Complements module 10's WGCNA coexpression modules: WGCNA gives
    co-expression clusters; this ranks DIRECTIONAL regulator-target edges.

    For each target gene:
      1. candidate regulators = TFs provided via --tfs (recommended: the
         PlantTFDB rice TF list shipped in assets/annotations/) minus the
         target itself; with --all-genes every gene is a candidate;
      2. fit GradientBoostingRegressor (depth-1 stumps, early stopping)
         predicting the target from candidate expression;
      3. record feature importances as edge scores.

    Output is a ranked edge list; --top-edges keeps the strongest K edges
    per target. Importance is an association score, NOT proof of direct
    binding or causality (same caveat as WGCNA/GENIE3/GRNBoost2).

Usage:
    python scripts/integration/07_grn_boost.py \
        --expr demo_data/transcriptomics/counts_norm.csv \
        --tfs assets/annotations/rice_tfs_planttfdb.csv \
        --top-edges 5 --outdir results/integration/07_grn
    # self-test with synthetic modules (no inputs)
    python scripts/integration/07_grn_boost.py --demo --outdir results/integration/07_grn_demo

Inputs:
    --expr   CSV, first column = sample_id, remaining columns = genes
    --tfs    optional CSV with a gene-id column (auto-detected name or
             --tf-column); plain TXT (one id per line) also accepted

Outputs:
    grn_edges.csv    regulator, target, importance (all edges, ranked)
    top_edges.csv    strongest --top-edges edges per target
    grn_summary.txt  parameters, run stats, demo self-test result

Notes:
    - Requires scikit-learn (system Python 3.14 env). Depth-1 stumps +
      early stopping keep this fast on ~1e4 genes x ~1e2 samples; for
      larger matrices subsample targets via --max-targets.
    - Thearboreto compatibility notes (PyPI 0.1.6 vs Dask) live in
      third_party/kdense-scientific-skills/skills/arboreto/ if you later
      want the upstream implementation.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor


# ------------------------------------------------------------------- io
def load_expr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    return df


def load_tfs(path: Path | None, tf_column: str | None) -> list[str] | None:
    if path is None:
        return None
    p = Path(path)
    if p.suffix.lower() in (".txt", ".tsv", ".list"):
        sep = "\t" if p.suffix.lower() == ".tsv" else None
        df = pd.read_csv(p, sep=sep or r"\s+", header=None, comment="#")
        return [str(x) for x in df.iloc[:, 0].tolist()]
    df = pd.read_csv(p)
    col = tf_column
    if col is None:
        for cand in ("gene_id", "Gene_id", "gene", "TF", "tf", "Gene_ID", "MSU_ID"):
            if cand in df.columns:
                col = cand
                break
        if col is None:
            col = df.columns[0]
    return [str(x) for x in df[col].dropna().astype(str).tolist()]


# ----------------------------------------------------------------- core
def infer_edges(expr: pd.DataFrame, tfs: list[str] | None, top_edges: int,
                max_targets: int | None, seed: int, verbose: bool = True) -> pd.DataFrame:
    genes = list(expr.columns)
    cand = sorted(set(tfs)) if tfs else genes
    targets = genes if max_targets is None else genes[:max_targets]
    edges: list[dict] = []
    for k, target in enumerate(targets, 1):
        regs = [g for g in cand if g != target]
        if len(regs) < 2:
            continue
        X = expr[regs].to_numpy(dtype=float)
        y = expr[target].to_numpy(dtype=float)
        gbr = GradientBoostingRegressor(
            n_estimators=250, learning_rate=0.1, max_depth=1,
            subsample=0.632, random_state=seed,
        )
        gbr.fit(X, y)
        imp = gbr.feature_importances_
        nz = [(regs[i], float(imp[i])) for i in np.argsort(-imp) if imp[i] > 0]
        for reg, score in nz[:top_edges]:
            edges.append({"regulator": reg, "target": target, "importance": score})
        if verbose and k % 10 == 0:
            print(f"[07] target {k}/{len(targets)} done", flush=True)
    return pd.DataFrame(edges, columns=["regulator", "target", "importance"])


# ----------------------------------------------------------------- demo
def demo_expr() -> pd.DataFrame:
    rng = np.random.default_rng(42)
    n = 60
    tf1 = rng.normal(size=n)
    tf2 = rng.normal(size=n)
    tf3 = rng.normal(size=n)
    data = {
        "TF1": tf1, "TF2": tf2, "TF3": tf3,
        "TF4": rng.normal(size=n),                       # decoy TF, no targets
        "TGT_A": 1.2 * tf1 + 0.3 * rng.normal(size=n),   # TF1 -> TGT_A
        "TGT_B": -0.9 * tf2 + 0.3 * rng.normal(size=n),  # TF2 -> TGT_B
        "TGT_C": 0.8 * tf3 + 0.3 * rng.normal(size=n),   # TF3 -> TGT_C
        "BG_1": rng.normal(size=n), "BG_2": rng.normal(size=n),
    }
    return pd.DataFrame(data, index=[f"S{i:02d}" for i in range(n)])


# ----------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--expr", default=None)
    ap.add_argument("--tfs", default=None)
    ap.add_argument("--tf-column", default=None)
    ap.add_argument("--all-genes", action="store_true", help="every gene as candidate regulator")
    ap.add_argument("--top-edges", type=int, default=5)
    ap.add_argument("--max-targets", type=int, default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="results/integration/07_grn")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        expr = demo_expr()
        tfs = ["TF1", "TF2", "TF3", "TF4"]
        source = "demo (synthetic modules)"
    else:
        if not args.expr:
            ap.error("--expr required (or --demo)")
        expr = load_expr(Path(args.expr))
        tfs = None if args.all_genes else load_tfs(Path(args.tfs) if args.tfs else None, args.tf_column)
        source = f"file:{args.expr}"

    edges = infer_edges(expr, tfs, args.top_edges, args.max_targets, args.seed)
    edges = edges.sort_values("importance", ascending=False).reset_index(drop=True)
    edges.to_csv(outdir / "grn_edges.csv", index=False)
    top = edges.groupby("target", as_index=False).head(args.top_edges)
    top = top.sort_values(["target", "importance"], ascending=[True, False])
    top.to_csv(outdir / "top_edges.csv", index=False)

    stats = [
        f"source: {source}",
        f"timestamp: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"samples: {expr.shape[0]}; genes: {expr.shape[1]}; candidate regulators: "
        f"{len(tfs) if tfs else expr.shape[1]}",
        f"edges kept (top {args.top_edges}/target): {len(edges)}",
    ]
    if len(edges):
        stats.append(f"strongest edge: {edges.iloc[0]['regulator']} -> "
                     f"{edges.iloc[0]['target']} (importance {edges.iloc[0]['importance']:.3f})")
    (outdir / "grn_summary.txt").write_text("\n".join(stats) + "\n", encoding="utf-8")
    print("\n".join("[07] " + s for s in stats))
    print(f"[07] wrote {outdir / 'grn_edges.csv'}")

    if args.demo:
        # per-target strongest edge (groupby head keeps 5 rows per target;
        # must take the argmax, not the last row)
        best = top.loc[top.groupby("target")["importance"].idxmax()]
        top_by_target = best.set_index("target")["regulator"].to_dict()
        assert top_by_target.get("TGT_A") == "TF1", f"TF1->TGT_A missing: {top_by_target}"
        assert top_by_target.get("TGT_B") == "TF2", f"TF2->TGT_B missing: {top_by_target}"
        assert top_by_target.get("TGT_C") == "TF3", f"TF3->TGT_C missing: {top_by_target}"
        assert not best.loc[best["target"].str.startswith("TGT"), "regulator"].eq("TF4").any(), \
            "decoy TF4 is the strongest regulator of a real target"
        print("[07] demo self-test PASSED (TF1/TF2/TF3 -> their targets; decoy TF4 absent)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
