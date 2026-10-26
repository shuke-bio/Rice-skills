#!/usr/bin/env python3
###############################################################################
# benchmark_suite.py -- run demo-based regression benchmarks
#
# 用途：发版前/CI 中验证关键脚本仍可运行且结果方向正确。
#   --quick : 3 个快速用例（脱靶、TWAS、MR）
#   --all   : 全部 6 个用例（含 mGWAS、k-mer GWAS、DTU）
#
# 用法：
#   python scripts/utils/benchmark_suite.py --quick
#   python scripts/utils/benchmark_suite.py --all --outdir results/benchmark
###############################################################################
import argparse, os, sys, csv, time, subprocess, shutil, datetime, glob

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PY = sys.executable
R = os.environ.get("RSCRIPT") or shutil.which("Rscript")
if not R:
    # generic Windows install location (no user-specific path)
    _cand = sorted(glob.glob("C:/Program Files/R/R-*/bin/Rscript.exe"))
    R = _cand[-1] if _cand else "Rscript"

def py(script, out):
    return [PY, os.path.join(ROOT, "scripts", script), "--demo", "--outdir", out]

def r(script, out):
    return [R, os.path.join(ROOT, "scripts", script), "--demo", "--outdir", out]

TESTS = [
    ("offtarget", lambda o: py("design/04_offtarget.py", o), ["offtarget_summary.csv"]),
    ("twas", lambda o: r("qtl/06_twas.R", o), ["twas_results.csv"]),
    ("mr", lambda o: r("qtl/07_mr.R", o), ["mr_results.csv"]),
    ("mgwas", lambda o: r("metabolomics/08_mgwas.R", o), ["mgwas_results.csv"]),
    ("kmer", lambda o: py("gwas/12_kmer_gwas.py", o), ["kmer_gwas_results.csv"]),
    ("dtu", lambda o: py("transcriptomics/09_dtu.py", o), ["dtu_results.csv"]),
]

def check(metric, out, files):
    for f in files:
        if not os.path.exists(os.path.join(out, f)):
            return False, f"missing {f}"
    try:
        import pandas as pd
        if metric == "twas":
            d = pd.read_csv(os.path.join(out, "twas_results.csv"))
            return (d.iloc[0]["gene_id"] == "Gene07"), "top gene " + str(d.iloc[0]["gene_id"])
        if metric == "mr":
            d = pd.read_csv(os.path.join(out, "mr_results.csv"))
            b = float(d.loc[d["method"] == "IVW_fixed", "beta"].iloc[0])
            return (0.25 <= b <= 0.45), f"IVW beta {b:.3f}"
        if metric == "dtu":
            d = pd.read_csv(os.path.join(out, "dtu_results.csv"))
            g2 = d[d["gene_id"] == "G2"]
            return (len(g2) and float(g2["padj"].iloc[0]) < 0.05), "G2 padj"
        if metric == "mgwas":
            d = pd.read_csv(os.path.join(out, "mgwas_results.csv"))
            return (int((d["padj"] < 0.05).sum()) >= 2), "FDR hits"
        if metric == "kmer":
            d = pd.read_csv(os.path.join(out, "kmer_gwas_results.csv"))
            return (int((d["padj"] < 0.05).sum()) >= 1), "FDR hits"
        if metric == "offtarget":
            d = pd.read_csv(os.path.join(out, "offtarget_summary.csv"))
            return (int(d["mm0"].sum()) >= 3), "on-target guides"
    except Exception as e:
        return False, f"metric error: {e}"
    return True, "ok"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true"); ap.add_argument("--all", action="store_true")
    ap.add_argument("--outdir", default="results/benchmark")
    a = ap.parse_args()
    if not a.quick and not a.all: a.quick = True
    tests = TESTS[:3] if a.quick else TESTS
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    for name, cmd, files in tests:
        out = os.path.join(a.outdir, name)
        if os.path.isdir(out): shutil.rmtree(out)
        os.makedirs(out, exist_ok=True)
        t0 = time.time()
        try:
            p = subprocess.run(cmd(out), cwd=ROOT, capture_output=True, text=True, timeout=300)
            ok, detail = check(name, out, files) if p.returncode == 0 else (False, f"exit {p.returncode}")
        except Exception as e:
            ok, detail = False, str(e)
        dt = time.time() - t0
        rows.append([name, "PASS" if ok else "FAIL", f"{dt:.1f}", detail])
        print(f"[bench] {name:10s} {'PASS' if ok else 'FAIL'} {dt:6.1f}s  {detail}")
    report = os.path.join(a.outdir, "benchmark_report.csv")
    with open(report, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["test", "status", "seconds", "detail"]); w.writerows(rows)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"mode={'quick' if a.quick else 'all'}\ntests={len(rows)}\n")
    nfail = sum(1 for r in rows if r[1] != "PASS")
    print(f"[bench] {len(rows)-nfail}/{len(rows)} passed -> {report}")
    sys.exit(0 if nfail == 0 else 1)

if __name__ == "__main__":
    main()
