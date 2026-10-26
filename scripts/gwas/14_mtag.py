#!/usr/bin/env python3
###############################################################################
# 14_mtag.py -- MTAG multi-trait GWAS plan/runner wrapper
#
# 说明：MTAG 需要 LD score 文件与各性状 sumstats（SNP, Z, N）；
#       本脚本校验输入并生成 mtag 命令，不重新实现 MTAG。
#
# 用法：
#   python scripts/gwas/14_mtag.py --sumstats t1.txt,t2.txt --ld eur_w_ld_chr/ --outdir results/gwas/14_mtag
#   python scripts/gwas/14_mtag.py --demo --outdir results/gwas/14_mtag_demo
###############################################################################
import argparse, os, sys, datetime

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sumstats"); ap.add_argument("--ld")
    ap.add_argument("--outdir", default="results/gwas/14_mtag")
    ap.add_argument("--n-min", type=int, default=1000)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.sumstats = os.path.join(a.outdir, "demo_t1.txt") + "," + os.path.join(a.outdir, "demo_t2.txt")
        a.ld = os.path.join(a.outdir, "ld/")
        for p in a.sumstats.split(","):
            open(p, "w").write("SNP\tZ\tN\nrs1\t1.2\t5000\nrs2\t-0.8\t5000\n")
        os.makedirs(a.ld, exist_ok=True)
    if not a.sumstats or not a.ld: sys.exit("Need --sumstats and --ld, or --demo")
    files = [x.strip() for x in a.sumstats.split(",") if x.strip()]
    for f in files:
        if not os.path.exists(f): sys.exit(f"sumstats not found: {f}")
    plan = [
        "# MTAG plan (requires mtag.py from https://github.com/omeed-maghzian/mtag)",
        "python mtag.py \\",
        f"  --sumstats {','.join(files)} \\",
        f"  --ld {a.ld} \\",
        f"  --n_min {a.n_min} \\",
        f"  --out {os.path.join(a.outdir, 'mtag_results')}",
        "# outputs: mtag_results_trait*.txt (contains MTAG Z/P)",
    ]
    open(os.path.join(a.outdir, "mtag_plan.sh"), "w").write("\n".join(plan) + "\n")
    open(os.path.join(a.outdir, "run_info.txt"), "w").write(
        f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n"
        f"traits={len(files)}\nn_min={a.n_min}\nmode=plan only\n")
    print(f"[mtag] traits={len(files)} -> {a.outdir}/mtag_plan.sh (plan only)")
    print("       sumstats must contain SNP, Z, N; LD score prefix required")

if __name__ == "__main__":
    main()
