#!/usr/bin/env python3
###############################################################################
# 06_menet_gs.py -- multi-environment genomic prediction (MeNet)
#
# Title:    Rice-skills multi-environment genomic selection (MeNet)
# Description:
#   Wraps MeNet (Plant Communications 2025; doi:10.1016/j.xplc.2025.101620),
#   a mixed-effects deep neural network that models phenotype-specific genetic
#   relatedness as random effects and cumulative variant effects as fixed
#   effects. Validated on 12 rice traits x 3 environments; outperformed 11
#   baselines in 34/36 evaluations.
#
#   Upstream: github.com/ganlab/MENET -- ACADEMIC-USE licence only (not
#   MIT/Apache). This wrapper never bundles upstream code; it prepares data in
#   MeNet's layout and calls the user's local clone.
#
# Usage:
#   python scripts/phenomics/06_menet_gs.py \
#       --genotype data/geno.csv --phenotype data/pheno.csv \
#       --menet-dir /path/to/MENET --trait trait1 \
#       --outdir results/phenomics/06_menet
#
#   python scripts/phenomics/06_menet_gs.py --demo
#
# Inputs:
#   --genotype   sample x marker dosage CSV (0/1/2, header = marker IDs)
#   --phenotype  sample x trait CSV (first column = sample_id)
#   --menet-dir  local clone of github.com/ganlab/MENET (academic licence)
#   --trait      trait column to model (default: first trait)
#   --windows    0/1, MeNet chromosome-window mechanism (default 1)
#   --device     cpu | cuda (default cpu)
#   --outdir     output directory
#   --demo       synthetic fixture + dependency check (SKIP when torch absent)
#
# Outputs:
#   menet_data/{gene/genotype.csv, phen/phenotype.csv, splits/}
#   menet.log, prediction.csv (if upstream writes predictions), run_info.txt
#
# Notes:
#   - Requires PyTorch (upstream pins torch 2.0 / CUDA 11.8; CPU works for
#     small data). Install: pip install torch --index-url .../cpu
#   - Academic-use licence: cite Li et al. 2025 Plant Communications.
###############################################################################
import argparse, csv, json, os, random, shutil, subprocess, sys

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--genotype"); p.add_argument("--phenotype")
    p.add_argument("--menet-dir", dest="menet_dir"); p.add_argument("--trait")
    p.add_argument("--windows", default="1"); p.add_argument("--device", default="cpu")
    p.add_argument("--outdir", default="results/phenomics/06_menet")
    p.add_argument("--demo", action="store_true")
    return p.parse_args()

def torch_ok():
    try:
        import torch  # noqa
        return torch.__version__
    except Exception:
        return None

def write_demo(outdir):
    fix = os.path.join(outdir, "_demo_fixture")
    os.makedirs(os.path.join(fix, "gene"), exist_ok=True)
    os.makedirs(os.path.join(fix, "phen"), exist_ok=True)
    random.seed(42)
    n, m = 60, 200
    ids = [f"Ind{i:03d}" for i in range(n)]
    markers = [f"S{j}" for j in range(m)]
    with open(os.path.join(fix, "gene", "genotype.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sample_id"] + markers)
        for i in range(n):
            w.writerow([ids[i]] + [random.choice([0, 1, 2]) for _ in range(m)])
    with open(os.path.join(fix, "phen", "phenotype.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sample_id", "trait1"])
        for i in range(n):
            w.writerow([ids[i], round(random.gauss(100, 10), 3)])
    return fix

def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    tver = torch_ok()
    report = [f"torch: {tver or 'MISSING'}"]
    if a.demo:
        fix = write_demo(a.outdir)
        report.append(f"fixture: {fix} (60 samples x 200 markers, 1 trait)")
        open(os.path.join(a.outdir, "run_info.txt"), "w").write(
            "\n".join(["status: DEMO", *report]) + "\n")
        print("== DEMO: synthetic MeNet data prepared ==")
        for r in report: print(" ", r)
        if not tver:
            print("SKIP: PyTorch not installed.")
            print("  pip install torch --index-url https://download.pytorch.org/whl/cpu")
            print("  (upstream pins torch 2.0/CUDA 11.8; CPU works for small data)")
        else:
            print("PyTorch present; run MeNet from your clone:")
            print(f"  python {a.menet_dir or '<MENET>'}/menet.py -p trait1 -d {a.device}")
        return
    if not (a.genotype and a.phenotype and a.menet_dir):
        sys.exit("Required: --genotype --phenotype --menet-dir (or --demo)")
    if not os.path.isdir(a.menet_dir):
        sys.exit(f"MeNet directory not found: {a.menet_dir}")
    if not tver:
        sys.exit("PyTorch is required: pip install torch (CPU wheel is fine)")
    data = os.path.join(a.outdir, "menet_data")
    os.makedirs(os.path.join(data, "gene"), exist_ok=True)
    os.makedirs(os.path.join(data, "phen"), exist_ok=True)
    shutil.copy(a.genotype, os.path.join(data, "gene", "genotype.csv"))
    shutil.copy(a.phenotype, os.path.join(data, "phen", "phenotype.csv"))
    os.makedirs(os.path.join(data, "splits"), exist_ok=True)
    log = open(os.path.join(a.outdir, "menet.log"), "w")
    try:
        subprocess.run([sys.executable, "utils/process_data/split.py"], cwd=a.menet_dir,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, "menet.py", "-p", a.trait or "trait1",
                        "-w", a.windows, "-d", a.device], cwd=a.menet_dir,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as e:
        log.close()
        sys.exit(f"MeNet failed (see menet.log): {e}")
    log.close()
    open(os.path.join(a.outdir, "run_info.txt"), "w").write("\n".join([
        "script: scripts/phenomics/06_menet_gs.py",
        f"menet_dir: {a.menet_dir}", f"trait: {a.trait}", f"device: {a.device}",
        f"torch: {tver}", "status: ok"]) + "\n")
    print(f"done -> {a.outdir}")

if __name__ == "__main__":
    main()
