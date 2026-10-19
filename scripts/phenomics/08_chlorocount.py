#!/usr/bin/env python3
###############################################################################
# 08_chlorocount.py -- chloroplast / bundle-sheath cell segmentation & counting
#
# Title:    Rice-skills chloroplast occupancy quantification (Chloro-Count)
# Description:
#   Wraps Chloro-Count (New Phytologist 2025, PMID 39668515), a deep-learning
#   tool that segments, counts and measures chloroplasts and bundle-sheath
#   cells from confocal (.lif / image) stacks. Used in rice hap3H mutants to
#   quantify chloroplast occupancy in bundle-sheath cells.
#
#   Upstream: github.com/pedropgusmao/chloro-count -- BSD-3-Clause. Not
#   vendored; point --repo at your clone. Pre-trained models are distributed
#   by the authors (Google Drive link in the upstream README).
#
# Usage:
#   python scripts/phenomics/08_chlorocount.py \
#       --repo /path/to/chloro-count --input sample.lif \
#       --models /path/to/models --channel 1 \
#       --outdir results/phenomics/08_chlorocount
#
#   python scripts/phenomics/08_chlorocount.py --demo
#
# Inputs:
#   --repo     local clone of chloro-count (BSD-3-Clause)
#   --input    .lif file or image stack
#   --models   directory containing the two pre-trained .pth models
#   --channel  image channel index (upstream note: adapt to your acquisition)
#   --outdir   output directory
#   --demo     fixture + dependency check (SKIP when torch/models absent)
#
# Outputs:
#   segmentations/, counts.csv, measurements.csv, run_info.txt
#
# Notes:
#   - Requires PyTorch + the authors' pre-trained models; the models are NOT
#     redistributed here.
#   - BSD-3-Clause: retain LICENSE + NOTICE if redistributing.
###############################################################################
import argparse, csv, os, subprocess, sys

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--repo"); p.add_argument("--input"); p.add_argument("--models")
    p.add_argument("--channel", default="1")
    p.add_argument("--outdir", default="results/phenomics/08_chlorocount")
    p.add_argument("--demo", action="store_true")
    return p.parse_args()

def torch_ok():
    try:
        import torch; return torch.__version__
    except Exception:
        return None

def main():
    a = parse_args(); os.makedirs(a.outdir, exist_ok=True)
    t = torch_ok()
    if a.demo:
        fix = os.path.join(a.outdir, "_demo_fixture"); os.makedirs(fix, exist_ok=True)
        # 2-channel synthetic stack (chloroplast signal + background)
        from PIL import Image
        import random
        random.seed(42)
        for ch in (0, 1):
            im = Image.new("L", (256, 256), 10 if ch == 0 else 40)
            px = im.load()
            for _ in range(120):
                x, y = random.randint(10, 245), random.randint(10, 245)
                r = random.randint(2, 5)
                for dx in range(-r, r + 1):
                    for dy in range(-r, r + 1):
                        if dx * dx + dy * dy <= r * r and 0 <= x+dx < 256 and 0 <= y+dy < 256:
                            px[x+dx, y+dy] = 230 if ch else 90
            im.save(os.path.join(fix, f"demo_ch{ch}.png"))
        print("== DEMO: synthetic chloroplast stack ==")
        print(f"  fixture: {fix} (2 channels, 120 bright organelles)")
        print(f"  torch: {t or 'MISSING'}")
        open(os.path.join(a.outdir, "run_info.txt"), "w").write(
            f"status: DEMO\ntorch: {t or 'MISSING'}\nfixture: {fix}\n")
        if not t:
            print("SKIP: PyTorch not installed.")
            print("  pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu")
            print("  Models: download the two .pth files from the upstream README into models/")
        else:
            print("PyTorch present; run with the authors' models:")
            print(f"  python {a.repo or '<clone>'}/segment.py --input <file> --models <dir>")
        return
    if not (a.repo and a.input and a.models):
        sys.exit("Required: --repo --input --models (or --demo)")
    if not os.path.isdir(a.repo): sys.exit(f"repo not found: {a.repo}")
    if not t: sys.exit("PyTorch is required (pip install torch)")
    models = [f for f in os.listdir(a.models) if f.endswith((".pth", ".pt"))]
    if not models:
        sys.exit(f"No .pth/.pt models found in {a.models} (upstream README has the download link)")
    seg = os.path.join(a.repo, "segment.py")
    if not os.path.isfile(seg): sys.exit(f"segment.py not found in {a.repo}")
    out_seg = os.path.join(a.outdir, "segmentations"); os.makedirs(out_seg, exist_ok=True)
    subprocess.run([sys.executable, seg, "--input", a.input, "--models", a.models,
                    "--output", out_seg, "--channel", str(a.channel)],
                   check=True, cwd=a.repo)
    assoc = os.path.join(a.repo, "associate_and_measure.py")
    if os.path.isfile(assoc):
        subprocess.run([sys.executable, assoc, "--input", out_seg,
                        "--output", a.outdir], check=True, cwd=a.repo)
    open(os.path.join(a.outdir, "run_info.txt"), "w").write("\n".join([
        "script: scripts/phenomics/08_chlorocount.py",
        f"repo: {a.repo}", f"input: {a.input}", f"models: {a.models}",
        f"channel: {a.channel}", "status: ok"]) + "\n")
    print(f"done -> {a.outdir}")

if __name__ == "__main__":
    main()
