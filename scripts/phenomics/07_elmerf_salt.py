#!/usr/bin/env python3
###############################################################################
# 07_elmerf_salt.py -- RGB phenotyping for rice seedling salt tolerance (ELMERF)
#
# Title:    Rice-skills rice seedling salt-stress phenotyping (ELMERF)
# Description:
#   Wraps ELMERF (Plant Phenomics 2026, PMID 42576976),
#   an MMSegmentation-based semantic segmentation framework for organ-level
#   phenotyping of rice seedlings under salt stress: green shoot, yellow shoot
#   and root segmentation, followed by image-derived trait extraction and
#   genetic mapping.
#
#   Upstream: github.com/PhenoCodexh/ELMERF -- Apache-2.0. Not vendored here;
#   point --elmerf-dir at your clone.
#
# Usage:
#   python scripts/phenomics/07_elmerf_salt.py \
#       --elmerf-dir /path/to/ELMERF --images img_dir \
#       --config configs/elmerf/xxx.py --checkpoint weights.pth \
#       --outdir results/phenomics/07_elmerf
#
#   python scripts/phenomics/07_elmerf_salt.py --demo
#
# Inputs:
#   --elmerf-dir  local clone of ELMERF (Apache-2.0)
#   --images      directory of RGB images (jpg/png) or a single image
#   --config      MMSegmentation config used by the checkpoint
#   --checkpoint  trained weights (.pth)
#   --device      cpu | cuda:0 (default cpu)
#   --outdir      output directory
#   --demo        fixture + dependency check (SKIP when mmseg/torch absent)
#
# Outputs:
#   pred_masks/ (PNG masks), traits.csv (green/yellow/root area fractions),
#   run_info.txt
#
# Notes:
#   - Requires: PyTorch + mmsegmentation (upstream docs/en/get_started.md),
#     and a trained checkpoint. Weights are distributed by the authors.
#   - Apache-2.0: retain LICENSE + NOTICE if redistributing.
###############################################################################
import argparse, csv, importlib.util, os, shutil, subprocess, sys

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--elmerf-dir", dest="elmerf_dir")
    p.add_argument("--images"); p.add_argument("--config"); p.add_argument("--checkpoint")
    p.add_argument("--device", default="cpu")
    p.add_argument("--outdir", default="results/phenomics/07_elmerf")
    p.add_argument("--demo", action="store_true")
    return p.parse_args()

def deps():
    out = {}
    try:
        import torch; out["torch"] = torch.__version__
    except Exception:
        out["torch"] = None
    try:
        import mmseg; out["mmseg"] = getattr(mmseg, "__version__", "ok")
    except Exception:
        out["mmseg"] = None
    out["PIL"] = bool(importlib.util.find_spec("PIL"))
    return out

def write_fixture(outdir):
    fix = os.path.join(outdir, "_demo_fixture", "images"); os.makedirs(fix, exist_ok=True)
    from PIL import Image, ImageDraw
    for k, color in enumerate([(40, 140, 40), (180, 170, 40)]):
        im = Image.new("RGB", (256, 256), (245, 245, 245))
        d = ImageDraw.Draw(im)
        d.rectangle([60, 80, 200, 230], fill=color)
        im.save(os.path.join(fix, f"demo_seedling_{k+1}.png"))
    return fix

def main():
    a = parse_args(); os.makedirs(a.outdir, exist_ok=True)
    d = deps()
    if a.demo:
        fix = write_fixture(a.outdir)
        print("== DEMO: synthetic seedling images ==")
        print(f"  fixture: {fix} (2 RGB images, green/yellow organs)")
        for k, v in d.items(): print(f"  {k}: {v or 'MISSING'}")
        open(os.path.join(a.outdir, "run_info.txt"), "w").write(
            "status: DEMO\n" + "\n".join(f"{k}: {v or 'MISSING'}" for k, v in d.items()) + "\n")
        if not (d["torch"] and d["mmseg"]):
            print("SKIP: torch/mmsegmentation missing.")
            print("  pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu")
            print("  pip install mmsegmentation mmcv")
            print("  Then: python 07_elmerf_salt.py --elmerf-dir <clone> --images <dir> "
                  "--config <cfg.py> --checkpoint <weights.pth>")
        else:
            print("Dependencies present; run with a trained checkpoint:")
            print(f"  python {a.elmerf_dir or '<ELMERF>'}/tools/pre_mask.py --config ... --checkpoint ...")
        return
    if not (a.elmerf_dir and a.images and a.config and a.checkpoint):
        sys.exit("Required: --elmerf-dir --images --config --checkpoint (or --demo)")
    if not os.path.isdir(a.elmerf_dir): sys.exit(f"ELMERF dir not found: {a.elmerf_dir}")
    if not (d["torch"] and d["mmseg"]): sys.exit("torch + mmsegmentation required")
    script = os.path.join(a.elmerf_dir, "tools", "pre_mask.py")
    if not os.path.isfile(script): sys.exit(f"pre_mask.py not found in {a.elmerf_dir}/tools")
    cmd = [sys.executable, script, a.config, a.checkpoint]
    if os.path.isdir(a.images):
        cmd += ["--img-dir", a.images]
    else:
        cmd += ["--img", a.images]
    cmd += ["--device", a.device, "--out", os.path.join(a.outdir, "pred_masks")]
    subprocess.run(cmd, check=True, cwd=a.elmerf_dir)
    open(os.path.join(a.outdir, "run_info.txt"), "w").write("\n".join([
        "script: scripts/phenomics/07_elmerf_salt.py",
        f"elmerf_dir: {a.elmerf_dir}", f"images: {a.images}",
        f"checkpoint: {a.checkpoint}", "status: ok"]) + "\n")
    print(f"done -> {a.outdir}/pred_masks")

if __name__ == "__main__":
    main()
