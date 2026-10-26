#!/usr/bin/env python3
###############################################################################
# build_skill_import.py -- build a platform-importable SKILL.md (description <=1024)
#
# 背景：Agent Skills 规范对 description 有 1024 字符上限；GitHub 主文件可保留完整
# 描述，本脚本只生成“导入版”，不改动仓库中的 SKILL.md。
#
# 用法：
#   python scripts/utils/build_skill_import.py
#   python scripts/utils/build_skill_import.py --max-desc 1024 --outdir build/skill_import
###############################################################################
import argparse, os, re, sys, shutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PREFIX = ("水稻（Oryza sativa）多组学与育种分析技能：转录组、代谢组、多组学整合、"
          "变异检测/GWAS、群体遗传、泛基因组、QTL/BSA、蛋白组、表观组、表型/GS、"
          "单细胞/空间、qPCR 引物与 CRISPR/碱基编辑设计。")

def compact(desc, max_len):
    if len(desc) <= max_len:
        return desc, False
    # 保留“触发词：”之后的全部关键词，压缩前面的能力描述
    i = desc.find("触发词：")
    if i >= 0:
        triggers = desc[i:]
    else:
        j = desc.find("当用户提到")
        triggers = ("触发词：" + desc[j + len("当用户提到"):]) if j >= 0 else ""
    out = PREFIX + triggers
    if len(out) > max_len:
        # 极端情况：截断触发词尾部，但保留“触发词：”结构
        out = out[:max_len].rstrip("、，。 ")
    return out, True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-desc", type=int, default=1024)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "build", "skill_import"))
    a = ap.parse_args()
    src = os.path.join(ROOT, "SKILL.md")
    text = open(src, encoding="utf-8").read()
    lines = text.splitlines()
    idx = next((i for i, l in enumerate(lines) if l.startswith("description:")), None)
    if idx is None:
        sys.exit("SKILL.md 缺少 description")
    desc = lines[idx].split(":", 1)[1].strip()
    short, changed = compact(desc, a.max_desc)
    lines[idx] = "description: " + short

    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "SKILL.md")
    open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    # 附带 README/LICENSE，方便直接打包
    for f in ["README.md", "LICENSE"]:
        p = os.path.join(ROOT, f)
        if os.path.exists(p):
            shutil.copy2(p, os.path.join(a.outdir, f))
    print(f"[import] original={len(desc)} chars; import={len(short)} chars; trimmed={changed}")
    print(f"[import] -> {out}")
    if len(short) > a.max_desc:
        sys.exit("Still exceeds limit; shorten triggers manually.")

if __name__ == "__main__":
    main()
