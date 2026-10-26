#!/usr/bin/env python3
###############################################################################
# 06_edit_efficiency.py -- heuristic base/prime editing guide ranking
#
# 说明：这是规则型打分（GC、homopolymer、编辑窗口位置），不是 BE-Hive/PRIDICT
# 的深度学习预测；正式实验请用外部模型复核（脚本会生成命令计划）。
#
# 输入 CSV：guide_id, sequence, edit_type（BE/PE）
#
# 用法：
#   python scripts/design/06_edit_efficiency.py --guides guides.csv --outdir results/design/06_edit
#   python scripts/design/06_edit_efficiency.py --demo --outdir results/design/06_edit_demo
###############################################################################
import argparse, os, sys, csv, datetime

def score(seq, edit):
    seq = seq.upper(); n = len(seq)
    gc = 100 * (seq.count("G") + seq.count("C")) / n
    s = 100.0
    if gc < 40 or gc > 70: s -= min(30, abs(gc - 55) * 1.5)
    # homopolymer
    run = 1
    for i in range(1, n):
        run = run + 1 if seq[i] == seq[i-1] else 1
        if run >= 5: s -= 10
    # editing window (BE C/A at 4-8; PE no strong rule)
    if edit.upper() == "BE":
        pos = next((i for i, b in enumerate(seq) if b in "CA"), None)
        if pos is None or not (3 <= pos <= 8): s -= 15
    return max(0, round(s, 1)), round(gc, 1)

def make_demo(outdir):
    rows = [
        ["g1", "ACGTACGTACGTACGTACGT", "BE"],
        ["g2", "GGGGGCCCCAAAATTTTACG", "BE"],
        ["g3", "ATATATATATATATATATAT", "PE"],
        ["g4", "GCGCGCGCGCGCGCGCGCGC", "BE"],
    ]
    p = os.path.join(outdir, "demo_guides.csv")
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["guide_id","sequence","edit_type"]); w.writerows(rows)
    return p

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--guides"); ap.add_argument("--outdir", default="results/design/06_edit")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    if a.demo: a.guides = make_demo(a.outdir)
    if not a.guides: sys.exit("Need --guides or --demo")
    out = []
    with open(a.guides, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            seq = r["sequence"].strip().upper()
            sc, gc = score(seq, r.get("edit_type", "BE"))
            out.append([r["guide_id"], seq, r.get("edit_type", "BE"), gc, sc])
    out.sort(key=lambda x: -x[4])
    with open(os.path.join(a.outdir, "edit_efficiency_ranking.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["guide_id","sequence","edit_type","gc_pct","heuristic_score"]); w.writerows(out)
    plan = [
        "# 外部效率模型复核（示例命令；需自行安装）",
        "# BE-Hive: https://github.com/maxwshen/be-hive",
        "# PRIDICT: https://github.com/uzh-dqbm-cmi/PRIDICT",
        "# 将 guides.csv 转为模型输入后运行；本脚本的 heuristic_score 仅用于初筛。",
    ]
    open(os.path.join(a.outdir, "external_plan.txt"), "w").write("\n".join(plan) + "\n")
    open(os.path.join(a.outdir, "run_info.txt"), "w").write(
        f"date={datetime.datetime.now().isoformat(timespec='seconds')}\nguides={len(out)}\n"
        "score=heuristic (GC/homopolymer/window); not BE-Hive/PRIDICT\n")
    print(f"[edit] guides={len(out)}; top={out[0][0]} score={out[0][4]}")
    print(f"[edit] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
