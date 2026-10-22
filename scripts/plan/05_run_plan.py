#!/usr/bin/env python3
###############################################################################
# 05_run_plan.py -- 方案确认闸门与执行
#
# 规则：先把「分析方案 + 图形方案」给用户看，用户明确确认后才允许执行。
#       没有 approval.json 或方案已变更（hash 不匹配）时，拒绝执行。
#
# 用法：
#   python scripts/plan/05_run_plan.py --plan results/plan/plan.json --tier minimal
#       # 只打印确认单（分析步骤 + 图形方案 + 阻断项），不执行
#   python scripts/plan/05_run_plan.py --plan ... --tier minimal --approve
#       # 用户确认后调用：写入 approval.json
#   python scripts/plan/05_run_plan.py --plan ... --tier minimal --execute
#       # 校验批准记录 → 依赖体检 → 逐个 --demo 执行
#
# 输出：approval.json、run.log
###############################################################################
import argparse, hashlib, json, os, re, subprocess, sys
from datetime import datetime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PY = sys.executable

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def load_tier(plan_path, tier):
    data = json.load(open(plan_path, encoding="utf-8"))
    want = {"minimal": "最小可行", "standard": "标准", "advanced": "进阶"}[tier]
    for t in data.get("tiers", []):
        if want in t.get("tier", ""):
            return data, t
    sys.exit(f"方案中没有档位：{tier}")

def sheet(data, t):
    d = data["input"]
    L = [f"# 方案确认单：{t['tier']}", "",
         f"- 数据类型：`{d['assay']}`｜目标：`{d['goal']}`｜设计：`{d['design']}`",
         f"- 样本量：{d['n_per_group']}｜重复：{d['replicates']}｜批次：{d['batch']}", ""]
    if data.get("blockers"):
        L += ["## ⛔ 阻断项（需先解决）", ""] + [f"- {b}" for b in data["blockers"]] + [""]
    if data.get("warnings"):
        L += ["## ⚠️ 风险提示", ""] + [f"- {w}" for w in data["warnings"]] + [""]
    L += ["## 分析步骤", ""] + [f"{i}. {s}" for i, s in enumerate(t["steps"], 1)] + [""]
    if t.get("figures"):
        L += ["## 图形方案", "", "| 图 | 类型 | 规格 | 脚本 |", "|---|---|---|---|"]
        L += [f"| {f['name']} | {f['type']} | {f['spec']} | `{f['script']}` |" for f in t["figures"]]
        L.append("")
    L += ["## 预期产出", "", "；".join(t.get("outputs", [])), "",
          "---", "请确认以上**分析步骤与图形方案**；确认后再执行。"]
    return "\n".join(L)

def resolve_scripts(t):
    idx = {}
    for dp, _, fns in os.walk(os.path.join(ROOT, "scripts")):
        for fn in fns:
            if fn.endswith((".R", ".py", ".sh")):
                rel = os.path.relpath(os.path.join(dp, fn), ROOT).replace("\\", "/")
                idx.setdefault(fn, []).append(rel)
    out, seen = [], set()
    for s in t["steps"]:
        for tok in re.findall(r"[A-Za-z0-9_./-]+\.(?:R|py|sh)", s):
            tok = tok.lstrip("./")
            if tok in seen: continue
            cand = tok if os.path.exists(os.path.join(ROOT, tok)) else None
            if not cand:
                hits = idx.get(os.path.basename(tok), [])
                cand = hits[0] if len(hits) == 1 else None
            if cand:
                out.append(cand); seen.add(cand)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", required=True)
    ap.add_argument("--tier", default="minimal", choices=["minimal", "standard", "advanced"])
    ap.add_argument("--approve", action="store_true")
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--skip-deps", action="store_true")
    a = ap.parse_args()
    plan = a.plan if os.path.isabs(a.plan) else os.path.join(ROOT, a.plan)
    if not os.path.exists(plan): sys.exit(f"方案不存在：{plan}")
    data, t = load_tier(plan, a.tier)
    phash = sha(plan)
    outdir = os.path.dirname(plan)
    approval_path = os.path.join(outdir, "approval.json")

    if a.approve:
        json.dump({"plan": os.path.abspath(plan), "plan_sha256": phash, "tier": a.tier,
                   "figures": len(t.get("figures", [])), "approved": True,
                   "approved_at": datetime.now().isoformat(timespec="seconds")},
                  open(approval_path, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
        print(f"已记录确认：{approval_path}")
        print(sheet(data, t))
        return

    if not a.execute:
        print(sheet(data, t))
        return

    # ---- execute: require matching approval ----
    if data.get("blockers"):
        sys.exit("存在阻断项，不能执行。请先解决后重新生成方案。")
    if not os.path.exists(approval_path):
        sys.exit("未找到用户确认记录 approval.json：请先展示方案，用户确认后运行 --approve。")
    ap_data = json.load(open(approval_path, encoding="utf-8"))
    if ap_data.get("plan_sha256") != phash or ap_data.get("tier") != a.tier:
        sys.exit("方案已变更或档位不一致：请重新展示方案并再次确认。")

    if not a.skip_deps:
        rc = subprocess.run([PY, os.path.join(ROOT, "scripts/plan/04_check_deps.py"),
                             "--plan", plan, "--tier", a.tier,
                             "--outdir", os.path.join(outdir, "deps")], cwd=ROOT)
        if rc.returncode != 0:
            sys.exit("依赖未就绪：请查看 deps/deps_report.md，用户同意后再安装。")

    scripts = resolve_scripts(t)
    log = open(os.path.join(outdir, "run.log"), "w", encoding="utf-8")
    for s in scripts:
        body = open(os.path.join(ROOT, s), encoding="utf-8", errors="ignore").read()
        if "--demo" not in body:
            log.write(f"# skip {s}: 不支持 --demo，需真实参数\n")
            print(f"  跳过 {s}（不支持 --demo，需按 docs/*_usage.md 传参）")
            continue
        ext = os.path.splitext(s)[1].lower()
        cmd = (["bash", s, "--demo"] if ext == ".sh" else
               [PY, s, "--demo"] if ext == ".py" else ["Rscript", s, "--demo"])
        log.write("$ " + " ".join(cmd) + "\n")
        print("$", " ".join(cmd))
        p = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, encoding="utf-8", errors="ignore")
        log.write(p.stdout + "\n")
        if p.returncode != 0:
            log.write(f"[exit {p.returncode}]\n")
            print(f"  脚本返回 {p.returncode}，停止。详见 run.log")
            break
    log.close()
    print(f"执行记录 -> {os.path.join(outdir, 'run.log')}")

if __name__ == "__main__":
    main()
