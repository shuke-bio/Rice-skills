#!/usr/bin/env python3
###############################################################################
# 03_validate_cards.py -- 校验方法卡
#
# 检查：JSON 可解析、必填字段、id 与文件名一致、脚本路径存在、
#       来源期刊在白名单内（或标记 tool_only）、DOI 格式、id 唯一。
#
# 用法：python scripts/plan/03_validate_cards.py [--cards-dir docs/decision/cards]
# 退出码：0 = 全部通过；1 = 有错误
###############################################################################
import argparse, json, os, re, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WHITELIST = {
    "Nature", "Science", "Cell", "Nature Biotechnology", "Nature Methods",
    "Nature Genetics", "Nature Plants", "Nature Communications", "Science Advances",
    "PNAS", "Molecular Plant", "Plant Cell", "Plant Communications", "Plant Physiology",
    "JIPB", "New Phytologist", "Plant Journal", "Briefings in Bioinformatics",
    "Nucleic Acids Research", "Genome Biology", "Genome Research", "Bioinformatics",
    "GigaScience", "Plant Phenomics", "Crop Journal", "Rice",
}
REQUIRED = ["id", "name", "category", "applies_to", "assumptions", "inputs",
            "outputs", "alternatives", "fail_modes", "scripts", "source", "tier"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cards-dir", default=os.path.join("docs", "decision", "cards"))
    a = ap.parse_args()
    d = a.cards_dir if os.path.isabs(a.cards_dir) else os.path.join(ROOT, a.cards_dir)
    if not os.path.isdir(d):
        sys.exit(f"目录不存在: {d}")

    errors, ids, n = [], set(), 0
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".json"):
            continue
        n += 1
        p = os.path.join(d, fn)
        try:
            c = json.load(open(p, encoding="utf-8"))
        except Exception as e:
            errors.append(f"{fn}: JSON 解析失败 {e}"); continue
        for k in REQUIRED:
            if k not in c:
                errors.append(f"{fn}: 缺少字段 {k}")
        cid = c.get("id", "")
        if cid != os.path.splitext(fn)[0]:
            errors.append(f"{fn}: id({cid}) 与文件名不一致")
        if cid in ids:
            errors.append(f"{fn}: id 重复 {cid}")
        ids.add(cid)
        for s in c.get("scripts", []):
            if not os.path.exists(os.path.join(ROOT, s)):
                errors.append(f"{fn}: 脚本不存在 {s}")
        src = c.get("source", {})
        if not src.get("tool_only"):
            j = src.get("journal")
            if j not in WHITELIST:
                errors.append(f"{fn}: 期刊不在白名单: {j}")
            doi = src.get("doi", "")
            if doi and not re.match(r"^10\.\d{4,9}/\S+$", doi):
                errors.append(f"{fn}: DOI 格式可疑: {doi}")
        if not c.get("fail_modes"):
            errors.append(f"{fn}: fail_modes 不能为空")
        if c.get("tier") not in ("S", "A", "B"):
            errors.append(f"{fn}: tier 必须是 S/A/B")
    print(f"cards={n} errors={len(errors)}")
    for e in errors:
        print(" -", e)
    sys.exit(0 if not errors else 1)

if __name__ == "__main__":
    main()
