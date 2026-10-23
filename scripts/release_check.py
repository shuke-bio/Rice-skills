#!/usr/bin/env python3
###############################################################################
# release_check.py -- 发版前一致性检查（零依赖）
#
# 检查：
#   1) SKILL.md 的 name 符合小写连字符规范，description <= 1024 字符；
#   2) SKILL.md / CITATION.cff / README 双语徽章 / npm/package.json /
#      package-lock.json / CHANGELOG 的版本号完全一致；
#   3) README 的 ZIP 链接 tag 与版本一致；
#   4) 关键触发词仍在 description 中（防止压缩描述时误删）。
#
# 用法：
#   python scripts/release_check.py
#   python scripts/release_check.py --expected 3.3.12
# 退出码：0 = 全部通过；1 = 有错误
###############################################################################
import argparse, json, os, re, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REQUIRED_TRIGGERS = [
    "水稻转录组", "水稻代谢组", "水稻GWAS", "多位点GWAS", "上位性", "QEI",
    "case-control GWAS", "MAGIC/NAM/ROAM", "群体遗传", "泛基因组", "eQTL",
    "蛋白质组", "表观组", "ATAC", "表型/BLUP", "单细胞", "scRNA",
    "空间转录组", "Stereo-seq", "引物设计", "CRISPR", "qPCR", "多组学",
    "WGCNA", "GSEA", "DESeq2", "OPLS-DA", "MOFA", "DIABLO", "差异表达",
    "富集分析", "CNS 级绘图", "发表图", "写论文", "论文润色", "查文献",
    "文献检索", "下载论文", "精读论文", "结果解读", "审稿回复", "论文转PPT",
    "开题报告",
]

def fail(errors, msg):
    errors.append(msg)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expected", help="期望版本号（可选），例如 3.3.12")
    a = ap.parse_args()
    errors, warnings = [], []
    versions = {}

    # ---- SKILL.md ----
    skill_path = os.path.join(ROOT, "SKILL.md")
    skill = open(skill_path, encoding="utf-8").read()
    m = re.search(r"^name:\s*(\S+)\s*$", skill, re.M)
    name = m.group(1) if m else ""
    if not re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", name):
        fail(errors, f"SKILL.md name 不合规: {name!r}（只允许小写字母/数字/连字符）")

    dm = re.search(r"^description:\s*(.+)$", skill, re.M)
    if not dm:
        fail(errors, "SKILL.md 缺少 description")
        desc = ""
    else:
        desc = dm.group(1).strip()
    if len(desc) > 1024:
        fail(errors, f"description 超长: {len(desc)} > 1024 字符")
    missing = [t for t in REQUIRED_TRIGGERS if t not in desc]
    if missing:
        warnings.append("description 缺少触发词: " + "、".join(missing))

    vm = re.search(r"^version:\s*([0-9]+\.[0-9]+\.[0-9]+)\s*$", skill, re.M)
    versions["SKILL.md"] = vm.group(1) if vm else ""

    # ---- CITATION.cff ----
    cit = open(os.path.join(ROOT, "CITATION.cff"), encoding="utf-8").read()
    m = re.search(r"^version:\s*([0-9]+\.[0-9]+\.[0-9]+)\s*$", cit, re.M)
    versions["CITATION.cff"] = m.group(1) if m else ""

    # ---- README badges / zip link ----
    for fn in ["README.md", "README.zh-CN.md"]:
        s = open(os.path.join(ROOT, fn), encoding="utf-8").read()
        m = re.search(r"badge/version-([0-9]+\.[0-9]+\.[0-9]+)-", s)
        versions[fn] = m.group(1) if m else ""
        zm = re.search(r"archive/refs/tags/v([0-9]+\.[0-9]+\.[0-9]+)\.zip", s)
        if not zm:
            fail(errors, f"{fn} 缺少 ZIP 版本链接")
        elif m and zm.group(1) != m.group(1):
            fail(errors, f"{fn} ZIP 链接 v{zm.group(1)} 与版本 {m.group(1)} 不一致")

    # ---- npm ----
    try:
        pkg = json.load(open(os.path.join(ROOT, "npm", "package.json"), encoding="utf-8"))
        versions["npm/package.json"] = pkg.get("version", "")
        if pkg.get("name") != "rice-skills":
            warnings.append(f"npm 包名是 {pkg.get('name')}，与技能名 rice-skills 不同")
        lock = json.load(open(os.path.join(ROOT, "npm", "package-lock.json"), encoding="utf-8"))
        versions["npm/package-lock.json"] = lock.get("packages", {}).get("", {}).get("version", "")
    except Exception as e:
        fail(errors, f"npm 元数据读取失败: {e}")

    # ---- CHANGELOG ----
    try:
        ch = open(os.path.join(ROOT, "CHANGELOG.md"), encoding="utf-8").read()
        if versions.get("SKILL.md"):
            if f"## [{versions['SKILL.md']}]" not in ch:
                fail(errors, f"CHANGELOG.md 缺少 ## [{versions['SKILL.md']}] 段落")
    except Exception as e:
        fail(errors, f"CHANGELOG 读取失败: {e}")

    # ---- 一致性 ----
    uniq = {v for v in versions.values() if v}
    if len(uniq) != 1:
        for k, v in versions.items():
            if v != versions.get("SKILL.md"):
                fail(errors, f"{k} 版本 {v or '(缺失)'} 与 SKILL.md {versions.get('SKILL.md')} 不一致")
        if "" in versions.values():
            fail(errors, "存在缺失的版本号")
    ver = versions.get("SKILL.md", "")
    if a.expected and ver != a.expected:
        fail(errors, f"版本 {ver} 与期望 {a.expected} 不一致")

    # ---- 输出 ----
    print(f"release_check: version={ver or '(missing)'} name={name or '(missing)'} description={len(desc)} chars")
    for k, v in versions.items():
        print(f"  {v or '(missing)':>8}  {k}")
    for w in warnings:
        print(f"  WARN  {w}")
    for e in errors:
        print(f"  ERROR {e}")
    print("PASS" if not errors else f"FAIL ({len(errors)} errors)")
    sys.exit(0 if not errors else 1)

if __name__ == "__main__":
    main()
