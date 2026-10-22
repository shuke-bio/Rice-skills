#!/usr/bin/env python3
###############################################################################
# 04_check_deps.py -- 分析前依赖体检 + 用户同意后安装
#
# 设计原则：
#   1) 默认只检测、不安装；安装必须显式 --install --yes（或交互确认）；
#   2) 只安装 R/Python 包，不安装系统级工具（PLINK/BLAST/TRN 等仅给指引）；
#   3) 安装前后各检测一次，失败如实报告，不假装成功；
#   4) 可接 --plan 读取方案生成的脚本清单，也可直接 --scripts 指定。
#
# 用法：
#   python scripts/plan/04_check_deps.py --plan results/plan/plan.json
#   python scripts/plan/04_check_deps.py --scripts scripts/transcriptomics/02_deseq2.R
#   python scripts/plan/04_check_deps.py --plan ... --install            # 只打印安装命令
#   python scripts/plan/04_check_deps.py --plan ... --install --yes      # 执行安装
#   python scripts/plan/04_check_deps.py --demo
#
# 输出：deps_report.md、deps.json、install.log（仅安装时）
###############################################################################
import argparse, json, os, re, shutil, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

PY_MAP = {"sklearn": "scikit-learn", "PIL": "pillow", "cv2": "opencv-python",
          "yaml": "pyyaml", "Bio": "biopython", "dateutil": "python-dateutil"}
BIOC = {"DESeq2", "edgeR", "limma", "clusterProfiler", "org.Osativa.eg.db",
        "ComplexHeatmap", "preprocessCore", "impute", "sva", "DMRcate",
        "ChIPseeker", "GenomicRanges", "IRanges", "rtracklayer", "MOFA2",
        "mixOmics", "ropls", "methylKit", "SummarizedExperiment", "S4Vectors",
        "BiocGenerics", "GenomeInfoDb", "GenomicFeatures", "apeglm", "fgsea",
        "GSEABase", "DOSE", "enrichplot", "STRINGdb", "HDF5Array"}
RECOMMENDED = {"MASS", "Matrix", "survival", "nlme", "stats", "utils", "methods"}
TOOLS = {"plink": "PLINK 1.9", "plink2": "PLINK 2", "bwa": "BWA", "samtools": "samtools",
         "bcftools": "bcftools", "gatk": "GATK", "minimap2": "minimap2", "syri": "SyRI",
         "snakemake": "Snakemake", "blastn": "BLAST+", "trf": "TRF", "node": "Node.js",
         "java": "Java", "fastp": "fastp", "hisat2": "HISAT2", "salmon": "salmon",
         "featureCounts": "Subread"}

def read_text(p):
    try:
        return open(p, encoding="utf-8", errors="ignore").read()
    except Exception:
        return ""

STDLIB = {"sys", "os", "re", "json", "csv", "argparse", "math", "random", "subprocess",
          "pathlib", "collections", "itertools", "datetime", "time", "shutil", "glob",
          "warnings", "functools", "typing", "dataclasses", "statistics", "io", "textwrap",
          "gzip", "zipfile", "urllib", "hashlib", "base64", "pickle", "multiprocessing",
          "concurrent", "logging", "unittest", "tempfile", "traceback", "importlib",
          "platform", "socket", "struct", "sqlite3", "xml", "html", "http", "email", "ssl", "uuid"}

def _local_py_modules():
    names = set()
    for dp, _, fns in os.walk(ROOT):
        if "third_party" in dp or ".git" in dp:
            continue
        for fn in fns:
            if fn.endswith(".py"):
                names.add(os.path.splitext(fn)[0])
    return names

def scan_script(path):
    s = read_text(path)
    ext = os.path.splitext(path)[1].lower()
    py, r, tools = set(), set(), set()
    if ext == ".py":
        local = _local_py_modules()
        for m in re.finditer(r"^\s*(?:import|from)\s+([A-Za-z0-9_\.]+)", s, re.M):
            mod = m.group(1).split(".")[0]
            if mod not in STDLIB and mod not in local:
                py.add(mod)
    elif ext == ".r":
        base_pkgs = {"base", "stats", "utils", "methods", "graphics", "grDevices",
                     "datasets", "tools", "parallel", "compiler"}
        for m in re.finditer(r"\b(?:library|require|requireNamespace)\(\s*['\"]?([A-Za-z0-9_.]+)", s):
            name = m.group(1)
            if name in base_pkgs or re.match(r"^[a-z][a-z0-9_]*$", name):
                continue
            r.add(name)
        for m in re.finditer(r"need_pkgs\(\s*c?\(([^)]*)\)", s):
            for q in re.findall(r"['\"]([A-Za-z0-9_.]+)['\"]", m.group(1)):
                r.add(q)
        for m in re.finditer(r"\b([A-Za-z][A-Za-z0-9.]*)::", s):
            if m.group(1) not in base_pkgs:
                r.add(m.group(1))
    elif ext == ".sh":
        for t in TOOLS:
            if re.search(r"(?<![A-Za-z0-9_-])" + re.escape(t) + r"(?![A-Za-z0-9_-])", s):
                tools.add(t)
    return py, r, tools

def _run_json(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", timeout=120)
        return json.loads(out.stdout.strip().splitlines()[-1])
    except Exception:
        return None

def check_python(pkgs, python):
    if not pkgs:
        return {}
    code = ("import importlib.util,json,sys;"
            "print(json.dumps({p: bool(importlib.util.find_spec(p)) for p in sys.argv[1:]}))")
    res = _run_json([python, "-c", code, *sorted(pkgs)])
    return res or {p: False for p in pkgs}

def check_r(pkgs, rscript):
    if not pkgs:
        return {}
    if not rscript or not os.path.exists(rscript):
        return {p: False for p in pkgs}
    vec = ",".join(f'"{p}"' for p in sorted(pkgs))
    code = (f'p<-c({vec}); ok<-vapply(p, function(x) requireNamespace(x, quietly=TRUE), logical(1));'
            'cat(paste(names(ok), ok, sep="="), sep="\\n")')
    try:
        out = subprocess.run([rscript, "-e", code], capture_output=True, text=True,
                             encoding="utf-8", timeout=300)
        lines = [l.strip() for l in out.stdout.splitlines()
                 if re.match(r"^[A-Za-z0-9_.]+=(TRUE|FALSE)$", l.strip())]
        return {k: v == "TRUE" for k, v in (l.split("=", 1) for l in lines)}
    except Exception:
        return {p: False for p in pkgs}

def find_rscript(explicit):
    if explicit and os.path.exists(explicit):
        return explicit
    w = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if w:
        return w
    for base in [r"C:\Program Files\R", r"C:\Program Files\R\R-4.6.0\bin"]:
        if os.path.isdir(base):
            for d in sorted(os.listdir(base), reverse=True):
                p = os.path.join(base, d, "bin", "Rscript.exe")
                if os.path.exists(p):
                    return p
    return ""

def collect_scripts(args):
    scripts = list(args.scripts or [])
    if args.plan:
        data = json.load(open(args.plan, encoding="utf-8"))
        want = {"minimal": "最小可行", "standard": "标准", "advanced": "进阶"}.get(args.tier)
        for t in data.get("tiers", []):
            if want and want not in t.get("tier", ""):
                continue
            for step in t.get("steps", []):
                for tok in re.findall(r"[A-Za-z0-9_./-]+\.(?:R|py|sh)", step):
                    scripts.append(tok)
    resolved, seen = [], set()
    index = {}
    for dp, _, fns in os.walk(os.path.join(ROOT, "scripts")):
        for fn in fns:
            if fn.endswith((".R", ".py", ".sh")):
                index.setdefault(fn, []).append(os.path.relpath(os.path.join(dp, fn), ROOT).replace("\\", "/"))
    for s in scripts:
        s = s.strip().lstrip("./")
        if s in seen: continue
        if os.path.exists(os.path.join(ROOT, s)):
            resolved.append(s); seen.add(s); continue
        hits = index.get(os.path.basename(s), [])
        if len(hits) == 1:
            resolved.append(hits[0]); seen.add(hits[0])
        elif hits:
            resolved.extend(h for h in hits if h not in seen)
    return resolved

def build_plan(scripts):
    py, r, tools = set(), set(), set()
    for s in scripts:
        a, b, c = scan_script(os.path.join(ROOT, s))
        py |= a; r |= b; tools |= c
    return py, r, tools

def install_commands(py_missing, r_missing, python):
    cmds = []
    if py_missing:
        pkgs = [PY_MAP.get(p, p) for p in sorted(py_missing)]
        cmds.append([python, "-m", "pip", "install", *pkgs])
    bioc = sorted(p for p in r_missing if p in BIOC)
    cran = sorted(p for p in r_missing if p not in BIOC and p not in RECOMMENDED)
    if bioc:
        cmds.append(["Rscript", "-e",
                     'if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org");'
                     f'BiocManager::install(c({",".join(repr(p) for p in bioc)}),ask=FALSE,update=FALSE)'])
    if cran:
        cmds.append(["Rscript", "-e",
                     f'install.packages(c({",".join(repr(p) for p in cran)}),repos="https://cloud.r-project.org")'])
    return cmds

def confirm(msg):
    try:
        if not sys.stdin.isatty():
            return False
        return input(msg + " [y/N] ").strip().lower() in ("y", "yes")
    except (EOFError, OSError, KeyboardInterrupt):
        return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan"); ap.add_argument("--scripts", nargs="*")
    ap.add_argument("--tier", default="minimal",
                    choices=["minimal", "standard", "advanced", "all"])
    ap.add_argument("--outdir", default="results/deps")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--rscript", default="")
    ap.add_argument("--install", action="store_true")
    ap.add_argument("--yes", action="store_true")
    ap.add_argument("--run", action="store_true", help="安装后以 --demo 运行脚本")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    outdir = a.outdir if os.path.isabs(a.outdir) else os.path.join(ROOT, a.outdir)
    os.makedirs(outdir, exist_ok=True)

    if a.demo:
        d = os.path.join(outdir, "_demo"); os.makedirs(d, exist_ok=True)
        p1 = os.path.join(d, "fake_a.py"); open(p1, "w").write("import fake_pkg_xyz\n")
        p2 = os.path.join(d, "fake_b.R"); open(p2, "w").write('library(fakePkgXYZ)\n')
        py, r, tools = set(), set(), set()
        for p in (p1, p2):
            x, y, z = scan_script(p); py |= x; r |= y; tools |= z
        cmds = install_commands(py, r, a.python)
        assert "fake_pkg_xyz" in py and "fakePkgXYZ" in r, "demo 扫描失败"
        assert cmds, "demo 未生成安装命令"
        print(f"DEMO OK: python={sorted(py)} r={sorted(r)} commands={len(cmds)}")
        return

    scripts = collect_scripts(a)
    if not scripts:
        sys.exit("未找到脚本：请给 --plan 或 --scripts")
    py, r, tools = build_plan(scripts)
    rscript = find_rscript(a.rscript)
    py_st = check_python(py, a.python)
    r_st = check_r(r, rscript)
    tool_st = {t: bool(shutil.which(t)) for t in tools}

    py_missing = sorted(p for p, ok in py_st.items() if not ok)
    r_missing = sorted(p for p, ok in r_st.items() if not ok)
    tool_missing = sorted(t for t, ok in tool_st.items() if not ok)
    cmds = install_commands(py_missing, r_missing, a.python)

    report = ["# 依赖体检报告", "",
              f"- 脚本数：{len(scripts)}", f"- Python：{a.python}",
              f"- Rscript：{rscript or '未找到'}", ""]
    report += ["## 缺失依赖", ""]
    if not (py_missing or r_missing or tool_missing):
        report += ["无缺失，可直接运行。", ""]
    else:
        report += ["| 类型 | 名称 | 安装方式 |", "|---|---|---|"]
        for p in py_missing: report.append(f"| Python | `{p}` | `pip install {PY_MAP.get(p,p)}` |")
        for p in r_missing:
            how = "BiocManager" if p in BIOC else "CRAN"
            report.append(f"| R | `{p}` | {how} |")
        for t in tool_missing:
            report.append(f"| 工具 | `{t}` | 需手动安装（{TOOLS.get(t,t)}），脚本不会自动装系统工具 |")
        report.append("")
    report += ["## 安装命令（需你确认后执行）", ""]
    if cmds:
        for c in cmds:
            report.append("```bash\n" + " ".join(c) + "\n```")
    else:
        report.append("无需要安装的 Python/R 包。")
    report += ["", "> 系统工具（PLINK/BLAST/TRF/Snakemake 等）不会被自动安装；"
               "请按仓库 `docs/installation.md` 或对应官网安装。", ""]
    open(os.path.join(outdir, "deps_report.md"), "w", encoding="utf-8").write("\n".join(report))
    json.dump({"scripts": scripts, "missing": {"python": py_missing, "r": r_missing, "tools": tool_missing},
               "commands": [" ".join(c) for c in cmds]},
              open(os.path.join(outdir, "deps.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    print(f"scripts={len(scripts)} missing_py={len(py_missing)} missing_r={len(r_missing)} missing_tools={len(tool_missing)}")
    for p in py_missing: print("  [py]", p)
    for p in r_missing: print("  [R ]", p)
    for t in tool_missing: print("  [tool]", t)

    if not a.install:
        print(f"报告 -> {os.path.join(outdir, 'deps_report.md')}（未安装任何东西）")
        return
    if not (py_missing or r_missing):
        print("无可安装的 Python/R 包。"); return
    if not a.yes and not confirm("是否安装以上 Python/R 包？"):
        print("已取消安装（未做任何改动）。"); sys.exit(2)

    log = open(os.path.join(outdir, "install.log"), "w", encoding="utf-8")
    ok = True
    for c in cmds:
        cmd = list(c)
        if cmd[0] == "Rscript": cmd[0] = rscript or "Rscript"
        log.write("$ " + " ".join(cmd) + "\n")
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           encoding="utf-8", errors="ignore")
        log.write(p.stdout + "\n")
        if p.returncode != 0:
            ok = False; print("安装失败：", " ".join(cmd))
    log.close()

    py2 = check_python(py, a.python); r2 = check_r(r, rscript)
    still = [p for p, v in py2.items() if not v] + [p for p, v in r2.items() if not v]
    print("安装后仍缺失：", still or "无")
    if still: ok = False

    if a.run and ok:
        print("开始以 --demo 运行脚本：")
        for s in scripts:
            ext = os.path.splitext(s)[1].lower()
            cmd = ["bash", s, "--demo"] if ext == ".sh" else (
                  [a.python, s, "--demo"] if ext == ".py" else [rscript, s, "--demo"])
            print(" $", " ".join(cmd))
            subprocess.run(cmd, cwd=ROOT)
    if not ok: sys.exit(1)

if __name__ == "__main__":
    main()
