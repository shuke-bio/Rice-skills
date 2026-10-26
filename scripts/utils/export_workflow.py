#!/usr/bin/env python3
###############################################################################
# export_workflow.py -- export a step list to Nextflow / Snakemake / Dockerfile
#
# 输入 steps JSON：
#   [{"name":"qc","command":"python scripts/...","cpus":2,"mem":"8 GB"}, ...]
#
# 输出：
#   workflow/main.nf, workflow/nextflow.config
#   workflow/Snakefile
#   workflow/Dockerfile
#   workflow/README.md
#
# 用法：
#   python scripts/utils/export_workflow.py --steps steps.json --outdir workflow
#   python scripts/utils/export_workflow.py --demo --outdir results/workflow_demo
###############################################################################
import argparse, json, os, sys, datetime

DEMO_STEPS = [
    {"name": "qc", "command": "python scripts/transcriptomics/01_qc_filter.py --counts counts.csv --metadata meta.csv", "cpus": 2, "mem": "8 GB"},
    {"name": "deseq2", "command": "Rscript scripts/transcriptomics/02_deseq2.R --counts counts_filtered.csv --metadata meta.csv --contrast salt_vs_control", "cpus": 2, "mem": "16 GB"},
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps"); ap.add_argument("--outdir", default="workflow")
    ap.add_argument("--version", default="3.7.0")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        sp = os.path.join(a.outdir, "demo_steps.json")
        json.dump(DEMO_STEPS, open(sp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
        a.steps = sp
    if not a.steps: sys.exit("Need --steps or --demo")
    steps = json.load(open(a.steps, encoding="utf-8"))
    if not isinstance(steps, list) or not steps: sys.exit("steps JSON must be a non-empty list")

    # Nextflow
    nf = ["nextflow.enable.dsl=2", ""]
    calls = []
    for s in steps:
        name = s["name"]
        nf += [f"process {name.upper()} {{",
               f"    container 'rice-skills:{a.version}'",
               f"    cpus {s.get('cpus', 2)}",
               f"    memory '{s.get('mem', '8 GB')}'",
               '    output: path "results/**"',
               '    script:',
               f'    """',
               f'    {s["command"]}',
               f'    """',
               "}", ""]
        calls.append(f"    {name.upper()}()")
    nf += ["workflow {"] + calls + ["}"]
    open(os.path.join(a.outdir, "main.nf"), "w", encoding="utf-8").write("\n".join(nf) + "\n")

    open(os.path.join(a.outdir, "nextflow.config"), "w", encoding="utf-8").write(
        f"docker.enabled = true\ndocker.runOptions = '--rm'\nparams.version = '{a.version}'\n")

    # Snakemake
    sk = ["SHELL := /bin/bash", ""]
    for s in steps:
        sk += [f"rule {s['name']}:",
               f"    threads: {s.get('cpus', 2)}",
               f"    resources: mem_mb={(int(str(s.get('mem','8 GB')).split()[0]) * 1024)}",
               f"    shell:",
               f"        \"{s['command']}\"", ""]
    sk += ["rule all:", "    input: " + ", ".join(f"'{s['name']}.done'" for s in steps)]
    open(os.path.join(a.outdir, "Snakefile"), "w", encoding="utf-8").write("\n".join(sk) + "\n")

    # Dockerfile
    open(os.path.join(a.outdir, "Dockerfile"), "w", encoding="utf-8").write(
        "FROM continuumio/miniconda3:latest\n"
        "WORKDIR /opt/rice-skills\n"
        "COPY environment.yml .\n"
        "RUN conda env create -f environment.yml && conda clean -afy\n"
        "SHELL [\"conda\", \"run\", \"-n\", \"rice-skills\", \"/bin/bash\", \"-c\"]\n"
        "COPY . .\n"
        "ENTRYPOINT [\"conda\", \"run\", \"-n\", \"rice-skills\"]\n"
        "CMD [\"bash\"]\n")

    readme = [f"# Workflow export ({datetime.date.today().isoformat()})", "",
              f"{len(steps)} steps exported.", "",
              "## Nextflow", "```bash", "nextflow run main.nf -with-docker", "```", "",
              "## Snakemake", "```bash", "snakemake -j 8 --use-conda", "```", "",
              "## Docker", "```bash", f"docker build -t rice-skills:{a.version} .", "```", "",
              "## Steps"]
    readme += [f"- `{s['name']}`: {s['command']}" for s in steps]
    open(os.path.join(a.outdir, "README.md"), "w", encoding="utf-8").write("\n".join(readme) + "\n")

    print(f"[workflow] exported {len(steps)} steps -> {a.outdir}")
    print("  main.nf, nextflow.config, Snakefile, Dockerfile, README.md")

if __name__ == "__main__":
    main()
