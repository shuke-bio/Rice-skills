#!/usr/bin/env python3
###############################################################################
# hpc_submit.py -- generate SLURM / PBS submission scripts
#
# 用法：
#   python scripts/utils/hpc_submit.py --cmd "Rscript scripts/gwas/03_association.R ..." \
#       --engine slurm --name gwas --cpus 16 --mem 64G --time 24:00:00
#   python scripts/utils/hpc_submit.py --demo --outdir results/hpc_demo
###############################################################################
import argparse, os, sys, datetime

def slurm(a):
    return f"""#!/usr/bin/env bash
#SBATCH --job-name={a.name}
#SBATCH --cpus-per-task={a.cpus}
#SBATCH --mem={a.mem}
#SBATCH --time={a.time}
#SBATCH --output={a.outdir}/%x.%j.out
#SBATCH --error={a.outdir}/%x.%j.err

set -euo pipefail
cd "{a.workdir}"
conda activate rice-skills
{a.cmd}
"""

def pbs(a):
    return f"""#!/usr/bin/env bash
#PBS -N {a.name}
#PBS -l nodes=1:ppn={a.cpus}
#PBS -l mem={a.mem}
#PBS -l walltime={a.time}
#PBS -o {a.outdir}/{a.name}.out
#PBS -e {a.outdir}/{a.name}.err

set -euo pipefail
cd "{a.workdir}"
conda activate rice-skills
{a.cmd}
"""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cmd"); ap.add_argument("--engine", choices=["slurm", "pbs"], default="slurm")
    ap.add_argument("--name", default="rice-skills")
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--mem", default="32G")
    ap.add_argument("--time", default="12:00:00")
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--outdir", default="results/hpc")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.cmd = "Rscript scripts/gwas/03_association.R --genotype demo_data/genomics/genotype.csv --snp-info demo_data/genomics/snp_info.csv --phenotype demo_data/genomics/phenotype.csv --trait salt_tolerance_index --engine rmvp --model MLM --outdir results/gwas/03_association"
    if not a.cmd: sys.exit("Need --cmd or --demo")
    text = slurm(a) if a.engine == "slurm" else pbs(a)
    path = os.path.join(a.outdir, f"submit_{a.name}.sh")
    open(path, "w", encoding="utf-8").write(text)
    os.chmod(path, 0o755)
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\nengine={a.engine}\n")
        f.write(f"cpus={a.cpus}\nmem={a.mem}\ntime={a.time}\n")
    print(f"[hpc] {a.engine} script -> {path}")
    print(f"  submit: {'sbatch' if a.engine=='slurm' else 'qsub'} {path}")

if __name__ == "__main__":
    main()
