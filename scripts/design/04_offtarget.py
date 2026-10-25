#!/usr/bin/env python3
###############################################################################
# 04_offtarget.py -- genome-wide CRISPR off-target search (SpCas9, NGG PAM)
#
# 说明：
#   - 内置扫描器适合 <= 50 Mb 的基因组（demo/叶绿体/小基因组/靶区）；
#   - 全基因组（水稻 ~400 Mb）请用 Cas-OFFinder/CRISPOR，本脚本会明确拒绝并提示；
#   - 输出的是错配位点与风险分级，不是 CFD 精确打分；需要 CFD 请用 CRISPOR。
#
# 用法：
#   python scripts/design/04_offtarget.py --sgrna sgrna.csv --genome genome.fa \
#       --outdir results/design/04_offtarget --max-mm 3
#   python scripts/design/04_offtarget.py --demo --outdir results/design/04_offtarget_demo
#
# 输入 sgrna.csv：sgrna_id,sequence（20 nt，可带 NGG）
# 输出：offtarget_hits.csv, offtarget_summary.csv, run_info.txt
###############################################################################
import argparse, csv, os, random, re, sys, json, datetime

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs

def read_sgrna(path):
    rows = []
    with open(path, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            sid = (r.get("sgrna_id") or r.get("id") or "").strip()
            seq = (r.get("sequence") or r.get("sgrna") or "").strip().upper()
            seq = re.sub(r"[^ACGTN]", "", seq)
            if len(seq) >= 23 and seq[20:23].endswith("GG"):
                seq = seq[:20]
            if not sid or len(seq) != 20:
                raise ValueError(f"Invalid sgRNA row: id={sid!r}, seq={seq!r} (need 20 nt)")
            rows.append((sid, seq))
    if not rows:
        raise ValueError("No sgRNA rows found")
    return rows

def pam_ok(tri: str) -> bool:
    return len(tri) == 3 and tri[1:] == "GG" and tri[0] in "ACGT"

def scan(seqs, sgrna, max_mm):
    """Return hits: (chrom, pos1, strand, protospacer, pam, mm, mismatches)"""
    hits = []
    L = len(sgrna)
    for chrom, g in seqs.items():
        # plus strand: PAM after protospacer
        for m in re.finditer(r"(?=([ACGT]GG))", g):
            pam_start = m.start()
            if pam_start < L:
                continue
            proto = g[pam_start - L:pam_start]
            if "N" in proto:
                continue
            mm = sum(1 for a, b in zip(proto, sgrna) if a != b)
            if mm <= max_mm:
                pos = pam_start - L + 1
                hits.append((chrom, pos, "+", proto, g[pam_start:pam_start+3], mm,
                             ",".join(str(i + 1) for i, (a, b) in enumerate(zip(proto, sgrna)) if a != b)))
        # minus strand: PAM before protospacer on reverse complement
        rc = revcomp(g)
        for m in re.finditer(r"(?=([ACGT]GG))", rc):
            pam_start = m.start()
            if pam_start < L:
                continue
            proto = rc[pam_start - L:pam_start]
            if "N" in proto:
                continue
            mm = sum(1 for a, b in zip(proto, sgrna) if a != b)
            if mm <= max_mm:
                # rc coordinate -> forward coordinate
                pos = len(g) - (pam_start - L + 1) + 1
                hits.append((chrom, pos, "-", proto, rc[pam_start:pam_start+3], mm,
                             ",".join(str(i + 1) for i, (a, b) in enumerate(zip(proto, sgrna)) if a != b)))
    return hits

def make_demo(outdir):
    random.seed(42)
    genome = list("".join(random.choice("ACGT") for _ in range(120_000)))
    sgrnas, targets = [], []
    for i in range(3):
        # place a PAM site
        pos = 5_000 + i * 30_000
        proto = "".join(random.choice("ACGT") for _ in range(20))
        genome[pos:pos+20] = list(proto)
        genome[pos+20:pos+23] = list("AGG")
        sgrnas.append((f"sg{i+1}", proto))
        targets.append((pos, proto))
    # plant 1-2 mismatch off-targets
    for i, (_, proto) in enumerate(sgrnas):
        for j, mm in enumerate((1, 2), start=1):
            pos = 70_000 + i * 10_000 + j * 1_000
            mut = list(proto)
            for k in range(mm):
                idx = 5 + k * 7
                mut[idx] = random.choice([b for b in "ACGT" if b != mut[idx]])
            genome[pos:pos+20] = mut
            genome[pos+20:pos+23] = list("TGG")
    os.makedirs(outdir, exist_ok=True)
    gpath = os.path.join(outdir, "demo_genome.fa")
    with open(gpath, "w", encoding="utf-8") as f:
        f.write(">chrDemo\n")
        seq = "".join(genome)
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
    spath = os.path.join(outdir, "demo_sgrna.csv")
    with open(spath, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f); w.writerow(["sgrna_id", "sequence"])
        for sid, seq in sgrnas: w.writerow([sid, seq])
    return gpath, spath

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sgrna"); ap.add_argument("--genome")
    ap.add_argument("--outdir", default="results/design/04_offtarget")
    ap.add_argument("--max-mm", type=int, default=3)
    ap.add_argument("--max-genome-mb", type=float, default=50.0)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        gpath, spath = make_demo(a.outdir)
        a.genome, a.sgrna = gpath, spath
    if not a.genome or not a.sgrna:
        sys.exit("Need --genome and --sgrna, or --demo")

    seqs = read_fasta(a.genome)
    total_mb = sum(len(s) for s in seqs.values()) / 1e6
    if total_mb > a.max_genome_mb:
        sys.exit(f"Genome is {total_mb:.1f} Mb > {a.max_genome_mb} Mb. "
                 "Use Cas-OFFinder/CRISPOR for genome-wide search, or raise --max-genome-mb "
                 "only for targeted regions.")

    sgrnas = read_sgrna(a.sgrna)
    all_hits, summary = [], []
    for sid, seq in sgrnas:
        hits = scan(seqs, seq, a.max_mm)
        hits.sort(key=lambda x: (x[5], x[0], x[1]))
        for h in hits:
            all_hits.append([sid, seq] + list(h))
        counts = {k: sum(1 for h in hits if h[5] == k) for k in range(a.max_mm + 1)}
        summary.append([sid, seq, len(hits), counts.get(0, 0), counts.get(1, 0),
                        counts.get(2, 0), sum(v for k, v in counts.items() if k >= 3)])

    hp = os.path.join(a.outdir, "offtarget_hits.csv")
    with open(hp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sgrna_id", "sgrna", "chrom", "pos", "strand", "protospacer",
                    "pam", "mismatches", "mismatch_positions"])
        w.writerows(all_hits)
    sp = os.path.join(a.outdir, "offtarget_summary.csv")
    with open(sp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sgrna_id", "sgrna", "total_hits", "mm0", "mm1", "mm2", "mm3plus"])
        w.writerows(summary)

    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"genome={a.genome}\nsgrna={a.sgrna}\nmax_mm={a.max_mm}\n")
        f.write(f"genome_mb={total_mb:.2f}\n")
        f.write("note=built-in scanner; not CFD; use CRISPOR for CFD scores\n")

    print(f"[offtarget] {len(sgrnas)} sgRNA; genome {total_mb:.2f} Mb; hits={len(all_hits)}")
    for row in summary:
        print(f"  {row[0]}: total={row[2]} mm0={row[3]} mm1={row[4]} mm2={row[5]}")
    print(f"[offtarget] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
