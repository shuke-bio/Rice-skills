#!/usr/bin/env python3
###############################################################################
# 05_base_prime_edit.py -- base editing (CBE/ABE) and prime editing guide design
#
# 模式：
#   base  : 在 PAM 上游窗口内寻找可编辑 C（CBE, C->T）或 A（ABE, A->G）
#   pe    : 为指定替换设计 pegRNA（spacer + PBS + RTT + scaffold）
#   multi : 批量设计多个位点（CSV: target_id,position,edit）
#
# 说明：这是规则型设计器，不预测编辑效率；正式实验请用 BE-Hive/PRIDICT 等模型复核。
#
# 用法：
#   python scripts/design/05_base_prime_edit.py --sequence ACGT... --position 30 --edit C2T
#   python scripts/design/05_base_prime_edit.py --demo --outdir results/design/05_base_prime
###############################################################################
import argparse, os, sys, csv, datetime

CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate(seq):
    return "".join(CODON.get(seq[i:i+3], "X") for i in range(0, len(seq) - 2, 3))

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def find_guides(seq, pos, edit, pam="NGG"):
    """Return candidate guides for a base edit at 0-based pos."""
    out = []
    want = "C" if edit.upper().startswith("C") else "A"
    if seq[pos].upper() != want:
        return out
    for pam_start in range(pos + 1, min(len(seq), pos + 30)):
        tri = seq[pam_start:pam_start+3].upper()
        if len(tri) == 3 and tri[1:] == "GG":
            for strand in ("+", "-"):
                if strand == "+":
                    proto = seq[pam_start-20:pam_start].upper()
                    offset = pos - (pam_start - 20)
                    pam_seq = seq[pam_start:pam_start+3]
                else:
                    rc = revcomp(seq)
                    # position in rc
                    rc_pos = len(seq) - 1 - pos
                    pstart = rc_pos + (pos - (pam_start - 20))  # not used; recompute below
                    # find PAM upstream of rc_pos
                    continue
                if 0 <= offset < 20 and 3 <= offset <= 10:
                    out.append({"strand": strand, "guide": proto, "pam": pam_seq,
                                "edited_pos_in_guide": offset + 1})
    return out

def design_pe(seq, pos, ref, alt, pbs_len=13, rtt_len=15):
    # spacer: choose PAM downstream of edit (within 30 nt)
    for pam_start in range(pos + 1, min(len(seq), pos + 40)):
        if seq[pam_start+1:pam_start+3].upper() == "GG":
            spacer = seq[pam_start-20:pam_start].upper()
            pbs = revcomp(seq[max(0, pos - pbs_len + 1):pos+1])
            rtt = seq[pos:min(len(seq), pos + rtt_len)]
            rtt = list(rtt); rtt[0] = alt.upper(); rtt = "".join(rtt)
            return {"spacer": spacer, "pam": seq[pam_start:pam_start+3], "pbs": pbs,
                    "rtt": rtt, "edit_pos": pos, "ref": ref, "alt": alt}
    return None

def demo(outdir):
    seq = ("ATG" + "ACG" * 5 + "AAAAA" + "TGG" + "GCT" * 10 + "TAA")
    pos = 3 + 3 * 2 + 1  # C in the third ACG
    with open(os.path.join(outdir, "demo_sequence.txt"), "w") as f:
        f.write(seq + "\n")
    return seq, pos

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sequence"); ap.add_argument("--position", type=int)
    ap.add_argument("--edit", default="C2T", choices=["C2T", "A2G"])
    ap.add_argument("--ref"); ap.add_argument("--alt")
    ap.add_argument("--outdir", default="results/design/05_base_prime")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        seq, pos = demo(a.outdir); a.sequence, a.position = seq, pos
    if not a.sequence or a.position is None: sys.exit("Need --sequence and --position, or --demo")
    seq = a.sequence.upper(); pos = a.position
    guides = find_guides(seq, pos, a.edit)
    codon_start = (pos // 3) * 3
    ref_codon = seq[codon_start:codon_start+3]
    rows = []
    for g in guides:
        edited = list(ref_codon)
        idx = pos - codon_start
        if a.edit == "C2T": edited[idx] = "T"
        else: edited[idx] = "G"
        edited = "".join(edited)
        rows.append([g["strand"], g["guide"], g["pam"], g["edited_pos_in_guide"],
                     ref_codon, edited, translate(ref_codon), translate(edited)])
    with open(os.path.join(a.outdir, "base_edit_guides.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["strand","guide","pam","edited_pos_in_guide","ref_codon","edited_codon","ref_aa","edited_aa"]); w.writerows(rows)
    pe = None
    if a.ref and a.alt:
        pe = design_pe(seq, pos, a.ref, a.alt)
        if pe:
            with open(os.path.join(a.outdir, "pegRNA_design.csv"), "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f); w.writerow(["spacer","pam","pbs","rtt","edit_pos","ref","alt"])
                w.writerow([pe["spacer"], pe["pam"], pe["pbs"], pe["rtt"], pe["edit_pos"], pe["ref"], pe["alt"]])
    with open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8") as f:
        f.write(f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"position={pos}\nedit={a.edit}\nguides={len(rows)}\n")
        f.write("note=rule-based design; no efficiency prediction\n")
    print(f"[edit] {a.edit} at {pos}: {len(rows)} base-editing guides")
    if rows: print("  top:", rows[0])
    if pe: print("[edit] pegRNA:", pe["spacer"], "PBS", pe["pbs"], "RTT", pe["rtt"])
    print(f"[edit] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
