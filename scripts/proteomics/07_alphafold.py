#!/usr/bin/env python3
###############################################################################
# 07_alphafold.py -- query AlphaFold DB for predicted structures
#
# 输入：--uniprot CSV（一列 uniprot_id）或逗号分隔 ID
# 输出：alphafold_models.csv（pLDDT、PDB/PAE 链接）；--download 可下载 PDB
#
# 说明：AlphaFold DB 对植物特有蛋白置信度可能偏低（pLDDT 仅作参考）；
#       本脚本只查询已收录模型，不运行 AlphaFold。
#
# 用法：
#   python scripts/proteomics/07_alphafold.py --uniprot ids.csv --download
#   python scripts/proteomics/07_alphafold.py --demo --outdir results/proteomics/07_alphafold_demo
###############################################################################
import argparse, os, sys, csv, json, urllib.request, datetime

API = "https://alphafold.ebi.ac.uk/api/prediction/{}"

def fetch(uid, timeout=20):
    url = API.format(uid)
    req = urllib.request.Request(url, headers={"User-Agent": "rice-skills/3.10"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = json.loads(r.read().decode("utf-8"))
    if isinstance(data, list) and data:
        d = data[0]
        return {
            "uniprot_id": uid,
            "entry_id": d.get("entryId", ""),
            "gene": d.get("gene", ""),
            "plddt": d.get("globalMetricValue", ""),
            "pdb_url": d.get("pdbUrl", ""),
            "pae_image_url": d.get("paeImageUrl", ""),
        }
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uniprot")
    ap.add_argument("--outdir", default="results/proteomics/07_alphafold")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        p = os.path.join(a.outdir, "demo_uniprot.txt")
        open(p, "w").write("P00533\nQ8W3K0\n")
        a.uniprot = p
        print("[alphafold] demo IDs ->", p, "(no network call in --demo)")
        open(os.path.join(a.outdir, "run_info.txt"), "w").write(
            f"date={datetime.datetime.now().isoformat(timespec='seconds')}\nmode=demo\n")
        return
    if not a.uniprot:
        sys.exit("Need --uniprot or --demo")
    if os.path.exists(a.uniprot):
        ids = [l.strip() for l in open(a.uniprot) if l.strip() and not l.startswith("#")]
    else:
        ids = [x.strip() for x in a.uniprot.split(",") if x.strip()]
    rows, failed = [], []
    for uid in ids:
        try:
            d = fetch(uid)
            if d: rows.append(d)
            else: failed.append(uid)
        except Exception as e:
            failed.append(f"{uid}: {e}")
    if rows:
        keys = list(rows[0].keys())
        with open(os.path.join(a.outdir, "alphafold_models.csv"), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    if a.download:
        for r in rows:
            if not r["pdb_url"]: continue
            try:
                out = os.path.join(a.outdir, f"{r['uniprot_id']}.pdb")
                urllib.request.urlretrieve(r["pdb_url"], out)
            except Exception as e:
                failed.append(f"download {r['uniprot_id']}: {e}")
    open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8").write(
        f"date={datetime.datetime.now().isoformat(timespec='seconds')}\nrequested={len(ids)}\n"
        f"found={len(rows)}\nfailed={'; '.join(failed) or 'none'}\n")
    print(f"[alphafold] requested={len(ids)} found={len(rows)} failed={len(failed)}")
    print(f"[alphafold] outputs -> {a.outdir}")

if __name__ == "__main__":
    main()
