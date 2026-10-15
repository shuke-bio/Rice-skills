"""fetch_rice_annotations.py -- download and cache official rice annotation resources.

Title:    Rice-skills annotation fetcher
Description:
    Downloads rice annotation resources from their OFFICIAL sources, caches
    them locally, and derives the files that the Rice-skills pipeline can
    consume directly:

      1. RAP-DB MSU <-> RAP-DB gene-ID conversion table
         -> msu_rapdb_map.csv          (msu,rapdb)
         -> rapdb_to_msu_map.csv       (from,to)  <-- pass to adapt_company_data.py --id-map
         -> msu_to_rapdb_map.csv       (from,to)
      2. PlantTFDB Oryza sativa Japonica transcription-factor list
         (already MSU-based, e.g. LOC_Os01g04750)
         -> rice_TF_PlantTFDB.csv      (tf_id,gene_id,family)
         -> rice_TF_PlantTFDB_msu.csv  (gene_id,family) <-- 06_tf_analysis.R --tf-annotation
         -> rice_TF_families.gmt       (GSEA gene sets per TF family) <-- 04_gsea.R --gmt

    Nothing is redistributed at download time: files are fetched from the
    official RAP-DB / PlantTFDB servers and cached under --outdir (git-ignored).
    A verified snapshot of the same resources also ships inside the repository
    (assets/annotations/) and is used automatically when the cache is empty,
    so the pipeline works offline out of the box; run with --force to refresh.

Usage:
    python scripts/utils/fetch_rice_annotations.py
    python scripts/utils/fetch_rice_annotations.py --outdir data/annotations --force
    # networks that intercept TLS (error: "self-signed certificate in chain"):
    python scripts/utils/fetch_rice_annotations.py --insecure
    # offline / mirror use (local file paths are accepted for both options):
    python scripts/utils/fetch_rice_annotations.py --rapdb-url /path/RAP-MSU_2018-03-29.txt.gz \
        --tf-url /path/Osj_TF_list.txt.gz

Notes:
    - Requires network access to rapdb.dna.affrc.go.jp and planttfdb.gao-lab.org.
    - Some networks receive HTTP 403 from the RAP-DB host for automated requests
      (browser access still works). In that case download the file in a browser
      and pass it via --rapdb-url <local path>, or rely on the bundled snapshot.
    - If those hosts are blocked, download the two files in a browser and pass
      the local paths via --rapdb-url / --tf-url; everything else is identical.
"""

from __future__ import annotations

import argparse
import gzip
import re
import ssl
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

from data_io import SEED, save_table, setup_logging, write_run_info  # noqa: E402

import pandas as pd  # noqa: E402

LOGGER = setup_logging("Rice-skills")

RAPDB_URL = "https://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2018-03-29.txt.gz"
PLANTTFDB_URL = "https://planttfdb.gao-lab.org/download/TF_list/Osj_TF_list.txt.gz"

MSU_FULL_RE = re.compile(r"^LOC_Os\d{2}g\d+$")
RAPDB_FULL_RE = re.compile(r"^Os\d{2}g\d{5,}$")
ANY_GENE_RE = re.compile(r"^(?:LOC_Os\d{2}g\d+|Os\d{2}g\d{5,})$")

# A browser-like User-Agent avoids HTTP 403 blocks from some hosts (the RAP-DB
# server rejects plain automated clients on several networks).
USER_AGENT = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
              "(KHTML, like Gecko) Chrome/126.0 Safari/537.36 Rice-skills/1.0")


def ssl_context(insecure: bool = False) -> ssl.SSLContext:
    """Build an SSL context (certifi bundle; unverified only when --insecure)."""
    if insecure:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        return ctx
    try:
        import certifi  # noqa: PLC0415 - optional dependency

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:  # noqa: BLE001 - fall back to system bundle
        return ssl.create_default_context()


def fetch_text(url: str, timeout: int = 120, insecure: bool = False) -> str:
    """Download a (possibly gzipped) text resource; local paths are supported."""
    local = Path(url)
    if local.is_file():
        LOGGER.info("Using local file: %s", url)
        data = local.read_bytes()
    else:
        LOGGER.info("Downloading: %s", url)
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout,
                                        context=ssl_context(insecure)) as response:
                data = response.read()
        except urllib.error.HTTPError as exc:
            raise SystemExit(
                f"HTTP error {exc.code} while downloading:\n  {url}\n"
                "Hint: the official archive may have moved; check the RAP-DB / PlantTFDB "
                "download pages, or download manually and pass the file path instead."
            ) from exc
        except (urllib.error.URLError, TimeoutError, ssl.SSLError) as exc:
            raise SystemExit(
                f"Download failed:\n  {url}\n  ({exc})\n"
                "Hints:\n"
                "  1. Check your connection/proxy;\n"
                "  2. If the error mentions a certificate problem (e.g. 'self-signed "
                "certificate in certificate chain'), your network intercepts TLS - retry "
                "with --insecure;\n"
                "  3. Or download the file in a browser and pass its local path "
                "(--rapdb-url / --tf-url)."
            ) from exc
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


def apply_bundled_snapshot(raw_path: Path, bundled_name: str) -> bool:
    """Seed the cache from the bundled snapshot (assets/annotations/) if present.

    The bundled files may be gzip-compressed; the cache always stores plain
    text, so the data is decompressed on the way in. Returns True when the
    cache was seeded.
    """
    if raw_path.is_file():
        return False
    bundled = Path(__file__).resolve().parents[2] / "assets" / "annotations" / bundled_name
    if bundled.is_file():
        data = bundled.read_bytes()
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        raw_path.write_bytes(data)
        LOGGER.info("Cache seeded from bundled snapshot: %s", bundled)
        return True
    return False


def parse_rapdb_table(text: str) -> pd.DataFrame:
    """Parse the RAP-DB conversion table ('rapdb / msu7' columns) into pairs."""
    pairs: list[tuple[str, str]] = []
    for line in text.splitlines():
        line = line.strip().lstrip("\ufeff")
        if not line or line.startswith("#"):
            continue
        fields = re.split(r"[\t ]+", line)
        rap = next((f for f in fields if RAPDB_FULL_RE.match(f)), None)
        msu = next((f for f in fields if MSU_FULL_RE.match(f)), None)
        if rap and msu:  # rows with 'None' in the msu column are skipped here
            pairs.append((msu, rap))
    if not pairs:
        raise SystemExit(
            "No MSU/RAP-DB ID pairs found in the conversion table.\n"
            "Hint: check that the file is the RAP-MSU_*.txt table (columns: rapdb, msu7)."
        )
    return pd.DataFrame(pairs, columns=["msu", "rapdb"]).drop_duplicates()


def parse_tf_list(text: str) -> pd.DataFrame:
    """Parse the PlantTFDB TF list (TF_ID / Gene_ID / Family) into rows."""
    rows: list[tuple[str, str, str]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = [f.strip() for f in line.split("\t") if f.strip()]
        if len(fields) < 2:
            fields = line.split()
        if len(fields) < 2:
            continue
        if not rows and any(f.lower() in ("family", "gene_id", "tf_id", "gene")
                            for f in fields[:4]):
            continue  # header line
        gene_idx = next((i for i, f in enumerate(fields) if ANY_GENE_RE.match(f)), None)
        if gene_idx is None:
            continue
        gene = fields[gene_idx]
        family = None
        for field in reversed(fields):  # family = last non-ID, non-TF_ID field
            if field == gene or field == fields[0] or ANY_GENE_RE.match(field):
                continue
            family = field
            break
        rows.append((fields[0], gene, family or "Unknown"))
    if not rows:
        raise SystemExit(
            "No transcription factors parsed from the PlantTFDB file.\n"
            "Hint: expected the Osj_TF_list.txt table (columns: TF_ID / Gene_ID / Family)."
        )
    return pd.DataFrame(rows, columns=["tf_id", "gene_id", "family"]).drop_duplicates("gene_id")


def build_gmt(tf_msu: pd.DataFrame, out_path: Path) -> int:
    """Write a GSEA GMT file with one gene set per TF family."""
    lines = []
    for family, group in tf_msu.groupby("family"):
        genes = sorted(set(group["gene_id"]))
        if len(genes) < 2:
            continue
        name = f"TF_{re.sub(r'[^A-Za-z0-9_]', '_', str(family))}"
        lines.append("\t".join([name, f"PlantTFDB TF family {family}", *genes]))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Download & cache official rice annotation resources "
                    "(RAP-DB ID conversion + PlantTFDB TF list + TF-family GMT)."
    )
    parser.add_argument("--outdir", default="data/annotations",
                        help="Cache/output directory (git-ignored by default)")
    parser.add_argument("--rapdb-url", default=RAPDB_URL,
                        help="RAP-DB conversion table URL or local path")
    parser.add_argument("--tf-url", default=PLANTTFDB_URL,
                        help="PlantTFDB Osj TF list URL or local path")
    parser.add_argument("--insecure", action="store_true",
                        help="Skip TLS verification (use only on networks that "
                             "intercept certificates; data are public files)")
    parser.add_argument("--force", action="store_true", help="Re-download even if cached")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    report: list[str] = ["Rice annotation fetch report", f"outdir: {outdir.resolve()}", ""]

    # --- 1. RAP-DB MSU <-> RAP-DB -------------------------------------------------
    raw_map = outdir / "raw_RAP-MSU_table.txt"
    map_seeded = False
    if args.force:
        raw_map.write_text(fetch_text(args.rapdb_url, insecure=args.insecure),
                           encoding="utf-8")
    elif not apply_bundled_snapshot(raw_map, "raw_RAP-MSU_2018-03-29.txt.gz"):
        raw_map.write_text(fetch_text(args.rapdb_url, insecure=args.insecure),
                           encoding="utf-8")
    else:
        map_seeded = True
    map_df = parse_rapdb_table(raw_map.read_text(encoding="utf-8"))
    save_table(map_df, outdir / "msu_rapdb_map.csv")
    save_table(map_df.rename(columns={"rapdb": "from", "msu": "to"})[["from", "to"]],
               outdir / "rapdb_to_msu_map.csv")
    save_table(map_df.rename(columns={"msu": "from", "rapdb": "to"})[["from", "to"]],
               outdir / "msu_to_rapdb_map.csv")
    map_src = "bundled snapshot" if map_seeded else "official download"
    report.append(f"[OK]   ID conversion table: {len(map_df)} MSU <-> RAP-DB pairs "
                  f"({map_src}) -> msu_rapdb_map.csv / rapdb_to_msu_map.csv / msu_to_rapdb_map.csv")
    LOGGER.info("Parsed %d ID pairs", len(map_df))

    # --- 2. PlantTFDB TF list -----------------------------------------------------
    raw_tf = outdir / "raw_Osj_TF_list.txt"
    tf_seeded = False
    if args.force:
        raw_tf.write_text(fetch_text(args.tf_url, insecure=args.insecure),
                          encoding="utf-8")
    elif not apply_bundled_snapshot(raw_tf, "raw_Osj_TF_list.txt.gz"):
        raw_tf.write_text(fetch_text(args.tf_url, insecure=args.insecure),
                          encoding="utf-8")
    else:
        tf_seeded = True
    tf_raw = parse_tf_list(raw_tf.read_text(encoding="utf-8"))
    save_table(tf_raw, outdir / "rice_TF_PlantTFDB.csv")
    tf_src = "bundled snapshot" if tf_seeded else "official download"
    report.append(f"[OK]   PlantTFDB rice TF list: {len(tf_raw)} TFs, "
                  f"{tf_raw['family'].nunique()} families ({tf_src}) -> rice_TF_PlantTFDB.csv")

    # --- 3. Normalise TF gene IDs to MSU -------------------------------------------
    lut = dict(zip(map_df["rapdb"], map_df["msu"]))
    mapped = tf_raw.copy()
    mapped["gene_id"] = mapped["gene_id"].map(
        lambda g: g if MSU_FULL_RE.match(g) else lut.get(g, g)
    )
    mapped = mapped[["gene_id", "family"]].drop_duplicates("gene_id")
    n_msu = int(mapped["gene_id"].str.startswith("LOC_Os").sum())
    save_table(mapped, outdir / "rice_TF_PlantTFDB_msu.csv")
    report.append(f"[OK]   TF list ready for 06_tf_analysis.R: {n_msu}/{len(mapped)} in "
                  "MSU format -> rice_TF_PlantTFDB_msu.csv")
    if n_msu < len(mapped):
        report.append(f"[WARN] {len(mapped) - n_msu} TF IDs are not MSU; keep the original "
                      "IDs when running 06_tf_analysis.R (it matches gene_id columns).")

    # --- 4. GMT for GSEA -----------------------------------------------------------
    n_sets = build_gmt(mapped, outdir / "rice_TF_families.gmt")
    report.append(f"[OK]   GSEA gene sets: {n_sets} TF families with >=2 genes -> "
                  "rice_TF_families.gmt")

    (outdir / "fetch_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    write_run_info(outdir, "fetch_rice_annotations.py",
                   params={"rapdb_url": args.rapdb_url, "tf_url": args.tf_url,
                           "force": bool(args.force), "insecure": bool(args.insecure)},
                   packages=("pandas", "numpy"))

    print("\n".join(report))
    print(f"\nCache directory: {outdir.resolve()}")
    print("Next steps:")
    print("  Rscript scripts/transcriptomics/06_tf_analysis.R --de-all <de_all.csv> \\")
    print(f"      --tf-annotation {outdir}/rice_TF_PlantTFDB_msu.csv")
    print("  Rscript scripts/transcriptomics/04_gsea.R --de-all <de_all.csv> \\")
    print(f"      --gmt {outdir}/rice_TF_families.gmt")
    print("  python scripts/utils/adapt_company_data.py --mode transcriptome \\")
    print(f"      --input <vendor.csv> --id-map {outdir}/rapdb_to_msu_map.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
