"""04_uniprot_mapping.py -- Automated protein_id -> gene_id mapping via UniProt.

Title:    Rice-skills proteome annotation mapping (step 04)
Description:
    Builds the protein_id -> MSU gene_id (LOC_Osxxgxxxxxx) map used by
    01_preprocess.py / 02_diff_protein.R automatically, from the UniProt
    reference proteome UP000059680 (Oryza sativa subsp. japonica, the same
    proteome behind the MSU/RGAP annotation). Replaces the previous workflow
    of asking the user to supply protein_gene_map.csv by hand.

    Data flow:
      1. fetch UniProtKB TSV for proteome UP000059680 (cached; --force to
         refresh). Fields: accession, id, protein_name, gene_names,
         gene_oligoid, length.
      2. extract the OrderedLocusNames token (LOC_Osxxgxxxxxx) per entry;
         fall back to any LOC_Os token found in gene_names.
      3. write the two-column map (protein_id,gene_id) consumed by the
         module-1 enrichment step (03_enrichment.R / 04_gsea.R via the
         protein->gene bridge documented in docs/proteomics_usage.md),
         plus unmatched-protein diagnostics.

    Network-restricted environments (TLS interception, offline): pass
    --tsv-file with a TSV previously downloaded in a browser from
    https://rest.uniprot.org/uniprotkb/stream?query=(proteome:UP000059680)&format=tsv
    or run with --insecure. --demo runs a built-in synthetic self-test.

Usage:
    python scripts/proteomics/04_uniprot_mapping.py \
        --outdir results/proteomics/04_mapping
    # refresh the cache
    python scripts/proteomics/04_uniprot_mapping.py --outdir ... --force
    # fully offline (browser-downloaded TSV)
    python scripts/proteomics/04_uniprot_mapping.py --outdir ... \
        --tsv-file data/annotations/uniprot_up000059680.tsv
    # self-test with synthetic data (no network)
    python scripts/proteomics/04_uniprot_mapping.py --demo --outdir results/proteomics/04_mapping_demo

Inputs:
    --tsv-file   Optional UniProtKB TSV (tab-separated, with UniProt headers).
    --cache-dir  Cache directory (default: data/annotations).

Outputs (written to --outdir):
    protein_gene_map.csv    two-column map: protein_id (UniProt accession),
                            gene_id (MSU LOC_Os...)
    mapping_stats.txt       counts, source, timestamp
    unmatched.csv           accessions with no LOC_Os token (diagnostic)

Notes:
    - MSU "LOC_" prefix is stripped/normalised consistently
      (both LOC_Os01g01010 and Os01g01010 map to Os01g01010 style output
      when --strip-loc is set, default on, matching demo_data convention).
    - One UniProt entry may map to several locus tokens (isoforms/gene
      models); duplicates are emitted as separate rows, flagged in stats.
    - Upstream data: UniProtKB (CC BY 4.0). We store only the fetched
      annotation snapshot, not UniProt source code.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.request
import urllib.error
import datetime
from pathlib import Path

UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=proteome:UP000059680&format=tsv"
    "&fields=accession,id,protein_name,gene_names,length"
)
CACHE_NAME = "uniprot_up000059680.tsv"

# MSU OrderedLocusNames style: LOC_Os04g55560 (7 digits after 'g');
# (?!\d) prevents truncating the 10-digit RAP-DB style.
LOC_MSU_RE = re.compile(r"LOC[_\s]?Os(\d{2})g(\d{5})(?!\d)", re.IGNORECASE)
# RAP-DB new style: Os04g0649100 (10 digits after 'g'); must NOT be
# truncated by the MSU pattern, hence a separate 7-digit-cap regex.
LOC_RAPDB_RE = re.compile(r"(?<!LOC)(?<![_A-Za-z])Os(\d{2})g(\d{7})(?:\.\d+)?", re.IGNORECASE)


# ----------------------------------------------------------------- helpers
def ssl_context(insecure: bool = False):
    """SSL context: certifi bundle by default; unverified with --insecure."""
    if insecure:
        import ssl
        return ssl._create_unverified_context()
    try:
        import certifi  # optional dependency
        import ssl
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        import ssl
        return ssl.create_default_context()


def fetch_tsv(url: str, insecure: bool, timeout: int) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Rice-skills/3.1"})
    with urllib.request.urlopen(request, timeout=timeout, context=ssl_context(insecure)) as resp:
        return resp.read().decode("utf-8", errors="replace")


def normalise_locus(token: str, strip_loc: bool) -> str:
    """Normalise a locus token: drop version suffix (Os02g0000200.1 ->
    Os02g0000200) and optionally the LOC_ prefix."""
    tok = token.strip().replace(" ", "_")
    tok = re.sub(r"\.\d+$", "", tok)  # strip version suffix
    if strip_loc and tok.upper().startswith("LOC_"):
        tok = tok[4:]
    return tok


def extract_loci(gene_names: str, oligoid: str) -> list[tuple[str, str]]:
    """Return (token, style) pairs; style is 'msu' (LOC_Os04g55560) or
    'rapdb' (Os04g0649100). MSU tokens take precedence downstream."""
    tokens: list[tuple[str, str]] = []
    for field in (oligoid, gene_names):
        if not field:
            continue
        for m in LOC_MSU_RE.finditer(field):
            raw = m.group(0).replace(" ", "_")
            if not raw.upper().startswith("LOC_"):
                raw = "LOC_" + raw[3:] if not raw.upper().startswith("LOC") else raw
                raw = raw if raw.upper().startswith("LOC_") else "LOC_" + raw
            tokens.append((raw, "msu"))
        for m in LOC_RAPDB_RE.finditer(field):
            tokens.append((m.group(0), "rapdb"))
    # dedupe on (upper token, style), keep order
    seen: set[tuple[str, str]] = set()
    out: list[tuple[str, str]] = []
    for t, s in tokens:
        k = (t.upper(), s)
        if k not in seen:
            seen.add(k)
            out.append((t, s))
    return out


def parse_uniprot_tsv(text: str, strip_loc: bool) -> tuple[list[dict], list[str], dict]:
    """Parse a UniProt TSV dump into mapping rows + unmatched + style stats."""
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    rows: list[dict] = []
    unmatched: list[str] = []
    style_count = {"msu": 0, "rapdb": 0}
    for rec in reader:
        acc = (rec.get("Entry") or rec.get("accession") or "").strip()
        if not acc:
            continue
        gene_names = rec.get("Gene Names") or rec.get("gene_names") or ""
        oligoid = rec.get("Gene Names (OLN)") or rec.get("gene_oligoid") or ""
        loci = extract_loci(gene_names, oligoid)
        if not loci:
            unmatched.append(acc)
            continue
        # MSU style wins; fall back to RAP-DB style when no MSU token exists.
        msu = [(t, s) for t, s in loci if s == "msu"]
        chosen = msu if msu else loci
        # normalise, then dedupe: OLN "LOC_Os01g00001" and Gene Names
        # "Os01g00001.1" are the same locus and must yield one row.
        seen: set[str] = set()
        for loc, style in chosen:
            nid = normalise_locus(loc, strip_loc)
            key = nid.upper()
            if key in seen:
                continue
            seen.add(key)
            rows.append({"protein_id": acc, "gene_id": nid, "id_style": style})
            style_count[style] += 1
    return rows, unmatched, style_count


# ------------------------------------------------------------------- demo
DEMO_TSV = "\t".join(["Entry", "Entry Name", "Protein names", "Gene Names", "Length"]) + "\n" + "\n".join(
    [
        # realistic row: MSU + RAP-DB styles coexist for one gene
        "\t".join(["OsProt00001", "P1_ORYSJ", "Synthetic protein 1", "P1 SHAT1 LOC_Os01g00001 Os01g0000100", "350"]),
        # RAP-DB style only
        "\t".join(["OsProt00002", "P2_ORYSJ", "Synthetic protein 2", "P2 Os02g0000200.1", "220"]),
        # MSU style only
        "\t".join(["OsProt00003", "P3_ORYSJ", "Synthetic protein 3", "P3 LOC_Os03g00003", "512"]),
        # no locus token at all
        "\t".join(["OsProt00004", "P4_ORYSJ", "No-locus protein", "SOMEGENE1", "180"]),
    ]
)


# ------------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="results/proteomics/04_mapping")
    ap.add_argument("--tsv-file", default=None, help="offline UniProt TSV (skips download)")
    ap.add_argument("--cache-dir", default="data/annotations")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    ap.add_argument("--insecure", action="store_true", help="skip TLS verification")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--keep-loc", dest="strip_loc", action="store_false",
                    help="keep the LOC_ prefix (default strips to Os01g01010)")
    ap.add_argument("--demo", action="store_true", help="run built-in synthetic self-test")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        tsv_text = DEMO_TSV
        source = "demo (synthetic)"
    elif args.tsv_file:
        tsv_text = Path(args.tsv_file).read_text(encoding="utf-8")
        source = f"file:{args.tsv_file}"
    else:
        cache = Path(args.cache_dir) / CACHE_NAME
        if cache.exists() and not args.force:
            tsv_text = cache.read_text(encoding="utf-8")
            source = f"cache:{cache}"
        else:
            print(f"[04] downloading UniProt proteome UP000059680 ...", flush=True)
            try:
                tsv_text = fetch_tsv(UNIPROT_URL, args.insecure, args.timeout)
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                print(
                    "[04] download failed:", exc,
                    "\n[04] fallback options:\n"
                    "     --insecure            TLS-intercepted networks\n"
                    "     --tsv-file <file>     browser-downloaded TSV\n"
                    "     --demo                synthetic self-test",
                    flush=True,
                )
                return 2
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(tsv_text, encoding="utf-8")
            source = f"download:{UNIPROT_URL}"

    rows, unmatched, style_count = parse_uniprot_tsv(tsv_text, args.strip_loc)

    map_path = outdir / "protein_gene_map.csv"
    with map_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["protein_id", "gene_id", "id_style"])
        w.writeheader()
        w.writerows(rows)

    un_path = outdir / "unmatched.csv"
    with un_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["protein_id"])
        w.writerows([[u] for u in unmatched])

    n_prot = len({r["protein_id"] for r in rows})
    n_dupes = len(rows) - len({(r["protein_id"], r["gene_id"].upper()) for r in rows})
    stats = [
        f"source: {source}",
        f"timestamp: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"UniProt entries parsed: {n_prot + len(unmatched)}",
        f"entries with >=1 locus: {n_prot}",
        f"mapping rows written: {len(rows)} (MSU style: {style_count['msu']}, RAP-DB style: {style_count['rapdb']})",
        f"unmatched entries: {len(unmatched)} (see unmatched.csv)",
        f"strip LOC_ prefix: {args.strip_loc}",
    ]
    (outdir / "mapping_stats.txt").write_text("\n".join(stats) + "\n", encoding="utf-8")
    print("\n".join("[04] " + s for s in stats))
    print(f"[04] wrote {map_path}")

    if args.demo:
        # self-test assertions
        got = {(r["protein_id"], r["gene_id"]) for r in rows}
        style = {r["protein_id"]: r["id_style"] for r in rows}
        assert ("OsProt00001", "Os01g00001") in got, f"MSU token extraction failed: {got}"
        assert ("OsProt00002", "Os02g0000200") in got, f"RAP-DB style extraction failed: {got}"
        assert ("OsProt00003", "Os03g00003") in got, f"MSU-only extraction failed: {got}"
        assert style["OsProt00001"] == "msu" and style["OsProt00002"] == "rapdb", style
        # MSU must win when both styles are present (no RAP-DB row for P1)
        assert ("OsProt00001", "Os01g0000100") not in got, "MSU preference violated"
        assert unmatched == ["OsProt00004"], f"unmatched detection failed: {unmatched}"
        assert n_dupes == 0
        print("[04] demo self-test PASSED (3 mapped, 1 unmatched, MSU-preference OK, "
              "LOC_ normalisation OK, no double-style duplicates)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
