"""adapt_company_data.py -- adapt vendor data to the Rice-skills format.

Title:    Rice-skills vendor-data adapter
Description:
    Sequencing companies (RNA-seq) and metabolomics platforms all deliver data
    in their own column-naming conventions. This utility converts those files
    into the strict Rice-skills input format:

        transcriptome : gene_id + sample columns (integer counts)
        metabolome    : metabolite_id + metabolite_name/class/kegg_id + samples

    It performs, with a full adaptation report:
      1. delimiter and encoding detection (CSV / TSV / TXT; XLSX if openpyxl);
      2. column-name normalisation (English and Chinese vendor headers);
      3. gene-ID system detection (MSU LOC_Os / RAP-DB Os..g / Ensembl Os..t)
         and optional conversion via a user-supplied ID map;
      4. count-value sanity checks (non-integer values, FPKM/TPM-like ranges);
      5. sample-name cleaning and optional metadata-overlap validation.

Usage:
    # transcriptome (vendor file with a 'Gene ID' column and RAP-DB identifiers)
    python scripts/utils/adapt_company_data.py --mode transcriptome \
        --input company_counts.xlsx --metadata metadata.csv \
        --id-map msu_rapdb_map.csv --outdir adapted/

    # metabolome (vendor table with Chinese annotation headers)
    python scripts/utils/adapt_company_data.py --mode metabolome \
        --input 代谢组数据.csv --outdir adapted/

Inputs:
    --input       vendor file (.csv/.tsv/.txt/.xlsx)
    --metadata    optional sample metadata for overlap validation
    --id-map      optional two-column CSV (from,to) for gene-ID conversion
Outputs:
    adapted/counts.csv or adapted/metabolites.csv + adapted/adaptation_report.txt
"""

from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

from data_io import RiceDataError, SEED, save_table, setup_logging, write_run_info  # noqa: E402

LOGGER = setup_logging("Rice-skills")

# ---------------------------------------------------------------------------
# Vendor-header dictionaries (English + Chinese)
# ---------------------------------------------------------------------------
GENE_ID_CANDIDATES = [
    "gene_id", "geneid", "gene id", "gene", "id", "gene_name", "gene name",
    "genes", "gene symbol", "symbol", "probe id", "probeid", "transcript_id",
    "基因id", "基因编号", "基因", "基因名称", "基因符号",
]
METAB_ID_CANDIDATES = [
    "metabolite_id", "metaboliteid", "metabolite id", "id", "compound_id",
    "compound id", "index", "编号", "代谢物id", "物质id",
]
METAB_NAME_CANDIDATES = [
    "metabolite_name", "name", "metabolite", "metabolites", "compound",
    "compound name", "compound_name", "物质", "代谢物", "代谢物名称",
    "化合物", "化合物名称", "名称",
]
METAB_CLASS_CANDIDATES = [
    "class", "chemical class", "class i", "category", "superclass",
    "分类", "化学分类", "物质类别", "类别", "化合物分类",
]
METAB_KEGG_CANDIDATES = [
    "kegg_id", "kegg", "kegg id", "kegg compound", "kegg compound id",
    "kegg pathway", "kegg编号",
]
METAB_EXTRA_ANNOT = ["cas", "cas号", "cas id", "formula", "molecular formula",
                     "分子式", "smiles", "pubchem", "hmdb", "inchikey", "inchi"]

MSU_RE = re.compile(r"^LOC_Os\d{2}g\d+", re.IGNORECASE)
RAPDB_RE = re.compile(r"^Os\d{2}g\d+", re.IGNORECASE)
ENSEMBL_RE = re.compile(r"^Os\d{2}t\d+", re.IGNORECASE)


def _norm(text: str) -> str:
    """Normalise a header for matching: lowercase, collapse spaces/underscores."""
    text = unicodedata.normalize("NFKC", str(text)).strip().lower()
    return re.sub(r"[\s_\-]+", " ", text)


def _find_column(columns: list[str], candidates: list[str]) -> str | None:
    """Return the first column whose normalised name matches a candidate."""
    norm_map = {_norm(c): c for c in columns}
    for cand in candidates:
        if _norm(cand) in norm_map:
            return norm_map[_norm(cand)]
    return None


def read_vendor_table(path: Path) -> pd.DataFrame:
    """Read a vendor table with delimiter/encoding detection.

    Supports .csv/.tsv/.txt and .xlsx (requires openpyxl).
    """
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        try:
            return pd.read_excel(path)
        except ImportError as exc:
            raise RiceDataError(
                f"Reading Excel files requires openpyxl ({exc}).\n"
                "Hint: pip install openpyxl, or export the sheet as CSV and re-run."
            ) from exc
        except Exception as exc:  # noqa: BLE001
            raise RiceDataError(
                f"Failed to read Excel file: {path}\n{exc}\n"
                "Hint: make sure the file is a valid .xlsx and the sheet name is default."
            ) from exc
    for sep in (",", "\t", ";"):
        try:
            df = pd.read_csv(path, sep=sep)
            if df.shape[1] > 1:
                LOGGER.info("Detected delimiter %r in %s", sep, path.name)
                return df
        except Exception:  # noqa: BLE001 - try next delimiter
            continue
    raise RiceDataError(
        f"Could not parse {path} as CSV/TSV.\n"
        "Hint: export the vendor table as comma-separated CSV (UTF-8)."
    )


def detect_id_system(values: pd.Series) -> str:
    """Detect the rice gene-ID system of a column."""
    sample = values.astype(str).head(200)
    msu = sample.str.match(MSU_RE).mean()
    rapdb = sample.str.match(RAPDB_RE).mean()
    ens = sample.str.match(ENSEMBL_RE).mean()
    if msu >= 0.5:
        return "MSU (LOC_Os...)"
    if rapdb >= 0.5:
        return "RAP-DB (Os01g...)"
    if ens >= 0.5:
        return "Ensembl (Os01t...)"
    return "unknown"


def apply_id_map(ids: pd.Series, map_path: Path) -> tuple[pd.Series, dict]:
    """Convert gene IDs using a two-column (from,to) mapping CSV."""
    mapping = pd.read_csv(map_path)
    if mapping.shape[1] < 2:
        raise RiceDataError("--id-map must be a CSV with two columns: from,to")
    frm, to = mapping.columns[0], mapping.columns[1]
    lut = dict(zip(mapping[frm].astype(str), mapping[to].astype(str)))
    converted = ids.astype(str).map(lut)
    n_unmapped = int(converted.isna().sum())
    unmapped_examples = ids[converted.isna()].head(5).tolist()
    converted = converted.fillna(ids.astype(str))
    stats = {"mapped": int(len(ids) - n_unmapped), "unmapped": n_unmapped,
             "unmapped_examples": unmapped_examples}
    return converted, stats


def adapt_transcriptome(df: pd.DataFrame, report: list[str],
                        id_map: Path | None) -> pd.DataFrame:
    """Normalise a vendor count matrix to gene_id + integer sample columns."""
    id_col = _find_column(list(df.columns), GENE_ID_CANDIDATES)
    if id_col is None:
        id_col = df.columns[0]
        report.append(f"[WARN] No recognised gene-ID column; using first column '{id_col}'.")
    else:
        report.append(f"[OK]   Gene-ID column: '{id_col}'.")

    df = df.rename(columns={id_col: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str).str.strip()

    id_system = detect_id_system(df["gene_id"])
    report.append(f"[INFO] Gene-ID system detected: {id_system}.")
    if id_system.startswith("MSU"):
        report.append("[OK]   MSU IDs are directly compatible with org.Osativa.eg.db.")
    else:
        report.append(
            "[WARN] IDs are not MSU (LOC_Os...). Enrichment steps expect SYMBOL=LOC_Os IDs.\n"
            "       Option A: convert with RAP-DB Gene ID Converter "
            "(https://rapdb.dna.affrc.go.jp/tools/converter);\n"
            "       Option B: pass a two-column mapping file via --id-map."
        )
    if id_map is not None:
        df["gene_id"], stats = apply_id_map(df["gene_id"], id_map)
        report.append(
            f"[OK]   Applied --id-map: {stats['mapped']} converted, "
            f"{stats['unmapped']} unmapped (kept original)"
            + (f", e.g. {stats['unmapped_examples']}" if stats["unmapped"] else "")
        )

    sample_cols = [c for c in df.columns if c != "gene_id"]
    cleaned = {c: str(c).strip() for c in sample_cols}
    df = df.rename(columns=cleaned)

    values = df[[c for c in df.columns if c != "gene_id"]]
    non_numeric = [c for c in values.columns
                   if not pd.api.types.is_numeric_dtype(values[c])]
    if non_numeric:
        raise RiceDataError(
            f"Non-numeric sample column(s) detected: {non_numeric[:5]}.\n"
            "Hint: vendor files sometimes append comments/units; remove them and re-run."
        )
    if values.isna().any().any():
        n_na = int(values.isna().sum().sum())
        report.append(f"[WARN] {n_na} missing values in the count matrix; "
                      "DESeq2 requires complete counts (impute or remove those genes).")

    rounded = values.round()
    diff = (values - rounded).abs()
    n_noninteger = int((diff > 1e-6).sum().sum())
    if n_noninteger:
        frac = n_noninteger / values.size
        report.append(
            f"[WARN] {n_noninteger} non-integer values ({frac:.1%} of cells). "
            "DESeq2/DESeqDataSetFromMatrix needs raw counts.\n"
            "       If this file stores rounded counts as floats, run again with --round-counts.\n"
            "       If it is TPM/FPKM/CleanData, request the raw count matrix from the vendor."
        )
    med = float(np.nanmedian(values.to_numpy(dtype=float)))
    if med < 10:
        report.append(f"[WARN] Median value is {med:.2f}; TPM/FPKM-like matrices "
                      "are not suitable for DESeq2 (request raw counts).")
    return df


def adapt_metabolome(df: pd.DataFrame, report: list[str]) -> pd.DataFrame:
    """Normalise a vendor metabolite table to the standard annotated format."""
    cols = list(df.columns)
    id_col = _find_column(cols, METAB_ID_CANDIDATES)
    name_col = _find_column(cols, METAB_NAME_CANDIDATES)
    class_col = _find_column(cols, METAB_CLASS_CANDIDATES)
    kegg_col = _find_column(cols, METAB_KEGG_CANDIDATES)
    extra_cols = [c for c in cols if _norm(c) in {_norm(e) for e in METAB_EXTRA_ANNOT}]

    annot = {}
    if id_col:
        annot[id_col] = "metabolite_id"; report.append(f"[OK]   ID column: '{id_col}'.")
    if name_col:
        annot[name_col] = "metabolite_name"; report.append(f"[OK]   Name column: '{name_col}'.")
    if class_col:
        annot[class_col] = "class"; report.append(f"[OK]   Class column: '{class_col}'.")
    if kegg_col:
        annot[kegg_col] = "kegg_id"; report.append(f"[OK]   KEGG column: '{kegg_col}'.")
    for c in extra_cols:
        if c not in annot:
            annot[c] = _norm(c).replace(" ", "_")
            report.append(f"[INFO] Extra annotation column kept: '{c}'.")

    df = df.rename(columns=annot)
    if "metabolite_id" not in df.columns:
        df.insert(0, "metabolite_id", [f"M{i + 1:04d}" for i in range(len(df))])
        report.append("[INFO] No metabolite ID column found; generated M0001... IDs.")
    if "metabolite_name" not in df.columns:
        report.append("[WARN] No metabolite name column recognised; check the report "
                      "and rename manually if a name column exists.")
    if "kegg_id" not in df.columns:
        report.append("[WARN] No KEGG-ID column found; pathway analysis (03/04) will skip "
                      "metabolites without annotations.")
    if "class" not in df.columns:
        report.append("[WARN] No chemical-class column found; classification step (04) "
                      "will be limited.")

    # Keep known annotations + extra annotation columns; the rest are samples.
    known_annot = [c for c in ("metabolite_id", "metabolite_name", "class", "kegg_id",
                               "cas", "cas_id", "formula", "molecular_formula", "smiles",
                               "pubchem", "hmdb", "inchikey", "inchi") if c in df.columns]
    sample_cols = [c for c in df.columns if c not in known_annot]
    cleaned = {c: str(c).strip() for c in sample_cols}
    df = df.rename(columns=cleaned)
    sample_cols = [str(c).strip() for c in sample_cols]

    values = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    n_na = int(values.isna().sum().sum())
    pre_missing = int(df[sample_cols].isna().sum().sum())
    if n_na - pre_missing > 0:
        report.append(f"[WARN] {n_na - pre_missing} non-numeric cells were converted to NA; "
                      "check for text markers such as '<LOD' in the vendor file.")
    frac = n_na / max(values.size, 1)
    if frac > 0.2:
        report.append(f"[WARN] {frac:.1%} missing values — consider stronger filtering "
                      "before analysis (keep metabolites missing in <20% of samples).")
    else:
        report.append(f"[OK]   Missing-value fraction: {frac:.1%} (handled by 01_preprocess.py).")
    df[values.columns] = values
    ordered = known_annot + sample_cols
    return df[ordered]


def check_metadata_overlap(df: pd.DataFrame, metadata_path: Path,
                           report: list[str]) -> None:
    """Validate that the adapted sample columns match the metadata sample IDs."""
    meta = pd.read_csv(metadata_path)
    meta_ids = set(meta.iloc[:, 0].astype(str).str.strip())
    sample_cols = [c for c in df.columns
                   if c not in ("gene_id", "metabolite_id", "metabolite_name",
                                "class", "kegg_id", "cas", "formula")]
    matrix_ids = set(map(str, sample_cols))
    shared = matrix_ids & meta_ids
    report.append(
        f"[INFO] Metadata overlap: {len(shared)}/{len(matrix_ids)} matrix samples found "
        f"in {metadata_path.name}."
    )
    only_matrix = sorted(matrix_ids - meta_ids)
    only_meta = sorted(meta_ids - matrix_ids)
    if only_matrix:
        report.append(f"[WARN] Samples only in the data file: {only_matrix[:8]}"
                      + (" ..." if len(only_matrix) > 8 else ""))
    if only_meta:
        report.append(f"[WARN] Samples only in metadata: {only_meta[:8]}"
                      + (" ..." if len(only_meta) > 8 else ""))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Adapt vendor (sequencing company / metabolomics platform) data "
                    "to the Rice-skills input format."
    )
    parser.add_argument("--mode", required=True, choices=["transcriptome", "metabolome"])
    parser.add_argument("--input", required=True, help="Vendor file (.csv/.tsv/.txt/.xlsx)")
    parser.add_argument("--metadata", default=None, help="Optional metadata.csv for overlap check")
    parser.add_argument("--id-map", dest="id_map", default=None,
                        help="Optional two-column CSV (from,to) for gene-ID conversion")
    parser.add_argument("--round-counts", dest="round_counts", action="store_true",
                        help="Round float counts to integers (only if they are true counts)")
    parser.add_argument("--outdir", default="adapted", help="Output directory")
    args = parser.parse_args(argv)

    input_path = Path(args.input)
    if not input_path.is_file():
        LOGGER.error("Input file not found: %s", args.input)
        print("Hint: check the path; company files are often still inside a ZIP archive.")
        return 2

    report: list[str] = [f"Vendor data adaptation report (mode = {args.mode})",
                         f"input  : {input_path}", f"outdir : {args.outdir}", ""]
    try:
        raw = read_vendor_table(input_path)
        report.append(f"[OK]   Loaded {raw.shape[0]} rows x {raw.shape[1]} columns.")
        report.append("")

        if args.mode == "transcriptome":
            out = adapt_transcriptome(
                raw, report, Path(args.id_map) if args.id_map else None
            )
            if args.round_counts:
                sample_cols = [c for c in out.columns if c != "gene_id"]
                out[sample_cols] = out[sample_cols].round().astype("int64")
                report.append("[OK]   Counts rounded to integers (--round-counts).")
            out_name = "counts.csv"
        else:
            out = adapt_metabolome(raw, report)
            out_name = "metabolites.csv"

        if args.metadata:
            check_metadata_overlap(out, Path(args.metadata), report)

        out_path = save_table(out, Path(args.outdir) / out_name)
        annot_names = {"gene_id", "metabolite_id", "metabolite_name", "class", "kegg_id",
                       "cas", "cas_id", "formula", "molecular_formula", "smiles",
                       "pubchem", "hmdb", "inchikey", "inchi"}
        n_samples = sum(1 for c in out.columns if c not in annot_names)
        report += ["", f"[OK]   Standardised file written: {out_path}",
                   f"[OK]   Rows: {out.shape[0]}, sample columns: {n_samples}"]
        report_path = Path(args.outdir) / "adaptation_report.txt"
        report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
        write_run_info(Path(args.outdir), "adapt_company_data.py",
                       params={"mode": args.mode, "input": str(input_path),
                               "round_counts": bool(args.round_counts)},
                       packages=("pandas", "numpy", "openpyxl"))

        LOGGER.info("Adaptation complete:")
        for line in report:
            print(line)
        print(f"\nReport saved to: {report_path}")
        print("\nNext step: run the pipeline, e.g.")
        if args.mode == "transcriptome":
            print(f"  python scripts/transcriptomics/01_qc_filter.py "
                  f"--counts {out_path} --metadata <metadata.csv>")
        else:
            print(f"  python scripts/metabolomics/01_preprocess.py "
                  f"--abundance {out_path} --metadata <metadata.csv>")
        return 0
    except RiceDataError as exc:
        LOGGER.error("Adaptation failed: %s", exc)
        print("\nThe file could not be adapted automatically. Common fixes:")
        print("  1. Export the vendor table as CSV (UTF-8) with a single header row;")
        print("  2. Remove merged cells, footnotes and unit rows below the header;")
        print("  3. Keep one ID column (gene or metabolite) and one column per sample.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
