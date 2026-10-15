"""data_io.py -- Rice-skills shared data I/O and reproducibility helpers.

Title:    Rice-skills shared data loading, validation and logging utilities
Description:
    Common helpers used by every Python script in the Rice-skills skill:
    standardised loaders for count matrices, metabolite abundance matrices and
    sample metadata, sample-alignment validation, table export, logging and
    run provenance (run_info.txt with versions, parameters and the random seed).

Data format conventions (see examples/README.md):
    - Transcriptome count matrix (CSV): first column ``gene_id``, remaining
      columns are samples; integer counts, rows = genes.
    - Metabolome abundance matrix (CSV): first column ``metabolite_id``,
      recommended annotation columns ``metabolite_name``, ``class``,
      ``kegg_id``; remaining columns are samples.
    - Sample metadata (CSV): first column ``sample_id`` (must match matrix
      column names exactly), required column ``group`` (e.g. control/treatment),
      optional columns ``batch``, ``variety`` and trait columns for WGCNA.

Usage:
    from data_io import load_expression_matrix, load_metadata, save_table
    counts = load_expression_matrix("examples/transcriptome_data/counts.csv")
    meta = load_metadata("examples/metadata.csv")
"""

from __future__ import annotations

import logging
import os
import platform
import random
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

SEED = 42

LOGGER_NAME = "Rice-skills"


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def set_seed(seed: int = SEED) -> int:
    """Seed Python and NumPy RNGs; returns the seed used."""
    random.seed(seed)
    np.random.seed(seed)
    return seed


def setup_logging(
    name: str = LOGGER_NAME,
    log_file: str | os.PathLike | None = None,
    level: int = logging.INFO,
) -> logging.Logger:
    """Configure and return a logger writing to stdout (and optionally a file).

    Args:
        name: Logger name.
        log_file: Optional path for an additional file handler.
        level: Logging level.

    Returns:
        Configured logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(formatter)
    logger.addHandler(stream)
    if log_file is not None:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    logger.propagate = False
    return logger


# ---------------------------------------------------------------------------
# Friendly errors
# ---------------------------------------------------------------------------

class RiceDataError(ValueError):
    """Raised when input data violate the Rice-skills format conventions."""


def _check_file(path: str | os.PathLike, what: str) -> Path:
    """Return the path if it exists, otherwise raise a friendly error."""
    p = Path(path)
    if not p.is_file():
        raise RiceDataError(
            f"{what} not found: {path}\n"
            "Hint: check the path, and make sure the previous pipeline step ran successfully."
        )
    return p


def _first_column_as_index(df: pd.DataFrame, id_col: str, what: str) -> pd.DataFrame:
    """Set the ID column as index; validate uniqueness and missing values."""
    if df.columns[0] != id_col:
        df = df.rename(columns={df.columns[0]: id_col})
    if df[id_col].duplicated().any():
        dupes = df.loc[df[id_col].duplicated(), id_col].head(5).tolist()
        raise RiceDataError(
            f"{what} contains duplicate '{id_col}' values, e.g. {dupes}. "
            "IDs must be unique."
        )
    if df[id_col].isna().any():
        raise RiceDataError(f"{what} contains missing '{id_col}' values.")
    return df.set_index(id_col)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_expression_matrix(
    path: str | os.PathLike,
    id_col: str = "gene_id",
    what: str = "Expression matrix",
) -> pd.DataFrame:
    """Load a gene expression matrix (rows = genes, columns = samples).

    Args:
        path: CSV path. First column is the gene ID (any name; renamed to
            ``id_col``), remaining columns are sample IDs.
        id_col: Name to assign to the ID column.
        what: Human-readable description for error messages.

    Returns:
        DataFrame indexed by gene ID with sample columns.

    Raises:
        RiceDataError: If the file is missing, unparsable or malformed.
    """
    p = _check_file(path, what)
    try:
        df = pd.read_csv(p)
    except Exception as exc:  # noqa: BLE001 - friendly re-raise
        raise RiceDataError(
            f"Failed to parse {what}: {path}\n{exc}\n"
            "Hint: expected a comma-separated file with a header row."
        ) from exc
    if df.shape[1] < 2:
        raise RiceDataError(
            f"{what} must contain an ID column and at least one sample column: {path}"
        )
    df = _first_column_as_index(df, id_col, what)
    sample_cols = [c for c in df.columns]
    non_numeric = [c for c in sample_cols if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        raise RiceDataError(
            f"{what} has non-numeric sample column(s): {non_numeric[:5]}. "
            "Expression values must be numeric (counts or normalised intensities)."
        )
    return df


def load_metabolome(
    path: str | os.PathLike,
    id_col: str = "metabolite_id",
    annotation_cols: Sequence[str] = ("metabolite_name", "class", "kegg_id"),
    what: str = "Metabolome matrix",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load a metabolite abundance matrix (rows = metabolites, columns = samples).

    Annotation columns (``metabolite_name``, ``class``, ``kegg_id`` ...) are
    separated from the numeric sample matrix.

    Args:
        path: CSV path.
        id_col: Name of the ID column.
        annotation_cols: Non-numeric annotation columns to split off.
        what: Human-readable description for error messages.

    Returns:
        Tuple ``(abundance, annotation)`` where ``abundance`` is indexed by
        metabolite ID with numeric sample columns, and ``annotation`` holds the
        annotation columns (may be empty if none were provided).

    Raises:
        RiceDataError: If the file is missing, unparsable or malformed.
    """
    p = _check_file(path, what)
    try:
        df = pd.read_csv(p)
    except Exception as exc:  # noqa: BLE001 - friendly re-raise
        raise RiceDataError(
            f"Failed to parse {what}: {path}\n{exc}\n"
            "Hint: expected a comma-separated file with a header row."
        ) from exc
    if df.shape[1] < 2:
        raise RiceDataError(
            f"{what} must contain an ID column and at least one sample column: {path}"
        )
    df = _first_column_as_index(df, id_col, what)
    present_annot = [c for c in annotation_cols if c in df.columns]
    annotation = df[present_annot].copy()
    abundance = df.drop(columns=present_annot)
    non_numeric = [c for c in abundance.columns if not pd.api.types.is_numeric_dtype(abundance[c])]
    if non_numeric:
        raise RiceDataError(
            f"{what} has non-numeric sample column(s): {non_numeric[:5]}. "
            "Abundance values must be numeric (missing values should be empty/NA)."
        )
    return abundance, annotation


def load_metadata(
    path: str | os.PathLike,
    id_col: str = "sample_id",
    required_cols: Sequence[str] = ("group",),
    what: str = "Sample metadata",
) -> pd.DataFrame:
    """Load sample metadata (first column = sample_id).

    Args:
        path: CSV path.
        id_col: Name of the sample ID column.
        required_cols: Columns that must be present (default: ``group``).
        what: Human-readable description for error messages.

    Returns:
        DataFrame indexed by sample ID.

    Raises:
        RiceDataError: If the file is missing or required columns are absent.
    """
    p = _check_file(path, what)
    try:
        df = pd.read_csv(p)
    except Exception as exc:  # noqa: BLE001 - friendly re-raise
        raise RiceDataError(
            f"Failed to parse {what}: {path}\n{exc}\n"
            "Hint: expected a comma-separated file with a header row."
        ) from exc
    df = _first_column_as_index(df, id_col, what)
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise RiceDataError(
            f"{what} is missing required column(s): {missing}.\n"
            f"Available columns: {list(df.columns)}\n"
            "Hint: the 'group' column defines the conditions used for comparisons."
        )
    return df


# ---------------------------------------------------------------------------
# Validation and alignment
# ---------------------------------------------------------------------------

def align_samples(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    matrix_name: str = "matrix",
    strict: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Align matrix sample columns with metadata rows.

    Args:
        matrix: DataFrame with sample columns.
        metadata: DataFrame indexed by sample_id.
        matrix_name: Name used in messages.
        strict: If True, raise on any mismatch; if False, only warn.

    Returns:
        Tuple ``(matrix_aligned, metadata_aligned)`` restricted to the shared
        samples, ordered by metadata.

    Raises:
        RiceDataError: In strict mode, when samples do not overlap at all or
            when only some samples overlap.
    """
    matrix_samples = list(matrix.columns)
    meta_samples = list(metadata.index)
    shared = [s for s in meta_samples if s in matrix_samples]
    only_matrix = [s for s in matrix_samples if s not in meta_samples]
    only_meta = [s for s in meta_samples if s not in matrix_samples]

    if not shared:
        raise RiceDataError(
            f"No overlapping samples between {matrix_name} and metadata.\n"
            f"matrix samples: {matrix_samples[:5]} ...\n"
            f"metadata samples: {meta_samples[:5]} ...\n"
            "Hint: sample IDs must match exactly (case-sensitive, no extra spaces)."
        )
    if strict and (only_matrix or only_meta):
        raise RiceDataError(
            f"Sample mismatch between {matrix_name} and metadata.\n"
            f"only in matrix: {only_matrix}\nonly in metadata: {only_meta}\n"
            "Hint: align_samples(strict=False) keeps the intersection with a warning."
        )
    if only_matrix or only_meta:
        logging.getLogger(LOGGER_NAME).warning(
            "Keeping %d shared samples; dropped %d matrix-only and %d metadata-only samples.",
            len(shared), len(only_matrix), len(only_meta),
        )
    return matrix[shared], metadata.loc[shared]


def check_group_sizes(metadata: pd.DataFrame, group_col: str = "group", min_size: int = 2) -> None:
    """Warn about groups with fewer than ``min_size`` replicates."""
    sizes = metadata[group_col].value_counts()
    small = sizes[sizes < min_size]
    if len(small):
        logging.getLogger(LOGGER_NAME).warning(
            "Group(s) with fewer than %d replicates: %s. "
            "Statistical tests need biological replicates; results may be unreliable.",
            min_size, small.to_dict(),
        )


# ---------------------------------------------------------------------------
# Export and provenance
# ---------------------------------------------------------------------------

def save_table(
    df: pd.DataFrame,
    path: str | os.PathLike,
    index: bool = False,
    **kwargs: Any,
) -> Path:
    """Save a DataFrame as CSV, creating parent directories as needed."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=index, **kwargs)
    return p


def _package_versions(packages: Iterable[str]) -> dict[str, str]:
    """Best-effort version lookup for provenance."""
    from importlib.metadata import PackageNotFoundError, version

    versions: dict[str, str] = {}
    for pkg in packages:
        try:
            versions[pkg] = version(pkg)
        except PackageNotFoundError:  # noqa: PERF203 - simple loop
            versions[pkg] = "NOT INSTALLED"
    return versions


def write_run_info(
    outdir: str | os.PathLike,
    script: str,
    params: dict[str, Any] | None = None,
    packages: Sequence[str] = (),
    extra: Sequence[str] = (),
) -> Path:
    """Write run_info.txt documenting parameters, versions and the seed.

    Args:
        outdir: Output directory (created if missing).
        script: Script file name for provenance.
        params: Parameters used for this run.
        packages: Key package names to version-check.
        extra: Extra note lines.

    Returns:
        Path to run_info.txt.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    lines = [
        f"script: {script}",
        f"run_time: {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}",
        f"python_version: {platform.python_version()}",
        f"platform: {platform.platform()}",
        f"seed: {SEED}",
        "parameters:",
    ]
    for key, value in (params or {}).items():
        lines.append(f"  {key} = {value}")
    versions = _package_versions(packages)
    if versions:
        lines.append("packages:")
        lines.extend(f"  {k} = {v}" for k, v in versions.items())
    if extra:
        lines.append("notes:")
        lines.extend(f"  {line}" for line in extra)
    target = out / "run_info.txt"
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target
