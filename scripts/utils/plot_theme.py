"""plot_theme.py -- Rice-skills shared Nature-style plotting theme (Python backend).

Title:    Rice-skills Nature-style matplotlib helpers
Description:
    Provides a low-saturation, colorblind-friendly palette and matplotlib
    rcParams that follow Nature-series illustration conventions (Arial 8 pt
    base, bold 9 pt titles, 0.5 pt lines, PDF/SVG/TIFF >=300 dpi output).
    The palette is loaded from config/nature_palette.json so that the Python
    (matplotlib) and R (ggplot2/ComplexHeatmap) backends stay visually
    identical. All code is implemented from scratch; no third-party theme
    source is copied.

Usage:
    from plot_theme import setup_nature_style, save_nature_figure
    setup_nature_style()
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ...
    save_nature_figure(fig, "results/fig1_pca", width_mm=89, height_mm=70)

Notes:
    - Requires matplotlib >= 3.7 and pillow (for TIFF output).
    - Arial is used when available; DejaVu Sans is the fallback.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")  # headless-safe default; scripts may override before import

import matplotlib.pyplot as plt  # noqa: E402

PALETTE_FILENAME = "nature_palette.json"
DEFAULT_FORMATS = ("pdf", "svg", "tiff")
_PALETTE_CACHE: dict | None = None


def find_palette(path: str | os.PathLike | None = None) -> Path:
    """Locate the shared palette JSON file.

    Search order: explicit path -> RICE_PALETTE env var -> config/ subfolder of
    the current working directory and its parents (up to 4 levels).

    Args:
        path: Optional explicit path to nature_palette.json.

    Returns:
        Absolute path to the palette file.

    Raises:
        FileNotFoundError: If the palette file cannot be located.
    """
    if path is not None:
        p = Path(path)
        if not p.is_file():
            raise FileNotFoundError(
                f"Palette file not found: {path}\n"
                "Hint: pass a valid path or leave None to auto-locate config/nature_palette.json."
            )
        return p.resolve()

    env = os.environ.get("RICE_PALETTE", "")
    if env and Path(env).is_file():
        return Path(env).resolve()

    current = Path.cwd()
    for _ in range(5):
        candidate = current / "config" / PALETTE_FILENAME
        if candidate.is_file():
            return candidate.resolve()
        if current.parent == current:
            break
        current = current.parent

    raise FileNotFoundError(
        "Could not locate config/nature_palette.json.\n"
        "Hint: run scripts from the Rice-skills repository root, or set the "
        "RICE_PALETTE environment variable to the absolute JSON path."
    )


def load_palette(path: str | os.PathLike | None = None) -> dict:
    """Load (and cache) the shared palette dictionary."""
    global _PALETTE_CACHE
    if path is None and _PALETTE_CACHE is not None:
        return _PALETTE_CACHE
    with open(find_palette(path), "r", encoding="utf-8") as handle:
        palette = json.load(handle)
    if path is None:
        _PALETTE_CACHE = palette
    return palette


def mm_to_inch(mm: float) -> float:
    """Convert millimetres to inches."""
    return mm / 25.4


def setup_nature_style() -> dict:
    """Apply Nature-style rcParams to matplotlib.

    Returns:
        The palette dictionary (for convenience).
    """
    palette = load_palette()
    typo = palette["typography"]
    sans = [typo["font_family"], typo.get("fallback_family", "Helvetica"), "DejaVu Sans"]
    rc = {
        "font.family": "sans-serif",
        "font.sans-serif": sans,
        "font.size": typo["base_size_pt"],
        "axes.titlesize": typo["title_size_pt"],
        "axes.titleweight": "bold",
        "axes.labelsize": typo["base_size_pt"],
        "axes.linewidth": typo["axis_line_width_pt"],
        "axes.edgecolor": "#4D4A45",
        "axes.labelcolor": "#33312E",
        "axes.grid": False,
        "xtick.labelsize": typo["tick_size_pt"],
        "ytick.labelsize": typo["tick_size_pt"],
        "xtick.color": "#33312E",
        "ytick.color": "#33312E",
        "xtick.major.width": typo["line_width_pt"],
        "ytick.major.width": typo["line_width_pt"],
        "xtick.major.size": 1.5,
        "ytick.major.size": 1.5,
        "legend.fontsize": typo["legend_size_pt"],
        "legend.frameon": False,
        "legend.title_fontsize": typo["legend_size_pt"],
        "lines.linewidth": 0.8,
        "lines.markersize": 3,
        "grid.linewidth": 0.25,
        "grid.color": "#E8E8E8",
        "figure.dpi": 300,
        "savefig.dpi": palette["figure"]["dpi_raster"],
        "savefig.bbox": "tight",
        "savefig.transparent": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,  # embed TrueType (editable text in PDF)
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # keep text as text in SVG
    }
    plt.rcParams.update(rc)
    return palette


def nature_palette(n: int | None = None) -> list[str]:
    """Return the first ``n`` categorical colours (all 12 when n is None)."""
    series = load_palette()["series12"]
    return series if n is None else series[: max(1, n)]


def group_colors(n: int | None = None) -> list[str]:
    """Return the group palette (5 colours) extended with series12 if needed."""
    palette = load_palette()
    groups = palette["groups5"]
    if n is None or n <= len(groups):
        return groups if n is None else groups[:n]
    extra = [c for c in palette["series12"] if c not in groups]
    return (groups + extra)[:n]


def semantic(name: str) -> str:
    """Return a semantic colour: up/down/not_significant/highlight/reference_line."""
    colors = load_palette()["semantic"]
    if name not in colors:
        raise KeyError(
            f"Unknown semantic colour '{name}'. Available: {', '.join(colors)}"
        )
    return colors[name]


def diverging() -> list[str]:
    """Diverging blue-white-red ramp for correlation heatmaps."""
    return load_palette()["diverging7"]


def sequential() -> list[str]:
    """Sequential single-hue blue ramp for intensity heatmaps."""
    return load_palette()["sequential6"]


def save_nature_figure(
    fig,
    file_stem: str | os.PathLike,
    width_mm: float | None = None,
    height_mm: float | None = None,
    formats: Sequence[str] | None = None,
    **kwargs,
) -> list[str]:
    """Save a matplotlib figure in Nature-compatible formats.

    Args:
        fig: Matplotlib figure.
        file_stem: Output path without extension.
        width_mm: Figure width in mm (default 89 mm = single column).
        height_mm: Figure height in mm (default 89 mm).
        formats: Output formats (default: pdf, svg, tiff).
        **kwargs: Extra arguments forwarded to ``Figure.savefig``.

    Returns:
        List of saved file paths.
    """
    palette = load_palette()
    figspec = palette["figure"]
    w_mm = figspec["single_column_mm"] if width_mm is None else width_mm
    h_mm = figspec["single_column_mm"] if height_mm is None else height_mm
    fmts = tuple(f.lower() for f in (formats or DEFAULT_FORMATS))
    stem = Path(file_stem)
    if stem.parent and not stem.parent.exists():
        stem.parent.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(mm_to_inch(w_mm), mm_to_inch(h_mm))
    saved: list[str] = []
    for fmt in fmts:
        out = stem.with_suffix(f".{fmt}")
        fig.savefig(out, format=fmt, dpi=figspec["dpi_raster"], **kwargs)
        saved.append(str(out))
    return saved
