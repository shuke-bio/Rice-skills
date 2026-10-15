"""matplotlib_templates.py -- Rice-skills publication figure templates (Python backend).

Title:    Rice-skills matplotlib figure template library
Description:
    Copy-and-adapt starting points for every publication figure of the
    Rice-skills pipeline. Each template_*() function returns a
    ``(fig, ax)`` tuple; the caller saves it with ``save_nature_figure()``
    from the shared theme (PDF + SVG + TIFF at 300 dpi). All colours come
    from ``config/nature_palette.json`` through the theme API
    (``setup_nature_style()``, ``nature_palette()``, ``semantic()``,
    ``diverging()``, ``sequential()``); no colour is hard-coded, so the
    Python (matplotlib) and R (ggplot2) backends stay visually identical.

    The finished scripts under ``scripts/`` are generated from these
    templates: copy the template function into the analysis script, wire it
    to the pipeline CSV outputs, and keep the theme calls untouched.

Usage:
    Option 1 (recommended): copy this file next to the shared theme, e.g.
        cp templates/figure_templates/matplotlib_templates.py scripts/utils/
        cd scripts/utils && python matplotlib_templates.py

    Option 2: run it in place from the repository root -- this file inserts
    the repository's ``scripts/utils`` directory into ``sys.path`` so that
    ``from plot_theme import ...`` resolves:
        python templates/figure_templates/matplotlib_templates.py

    Option 3: import it from your own script after extending sys.path:
        import sys
        sys.path.insert(0, "scripts/utils")
        from matplotlib_templates import template_volcano

Dependencies: numpy, pandas, matplotlib, seaborn (for template_heatmap).
"""

from __future__ import annotations

import os
import sys

# ---------------------------------------------------------------------------
# Locate the shared theme module (plot_theme.py lives in scripts/utils/, while
# this template lives in templates/figure_templates/). When run in place the
# repository's scripts/utils directory is added to sys.path; when this file
# has been copied next to plot_theme.py its own directory is used.
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
for _candidate in (
    os.path.normpath(os.path.join(_HERE, "..", "..", "scripts", "utils")),
    _HERE,
):
    if os.path.isfile(os.path.join(_candidate, "plot_theme.py")):
        if _candidate not in sys.path:
            sys.path.insert(0, _candidate)
        break

# When the current working directory is outside the repository, point the
# palette loader at the repository config so figures stay reproducible.
_CFG = os.path.normpath(os.path.join(_HERE, "..", "..", "config", "nature_palette.json"))
if os.path.isfile(_CFG) and not os.environ.get("RICE_PALETTE"):
    os.environ["RICE_PALETTE"] = _CFG

from plot_theme import (  # noqa: E402
    diverging,
    group_colors,
    mm_to_inch,
    nature_palette,
    save_nature_figure,
    semantic,
    sequential,
    setup_nature_style,
)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Ellipse, Patch  # noqa: E402

# Fixed figure geometry: single column (89 mm) with sensible heights.
_W_MM, _H_MM = 89.0, 80.0
# chi-square 0.95 quantile with 2 degrees of freedom (scipy.stats.chi2.ppf(0.95, 2));
# hard-coded so the templates do not need scipy.
_CHI2_95_DF2 = 5.991


# ---------------------------------------------------------------------------
# Internal helpers (not part of the public template API)
# ---------------------------------------------------------------------------

def _normalise_sig(values) -> np.ndarray:
    """Normalise a ``significant`` column to the labels Up / DOWN / NS.

    Accepts the common spellings produced by DESeq2 wrappers and by the
    pipeline scripts: "Up", "UP", "up", "DOWN", "Down", "down", "NS", "ns",
    "NOT_SIGNIFICANT", ... Anything that is not up/down is treated as NS.
    """
    out = []
    for value in pd.Series(values).astype(str):
        token = value.strip().upper()
        if token.startswith("UP"):
            out.append("Up")
        elif token.startswith("DOWN"):
            out.append("DOWN")
        else:
            out.append("NS")
    return np.asarray(out, dtype=object)


def _sig_colors() -> dict:
    """Semantic colours keyed by significance label."""
    return {
        "Up": semantic("up"),
        "DOWN": semantic("down"),
        "NS": semantic("not_significant"),
    }


def _ink() -> str:
    """Theme ink colour for text labels (taken from the active rcParams)."""
    return plt.rcParams.get("axes.labelcolor", "black")


def _seq_cmap():
    """Continuous colormap built from the sequential palette."""
    return LinearSegmentedColormap.from_list("nature_sequential", sequential())


def _div_cmap():
    """Continuous colormap built from the diverging palette."""
    return LinearSegmentedColormap.from_list("nature_diverging", diverging())


def _confidence_ellipse(ax, x, y, color, level: float = 0.95) -> None:
    """Draw a normal-theory confidence ellipse on ``ax``.

    The ellipse is derived from the sample covariance and the chi-square
    quantile for ``level``; groups with fewer than three points (or a
    singular covariance) are skipped silently.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    if vals[0] <= 0 or not np.all(np.isfinite(vals)):
        return
    scale = np.sqrt(_CHI2_95_DF2) if abs(level - 0.95) < 1e-9 else None
    if scale is None:  # generic level: chi-square quantile via Wilson-Hilferty
        z = {0.90: 2.706, 0.95: 3.841, 0.99: 6.635}[round(level, 2)]
        scale = np.sqrt(1.0 - 2.0 / (9.0 * 2.0) + z * np.sqrt(2.0 / (9.0 * 2.0))) * 2.0
    theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width = 2.0 * scale * np.sqrt(vals[0])
    height = 2.0 * scale * np.sqrt(vals[1])
    ax.add_patch(
        Ellipse(
            (float(np.mean(x)), float(np.mean(y))),
            width,
            height,
            angle=theta,
            fill=False,
            edgecolor=color,
            linewidth=0.5,
            linestyle="--",
            alpha=0.8,
        )
    )


def _annotation_color_maps(annotation_df: pd.DataFrame) -> dict:
    """Map each annotation level to a group colour (groups5 first, then series12)."""
    palette = group_colors(12)
    maps = {}
    for col in annotation_df.columns:
        levels = list(dict.fromkeys(annotation_df[col].astype(str)))
        maps[col] = {lev: palette[i % len(palette)] for i, lev in enumerate(levels)}
    return maps


def _strip_rgb(annotation_df: pd.DataFrame, maps: dict) -> np.ndarray:
    """RGB matrix (n_annotation_columns x n_samples x 3) for a colour strip."""
    rgb = np.zeros((annotation_df.shape[1], annotation_df.shape[0], 3))
    for r, col in enumerate(annotation_df.columns):
        for c, value in enumerate(annotation_df[col].astype(str)):
            rgb[r, c] = to_rgb(maps[col][value])
    return rgb


def _neglog10(p_values) -> np.ndarray:
    """-log10(p) with a guard against p = 0 (mapped to -log10(1e-300))."""
    return -np.log10(np.maximum(pd.to_numeric(pd.Series(p_values), errors="coerce").fillna(1.0).to_numpy(dtype=float), 1e-300))


# ---------------------------------------------------------------------------
# Public templates -- each returns a (fig, ax) tuple; save with
# save_nature_figure(fig, "results/fig1", width_mm=89, height_mm=80).
# ---------------------------------------------------------------------------

def template_volcano(
    df: pd.DataFrame,
    log2fc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    id_col: str = "gene_id",
    sig_col: str = "significant",
    log2fc_cutoff: float = 1.0,
    padj_cutoff: float = 0.05,
    label_top: int = 5,
    title: str = "Volcano plot",
    ax=None,
):
    """Volcano plot of differential expression.

    Points are coloured by the semantic significance colours (up = red,
    down = blue, not significant = grey). Dashed grey reference lines mark
    the log2FC and adjusted-p thresholds. The top ``label_top`` significant
    genes are ranked by |log2FC| (ties broken by padj) and labelled.

    Args:
        df: Data frame with columns ``log2FoldChange``, ``padj``, ``gene_id``
            and (optionally) ``significant`` ("Up"/"DOWN"/"NS"); when
            ``significant`` is absent it is derived from the cutoffs.
        log2fc_col: Name of the log2 fold-change column.
        padj_col: Name of the adjusted p-value column.
        id_col: Name of the gene-ID column (used for labels).
        sig_col: Name of the significance column.
        log2fc_cutoff: |log2FC| threshold for the vertical lines.
        padj_cutoff: Adjusted p-value threshold for the horizontal line.
        label_top: Number of genes to label (0 disables labelling).
        title: Plot title.
        ax: Optional existing matplotlib axis.

    Returns:
        Tuple ``(fig, ax)``.

    Example:
        >>> de = pd.read_csv("results/transcriptome/02_de/de_all_salt_vs_control.csv")
        >>> fig, ax = template_volcano(de)
        >>> save_nature_figure(fig, "results/transcriptome/07_figures/fig_volcano",
        ...                    width_mm=89, height_mm=80)
    """
    setup_nature_style()
    for col in (log2fc_col, padj_col, id_col):
        if col not in df.columns:
            raise KeyError(f"template_volcano: missing column '{col}'")
    d = df
    if sig_col in d.columns:
        sig = _normalise_sig(d[sig_col])
    else:
        l2 = pd.to_numeric(d[log2fc_col], errors="coerce").to_numpy(dtype=float)
        p = pd.to_numeric(d[padj_col], errors="coerce").to_numpy(dtype=float)
        sig = np.where(
            p < padj_cutoff,
            np.where(l2 >= log2fc_cutoff, "Up", np.where(l2 <= -log2fc_cutoff, "DOWN", "NS")),
            "NS",
        ).astype(object)
    x = pd.to_numeric(d[log2fc_col], errors="coerce").to_numpy(dtype=float)
    y = _neglog10(d[padj_col])

    if ax is None:
        fig, ax = plt.subplots(figsize=(mm_to_inch(_W_MM), mm_to_inch(_H_MM)))
    else:
        fig = ax.figure
    colors = _sig_colors()
    for name in ("Up", "DOWN", "NS"):
        mask = sig == name
        if mask.any():
            ax.scatter(x[mask], y[mask], s=9, c=colors[name], label=name,
                       alpha=0.85, linewidths=0, zorder=2)
    ax.axvline(log2fc_cutoff, ls="--", lw=0.5, color=semantic("reference_line"), zorder=1)
    ax.axvline(-log2fc_cutoff, ls="--", lw=0.5, color=semantic("reference_line"), zorder=1)
    ax.axhline(-np.log10(padj_cutoff), ls="--", lw=0.5, color=semantic("reference_line"), zorder=1)
    ax.set_title(title)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 (adjusted P)")
    ax.legend(frameon=False, fontsize=6, loc="upper right", handletextpad=0.2,
              borderaxespad=0.2, labelspacing=0.3)

    if label_top and label_top > 0:
        sig_idx = np.flatnonzero(sig != "NS")
        if sig_idx.size:
            sub = pd.DataFrame(
                {"id": d[id_col].to_numpy()[sig_idx].astype(str),
                 "x": x[sig_idx], "y": y[sig_idx],
                 "p": pd.to_numeric(d[padj_col], errors="coerce").to_numpy(dtype=float)[sig_idx]}
            )
            sub["absx"] = sub["x"].abs()
            sub = sub.sort_values(["absx", "p"], ascending=[False, True]).head(label_top)
            for _, row in sub.iterrows():
                ax.annotate(
                    row["id"],
                    (row["x"], row["y"]),
                    fontsize=5.5,
                    color=_ink(),
                    xytext=(2, 2),
                    textcoords="offset points",
                )
    return fig, ax


def template_ma(
    df: pd.DataFrame,
    basemean_col: str = "baseMean",
    log2fc_col: str = "log2FoldChange",
    id_col: str = "gene_id",
    sig_col: str = "significant",
    padj_col: str = "padj",
    log2fc_cutoff: float = 1.0,
    padj_cutoff: float = 0.05,
    title: str = "MA plot",
    ax=None,
):
    """MA plot (log2 fold change vs mean expression).

    The x-axis is log10(baseMean); genes with baseMean == 0 (removed by the
    QC filter in practice) are not plotted. Colours follow the semantic
    palette; a dashed grey line marks y = 0 and the |log2FC| cutoffs.

    Args:
        df: Data frame with columns ``baseMean``, ``log2FoldChange``,
            ``gene_id`` and (optionally) ``significant`` or ``padj``.
        basemean_col: Name of the base-mean column.
        log2fc_col: Name of the log2 fold-change column.
        id_col: Name of the gene-ID column.
        sig_col: Name of the significance column.
        padj_col: Name of the adjusted p-value column (used when deriving).
        log2fc_cutoff: |log2FC| threshold for the vertical lines.
        padj_cutoff: Adjusted p-value threshold (used when deriving labels).
        title: Plot title.
        ax: Optional existing matplotlib axis.

    Returns:
        Tuple ``(fig, ax)``.

    Example:
        >>> de = pd.read_csv("results/transcriptome/02_de/de_all_salt_vs_control.csv")
        >>> fig, ax = template_ma(de)
        >>> save_nature_figure(fig, "results/transcriptome/07_figures/fig_ma",
        ...                    width_mm=89, height_mm=80)
    """
    setup_nature_style()
    for col in (basemean_col, log2fc_col, id_col):
        if col not in df.columns:
            raise KeyError(f"template_ma: missing column '{col}'")
    d = df
    if sig_col in d.columns:
        sig = _normalise_sig(d[sig_col])
    elif padj_col in d.columns:
        l2 = pd.to_numeric(d[log2fc_col], errors="coerce").to_numpy(dtype=float)
        p = pd.to_numeric(d[padj_col], errors="coerce").to_numpy(dtype=float)
        sig = np.where(
            p < padj_cutoff,
            np.where(l2 >= log2fc_cutoff, "Up", np.where(l2 <= -log2fc_cutoff, "DOWN", "NS")),
            "NS",
        ).astype(object)
    else:
        raise KeyError(
            "template_ma: the input has neither a 'significant' column nor a "
            f"'{padj_col}' column to derive it from."
        )
    base = pd.to_numeric(d[basemean_col], errors="coerce").to_numpy(dtype=float)
    x = np.where(base > 0, np.log10(np.where(base > 0, base, 1.0)), np.nan)
    y = pd.to_numeric(d[log2fc_col], errors="coerce").to_numpy(dtype=float)

    if ax is None:
        fig, ax = plt.subplots(figsize=(mm_to_inch(_W_MM), mm_to_inch(_H_MM)))
    else:
        fig = ax.figure
    colors = _sig_colors()
    for name in ("Up", "DOWN", "NS"):
        mask = (sig == name) & np.isfinite(x)
        if mask.any():
            ax.scatter(x[mask], y[mask], s=9, c=colors[name], label=name,
                       alpha=0.85, linewidths=0, zorder=2)
    ax.axhline(0, ls="--", lw=0.5, color=semantic("reference_line"), zorder=1)
    ax.axvline(log2fc_cutoff, ls="--", lw=0.5, color=semantic("reference_line"), zorder=1)
    ax.axvline(-log2fc_cutoff, ls="--", lw=0.5, color=semantic("reference_line"), zorder=1)
    ax.set_title(title)
    ax.set_xlabel("log10 base mean")
    ax.set_ylabel("log2 fold change")
    ax.legend(frameon=False, fontsize=6, loc="upper right", handletextpad=0.2,
              borderaxespad=0.2, labelspacing=0.3)
    return fig, ax


def template_pca(
    scores_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    color_col: str,
    shape_col: str | None = None,
    ellipse: bool = True,
    title: str = "PCA score plot",
    ax=None,
):
    """PCA score plot with per-group 95% confidence ellipses.

    ``meta_df`` is indexed by ``sample_id`` (or carries a ``sample_id``
    column). Point colours use the categorical palette (whose first five
    entries are the group colours); ``shape_col`` optionally maps a second
    metadata column to marker shapes. The dashed ellipse is a normal-theory
    95% confidence ellipse computed from the sample covariance.

    Args:
        scores_df: Data frame with columns ``sample_id``, ``PC1``, ``PC2``.
        meta_df: Metadata data frame indexed by sample_id.
        color_col: Metadata column used for point colours.
        shape_col: Optional metadata column used for marker shapes.
        ellipse: Draw 95% confidence ellipses (default True).
        title: Plot title.
        ax: Optional existing matplotlib axis.

    Returns:
        Tuple ``(fig, ax)``.

    Example:
        >>> scores = pd.read_csv("results/transcriptome/01_qc/pca_scores.csv")
        >>> meta = pd.read_csv("metadata.csv").set_index("sample_id")
        >>> fig, ax = template_pca(scores, meta, color_col="group", shape_col="batch")
        >>> save_nature_figure(fig, "results/transcriptome/07_figures/fig_pca_qc",
        ...                    width_mm=89, height_mm=80)
    """
    setup_nature_style()
    for col in ("sample_id", "PC1", "PC2"):
        if col not in scores_df.columns:
            raise KeyError(f"template_pca: scores_df is missing column '{col}'")
    if color_col not in meta_df.columns:
        raise KeyError(
            f"template_pca: meta_df is missing column '{color_col}'. "
            f"Available: {list(meta_df.columns)}"
        )
    if shape_col is not None and shape_col not in meta_df.columns:
        raise KeyError(
            f"template_pca: meta_df is missing column '{shape_col}'. "
            f"Available: {list(meta_df.columns)}"
        )
    meta = meta_df.copy()
    if "sample_id" not in meta.columns:
        # meta_df is indexed by sample_id: move the index into a column
        # (handles a named, unnamed or differently named index level).
        meta = meta.reset_index()
        if meta.columns[0] != "sample_id":
            meta = meta.rename(columns={meta.columns[0]: "sample_id"})
    d = scores_df.merge(meta, on="sample_id", how="left")
    # Discrete colours/shapes need discrete values (read_csv gives batch as int).
    d[color_col] = d[color_col].astype(str)
    if shape_col is not None:
        d[shape_col] = d[shape_col].astype(str)

    if ax is None:
        fig, ax = plt.subplots(figsize=(mm_to_inch(_W_MM), mm_to_inch(_H_MM)))
    else:
        fig = ax.figure
    palette = nature_palette()
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    groups = list(dict.fromkeys(d[color_col]))
    shapes = list(dict.fromkeys(d[shape_col])) if shape_col is not None else [None]
    for gi, group in enumerate(groups):
        color = palette[gi % len(palette)]
        for si, shape in enumerate(shapes):
            mask = (d[color_col] == group)
            if shape is not None:
                mask &= d[shape_col] == shape
            sub = d[mask]
            if sub.empty:
                continue
            marker = "o" if shape is None else markers[si % len(markers)]
            label = group if shape is None else f"{group}, {shape}"
            ax.scatter(
                sub["PC1"], sub["PC2"], s=14, c=color, marker=marker,
                label=label, edgecolors="none", alpha=0.9, zorder=2,
            )
            if ellipse and len(sub) >= 3:
                _confidence_ellipse(ax, sub["PC1"].to_numpy(dtype=float),
                                    sub["PC2"].to_numpy(dtype=float), color)
    ax.set_title(title)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(frameon=False, fontsize=6, loc="best", handletextpad=0.2,
              borderaxespad=0.2, labelspacing=0.3)
    return fig, ax


def template_dotplot(
    enrich_df: pd.DataFrame,
    top_n: int | None = 20,
    title: str = "Enrichment dot plot",
    ax=None,
):
    """Enrichment dot plot (GO / KEGG / Reactome style).

    x = GeneRatio (or Count when GeneRatio is absent), y = Description,
    size = Count, colour = -log10(p.adjust) on the sequential ramp. When the
    table has a ``Category`` column the ``top_n`` most significant terms of
    each category are shown.

    Args:
        enrich_df: Enrichment result table (clusterProfiler-style columns:
            ``Description``, ``p.adjust``, plus ``GeneRatio``/``Count`` and
            optional ``Category``).
        top_n: Maximum number of terms per category (None keeps all).
        title: Plot title.
        ax: Optional existing matplotlib axis.

    Returns:
        Tuple ``(fig, ax)``.

    Example:
        >>> ego = pd.read_csv("results/transcriptome/03_enrichment/go_bp.csv")
        >>> fig, ax = template_dotplot(ego, top_n=15)
        >>> save_nature_figure(fig, "results/transcriptome/07_figures/fig_go_bp",
        ...                    width_mm=89, height_mm=120)
    """
    setup_nature_style()
    for col in ("Description", "p.adjust"):
        if col not in enrich_df.columns:
            raise KeyError(f"template_dotplot: missing column '{col}'")
    d = enrich_df.copy()
    xvar = "GeneRatio" if "GeneRatio" in d.columns else ("Count" if "Count" in d.columns else None)
    if xvar is None:
        raise KeyError("template_dotplot: needs a 'GeneRatio' or 'Count' column for the x-axis.")
    svar = "Count" if "Count" in d.columns else xvar
    d["_neglog10p"] = _neglog10(d["p.adjust"])
    if top_n and len(d) > top_n:
        if "Category" in d.columns:
            d = d.sort_values("_neglog10p", ascending=False).groupby("Category", sort=False).head(top_n)
        else:
            d = d.sort_values("_neglog10p", ascending=False).head(top_n)
    d = d.sort_values(xvar)  # smallest value first; the y-axis is inverted below

    if ax is None:
        fig, ax = plt.subplots(figsize=(mm_to_inch(_W_MM), mm_to_inch(max(_H_MM, 40 + 8 * len(d)))))
    else:
        fig = ax.figure

    v = pd.to_numeric(d[svar], errors="coerce").fillna(0.0).to_numpy(dtype=float)
    lo, hi = float(v.min()), float(v.max())
    denom = (hi - lo) if hi > lo else 1.0
    sizes = 8.0 + 34.0 * (v - lo) / denom
    sc = ax.scatter(
        pd.to_numeric(d[xvar], errors="coerce"),
        np.arange(len(d)),
        s=sizes,
        c=d["_neglog10p"],
        cmap=_seq_cmap(),
        edgecolors="none",
        alpha=0.9,
        zorder=2,
    )
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["Description"].astype(str), fontsize=6)
    ax.invert_yaxis()  # largest x value on top
    ax.set_title(title)
    ax.set_xlabel("Gene ratio" if xvar == "GeneRatio" else xvar)
    ax.set_ylabel(None)

    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("-log10 (adjusted P)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.25)

    # Size legend with three proxy markers (small / medium / large).
    ref_vals = np.linspace(lo, hi, 3) if hi > lo else np.array([lo])
    handles = [
        ax.scatter([], [], s=8.0 + 34.0 * (r - lo) / denom, c=sequential()[3], alpha=0.9)
        for r in ref_vals
    ]
    ax.legend(
        handles=handles,
        labels=[f"{r:g}" for r in ref_vals],
        title=svar,
        frameon=False,
        fontsize=6,
        title_fontsize=6,
        loc="lower right",
        handletextpad=0.2,
        borderaxespad=0.2,
        labelspacing=0.3,
    )
    return fig, ax


def template_bubble(
    enrich_df: pd.DataFrame,
    top_n: int | None = 20,
    title: str = "Enrichment bubble plot",
    ax=None,
):
    """Enrichment bubble plot (GSEA style).

    x = NES (or a fold-enrichment / ratio column), y = Description,
    size = gene-set size (``Count`` / ``set_size`` when available),
    colour = -log10(p.adjust) on the sequential ramp.

    Args:
        enrich_df: Enrichment result table; the x column is auto-detected
            from ``("NES", "fold_enrichment", "FoldEnrichment",
            "foldEnrichment", "GeneRatio", "Count")``, the size column from
            ``("set_size", "setSize", "set size", "Count", "GeneRatio")`` and
            the p-value column from ``("p.adjust", "padj", "pvalue",
            "p_value", "Pvalue", "P.adjust")``.
        top_n: Maximum number of terms per category (None keeps all).
        title: Plot title.
        ax: Optional existing matplotlib axis.

    Returns:
        Tuple ``(fig, ax)``.

    Example:
        >>> gsea = pd.read_csv("results/transcriptome/04_gsea/gsea_results.csv")
        >>> fig, ax = template_bubble(gsea, top_n=15)
        >>> save_nature_figure(fig, "results/transcriptome/07_figures/fig_gsea",
        ...                    width_mm=89, height_mm=120)
    """
    setup_nature_style()
    if "Description" not in enrich_df.columns:
        raise KeyError("template_bubble: missing column 'Description'")

    def _pick(candidates):
        for cand in candidates:
            if cand in enrich_df.columns:
                return cand
        return None

    pcol = _pick(("p.adjust", "padj", "pvalue", "p_value", "Pvalue", "P.adjust"))
    if pcol is None:
        raise KeyError("template_bubble: needs a p-value column (e.g. 'p.adjust').")
    xvar = _pick(("NES", "fold_enrichment", "FoldEnrichment", "foldEnrichment",
                  "GeneRatio", "Count"))
    if xvar is None:
        raise KeyError("template_bubble: needs an x-axis column (e.g. 'NES').")
    svar = _pick(("set_size", "setSize", "set size", "Count", "GeneRatio"))
    xlab = {
        "NES": "Normalised enrichment score",
        "fold_enrichment": "Fold enrichment",
        "FoldEnrichment": "Fold enrichment",
        "foldEnrichment": "Fold enrichment",
        "GeneRatio": "Gene ratio",
        "Count": "Count",
    }[xvar]

    d = enrich_df.copy()
    d["_neglog10p"] = _neglog10(d[pcol])
    if top_n and len(d) > top_n:
        if "Category" in d.columns:
            d = d.sort_values("_neglog10p", ascending=False).groupby("Category", sort=False).head(top_n)
        else:
            d = d.sort_values("_neglog10p", ascending=False).head(top_n)
    d = d.sort_values(xvar)

    if ax is None:
        fig, ax = plt.subplots(figsize=(mm_to_inch(_W_MM), mm_to_inch(max(_H_MM, 40 + 8 * len(d)))))
    else:
        fig = ax.figure

    sc = ax.scatter(
        pd.to_numeric(d[xvar], errors="coerce"),
        np.arange(len(d)),
        s=pd.to_numeric(d[svar], errors="coerce").fillna(0.0).to_numpy(dtype=float) * 0.35 + 8.0
        if svar is not None
        else 20.0,
        c=d["_neglog10p"],
        cmap=_seq_cmap(),
        edgecolors="none",
        alpha=0.9,
        zorder=2,
    )
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["Description"].astype(str), fontsize=6)
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(None)

    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("-log10 (adjusted P)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.25)

    if svar is not None:
        v = pd.to_numeric(d[svar], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        lo, hi = float(v.min()), float(v.max())
        ref_vals = np.linspace(lo, hi, 3) if hi > lo else np.array([lo])
        handles = [
            ax.scatter([], [], s=r * 0.35 + 8.0, c=sequential()[3], alpha=0.9) for r in ref_vals
        ]
        ax.legend(
            handles=handles,
            labels=[f"{r:g}" for r in ref_vals],
            title=svar,
            frameon=False,
            fontsize=6,
            title_fontsize=6,
            loc="lower right",
            handletextpad=0.2,
            borderaxespad=0.2,
            labelspacing=0.3,
        )
    return fig, ax


def template_heatmap(
    mat,
    annotation_df: pd.DataFrame | None = None,
    type: str = "diverging",
    figsize=None,
    **kwargs,
):
    """Heatmap template (seaborn heatmap with optional annotation strips).

    ``type = "diverging"`` uses the blue-white-red ramp centred at zero
    (correlation matrices); ``type = "sequential"`` uses the single-hue blue
    ramp (expression / abundance matrices). ``annotation_df`` (indexed by
    the matrix columns) draws colour strips above the heatmap; for a square,
    symmetric matrix the strips are drawn on the rows as well. Annotation
    colours come from the group palette (groups5, extended by series12).

    Note:
        ``sns.heatmap`` plots rows in the given order; sort ``mat`` before
        calling (e.g. by module or log2FC). For hierarchical clustering with
        coloured annotations you can switch to ``sns.clustermap(...,
        col_colors=..., row_colors=...)`` and save the returned grid figure.

    Args:
        mat: Numeric matrix (DataFrame or array-like) with row/column names.
        annotation_df: Optional annotation data frame indexed like the
            matrix columns (and rows, for symmetric matrices).
        type: "diverging" (default) or "sequential".
        figsize: Optional matplotlib figure size in inches.
        **kwargs: Extra arguments forwarded to ``seaborn.heatmap``.

    Returns:
        Tuple ``(fig, ax)`` of the figure and the heatmap axis.

    Raises:
        ValueError: If ``type`` is invalid or annotation samples are missing.

    Example:
        >>> mat = pd.read_csv("results/transcriptome/02_de/vst_matrix.csv", index_col=0)
        >>> meta = pd.read_csv("metadata.csv").set_index("sample_id")
        >>> fig, ax = template_heatmap(mat, annotation_df=meta[["group"]], type="sequential")
        >>> save_nature_figure(fig, "results/transcriptome/07_figures/fig_heatmap_de",
        ...                    width_mm=120, height_mm=180)
    """
    setup_nature_style()
    if type not in ("diverging", "sequential"):
        raise ValueError(f"template_heatmap: type must be 'diverging' or 'sequential', got '{type}'")
    data = mat.copy() if isinstance(mat, pd.DataFrame) else pd.DataFrame(np.asarray(mat))
    if data.shape[0] < 1 or data.shape[1] < 1:
        raise ValueError("template_heatmap: the matrix must have at least one row and column")

    if type == "diverging":
        vmax = float(np.nanmax(np.abs(data.to_numpy(dtype=float))))
        heat_kwargs = dict(cmap=_div_cmap(), center=0.0, vmin=-vmax, vmax=vmax)
    else:
        heat_kwargs = dict(cmap=_seq_cmap())

    col_ann = row_ann = None
    if annotation_df is not None:
        ann = annotation_df.copy()
        symmetric = (
            data.shape[0] == data.shape[1]
            and [str(i) for i in data.index] == [str(c) for c in data.columns]
        )
        missing_cols = [s for s in data.columns if s not in ann.index]
        if missing_cols:
            raise ValueError(
                "template_heatmap: annotation_df is missing samples: "
                f"{missing_cols[:5]} (matrix columns must be annotated)."
            )
        ann = ann.reindex(data.columns)
        col_ann = ann
        if symmetric:
            row_ann = ann

    n_col_ann = 0 if col_ann is None else col_ann.shape[1]
    n_row_ann = 0 if row_ann is None else row_ann.shape[1]

    if n_col_ann or n_row_ann:
        fig = plt.figure(figsize=figsize or (3.6, 3.4))
        gs = fig.add_gridspec(
            2,
            2,
            width_ratios=[max(0.10 * n_row_ann, 0.04), 1.0],
            height_ratios=[max(0.10 * n_col_ann, 0.04), 1.0],
            hspace=0.03,
            wspace=0.03,
        )
        ax = fig.add_subplot(gs[1, 1])
        sns.heatmap(data, ax=ax, **heat_kwargs, **kwargs)

        col_maps = _annotation_color_maps(col_ann)
        if n_col_ann:
            ax_col = fig.add_subplot(gs[0, 1], sharex=ax)
            rgb = _strip_rgb(col_ann, col_maps)
            ax_col.imshow(rgb, aspect="auto", interpolation="none")
            ax_col.set_yticks(range(n_col_ann))
            ax_col.set_yticklabels(col_ann.columns.astype(str), fontsize=6)
            ax_col.set_xticks([])
            ax_col.tick_params(length=0)
            for spine in ax_col.spines.values():
                spine.set_visible(False)
        if n_row_ann:
            ax_row = fig.add_subplot(gs[1, 0], sharey=ax)
            rgb_t = np.transpose(_strip_rgb(row_ann, col_maps), (1, 0, 2))
            ax_row.imshow(rgb_t, aspect="auto", interpolation="none")
            ax_row.set_xticks(range(n_row_ann))
            ax_row.set_xticklabels(row_ann.columns.astype(str), fontsize=6, rotation=90)
            ax_row.set_yticks([])
            ax_row.tick_params(length=0)
            for spine in ax_row.spines.values():
                spine.set_visible(False)
        # Combined legend for the annotation strips.
        handles, labels = [], []
        for col, mapping in col_maps.items():
            for lev, colr in mapping.items():
                handles.append(Patch(facecolor=colr, edgecolor="none"))
                labels.append(f"{col}: {lev}")
        ax.legend(
            handles=handles,
            labels=labels,
            frameon=False,
            fontsize=6,
            loc="upper left",
            bbox_to_anchor=(1.15, 1.0),
            handlelength=1.0,
            handletextpad=0.4,
            labelspacing=0.3,
        )
    else:
        fig, ax = plt.subplots(figsize=figsize or (3.5, 3.2))
        sns.heatmap(data, ax=ax, **heat_kwargs, **kwargs)
    return fig, ax


# ---------------------------------------------------------------------------
# Self-contained demo: synthetic data (numpy seed 42) rendered into a system
# temporary directory (never inside the repository).
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import tempfile

    outdir = tempfile.mkdtemp(prefix="rice_figure_templates_")
    rng = np.random.default_rng(42)

    # -- synthetic differential-expression table (volcano / MA) ---------------
    n_genes = 400
    gene_ids = [f"LOC_Os{1 + (i % 12):02d}g{10000 + i:05d}" for i in range(n_genes)]
    log2fc = rng.normal(0.0, 1.2, n_genes)
    padj = 10.0 ** rng.uniform(-8.0, 0.0, n_genes)
    significant = np.where(
        padj < 0.05,
        np.where(log2fc >= 1.0, "Up", np.where(log2fc <= -1.0, "DOWN", "NS")),
        "NS",
    ).astype(object)
    de_df = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "baseMean": np.abs(rng.normal(5000, 3000, n_genes)) + 10,
            "log2FoldChange": log2fc,
            "padj": padj,
            "significant": significant,
        }
    )

    # -- synthetic PCA (two groups, two batches) ------------------------------
    samples = [f"{g}_{r}" for g in ("CTRL", "SALT") for r in (1, 2, 3)]
    pc1 = np.concatenate([rng.normal(0.0, 0.6, 3), rng.normal(3.0, 0.6, 3)])
    pc2 = rng.normal(0.0, 0.6, 6)
    scores_df = pd.DataFrame({"sample_id": samples, "PC1": pc1, "PC2": pc2})
    meta_df = pd.DataFrame(
        {
            "sample_id": samples,
            "group": ["control"] * 3 + ["salt"] * 3,
            "batch": [1, 2, 1, 2, 1, 2],
        }
    ).set_index("sample_id")

    # -- synthetic enrichment tables ------------------------------------------
    categories = ["BP"] * 8 + ["MF"] * 8
    enrich_df = pd.DataFrame(
        {
            "Category": categories,
            "Description": [f"term_{i + 1}" for i in range(16)],
            "GeneRatio": rng.uniform(0.02, 0.30, 16),
            "Count": rng.integers(5, 60, 16),
            "p.adjust": np.sort(10.0 ** rng.uniform(-8, -1.3, 16)),
        }
    )
    gsea_df = pd.DataFrame(
        {
            "Description": [f"pathway_{i + 1}" for i in range(10)],
            "NES": rng.normal(0.0, 1.5, 10),
            "set_size": rng.integers(10, 200, 10),
            "p.adjust": np.sort(10.0 ** rng.uniform(-6, -1.3, 10)),
        }
    )

    # -- synthetic expression / correlation matrices ---------------------------
    expr_mat = pd.DataFrame(
        rng.normal(0.0, 1.0, (8, 6)),
        index=[f"LOC_Os02g{10000 + i:05d}" for i in range(8)],
        columns=samples,
    )
    corr_mat = pd.DataFrame(
        np.corrcoef(rng.normal(0.0, 1.0, (8, 6)), rowvar=False),
        index=samples,
        columns=samples,
    )

    jobs = [
        ("template_volcano", template_volcano(de_df)),
        ("template_ma", template_ma(de_df)),
        ("template_pca", template_pca(scores_df, meta_df, color_col="group", shape_col="batch")),
        ("template_dotplot", template_dotplot(enrich_df, top_n=5)),
        ("template_bubble", template_bubble(gsea_df, top_n=8)),
        ("template_heatmap_sequential",
         template_heatmap(expr_mat, annotation_df=meta_df[["group", "batch"]], type="sequential")),
        ("template_heatmap_diverging",
         template_heatmap(corr_mat, annotation_df=meta_df[["group"]], type="diverging")),
    ]
    for name, (fig, _ax) in jobs:
        saved = save_nature_figure(fig, os.path.join(outdir, name),
                                   width_mm=89 if "heatmap" not in name else 120,
                                   height_mm=80 if "heatmap" not in name else 120)
        plt.close(fig)
        print(f"{name}: {len(saved)} files -> {os.path.basename(saved[0])}")

    print(f"\nDemo figures written to: {outdir}")
    print("NOTE: synthetic data are for template testing only, not for biology.")
