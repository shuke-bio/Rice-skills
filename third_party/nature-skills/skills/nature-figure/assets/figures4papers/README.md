# figures4papers assets — intentionally not bundled

The original `figures4papers` demo scripts and preview images are **not included**
in this copy. This is intentional:

1. **Size** — the full asset tree is about 30 MB, mostly preview images, and
   contributes nothing to running the RiceMultiOmics pipeline.
2. **Licensing** — `figures4papers` is a separate upstream project
   (<https://github.com/ChenLiu-1996/figures4papers>) whose asset licensing
   differs from the MIT license of nature-skills; we do not redistribute it.

## What is still available in this copy

The **distilled rules** that nature-figure derived from figures4papers are all
present and are what the skill actually enforces at drawing time:

- `../../references/api.md` — palette, helper signatures, export conventions
- `../../references/design-theory.md` — typography, colour, composition rationale
- `../../references/common-patterns.md` — reusable layout patterns
- `../../references/chart-types.md` — chart-family atlas
- `../../references/tutorials.md` — end-to-end scaffold examples
- `../../references/demos.md` — routing map (points here for raw demo scripts)

`demos.md` refers to raw `plot_*.py` demo scripts under this directory. If you
need those concrete examples (method-comparison bars, radar/polar charts,
trajectory illustrations, dense multi-panel GridSpec layouts), obtain them from
the upstream repository above. For rice multi-omics figures (volcano, MA, PCA,
OPLS-DA, VIP, ROC, heatmaps, networks, dotplots, sankey, circos, UpSet), use the
templates in `templates/figure_templates/` and the shared theme in
`scripts/utils/plot_theme.R` / `plot_theme.py` instead.
