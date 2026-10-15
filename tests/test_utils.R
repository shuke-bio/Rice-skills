###############################################################################
# test_utils.R -- Smoke test for the Rice-skills shared R utilities.
#
# Validates palette loading, the Nature-style ggplot2 theme, figure export
# (PDF/SVG/TIFF) and the run-provenance helper. Run from the repository root:
#
#   Rscript tests/test_utils.R
#
# Exit code 0 = all checks passed. Outputs go to a temporary directory.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
utils_file <- file.path(dirname(script_path), "..", "scripts", "utils", "plot_theme.R")

# Point the theme at the repository palette via the environment variable.
Sys.setenv(RICE_PALETTE = file.path(dirname(script_path), "..", "config",
                                    "nature_palette.json"))
source(utils_file)

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

# --- palette -----------------------------------------------------------------
pal <- rice_palette()
if (pal$typography$base_size_pt != 8) fail("unexpected base font size")
if (length(pal$series12) != 12) fail("expected 12 categorical colours")
if (length(pal$groups5) != 5) fail("expected 5 group colours")
if (nature_semantic("up") != "#CE6A6B") fail("wrong 'up' semantic colour")
if (nature_semantic("down") != "#5B84AE") fail("wrong 'down' semantic colour")
if (length(diverging_palette()) != 7) fail("expected 7 diverging colours")
if (length(sequential_palette()) != 6) fail("expected 6 sequential colours")
if (abs(mm_to_inch(89) - 89 / 25.4) > 1e-9) fail("mm_to_inch conversion wrong")
cat("PASS: palette checks\n")

# --- theme and figure export -------------------------------------------------
set.seed(RICE_SEED)
df <- data.frame(
  x = rnorm(60), y = rnorm(60),
  g = rep(c("control", "treatment"), each = 30)
)
p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = g)) +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::labs(title = "Self-test scatter", x = "PC1", y = "PC2") +
  scale_colour_groups() +
  theme_nature()

outdir <- tempfile("rice_selftest")
stem <- file.path(outdir, "fig_selftest")
saved <- save_nature_figure(p, stem, width_mm = 89, height_mm = 70,
                            formats = c("pdf", "svg", "tiff"))
if (length(saved) != 3) fail("expected 3 saved files")
for (f in saved) {
  if (!file.exists(f) || file.info(f)$size <= 0) fail(paste("empty or missing output:", f))
}
cat("PASS: theme and figure export (", paste(basename(saved), collapse = ", "), ")\n", sep = "")

# --- provenance --------------------------------------------------------------
info <- write_run_info_r(outdir, "test_utils.R",
                         params = list(n_points = 60, seed = RICE_SEED),
                         packages = c("ggplot2", "jsonlite"))
if (!file.exists(info)) fail("run_info.txt was not written")
cat("PASS: run_info provenance\n")

cat("ALL R UTILITY TESTS PASSED\n")
