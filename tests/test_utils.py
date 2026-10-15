"""Smoke test for the Rice-skills shared Python utilities.

Validates palette loading, the Nature-style matplotlib theme, figure export
(PDF/SVG/TIFF) and the data_io loaders/validators. Run from the repository
root:

    python tests/test_utils.py

Exit code 0 = all checks passed.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "utils"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import data_io  # noqa: E402
import plot_theme  # noqa: E402


def test_palette() -> None:
    pal = plot_theme.load_palette()
    assert pal["typography"]["base_size_pt"] == 8, "unexpected base font size"
    assert len(pal["series12"]) == 12, "expected 12 categorical colours"
    assert len(pal["groups5"]) == 5, "expected 5 group colours"
    assert plot_theme.semantic("up") == "#CE6A6B"
    assert plot_theme.semantic("down") == "#5B84AE"
    assert len(plot_theme.diverging()) == 7
    assert len(plot_theme.sequential()) == 6
    assert len(plot_theme.nature_palette(3)) == 3


def test_theme_and_export() -> None:
    plot_theme.setup_nature_style()
    fig, ax = plt.subplots()
    rng = np.random.default_rng(data_io.SEED)
    ax.scatter(rng.normal(size=30), rng.normal(size=30),
               c=plot_theme.group_colors(2) * 15)
    ax.set_title("Self-test scatter")
    with tempfile.TemporaryDirectory() as tmp:
        saved = plot_theme.save_nature_figure(
            fig, str(Path(tmp) / "fig_selftest"), width_mm=89, height_mm=70
        )
        assert len(saved) == 3, f"expected 3 files, got {saved}"
        for path in saved:
            assert Path(path).stat().st_size > 0, f"empty output: {path}"
        plt.close(fig)


def test_data_io() -> None:
    rng = np.random.default_rng(data_io.SEED)
    genes = [f"LOC_Os01g{10000 + i:05d}" for i in range(20)]
    samples = [f"s{i:02d}" for i in range(6)]
    counts = pd.DataFrame(
        rng.poisson(lam=100, size=(20, 6)).astype(int),
        index=genes, columns=samples,
    )
    counts.index.name = "gene_id"
    meta = pd.DataFrame(
        {"sample_id": samples,
         "group": ["control"] * 3 + ["treatment"] * 3,
         "variety": ["Nipponbare"] * 6},
    ).set_index("sample_id")

    with tempfile.TemporaryDirectory() as tmp:
        cpath = Path(tmp) / "counts.csv"
        mpath = Path(tmp) / "meta.csv"
        counts.to_csv(cpath)
        meta.to_csv(mpath)

        loaded = data_io.load_expression_matrix(cpath)
        loaded_meta = data_io.load_metadata(mpath)
        aligned_matrix, aligned_meta = data_io.align_samples(loaded, loaded_meta)
        assert aligned_matrix.shape == (20, 6)
        data_io.check_group_sizes(aligned_meta)

        # Friendly error on missing file
        try:
            data_io.load_expression_matrix(Path(tmp) / "missing.csv")
        except data_io.RiceDataError:
            pass
        else:  # pragma: no cover - failure path
            raise AssertionError("expected RiceDataError for missing file")

        # Friendly error on sample mismatch
        bad_meta = loaded_meta.copy()
        bad_meta.index = ["wrong_" + s for s in bad_meta.index]
        try:
            data_io.align_samples(loaded, bad_meta)
        except data_io.RiceDataError:
            pass
        else:  # pragma: no cover - failure path
            raise AssertionError("expected RiceDataError for sample mismatch")

        # Provenance file
        info = data_io.write_run_info(
            Path(tmp) / "out", "test_utils.py",
            params={"n_genes": 20}, packages=("pandas", "numpy"),
        )
        assert info.is_file()


def main() -> int:
    data_io.setup_logging(log_file=None)
    for test in (test_palette, test_theme_and_export, test_data_io):
        test()
        print(f"PASS: {test.__name__}")
    print("ALL PYTHON UTILITY TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
