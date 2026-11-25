"""End-to-end tests on synthetic images with a known ratio.

Run with:  python -m pytest tests/
"""

import numpy as np
import pandas as pd

from biosensor_quant import build_settings, run_acquisition
from biosensor_quant.demo import make_demo_dataset, make_fov
from biosensor_quant.pipeline import erode_cell

BASE = "SD--W00001--P00001--Z00000--T00000--"


def settings(**overrides):
    return build_settings("queen", overrides={"max_workers": 1, "output_format": "csv", **overrides})


def test_known_ratio(tmp_path):
    make_fov(tmp_path, BASE, ratio=1.8)
    cells = run_acquisition(tmp_path, tmp_path / "out", settings())
    assert len(cells) == 12
    assert np.allclose(cells["ratio"], 1.8, atol=0.02)
    assert cells["well"].unique().tolist() == ["SD"]
    assert (tmp_path / "out" / "summary_all_cells.csv").exists()
    assert (tmp_path / "out" / "settings.json").exists()


def test_median_and_no_erosion_options(tmp_path):
    make_fov(tmp_path, BASE, ratio=1.5)
    eroded = run_acquisition(tmp_path, tmp_path / "a", settings())
    plain = run_acquisition(tmp_path, tmp_path / "b", settings(statistic="median", erosion_px=0))
    assert np.allclose(plain["ratio"], 1.5, atol=0.02)
    assert (plain["roi_area"].values > eroded["roi_area"].values).all()


def test_entropy_filter_skips_image(tmp_path):
    make_fov(tmp_path, BASE, ratio=1.5)
    cells = run_acquisition(tmp_path, tmp_path / "out", settings(entropy_max=0.1))
    assert cells is None
    assert (tmp_path / "out" / "skipped_entropy.txt").read_text().startswith(BASE)


def test_erode_cell():
    mask = np.zeros((20, 20), bool)
    mask[5:15, 5:15] = True
    eroded = erode_cell(mask, (5, 5, 15, 15), 2)
    assert eroded.sum() == 6 * 6


def test_demo_dataset(tmp_path):
    folders = make_demo_dataset(tmp_path, replicates=1, fovs_per_well=2)
    run_acquisition(folders[0], tmp_path / "out", settings())
    cells = pd.read_csv(tmp_path / "out" / "summary_all_cells.csv")
    med = cells.groupby("well")["ratio"].median()
    assert med["SGly"] > med["SD"]


def test_min_seg_fraction_filter(tmp_path):
    make_fov(tmp_path, BASE, ratio=1.5)
    kept = run_acquisition(tmp_path, tmp_path / "a", settings())
    assert (kept["segmented_area"] / kept["roi_area"] > 0.5).all()
    assert run_acquisition(tmp_path, tmp_path / "b", settings(min_seg_fraction=1.0)) is None
