"""Synthetic demo dataset to check the pipeline without real microscopy images.

Creates cellSens-style filenames (brightfield, Cellpose-like mask, two channels) for several
replicates and conditions, with known per-condition biosensor ratios.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from imageio.v3 import imwrite
from skimage.draw import disk

DEMO_RATIOS = {"SD": 1.30, "SGly": 1.45}


def make_fov(
    folder: Path,
    base: str,
    ratio: float,
    numerator_suffix: str = "405.tif",
    denominator_suffix: str = "488.tif",
    shape: tuple[int, int] = (256, 256),
    n_cells: int = 12,
    ratio_sd: float = 0.0,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Write one field of view with round non-overlapping 'cells'. Returns the label image."""
    rng = rng or np.random.default_rng(0)
    labels = np.zeros(shape, np.uint16)
    lab = 0
    for _ in range(n_cells * 20):  # rejection sampling of non-overlapping disks
        if lab == n_cells:
            break
        r = int(rng.integers(12, 19))
        cy, cx = rng.integers(r + 2, shape[0] - r - 2), rng.integers(r + 2, shape[1] - r - 2)
        rr, cc = disk((cy, cx), r + 2, shape=shape)
        if labels[rr, cc].any():
            continue
        lab += 1
        labels[disk((cy, cx), r, shape=shape)] = lab

    background, den_signal = 100.0, 200.0
    cell_ratio = np.zeros(lab + 1)
    cell_ratio[1:] = ratio + rng.normal(0, ratio_sd, lab)
    den = background + den_signal * (labels > 0)
    num = background + den_signal * cell_ratio[labels]

    def noisy(img, sd=3.0):
        return np.clip(img + rng.normal(0, sd, shape), 0, None).astype(np.uint16)

    folder.mkdir(parents=True, exist_ok=True)
    imwrite(folder / f"{base}{numerator_suffix}", noisy(num))
    imwrite(folder / f"{base}{denominator_suffix}", noisy(den))
    imwrite(folder / f"{base}Trans (CSU).tif", noisy(np.full(shape, 1000.0), 30))
    imwrite(folder / f"{base}Trans (CSU)_cp_masks.png", labels)
    return labels


def make_demo_dataset(
    out_dir: str | Path, replicates: int = 3, fovs_per_well: int = 4, seed: int = 1
) -> list[Path]:
    """Write ``replicates`` folders (rep1, rep2, ...) with SD and SGly wells.

    SGly ratios are ~11% higher than SD, and each replicate gets a small offset, so the
    paired t-test and bootstrap CIs in R have something realistic to work on.
    """
    rng = np.random.default_rng(seed)
    folders = []
    for r in range(1, replicates + 1):
        folder = Path(out_dir) / f"rep{r}"
        offset = rng.normal(0, 0.02)
        for w, (well, ratio) in enumerate(DEMO_RATIOS.items(), start=1):
            for p in range(1, fovs_per_well + 1):
                base = f"{well}--W{w:05d}--P{p:05d}--Z00000--T00000--"
                make_fov(folder, base, ratio + offset, ratio_sd=0.08, rng=rng)
        folders.append(folder)
    return folders
