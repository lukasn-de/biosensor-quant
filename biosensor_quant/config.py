"""All user-adjustable settings, their defaults, and the QUEEN / GEQO channel presets.

Order of precedence (later wins):  defaults  <  --preset  <  --config file  <  command-line options.
Defaults follow the method described in the thesis (Section 2.4).
"""

import tomllib

# name: (default, type, help)
PARAMS = {
    # --- files ---
    "numerator_suffix": (None, str, "filename ending of the ratio numerator channel, e.g. '405.tif'"),
    "denominator_suffix": (None, str, "filename ending of the ratio denominator channel, e.g. '488.tif'"),
    "numerator_name": ("num", str, "label for the numerator channel in output columns"),
    "denominator_name": ("den", str, "label for the denominator channel in output columns"),
    "mask_suffix": ("Trans (CSU)_cp_masks.png", str, "filename ending of the Cellpose masks"),
    "brightfield_suffix": ("Trans (CSU).tif", str, "filename ending of the brightfield images"),
    "filename_pattern": (
        r"^(?P<well>[^-]+)--W(?P<wellnum>\d+)--P(?P<pos>\d+)--Z\d+--T(?P<time>\d+)--",
        str, "regex; named groups become output columns",
    ),
    "target_size": (None, "size", "resize all images to H,W (e.g. 2048,2048); default: no resizing"),
    # --- preprocessing ---
    "blur": (True, bool, "Gaussian blur on both channels"),
    "sigma": (2.0, float, "sigma of the channel blur"),
    "background_correction": (True, bool, "subtract median signal outside all ROIs"),
    "clip_negative": (True, bool, "set negative values to 0 after background subtraction"),
    # --- focus filter ---
    "entropy_source": ("brightfield", ["brightfield", "fluorescence"], "image(s) used for the entropy focus filter"),
    "entropy_min": (None, float, "skip images with entropy below this"),
    "entropy_max": (None, float, "skip images with entropy above this"),
    # --- fluorescence segmentation (average of both channels -> Otsu) ---
    "otsu_sigma": (2.0, float, "blur of the average image before Otsu (0 = none)"),
    "tophat_radius": (0, int, "white top-hat radius before Otsu (0 = off; needs pyclesperanto-prototype)"),
    # --- per-cell quantification ---
    "erosion_px": (2, int, "ROI erosion in pixels (0 = none)"),
    "statistic": ("mean", ["mean", "median"], "per-cell statistic used for the ratio"),
    "ratio_region": ("roi", ["roi", "segmented"], "measure in the whole ROI or only its Otsu-segmented part"),
    "min_seg_fraction": (0.5, float, "keep cells with segmented_area / roi_area above this; "
                         "removes out-of-focus or dead cells (0 = keep all)"),
    # --- output / performance ---
    "output_name": ("summary_all_cells", str, "name of the output table (without extension)"),
    "output_format": ("xlsx", ["xlsx", "csv"], "format of the output table"),
    "max_workers": (4, int, "parallel worker processes"),
    "max_in_flight": (8, int, "images held in memory at once"),
}

DEFAULTS = {name: spec[0] for name, spec in PARAMS.items()}

# Ratios (subscript = excitation, superscript = emission, nm):
#   QUEEN = Mean(ex405/em561) / Mean(ex488/em561)
#   GEQO  = Mean(ex488/em561) / Mean(ex488/em640)   (SFP / LSSFP)
PRESETS = {
    "queen": {
        "numerator_suffix": "405.tif",
        "denominator_suffix": "488.tif",
        "numerator_name": "405",
        "denominator_name": "488",
    },
    "geqo": {
        "numerator_suffix": "SFP.tif",
        "denominator_suffix": "LSSFP.tif",
        "numerator_name": "SFP",
        "denominator_name": "LSSFP",
    },
}


def load_config_file(path):
    """Read settings from a TOML file (keys as in PARAMS)."""
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    unknown = set(cfg) - set(PARAMS)
    if unknown:
        raise ValueError(f"unknown settings in {path}: {', '.join(sorted(unknown))}")
    return cfg


def build_settings(preset=None, config_file=None, overrides=None):
    settings = dict(DEFAULTS)
    if preset:
        settings.update(PRESETS[preset])
    if config_file:
        settings.update(load_config_file(config_file))
    if overrides:
        settings.update({k: v for k, v in overrides.items() if v is not None})
    if not settings["numerator_suffix"] or not settings["denominator_suffix"]:
        raise ValueError("set the channels with --preset or numerator_suffix / denominator_suffix")
    return settings
