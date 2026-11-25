# biosensor-quant

Single-cell quantification of ratiometric fluorescent biosensors in budding yeast microscopy
images, starting from Cellpose segmentations. This repository holds the image analysis
scripts from my Master's thesis, restructured into a reusable, configurable tool:

> Lukas Nutz (2026). *Quantification of ATP Dynamics in S. cerevisiae and Preparation of a
> Genetic Screen Using Fluorescent Biosensors.* Master's thesis, RPTU Kaiserslautern-Landau.

The code was developed during my Master's thesis in the Quantitative Cell Biology group
(Jun.-Prof. Yury Bykov) at RPTU Kaiserslautern-Landau.

It has presets for the ATP biosensors **QUEEN** and **GEQO**. Every processing parameter can be
changed, so the tool can be adapted to other ratiometric sensors, microscopes and acquisition
settings. The defaults follow the method described in the thesis (Section 2.4).

## Pipeline

| Step | What happens | Settings |
|---|---|---|
| Segmentation | Cellpose (v4, `cpsam` model) on brightfield (Trans CSU) images, run separately | – |
| Preprocessing | Gaussian blur of both channels, then flat background correction (median signal outside all ROIs subtracted) | `blur`, `sigma`, `background_correction`, `clip_negative` |
| Focus filter | Shannon entropy of the brightfield (or fluorescence) images; images beyond a manually chosen threshold are skipped | `entropy_source`, `entropy_min`, `entropy_max` |
| Fluorescence segmentation | Average of both channels, blurred (optionally white top-hat), then Otsu threshold, giving the segmented area per cell | `otsu_sigma`, `tophat_radius` |
| Cell QC | Cells whose ROI is more than 2× their segmented area (out of focus or dead) are omitted | `min_seg_fraction` |
| Readout | ROI eroded; mean intensity per channel; ratio = numerator / denominator | `erosion_px`, `statistic`, `ratio_region` |

Ratios (subscript = excitation, superscript = emission, nm):

- **QUEEN** = Mean(ex405/em561) / Mean(ex488/em561)
- **GEQO** = Mean(ex488/em561) / Mean(ex488/em640), i.e. SFP / LSSFP

## Repository layout

```
biosensor_quant/
  pipeline.py          image loading, preprocessing and per-cell quantification
  config.py            all settings with defaults, QUEEN/GEQO presets
  cli.py               `biosensor-quant` command
  demo.py              synthetic demo dataset
config/example.toml    every setting with its default, as a template
tests/                 tests on synthetic images (python -m pytest)
```

## Installation

```bash
conda env create -f environment.yml
conda activate biosensor-quant
pip install -e .
```

The white top-hat filter (off by default) runs on the GPU and needs
`pip install pyclesperanto-prototype`. Cellpose is only needed for segmentation; installing it in a
separate environment (`pip install cellpose`) keeps this one light.

## Try it without real data

`biosensor-quant demo` writes a synthetic dataset: 3 replicates, each with SD and SGly wells
(4 fields of view per well). SGly ratios are set about 11% higher than SD.

```bash
biosensor-quant demo                                      # writes demo/rep1 ... demo/rep3
for r in rep1 rep2 rep3; do
  biosensor-quant run demo/$r --preset queen -o demo/$r/results
done
```

The median ratios should come out at about 1.3 (SD) and 1.45 (SGly).

## Usage

### 1. Image naming

All images of one field of view share a common prefix and differ only in the ending:

```
SD--W00001--P00001--Z00000--T00000--Trans (CSU).tif              brightfield
SD--W00001--P00001--Z00000--T00000--Trans (CSU)_cp_masks.png     Cellpose mask
SD--W00001--P00001--Z00000--T00000--405.tif                      channel 1
SD--W00001--P00001--Z00000--T00000--488.tif                      channel 2
```

This is the Olympus cellSens export layout; `well`, `wellnum`, `pos` and `time` are parsed
from the prefix. Other naming schemes work by changing `filename_pattern` (a regex with named
groups) and the `*_suffix` settings.

### 2. Segment cells with Cellpose

```bash
python -m cellpose --dir /abs/path/to/data --img_filter "Trans (CSU)" \
    --pretrained_model cpsam --save_png --no_npy --use_gpu
```

### 3. Choose a focus threshold (optional)

```bash
biosensor-quant entropy data/ -o entropy.csv
```

Look at the images with the most extreme entropy values, then set `entropy_max` and/or
`entropy_min`.

### 4. Quantify

```bash
biosensor-quant run data/ --preset queen -o results/
biosensor-quant run data/ --preset geqo -o results/ --entropy-max 6.5

# any setting can be changed on the command line ...
biosensor-quant run data/ --preset queen --sigma 3 --statistic median --erosion-px 0

# ... or collected in a settings file (copy and edit config/example.toml)
biosensor-quant run data/ --config my_settings.toml
```

`biosensor-quant run -h` lists every setting with its default. Precedence: defaults < `--preset`
< `--config` < command-line options.

Outputs in the results folder:

- `summary_all_cells.xlsx` (or `.csv`): one row per cell, with the filename metadata,
  `cell_id`, `roi_area` (eroded ROI), `segmented_area`, `ratio`, `mean_*`, `median_*`,
  `circularity`, `eccentricity`, `solidity` and `extent`
- `settings.json`: the exact settings of the run (can be converted into a `--config` file)
- `skipped_entropy.txt`: images removed by the focus filter, if any

## Testing

```bash
python -m pytest tests/
```

## License

MIT, see [LICENSE](LICENSE). If you use this code, please cite the thesis (see
[CITATION.cff](CITATION.cff)).
