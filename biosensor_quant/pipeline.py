# Author: Lukas Nutz
#
# Single-cell quantification of ratiometric biosensors from Cellpose masks.
# Restructured from amd_Queen_v1_2.py and amd_GEQO_v1.py (09/2025), with the options of
# amd_v1.py (blur, top-hat, ROI erosion, mean ratio). The QUEEN and GEQO scripts were identical
# apart from channel names, so they are merged here; channels are set via the config.
#
# Workflow:
#   1) load Cellpose mask + two fluorescence channels (+ brightfield for the focus filter)
#   2) Gaussian blur, flat background correction (median outside all ROIs)
#   3) entropy-based focus filter
#   4) average image of both channels -> blur (+ optional white top-hat) -> Otsu threshold
#   5) per cell: erode ROI, QC filter, measure both channels, compute ratio, shape descriptors

import json
import os
import re
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
from imageio.v3 import imread
from skimage import filters
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops, shannon_entropy
from skimage.morphology import disk, erosion
from skimage.transform import resize


def run_acquisition(input_folder, output_folder, settings):
    """Quantify all images in input_folder and write one table with all cells."""
    p = settings
    os.makedirs(output_folder, exist_ok=True)
    masks = sorted(f for f in os.listdir(input_folder) if f.endswith(p["mask_suffix"]))
    results = []

    def loader():
        for mask in masks:
            images = load_and_preprocess(mask, input_folder, output_folder, p)
            if images:
                yield images

    print(f"Starting analysis of {len(masks)} files...")
    with ProcessPoolExecutor(p["max_workers"]) as executor:
        futures = set()
        images_iter = loader()
        for _ in range(p["max_in_flight"]):
            try:
                futures.add(executor.submit(analyze_cells, next(images_iter), p))
            except StopIteration:
                break
        completed = 0
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for fut in done:
                futures.remove(fut)
                try:
                    res = fut.result()
                    if res:
                        results.extend(res)
                except Exception as e:
                    print(f"⚠️ Task failed: {e}")
                completed += 1
                print(f"\rProcessed {completed}/{len(masks)}", end="")
                try:
                    futures.add(executor.submit(analyze_cells, next(images_iter), p))
                except StopIteration:
                    pass

    # keep a record of the settings used for this run
    with open(os.path.join(output_folder, "settings.json"), "w") as f:
        json.dump(p, f, indent=2)

    if results:
        cell_data = pd.DataFrame(results)
        output_file = os.path.join(output_folder, f"{p['output_name']}.{p['output_format']}")
        if p["output_format"] == "csv":
            cell_data.to_csv(output_file, index=False)
        else:
            cell_data.to_excel(output_file, index=False)
        print(f"\nSummary saved: {output_file}")
    else:
        print("\nNo results found.")
        return None

    return cell_data


def load_and_preprocess(mask_name, input_folder, output_folder, settings):
    """
    Steps:
    1) Load mask and channels (resize to target_size if set)
    2) Create inverse (background) mask from the labeled mask (label==0)
    3) Apply Gaussian blur to each channel
    4) Flat background correction: subtract median of background pixels from each blurred channel
    5) Entropy-based focus filter (brightfield or fluorescence channels)
    6) Build an average image from corrected channels, optionally white-tophat, Otsu threshold
    """
    p = settings
    m = re.search(p["filename_pattern"], mask_name, re.IGNORECASE)
    if not m:
        return None
    metadata = {k: int(v) if v is not None and v.isdigit() else v for k, v in m.groupdict().items()}

    loaded = load_images(mask_name, input_folder, p)
    if loaded is None:
        return None
    img_mask, img_num_raw, img_den_raw, img_bf = loaded

    # Inverse background mask from labeled mask
    bg_mask = img_mask == 0
    if not np.any(bg_mask):
        # Fallback: if mask covers the whole image for some reason, use everything as background
        bg_mask = np.ones_like(img_mask, dtype=bool)

    # Gaussian blur on both channels
    if p["blur"]:
        img_num_blur = filters.gaussian(img_num_raw, sigma=p["sigma"], preserve_range=True)
        img_den_blur = filters.gaussian(img_den_raw, sigma=p["sigma"], preserve_range=True)
    else:
        img_num_blur = img_num_raw
        img_den_blur = img_den_raw

    # Flat background correction from background median (inverse of mask)
    if p["background_correction"]:
        img_num = (img_num_blur - np.median(img_num_blur[bg_mask])).astype(np.float32)
        img_den = (img_den_blur - np.median(img_den_blur[bg_mask])).astype(np.float32)
    else:
        img_num = img_num_blur.astype(np.float32)
        img_den = img_den_blur.astype(np.float32)

    # keep non-negative after subtraction (flat-field like offset correction)
    if p["clip_negative"]:
        np.maximum(img_num, 0, out=img_num)
        np.maximum(img_den, 0, out=img_den)

    # Entropy check: skip out-of-focus images
    if p["entropy_min"] is not None or p["entropy_max"] is not None:
        if p["entropy_source"] == "brightfield":
            entropies = {"brightfield": shannon_entropy(img_bf)}
        else:
            entropies = {p["numerator_name"]: shannon_entropy(img_num),
                         p["denominator_name"]: shannon_entropy(img_den)}
        too_low = p["entropy_min"] is not None and min(entropies.values()) < p["entropy_min"]
        too_high = p["entropy_max"] is not None and max(entropies.values()) > p["entropy_max"]
        if too_low or too_high:
            values = ", ".join(f"entropy_{k}={v:.2f}" for k, v in entropies.items())
            with open(os.path.join(output_folder, "skipped_entropy.txt"), "a") as f:
                f.write(f"{mask_name}, {values}\n")
            return None

    # Average image and thresholding on corrected channels
    img_average = (img_num + img_den) / 2.0
    filtered_average = filters.gaussian(img_average, p["otsu_sigma"]) if p["otsu_sigma"] > 0 else img_average
    if p["tophat_radius"] != 0:
        filtered_average = white_tophat(filtered_average, p["tophat_radius"])

    threshold_average = threshold_otsu(filtered_average)

    # Return keys that analyze_cells expects
    return {
        "mask_name": mask_name,
        "metadata": metadata,
        "mask": img_mask,
        "img_num": img_num,  # corrected & blurred
        "img_den": img_den,  # corrected & blurred
        "filtered_average": filtered_average,
        "threshold_average": threshold_average,
    }


def load_images(mask_name, input_folder, settings):
    p = settings
    path_mask = os.path.join(input_folder, mask_name)
    if not os.path.exists(path_mask):
        return None

    # strip the mask suffix (e.g. "Trans (CSU)_cp_masks.png") and append the channel endings
    base = mask_name[: -len(p["mask_suffix"])]
    path_num = os.path.join(input_folder, f"{base}{p['numerator_suffix']}")
    path_den = os.path.join(input_folder, f"{base}{p['denominator_suffix']}")
    path_bf = os.path.join(input_folder, f"{base}{p['brightfield_suffix']}")

    if not (os.path.exists(path_num) and os.path.exists(path_den)):
        return None

    # channels -> float32
    img_num = imread(path_num).astype(np.float32)
    img_den = imread(path_den).astype(np.float32)
    img_bf = None
    if p["entropy_source"] == "brightfield" and \
            (p["entropy_min"] is not None or p["entropy_max"] is not None):
        if not os.path.exists(path_bf):
            print(f"\n⚠️ No brightfield image for {mask_name}, skipped")
            return None
        img_bf = imread(path_bf).astype(np.float32)

    target_size = tuple(p["target_size"]) if p["target_size"] else img_num.shape

    # mask -> label and resize with nearest (order=0)
    mask = imread(path_mask).astype(int)
    img_mask = label(mask)
    if img_mask.shape != target_size:
        img_mask = resize(img_mask, target_size, order=0, preserve_range=True, anti_aliasing=False).astype(int)

    # channels -> resize to target_size (bilinear)
    if img_num.shape != target_size:
        img_num = resize(img_num, target_size, order=1, preserve_range=True, anti_aliasing=True).astype(np.float32)
    if img_den.shape != target_size:
        img_den = resize(img_den, target_size, order=1, preserve_range=True, anti_aliasing=True).astype(np.float32)

    return img_mask, img_num, img_den, img_bf


def analyze_cells(images, settings):
    p = settings
    num, den = p["numerator_name"], p["denominator_name"]
    statistic = np.mean if p["statistic"] == "mean" else np.median
    results = []

    for cell in regionprops(images["mask"]):
        cell_id = cell.label
        cell_mask = np.zeros_like(images["mask"], dtype=bool)
        cell_mask[tuple(cell.coords.T)] = True
        if p["erosion_px"] > 0:
            cell_mask = erode_cell(cell_mask, cell.bbox, p["erosion_px"])
        if not cell_mask.any():
            continue  # cell vanished after erosion

        threshold_average = images["threshold_average"]
        segment_average = (images["filtered_average"] > threshold_average) & cell_mask

        # cell QC: omit cells whose ROI is much larger than their segmented signal
        # (out of focus or dead); 0.5 = ROI more than 2x the segmented area
        if segment_average.sum() / cell_mask.sum() <= p["min_seg_fraction"]:
            continue

        # measure in the whole (eroded) ROI or only in its Otsu-segmented part
        region = segment_average if p["ratio_region"] == "segmented" else cell_mask
        values_num = images["img_num"][region]
        values_den = images["img_den"][region]

        if values_num.size and statistic(values_den) > 0:
            values_ratio = statistic(values_num) / statistic(values_den)
        else:
            values_ratio = np.nan

        perim = cell.perimeter
        area = cell_mask.sum()
        circularity = 4 * np.pi * area / (perim**2 + 1e-6)

        results.append({
            **images["metadata"],
            "cell_id": cell_id,
            "roi_area": area,
            "segmented_area": segment_average.sum(),
            "ratio": values_ratio,
            f"mean_{num}": values_num.mean() if values_num.size else 0,
            f"mean_{den}": values_den.mean() if values_den.size else 0,
            f"median_{num}": np.median(values_num) if values_num.size else 0,
            f"median_{den}": np.median(values_den) if values_den.size else 0,
            "circularity": circularity,
            "eccentricity": cell.eccentricity,
            "solidity": cell.solidity,
            "extent": cell.extent,
        })

    return results


def erode_cell(cell_mask, bbox, pixels):
    """Erode one cell mask by `pixels`, working on its bounding box only (faster)."""
    pad = pixels + 1
    r0, c0 = max(bbox[0] - pad, 0), max(bbox[1] - pad, 0)
    r1, c1 = min(bbox[2] + pad, cell_mask.shape[0]), min(bbox[3] + pad, cell_mask.shape[1])
    eroded = np.zeros_like(cell_mask)
    eroded[r0:r1, c0:c1] = erosion(cell_mask[r0:r1, c0:c1], disk(pixels))
    return eroded


def white_tophat(img, radius):
    import pyclesperanto_prototype as cle  # GPU; only needed when tophat_radius != 0

    img_gpu = cle.push(img.astype(np.float32))
    eroded = cle.minimum_sphere(img_gpu, None, radius_x=radius, radius_y=radius)
    opened = cle.maximum_sphere(eroded, None, radius_x=radius, radius_y=radius)
    wth = cle.subtract_images(img_gpu, opened)
    return cle.pull(wth)


def entropy_table(input_folder, settings):
    """Shannon entropy of every brightfield image, to choose a focus threshold by inspection."""
    p = settings
    rows = []
    for f in sorted(os.listdir(input_folder)):
        if f.endswith(p["brightfield_suffix"]):
            img = imread(os.path.join(input_folder, f)).astype(np.float32)
            rows.append({"image": f, "entropy": shannon_entropy(img)})
    return pd.DataFrame(rows)
