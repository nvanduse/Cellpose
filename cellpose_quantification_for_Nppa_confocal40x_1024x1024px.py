"""
cellpose_quantify_Nppa_v1_requant.py

Quantify cell area and cytoplasm area (cell minus nuclei) using edited masks.
Uses per-channel TIFF naming to find sample bases, but metrics are derived
only from the masks.

Usage:
# put the script in the working directory, which should contain subfolders that have masks and images for individual channels
# open command prompt and navigate to the working directory
# Typical run:
# python cellpose_quantification_for_Nppa_confocal40x_1024x1024px.py --image-dir <image_dir> --output-dir <output_dir>
# Optional:
# --mask-cell-dir <path>   Override output_dir/masks_cell
# --mask-nuc-dir <path>    Override output_dir/masks_nuc
# --channel-tags c1 c2 c3 c4  Custom tags for far_red/red/green/blue channels
#
# Example:
# python cellpose_quantification_for_Nppa_confocal40x_1024x1024px --image-dir all_images --mask-cell-dir masks/masks_cell --mask-nuc-dir masks/masks_nuc --output-dir quantification_data

Expected inputs (per sample):
  <output_dir>/masks_cell/<base>_cell_seg.npy
  <output_dir>/masks_nuc/<base>_nuc_seg.npy
  <image_dir>/*_c1.tif ... *_c4.tif (used to discover sample bases)

Expected outputs:
    Two files
    - object_data_run_log.txt in the same --output-dir (run log with date/time, script path, and input/output
    directories)
    - object_data.csv in the specified --output-dir (cell-level measurements)
    - object_data_nuclei.csv in the specified --output-dir (nucleus-level measurements)
        The csv contains the following:
            - Image
            - Cell label and location
            - Areas, length, width, circularity
            - Nuclei count per cell
            - Mean and p99 intensity for each channel in cell and cytoplasm


"""

import argparse
import glob
import os
import warnings
import datetime
import sys
import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from skimage.io import imread
from skimage.measure import regionprops


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Quantify cell/cytoplasm areas from edited masks.")
    p.add_argument("--image-dir", required=True, help="Directory containing per-channel TIFFs (c1–c4).")
    p.add_argument("--output-dir", required=True, help="Directory with masks_cell/ and masks_nuc/.")
    p.add_argument("--mask-cell-dir", default=None, help="Override cell mask directory.")
    p.add_argument("--mask-nuc-dir", default=None, help="Override nucleus mask directory.")
    p.add_argument("--channel-tags", nargs=4, default=["c1", "c2", "c3", "c4"],
                   help="Filename tags for far-red/red/green/blue channels (default: c1 c2 c3 c4).")
    return p.parse_args()


def load_mask(mask_path: str) -> np.ndarray:
    arr = np.load(mask_path, allow_pickle=True)
    try:
        obj = arr.item()
        if isinstance(obj, dict) and "masks" in obj:
            return np.asarray(obj["masks"])
    except Exception:
        pass
    return np.asarray(arr)

def find_channel_file(base: str, image_dir: str, tag: str) -> Optional[str]:
    patterns = [
        os.path.join(image_dir, f"{base}*_{tag}.tif"),
        os.path.join(image_dir, f"{base}*_{tag}.tiff"),
    ]
    matches: List[str] = []
    for pat in patterns:
        matches.extend(sorted(glob.glob(pat)))
    return matches[0] if matches else None

def to_grayscale(img: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if img is None:
        return None
    if img.ndim == 2:
        return img.astype(np.float32)
    if img.ndim == 3:
        return img.mean(axis=2).astype(np.float32)
    return None

def measure_intensity(gray: Optional[np.ndarray], mask: np.ndarray) -> Tuple[float, float]:
    if gray is None or not np.any(mask):
        return math.nan, math.nan
    pix = gray[mask].astype(np.float32)
    return float(pix.mean()), float(np.percentile(pix, 99))

def fmt(val: float, ndigits: int = 3):
    return "" if np.isnan(val) else round(float(val), ndigits)


def collect_bases(image_dir: str, channel_tags: List[str]) -> List[str]:
    bases = set()
    for tag in channel_tags:
        for ext in ("tif", "tiff"):
            for path in glob.glob(os.path.join(image_dir, f"*_{tag}.{ext}")):
                stem = os.path.splitext(os.path.basename(path))[0]
                suffix = f"_{tag}"
                if stem.endswith(suffix):
                    bases.add(stem[:-len(suffix)])
    return sorted(bases)


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    mask_cell_dir = args.mask_cell_dir or os.path.join(args.output_dir, "masks_cell")
    mask_nuc_dir = args.mask_nuc_dir or os.path.join(args.output_dir, "masks_nuc")
    channel_order = ["far_red", "red", "green", "blue"]
    channel_tag_map = dict(zip(channel_order, args.channel_tags))

    bases = collect_bases(args.image_dir, args.channel_tags)
    if not bases:
        print("[WARN] No per-channel TIFFs found. Check --channel-tags and image_dir.")
        return

    cell_rows = []
    nuc_rows = []
    for base in bases:
        cell_path = os.path.join(mask_cell_dir, f"{base}_cell_seg.npy")
        nuc_path = os.path.join(mask_nuc_dir, f"{base}_nuc_seg.npy")
        if not (os.path.exists(cell_path) and os.path.exists(nuc_path)):
            warnings.warn(f"Missing mask(s) for {base}; skipping.")
            continue

        cell_masks = load_mask(cell_path)
        nuc_masks = load_mask(nuc_path)
        channel_images: Dict[str, Optional[np.ndarray]] = {}
        for chan in channel_order:
            tag = channel_tag_map[chan]
            img_path = find_channel_file(base, args.image_dir, tag)
            if not img_path:
                warnings.warn(f"Missing channel image for {base} ({tag}); intensities will be blank.")
                channel_images[chan] = None
                continue
            channel_images[chan] = to_grayscale(imread(img_path))

        nuc_props: Dict[int, dict] = {}
        nuc_to_cell: Dict[int, int] = {}
        cell_to_nucs: Dict[int, List[int]] = {}
        for rp in regionprops(nuc_masks):
            label = int(rp.label)
            coords = rp.coords
            cell_labels = cell_masks[coords[:, 0], coords[:, 1]]
            nonzero = cell_labels[cell_labels != 0]
            parent = 0
            if nonzero.size == cell_labels.size:
                unique = np.unique(nonzero)
                if unique.size == 1:
                    parent = int(unique[0])
                    cell_to_nucs.setdefault(parent, []).append(label)
            nuc_to_cell[label] = parent

            major = float(rp.major_axis_length)
            minor = float(rp.minor_axis_length) if rp.minor_axis_length else math.nan
            perimeter = float(getattr(rp, "perimeter_crofton", rp.perimeter))
            circularity = (4 * math.pi * rp.area) / (perimeter ** 2) if perimeter else math.nan

            nuc_props[label] = {
                "area": int(rp.area),
                "major": major,
                "minor": minor,
                "circularity": circularity,
                "centroid": (int(rp.centroid[1]), int(rp.centroid[0])),
            }

        for rp in regionprops(cell_masks):
            lbl = int(rp.label)
            bin_cell = (cell_masks == lbl)
            nuc_labels = cell_to_nucs.get(lbl, [])
            nuc_in_cell_mask = np.isin(nuc_masks, nuc_labels)
            cyto_mask = bin_cell & (~nuc_in_cell_mask)

            cell_area = int(rp.area)
            nuc_in_cell = int(np.count_nonzero(nuc_in_cell_mask))
            cyto_area = int(np.count_nonzero(cyto_mask))
            cyto_frac = (cyto_area / float(cell_area)) if cell_area > 0 else np.nan

            major = float(rp.major_axis_length)
            minor = float(rp.minor_axis_length) if rp.minor_axis_length else math.nan
            perimeter = float(getattr(rp, "perimeter_crofton", rp.perimeter))
            circularity = (4 * math.pi * rp.area) / (perimeter ** 2) if perimeter else math.nan
            cx, cy = int(rp.centroid[1]), int(rp.centroid[0])

            row = [
                base,
                lbl,
                f"{cx},{cy}",
                cell_area,
                nuc_in_cell,
                cyto_area,
                fmt(cyto_frac, 4),
                fmt(major, 2),
                fmt(minor, 2),
                fmt(circularity, 3),
                len(nuc_labels),
            ]
            for chan in channel_order:
                mean_cell, p99_cell = measure_intensity(channel_images[chan], bin_cell)
                mean_cyto, p99_cyto = measure_intensity(channel_images[chan], cyto_mask)
                row.extend([
                    fmt(mean_cell, 3),
                    fmt(p99_cell, 3),
                    fmt(mean_cyto, 3),
                    fmt(p99_cyto, 3),
                ])
            cell_rows.append(row)

        for lbl, props in nuc_props.items():
            mask = (nuc_masks == lbl)
            cx, cy = props["centroid"]
            row = [
                base,
                lbl,
                nuc_to_cell.get(lbl, 0),
                f"{cx},{cy}",
                props["area"],
                fmt(props["major"], 2),
                fmt(props["minor"], 2),
                fmt(props["circularity"], 3),
            ]
            for chan in channel_order:
                mean_nuc, p99_nuc = measure_intensity(channel_images[chan], mask)
                row.extend([
                    fmt(mean_nuc, 3),
                    fmt(p99_nuc, 3),
                ])
            nuc_rows.append(row)

    out_csv = os.path.join(args.output_dir, "object_data.csv")
    out_nuc_csv = os.path.join(args.output_dir, "object_data_nuclei.csv")
    log_path = os.path.join(args.output_dir, "object_data_run_log.txt")
    log_lines = [
        f"===== Quantification run log: {datetime.datetime.now().isoformat()} =====",
        f"Script: {os.path.abspath(sys.argv[0])}",
        f"Image dir: {os.path.abspath(args.image_dir)}",
        f"Output dir: {os.path.abspath(args.output_dir)}",
        f"Mask cell dir: {os.path.abspath(mask_cell_dir)}",
        f"Mask nuc dir: {os.path.abspath(mask_nuc_dir)}",
        f"Channel tags: {', '.join(args.channel_tags)}",
        f"Bases detected: {len(bases)}",
        f"Cells quantified: {len(cell_rows)}",
        f"Nuclei quantified: {len(nuc_rows)}",
    ]
    with open(log_path, "w", newline="") as f:
        f.write("\n".join(log_lines) + "\n")
    header = [
        "Image",
        "Cell Label",
        "Cell Location",
        "Cell Area (total)",
        "Nucleus Area (in cell)",
        "Cytoplasm Area",
        "Cytoplasm Fraction",
        "Cell Length",
        "Cell Width",
        "Cell Circularity",
        "Nuclei Count",
    ]
    for chan in channel_order:
        label = chan.replace("_", " ").title().replace(" ", "")
        header.extend([
            f"{label} Cell Mean",
            f"{label} Cell p99",
            f"{label} Cytoplasm Mean",
            f"{label} Cytoplasm p99",
        ])

    nuc_header = [
        "Image",
        "Nucleus Label",
        "Parent Cell Label",
        "Nucleus Location",
        "Nucleus Area",
        "Nucleus Length",
        "Nucleus Width",
        "Nucleus Circularity",
    ]
    for chan in channel_order:
        label = chan.replace("_", " ").title().replace(" ", "")
        nuc_header.extend([
            f"{label} Nucleus Mean",
            f"{label} Nucleus p99",
        ])

    with open(out_csv, "w", newline="") as f:
        f.write(",".join(header) + "\n")
        for row in cell_rows:
            f.write(",".join(map(str, row)) + "\n")

    with open(out_nuc_csv, "w", newline="") as f:
        f.write(",".join(nuc_header) + "\n")
        for row in nuc_rows:
            f.write(",".join(map(str, row)) + "\n")

    print(f"[DONE] Cell table -> {out_csv}")
    print(f"[DONE] Nucleus table -> {out_nuc_csv}")


if __name__ == "__main__":
    main()
