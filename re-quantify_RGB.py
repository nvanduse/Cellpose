"""
Navigate to where this script is stored in CMD. Run this script directly from command line, without a .bat launcher.
Specify paths to input RGB overlay image folder, input mask folder, and output folder. Specify the other options.

Re-quantify curated Cellpose masks against RGB overlay images.

Expected layout:
    <output_dir>/masks/*.npy      # curated *_seg.npy files
    <image_dir>/*.{tif,tiff,png,jpg,jpeg}  # RGB overlays

Outputs:
    <output_dir>/object_data_requant_rgb.csv   # per-object intensities for RGB channels
    <output_dir>/masks/<basename>_overlay.png  # mask outlines for QC

Usage:
    python re-quantify_RGB \
        --image-dir /path/to/images \
        --output-dir /path/to/output \
        [--mask-dir /path/to/masks] \
        [--pixel-size 0.5]
"""

import argparse
import glob
import math
import os
import sys
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
from skimage import img_as_ubyte
from skimage.io import imread, imsave
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries


# ------------------------- CLI -------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Re-quantify Cellpose masks against RGB overlays.")
    p.add_argument("--image-dir", required=True, help="Directory containing RGB overlay images.")
    p.add_argument("--output-dir", required=True, help="Directory that contains masks/ and where outputs will be written.")
    p.add_argument("--mask-dir", default=None, help="Override mask directory (defaults to <output-dir>/masks).")
    p.add_argument("--pixel-size", type=float, default=None,
                   help="Microns per pixel; if provided, length/area are reported in microns.")
    return p.parse_args()


# ------------------------- Helpers -------------------------
def load_mask(mask_path: str) -> np.ndarray:
    """Load *_seg.npy saved by Cellpose (dict or ndarray)."""
    arr = np.load(mask_path, allow_pickle=True)
    try:
        obj = arr.item()
        if isinstance(obj, dict) and "masks" in obj:
            return np.asarray(obj["masks"])
    except Exception:
        pass
    return np.asarray(arr)


def find_image_file(base: str, image_dir: str) -> Optional[str]:
    """Return first matching overlay image for a given base."""
    patterns = [
        os.path.join(image_dir, f"{base}.tif"),
        os.path.join(image_dir, f"{base}.tiff"),
        os.path.join(image_dir, f"{base}.png"),
        os.path.join(image_dir, f"{base}.jpg"),
        os.path.join(image_dir, f"{base}.jpeg"),
        os.path.join(image_dir, f"{base}*.tif"),
        os.path.join(image_dir, f"{base}*.tiff"),
        os.path.join(image_dir, f"{base}*.png"),
        os.path.join(image_dir, f"{base}*.jpg"),
        os.path.join(image_dir, f"{base}*.jpeg"),
    ]
    matches: List[str] = []
    for pat in patterns:
        matches.extend(sorted(glob.glob(pat)))
    return matches[0] if matches else None


def extract_rgb(img: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return float RGB channels; grayscale is replicated, alpha is ignored."""
    if img.ndim == 2:
        gray = img.astype(np.float32)
        return gray, gray, gray
    if img.ndim == 3 and img.shape[2] >= 3:
        r = img[..., 0].astype(np.float32)
        g = img[..., 1].astype(np.float32)
        b = img[..., 2].astype(np.float32)
        return r, g, b
    raise ValueError(f"Unsupported image shape {img.shape}")


def compute_shape_props(mask: np.ndarray, pixel_size: Optional[float]) -> Dict[int, dict]:
    """Return per-label geometric properties."""
    props = {}
    for rp in regionprops(mask):
        major = float(rp.major_axis_length)
        minor = float(rp.minor_axis_length) if rp.minor_axis_length else math.nan
        lw_ratio = major / minor if minor and not math.isnan(minor) else math.nan
        area_px = int(rp.area)

        perimeter = float(getattr(rp, "perimeter_crofton", rp.perimeter))
        circularity = (4 * math.pi * rp.area) / (perimeter ** 2) if perimeter else math.nan

        entry = {
            "centroid": (int(rp.centroid[1]), int(rp.centroid[0])),  # x, y
            "coords": rp.coords,
            "area_px": area_px,
            "major_px": round(major, 3),
            "minor_px": round(minor, 3) if not math.isnan(minor) else math.nan,
            "lw_ratio": round(lw_ratio, 4) if not math.isnan(lw_ratio) else math.nan,
            "circularity": round(circularity, 4) if not math.isnan(circularity) else math.nan,
        }

        if pixel_size:
            entry.update({
                "area_um2": round(area_px * (pixel_size ** 2), 3),
                "major_um": round(major * pixel_size, 3),
                "minor_um": round(minor * pixel_size, 3) if not math.isnan(minor) else math.nan,
            })
        props[rp.label] = entry
    return props


def measure_intensity(channel: np.ndarray, coords: np.ndarray) -> Tuple[float, float]:
    """Mean and p99 for pixels at coords (row, col)."""
    if coords.size == 0:
        return math.nan, math.nan
    pix = channel[coords[:, 0], coords[:, 1]].astype(np.float32)
    return float(pix.mean()), float(np.percentile(pix, 99))


def write_overlay(mask: np.ndarray, ref_image: np.ndarray, out_path: str) -> None:
    """Save an outline overlay for quick QC."""
    outline = find_boundaries(mask, mode="outer")
    overlay = ref_image if ref_image.ndim == 3 else np.repeat(ref_image[..., None], 3, axis=2)
    overlay_uint8 = img_as_ubyte(overlay / (overlay.max() or 1))
    overlay_uint8[outline] = [255, 255, 255]
    imsave(out_path, overlay_uint8)


# ------------------------- Main -------------------------
def main():
    args = parse_args()

    image_dir = args.image_dir
    output_dir = args.output_dir
    mask_dir = args.mask_dir or os.path.join(output_dir, "masks")
    overlay_dir = os.path.join(output_dir, "masks")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(overlay_dir, exist_ok=True)

    mask_files = sorted(f for f in os.listdir(mask_dir) if f.endswith("_seg.npy"))
    if not mask_files:
        print(f"[ERR] No *_seg.npy files found in {mask_dir}")
        sys.exit(1)

    csv_rows: List[List[object]] = []

    for mask_name in mask_files:
        base = mask_name[:-8]  # drop _seg.npy
        mask_path = os.path.join(mask_dir, mask_name)
        mask = load_mask(mask_path)

        image_path = find_image_file(base, image_dir)
        if not image_path:
            warnings.warn(f"No RGB overlay image found for {base}; skipping.")
            continue

        image = imread(image_path)
        r, g, b = extract_rgb(image)

        shape_props = compute_shape_props(mask, args.pixel_size)
        labels = sorted(shape_props.keys())

        intensity_r: Dict[int, Tuple[float, float]] = {}
        intensity_g: Dict[int, Tuple[float, float]] = {}
        intensity_b: Dict[int, Tuple[float, float]] = {}
        for label in labels:
            coords = shape_props[label]["coords"]
            intensity_r[label] = measure_intensity(r, coords)
            intensity_g[label] = measure_intensity(g, coords)
            intensity_b[label] = measure_intensity(b, coords)

        for label in labels:
            props = shape_props[label]
            row = [
                base,
                int(label),
                f"{props['centroid'][0]},{props['centroid'][1]}",
                props["area_px"],
                round(intensity_r[label][0], 3),
                round(intensity_g[label][0], 3),
                round(intensity_b[label][0], 3),
                round(intensity_r[label][1], 3),
                round(intensity_g[label][1], 3),
                round(intensity_b[label][1], 3),
                props["major_px"],
                props["minor_px"] if not math.isnan(props["minor_px"]) else "",
                props["lw_ratio"] if not math.isnan(props["lw_ratio"]) else "",
                props["circularity"] if not math.isnan(props["circularity"]) else "",
            ]

            if args.pixel_size:
                row.extend([
                    props.get("major_um", ""),
                    props.get("minor_um", ""),
                    props.get("area_um2", ""),
                ])

            csv_rows.append(row)

        overlay_path = os.path.join(overlay_dir, f"{base}_overlay.png")
        try:
            write_overlay(mask, image, overlay_path)
        except Exception as exc:  # pragma: no cover
            warnings.warn(f"Could not write overlay for {base}: {exc}")

    # ------------------------- Write CSV -------------------------
    header = [
        "Image",
        "Object Number",
        "Object Location (x,y)",
        "Area (px)",
        "Red Mean",
        "Green Mean",
        "Blue Mean",
        "Red p99",
        "Green p99",
        "Blue p99",
        "Length (px)",
        "Width (px)",
        "L/W Ratio",
        "Circularity",
    ]
    if args.pixel_size:
        header.extend(["Length (um)", "Width (um)", "Area (um^2)"])

    out_csv = os.path.join(output_dir, "object_data_requant_rgb.csv")
    import csv  # local import to keep header nearby

    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(csv_rows)

    print(f"[INFO] Saved quantification -> {out_csv} ({len(csv_rows)} objects)")


if __name__ == "__main__":
    main()
