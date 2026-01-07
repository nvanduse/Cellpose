"""
Navigate to where this script is stored in CMD. Run this script directly from command line, without a .bat launcher.
Specify paths to input image folder, input mask folder, and output folder. Specify the other options.
Note this script is for four color images - it will look for and quantify from the individual channel images, not overlays.

Re-quantify curated Cellpose masks against per-channel TIFFs and optionally
run SarcAsM on red-channel ROIs.

Expected layout (same as prior cellpose scripts):
    <output_dir>/masks/*.npy      # curated *_seg.npy files
    <image_dir>/*.tif[f]          # per-channel images, e.g. *_c1.tif ... *_c4.tif

Channel mapping (can be changed via CLI flags):
    c1 -> far_red, c2 -> red, c3 -> green, c4 -> blue

Outputs:
    <output_dir>/object_data_requant.csv   # per-object intensities for all channels
    <output_dir>/sarcasm_inputs/*.tif      # red-channel masked ROIs (one per object)
    SarcAsM outputs alongside each ROI if SarcAsM is installed and not skipped

Usage:
    python re-quantify_and_sarcasm.py \
        --image-dir /path/to/images \
        --output-dir /path/to/output \
        [--mask-dir /path/to/masks] \
        [--pixel-size 0.5] \
        [--roi-padding 5] \
        [--sarcasm-workers 3] \
        [--skip-sarcasm]
"""

import argparse
import glob
import math
import os
import sys
import warnings
from multiprocessing import Pool
from typing import Dict, List, Optional, Tuple

import numpy as np
from skimage import img_as_ubyte
from skimage.io import imread, imsave
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries


# ------------------------- CLI -------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Re-quantify Cellpose masks and run SarcAsM on red-channel ROIs.")
    p.add_argument("--image-dir", required=True, help="Directory containing per-channel TIFFs (c1–c4).")
    p.add_argument("--output-dir", required=True, help="Directory that contains masks/ and where outputs will be written.")
    p.add_argument("--mask-dir", default=None, help="Override mask directory (defaults to <output-dir>/masks).")
    p.add_argument("--pixel-size", type=float, default=None,
                   help="Microns per pixel; if provided, length/area are reported in microns.")
    p.add_argument("--roi-padding", type=int, default=5, help="Padding (pixels) around each ROI when writing SarcAsM inputs.")
    p.add_argument("--sarcasm-workers", type=int, default=1, help="Parallel workers for SarcAsM.")
    p.add_argument("--skip-sarcasm", action="store_true", help="Only write ROI TIFFs; do not run SarcAsM.")
    p.add_argument("--sarcasm-channel", default="red",
                   choices=["far_red", "red", "green", "blue"],
                   help="Which channel to use for SarcAsM ROIs (default: red).")
    p.add_argument("--channel-tags", nargs=4, default=["c1", "c2", "c3", "c4"],
                   help="Filename tags for far-red/red/green/blue channels (default: c1 c2 c3 c4).")
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


def find_channel_file(base: str, image_dir: str, tag: str) -> Optional[str]:
    """Return first matching tif/tiff for a given base and channel tag."""
    patterns = [
        os.path.join(image_dir, f"{base}*_{tag}.tif"),
        os.path.join(image_dir, f"{base}*_{tag}.tiff"),
    ]
    matches: List[str] = []
    for pat in patterns:
        matches.extend(sorted(glob.glob(pat)))
    return matches[0] if matches else None


def find_overlay_file(base: str, image_dir: str) -> Optional[str]:
    """Return the c1-4 overlay image if present."""
    patterns = [
        os.path.join(image_dir, f"{base}*_c1-4.tif"),
        os.path.join(image_dir, f"{base}*_c1-4.tiff"),
    ]
    matches: List[str] = []
    for pat in patterns:
        matches.extend(sorted(glob.glob(pat)))
    return matches[0] if matches else None


def to_grayscale(img: np.ndarray) -> np.ndarray:
    """Convert multi-channel images to grayscale float for quantification."""
    if img.ndim == 2:
        return img.astype(np.float32)
    if img.ndim == 3:
        return img.mean(axis=2).astype(np.float32)
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


def measure_intensity(gray: np.ndarray, coords: np.ndarray) -> Tuple[float, float]:
    """Mean and p99 for pixels at coords (row, col)."""
    if coords.size == 0:
        return math.nan, math.nan
    pix = gray[coords[:, 0], coords[:, 1]].astype(np.float32)
    return float(pix.mean()), float(np.percentile(pix, 99))


def write_overlay(mask: np.ndarray, ref_image: np.ndarray, out_path: str) -> None:
    """Save an outline overlay for quick QC."""
    outline = find_boundaries(mask, mode="outer")
    overlay = ref_image if ref_image.ndim == 3 else np.repeat(ref_image[..., None], 3, axis=2)
    overlay_uint8 = img_as_ubyte(overlay / (overlay.max() or 1))
    overlay_uint8[outline] = [255, 255, 255]
    imsave(out_path, overlay_uint8)


def crop_roi(img: np.ndarray, mask: np.ndarray, label: int, padding: int) -> np.ndarray:
    """Return cropped ROI image (masked) for a label."""
    rp = regionprops((mask == label).astype(np.uint8))[0]
    minr, minc, maxr, maxc = rp.bbox
    minr = max(minr - padding, 0)
    minc = max(minc - padding, 0)
    maxr = min(maxr + padding, mask.shape[0])
    maxc = min(maxc + padding, mask.shape[1])
    roi = np.zeros((maxr - minr, maxc - minc), dtype=img.dtype)
    submask = (mask[minr:maxr, minc:maxc] == label)
    roi[submask] = img[minr:maxr, minc:maxc][submask]
    return roi


def run_sarcasm_on_file(
    tif_path: str,
    max_patch_size=(2048, 2048),
    pixel_size: Optional[float] = None,
) -> Tuple[str, bool, str]:
    """Wrapper to execute SarcAsM on a single TIFF; returns (path, success, message)."""
    try:
        from sarcasm import Structure  # type: ignore
    except ImportError:
        return tif_path, False, "sarcasm package not installed"

    try:
        if pixel_size:
            sarc = Structure(tif_path, pixelsize=pixel_size)
        else:
            sarc = Structure(tif_path)

        # pass pixel size into detect_sarcomeres if supported
        detect_kwargs = {"max_patch_size": max_patch_size}
        try:
            import inspect
            params = inspect.signature(sarc.detect_sarcomeres).parameters
            if pixel_size:
                if "pixel_size" in params:
                    detect_kwargs["pixel_size"] = pixel_size
                elif "pixel_size_um" in params:
                    detect_kwargs["pixel_size_um"] = pixel_size
                elif "px_size" in params:
                    detect_kwargs["px_size"] = pixel_size
        except Exception:
            pass

        sarc.detect_sarcomeres(**detect_kwargs)
        sarc.full_analysis_structure(frames="all")
        return tif_path, True, "ok"
    except Exception as exc:  # pragma: no cover - run-time protection
        return tif_path, False, str(exc)


# ------------------------- Main -------------------------
def main():
    args = parse_args()

    image_dir = args.image_dir
    output_dir = args.output_dir
    mask_dir = args.mask_dir or os.path.join(output_dir, "masks")

    sarcasm_input_dir = os.path.join(output_dir, "sarcasm_inputs")
    os.makedirs(sarcasm_input_dir, exist_ok=True)

    mask_files = sorted(f for f in os.listdir(mask_dir) if f.endswith("_seg.npy"))
    if not mask_files:
        print(f"[ERR] No *_seg.npy files found in {mask_dir}")
        sys.exit(1)

    channel_order = ["far_red", "red", "green", "blue"]
    channel_tag_map = dict(zip(channel_order, args.channel_tags))

    csv_rows: List[List[object]] = []
    sarcasm_tifs: List[str] = []

    for mask_name in mask_files:
        base = mask_name[:-8]  # drop _seg.npy
        mask_path = os.path.join(mask_dir, mask_name)
        mask = load_mask(mask_path)

        # locate channel images
        channel_paths: Dict[str, Optional[str]] = {}
        for chan in channel_order:
            tag = channel_tag_map[chan]
            channel_paths[chan] = find_channel_file(base, image_dir, tag)

        if not any(channel_paths.values()):
            warnings.warn(f"No channel images found for {base}; skipping.")
            continue

        # pick a reference image for overlays (prefer c1-4 overlay)
        ref_path = find_overlay_file(base, image_dir)
        if not ref_path:
            ref_path = next((p for p in channel_paths.values() if p), None)
        ref_image = imread(ref_path) if ref_path else mask.astype(np.uint8)

        shape_props = compute_shape_props(mask, args.pixel_size)
        labels = sorted(shape_props.keys())

        # load grayscale images per available channel once
        gray_cache: Dict[str, np.ndarray] = {}
        for chan, path in channel_paths.items():
            if path:
                gray_cache[chan] = to_grayscale(imread(path))

        # intensity per label per channel
        intensity: Dict[str, Dict[int, Tuple[float, float]]] = {chan: {} for chan in channel_order}
        for chan, gray in gray_cache.items():
            for label in labels:
                mean_val, p99_val = measure_intensity(gray, shape_props[label]["coords"])
                intensity[chan][label] = (mean_val, p99_val)

        # rows
        for label in labels:
            props = shape_props[label]
            row = [
                base,
                int(label),
                f"{props['centroid'][0]},{props['centroid'][1]}",
                props["area_px"],
            ]
            # intensities
            for chan in channel_order:
                mean_val = intensity.get(chan, {}).get(label, (math.nan, math.nan))[0]
                row.append(round(mean_val, 3) if not math.isnan(mean_val) else "")
            for chan in channel_order:
                p99_val = intensity.get(chan, {}).get(label, (math.nan, math.nan))[1]
                row.append(round(p99_val, 3) if not math.isnan(p99_val) else "")

            row.extend([
                props["major_px"],
                props["minor_px"] if not math.isnan(props["minor_px"]) else "",
                props["lw_ratio"] if not math.isnan(props["lw_ratio"]) else "",
                props["circularity"] if not math.isnan(props["circularity"]) else "",
            ])

            if args.pixel_size:
                row.extend([
                    props.get("major_um", ""),
                    props.get("minor_um", ""),
                    props.get("area_um2", ""),
                ])

            csv_rows.append(row)

        # QC overlay saved next to masks
        overlay_path = os.path.join(mask_dir, f"{base}_overlay.png")
        try:
            write_overlay(mask, ref_image, overlay_path)
        except Exception as exc:  # pragma: no cover
            warnings.warn(f"Could not write overlay for {base}: {exc}")

        # Prepare SarcAsM ROIs (user-selected channel)
        sarc_channel = args.sarcasm_channel
        sarc_path = channel_paths[sarc_channel]
        if sarc_path:
            sarc_img_full = imread(sarc_path)
            sarc_gray = to_grayscale(sarc_img_full)
            for label in labels:
                roi = crop_roi(sarc_gray, mask, label, padding=args.roi_padding)
                roi_fname = f"{base}_obj{label:03d}_{sarc_channel}ROI.tif"
                roi_path = os.path.join(sarcasm_input_dir, roi_fname)
                imsave(roi_path, roi.astype(sarc_img_full.dtype))
                sarcasm_tifs.append(roi_path)
        else:
            warnings.warn(f"No {sarc_channel} channel for {base}; SarcAsM skipped for this image.")

    # ------------------------- Write CSV -------------------------
    header = [
        "Image",
        "Object Number",
        "Object Location (x,y)",
        "Area (px)",
        "FarRed Mean",
        "Red Mean",
        "Green Mean",
        "Blue Mean",
        "FarRed p99",
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

    out_csv = os.path.join(output_dir, "object_data_requant.csv")
    import csv  # local import to keep header nearby

    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(csv_rows)

    print(f"[INFO] Saved quantification -> {out_csv} ({len(csv_rows)} objects)")
    print(f"[INFO] SarcAsM ROIs -> {sarcasm_input_dir} ({len(sarcasm_tifs)} files)")

    # ------------------------- Run SarcAsM -------------------------
    if sarcasm_tifs and not args.skip_sarcasm:
        print("[INFO] Running SarcAsM on red-channel ROIs...")
        if args.sarcasm_workers > 1:
            with Pool(args.sarcasm_workers) as pool:
                results = pool.starmap(
                    run_sarcasm_on_file,
                    [(p, (2048, 2048), args.pixel_size) for p in sarcasm_tifs],
                )
        else:
            results = [run_sarcasm_on_file(p, (2048, 2048), args.pixel_size) for p in sarcasm_tifs]

        failures = [r for r in results if not r[1]]
        if failures:
            print("[WARN] Some SarcAsM runs failed:")
            for path, _, msg in failures:
                print(f"    {path}: {msg}")
        else:
            print("[INFO] SarcAsM completed for all ROIs.")
    else:
        print("[INFO] SarcAsM run skipped.")


if __name__ == "__main__":
    main()
