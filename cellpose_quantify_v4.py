# cellpose_quantify_v4.py
"""
Quantify Cellpose masks, optionally filtering objects by Length/Width ratio
and generating overlays that show only retained objects. Outputs per-object
channel averages and p99 intensity (99th percentile) per channel.

CLI usage (args in order):
    python cellpose_quantify.py image_dir output_dir [lwr_filter] [lwr_cutoff]
                                 [channel] [model] [diam_mode] [diameter]

- lwr_filter: 1 to enable filtering; 0 to disable (default 0).
- lwr_cutoff: numeric; objects must have (major_axis_length / minor_axis_length) > cutoff.
              Ignored if lwr_filter=0. Default 0 (no filter).
- Remaining args are optional; used only for recording run parameters in CSV header log.

Outputs:
- object_data.csv in output_dir
- *_overlay.png files in output_dir/masks (filtered outlines if filter active; else all)
"""

import sys
import os
import csv
import numpy as np
from skimage.io import imread, imsave
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries
from skimage.color import gray2rgb
from skimage import img_as_ubyte

# -------------------- Parse args --------------------
# Required
try:
    image_dir = sys.argv[1]
    output_dir = sys.argv[2]
except IndexError:
    print("Usage: python cellpose_quantify.py image_dir output_dir [lwr_filter] [lwr_cutoff] [channel] [model] [diam_mode] [diameter]")
    sys.exit(1)

def _get_arg(idx, default=None):
    return sys.argv[idx] if len(sys.argv) > idx else default

lwr_filter_arg = _get_arg(3, "0")
lwr_cutoff_arg = _get_arg(4, "0")
channel_arg    = _get_arg(5, "NA")
model_arg      = _get_arg(6, "NA")
diam_mode_arg  = _get_arg(7, "NA")  # "auto" or "manual"
diameter_arg   = _get_arg(8, "NA")

# normalize
try:
    LWR_FILTER = bool(int(lwr_filter_arg))
except (TypeError, ValueError):
    LWR_FILTER = False

try:
    LWR_CUTOFF = float(lwr_cutoff_arg)
except (TypeError, ValueError):
    LWR_CUTOFF = 0.0

# -------------------- Prep output --------------------
mask_dir = os.path.join(output_dir, "masks")
output_csv = os.path.join(output_dir, "object_data.csv")

# For log column header
prompt_log = (
    f"Prompt_Log: channel={channel_arg}; model={model_arg}; "
    f"diam_mode={diam_mode_arg}; diameter={diameter_arg}; "
    f"L/W_filter={int(LWR_FILTER)}; cutoff={LWR_CUTOFF}"
)

# -------------------- Data accumulation --------------------
data_rows = []  # each row = [..., length, width, lw_ratio]

# -------------------- Process masks --------------------
# We will generate filtered overlays; to do that we need to know which objects pass.
for fname in os.listdir(mask_dir):
    if not fname.endswith("_seg.npy"):
        continue

    base = fname[:-8]  # strip _seg.npy
    mask_path = os.path.join(mask_dir, fname)

    # Look for corresponding image
    image_path = None
    for ext in (".png", ".tif", ".tiff", ".jpg", ".jpeg"):
        candidate = os.path.join(image_dir, base + ext)
        if os.path.exists(candidate):
            image_path = candidate
            break
    if image_path is None:
        print(f"No image found for mask {fname}")
        continue

    # Load image & mask
    mask_data = np.load(mask_path, allow_pickle=True).item()
    mask = np.asarray(mask_data["masks"])
    image = imread(image_path)

    labels = np.unique(mask)
    labels = labels[labels != 0]  # exclude background

    kept_labels = []

    for label in labels:
        binary = (mask == label).astype(np.uint8)
        # regionprops: intensity_image optional; we supply original image
        props_list = regionprops(binary, intensity_image=image)
        if not props_list:
            continue
        p = props_list[0]

        centroid = p.centroid
        area = p.area
        major = float(p.major_axis_length)
        minor = float(p.minor_axis_length) if p.minor_axis_length != 0 else np.nan

        # compute mean and p99 per channel (robust to grayscale)
        mean_rgb = [np.nan, np.nan, np.nan]
        p99_rgb = [np.nan, np.nan, np.nan]
        coords = p.coords
        pix = image[coords[:, 0], coords[:, 1]]
        if pix.ndim == 1:  # grayscale -> shape (N,)
            pix = pix[:, None]

        if pix.shape[1] >= 3:
            first_three = pix[:, :3]
            mean_rgb = first_three.mean(axis=0)
            p99_rgb = np.percentile(first_three, 99, axis=0)
        else:
            mean_val = pix.mean()
            p99_val = np.percentile(pix, 99)
            mean_rgb = [mean_val] * 3
            p99_rgb = [p99_val] * 3

        # L/W ratio
        if np.isnan(minor) or minor == 0:
            lw_ratio = np.nan
        else:
            lw_ratio = major / minor

        # filter?
        keep = True
        if LWR_FILTER and LWR_CUTOFF > 0 and not np.isnan(lw_ratio):
            keep = lw_ratio > LWR_CUTOFF

        if keep:
            kept_labels.append(label)
            # Store original label for filtering, but will renumber later if needed
            data_rows.append([
                base,
                int(label),  # Will be renumbered later if filtering is used
                f"{int(centroid[1])},{int(centroid[0])}",  # x,y
                int(area),
                round(float(mean_rgb[0]), 2),
                round(float(mean_rgb[1]), 2),
                round(float(mean_rgb[2]), 2),
                round(float(p99_rgb[0]), 2),
                round(float(p99_rgb[1]), 2),
                round(float(p99_rgb[2]), 2),
                round(major, 2),
                round(minor, 2),
                round(float(lw_ratio), 3) if not np.isnan(lw_ratio) else ""
            ])

    # ----- Overlay generation -----
    if LWR_FILTER and len(kept_labels) > 0:
        # make a filtered mask
        filtered_mask = np.zeros_like(mask, dtype=mask.dtype)
        for k in kept_labels:
            filtered_mask[mask == k] = k
        outline = find_boundaries(filtered_mask, mode="outer")
    else:
        # Show all outlines (either no filter, or filter but nothing passed -> show nothing?)
        # Choice: show all if none kept, so user sees that segmentation ran. Adjust if you prefer.
        if LWR_FILTER and len(kept_labels) == 0:
            outline = np.zeros_like(mask, dtype=bool)
        else:
            outline = find_boundaries(mask, mode="outer")

    overlay = np.copy(image)
    if overlay.ndim == 2 or (overlay.ndim == 3 and overlay.shape[2] == 1):
        overlay = gray2rgb(overlay)

    # normalize to uint8
    # guard for all-zero images
    maxval = overlay.max()
    if maxval > 0:
        overlay_uint8 = img_as_ubyte(overlay / maxval)
    else:
        overlay_uint8 = np.zeros_like(overlay, dtype=np.uint8)

    # draw white outlines
    overlay_uint8[outline] = [255, 255, 255]

    overlay_path = os.path.join(mask_dir, f"{base}_overlay.png")
    imsave(overlay_path, overlay_uint8)

# -------------------- Write CSV --------------------
# Renumber objects if filtering was used
if LWR_FILTER:
    # Group data by image and renumber within each image
    from collections import defaultdict
    image_data = defaultdict(list)
    
    # Group rows by image
    for row in data_rows:
        image_name = row[0]
        image_data[image_name].append(row)
    
    # Renumber within each image and rebuild data_rows
    renumbered_data = []
    for image_name, rows in image_data.items():
        for new_num, row in enumerate(rows, 1):
            row[1] = new_num  # Replace original object number with new sequential number
            renumbered_data.append(row)
    
    data_rows = renumbered_data

header = [
    'Image', 'Object Number', 'Object Location', 'Area',
    'Red Avg', 'Green Avg', 'Blue Avg',
    'Red p99', 'Green p99', 'Blue p99',
    'Length', 'Width', 'L/W Ratio',
    prompt_log  # appended log column
]

with open(output_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    # For each data row, add empty cell for prompt_log column
    for row in data_rows:
        writer.writerow(row + [""])

print(f"Saved results to {output_csv}")
print(f"Objects retained: {len(data_rows)}; L/W filter={'ON' if LWR_FILTER else 'OFF'}; cutoff={LWR_CUTOFF}")
