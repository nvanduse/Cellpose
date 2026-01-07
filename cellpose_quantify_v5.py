"""
cellpose_quantify_v5.py
--------------------------------------------------------------
• Reads Cellpose *_seg.npy masks in <output>/masks
• Measures per-object features (area, length, width, L/W ratio, mean RGB)
• Writes   <output>/object_data.csv
• Makes    <output>/masks/<basename>_overlay.png       (white outlines)
• Makes    <output>/masks_for_sarcasm/<basename>_centreROI_ch#.png
           -– grayscale image keeping only the mask nearest the centre
--------------------------------------------------------------
Usage (called by the batch file):

    python cellpose_quantify_v5.py image_dir output_dir \
           [lwr_filter] [lwr_cutoff] [channel] [model] [diam_mode] [diameter]
"""

import sys
import os
import csv
import math
import numpy as np

from skimage.io import imread, imsave
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries
from skimage.color import gray2rgb
from skimage import img_as_ubyte


# ────────────────────── Parse CLI arguments ──────────────────────
try:
    image_dir  = sys.argv[1]
    output_dir = sys.argv[2]
except IndexError:
    print("Usage: python cellpose_quantify_v5.py image_dir output_dir [options]")
    sys.exit(1)

def _arg(idx, default=None):
    return sys.argv[idx] if len(sys.argv) > idx else default

lwr_filter_arg = _arg(3, "0")
lwr_cutoff_arg = _arg(4, "0")
channel_arg    = _arg(5, "NA")     # 0–3 or "NA"
model_arg      = _arg(6, "NA")
diam_mode_arg  = _arg(7, "NA")
diameter_arg   = _arg(8, "NA")

LWR_FILTER = bool(int(lwr_filter_arg))
LWR_CUTOFF = float(lwr_cutoff_arg)

# ────────────────────── Paths ──────────────────────
mask_dir    = os.path.join(output_dir, "masks")
sarcasm_dir = os.path.join(output_dir, "masks_for_sarcasm")
os.makedirs(sarcasm_dir, exist_ok=True)

output_csv  = os.path.join(output_dir, "object_data.csv")
prompt_log  = (f"Prompt_Log: channel={channel_arg}; model={model_arg}; "
               f"diam_mode={diam_mode_arg}; diameter={diameter_arg}; "
               f"L/W_filter={int(LWR_FILTER)}; cutoff={LWR_CUTOFF}")

# ────────────────────── Helper: single-channel extraction ──────────────────────
def extract_channel(img, idx: int) -> np.ndarray:
    """
    Return uint8 grayscale image for channel *idx* (0-2).
    If img is already grayscale, just rescale to uint8.
    """
    if img.ndim == 2:                      # grayscale
        return img_as_ubyte(img / (img.max() or 1))
    if idx not in (0, 1, 2):
        raise ValueError("Channel index must be 0, 1 or 2")
    ch = img[..., idx].astype(np.float32)
    return img_as_ubyte(ch / (ch.max() or 1))

# ────────────────────── Main loop over each *_seg.npy ──────────────────────
data_rows = []

for fname in os.listdir(mask_dir):
    if not fname.endswith("_seg.npy"):
        continue

    base      = fname[:-8]                     # strip "_seg.npy"
    mask_path = os.path.join(mask_dir, fname)

    # Find corresponding raw image (any common extension)
    image_path = next(
        (os.path.join(image_dir, base + ext)
         for ext in (".png", ".tif", ".tiff", ".jpg", ".jpeg")
         if os.path.exists(os.path.join(image_dir, base + ext))),
        None)
    if image_path is None:
        print(f"[WARN] No raw image found for {fname}")
        continue

    # ── Load mask & image ──
    mask_data = np.load(mask_path, allow_pickle=True).item()
    mask      = np.asarray(mask_data["masks"])
    image     = imread(image_path)
    H, W      = mask.shape

    # ── Identify label whose centroid is nearest image centre ──
    labels = np.unique(mask)
    labels = labels[labels != 0]
    if labels.size == 0:
        print(f"[WARN] Empty mask in {fname}")
        continue

    centroids = {}
    for lab in labels:
        rp = regionprops((mask == lab).astype(np.uint8))
        if rp:
            cy, cx = rp[0].centroid        # row, col
            centroids[lab] = (cx, cy)

    cx_mid, cy_mid = W / 2.0, H / 2.0
    centre_label   = min(centroids, key=lambda k: math.hypot(
                         centroids[k][0] - cx_mid, centroids[k][1] - cy_mid))

    # ── Measurements per object ──
    kept_labels = []
    for lab in labels:
        bin_mask = (mask == lab).astype(np.uint8)
        # supply only first channel for intensity_image to keep scikit-image happy
        rp = regionprops(bin_mask, intensity_image=image[..., 0] if image.ndim==3 else image)[0]

        major = float(rp.major_axis_length)
        minor = float(rp.minor_axis_length) or np.nan
        lw_ratio = major / minor if not np.isnan(minor) else np.nan

        keep = True
        if LWR_FILTER and not np.isnan(lw_ratio):
            keep = lw_ratio > LWR_CUTOFF
        if keep:
            kept_labels.append(lab)

        # mean RGB (or grey) values
        if image.ndim == 2:
            mean_rgb = [float(rp.mean_intensity)] * 3
        else:
            pix = image[rp.coords[:, 0], rp.coords[:, 1]]
            mean_rgb = pix[:, :3].mean(axis=0)

        data_rows.append([
            base, int(lab),
            f"{int(rp.centroid[1])},{int(rp.centroid[0])}",
            int(rp.area),
            round(float(mean_rgb[0]), 2),
            round(float(mean_rgb[1]), 2),
            round(float(mean_rgb[2]), 2),
            round(major, 2),
            round(minor, 2),
            round(lw_ratio, 3) if not np.isnan(lw_ratio) else ""
        ])

    # ── Overlay with white outlines ──
    outline_mask = mask if not LWR_FILTER else np.where(
        np.isin(mask, kept_labels), mask, 0)

    outline = find_boundaries(outline_mask, mode="outer")

    overlay = gray2rgb(image) if image.ndim == 2 else np.copy(image)
    overlay_uint8 = img_as_ubyte(overlay / (overlay.max() or 1))
    overlay_uint8[outline] = [255, 255, 255]

    imsave(os.path.join(mask_dir, f"{base}_overlay.png"), overlay_uint8)

    # ── Centre-ROI grayscale image ──
    ch = 0  # default
    if image.ndim == 3:
        try:
            ch = int(channel_arg)
        except ValueError:
            ch = 0
        ch = max(0, min(ch, 2))  # clamp 0-2

    single_channel = extract_channel(image, ch)
    centre_img = np.zeros_like(single_channel)
    centre_img[mask == centre_label] = single_channel[mask == centre_label]

    imsave(os.path.join(
        sarcasm_dir, f"{base}_centreROI_ch{ch}.png"), centre_img)

# ────────────────────── Write CSV (renumber if filter on) ──────────────────────
if LWR_FILTER:
    from collections import defaultdict
    grouped = defaultdict(list)
    for row in data_rows:
        grouped[row[0]].append(row)
    renum = []
    for rows in grouped.values():
        for i, row in enumerate(rows, 1):
            row[1] = i
            renum.append(row)
    data_rows = renum

header = ['Image', 'Object Number', 'Object Location', 'Area',
          'Red Avg', 'Green Avg', 'Blue Avg',
          'Length', 'Width', 'L/W Ratio', prompt_log]

with open(output_csv, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(header)
    for r in data_rows:
        writer.writerow(r + [""])           # blank cell for prompt_log

print(f"[INFO] Saved {len(data_rows)} objects --> {output_csv}")
print(f"[INFO] Overlays  --> {mask_dir}")
print(f"[INFO] Sarcasm   --> {sarcasm_dir}")
