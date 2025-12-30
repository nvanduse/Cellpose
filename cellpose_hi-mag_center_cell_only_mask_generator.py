"""
cellpose_hi-mag_center_cell_only_mask_generator.py
---------------------------------------------------
Add this script to the working directory, containing input images. It will run Cellpose on every RGB image 
in the working directory and write the outputs shown below to a newly created "preliminary_masks" subfolder:
  preliminary_masks/<basename>_seg.npy
  preliminary_masks/<basename>_overlay.png

Usage:
  Put this script in a folder of images and run:
    python cellpose_hi-mag_center_cell_only_mask_generator.py
"""

import os
import sys
import numpy as np

from skimage.io import imread, imsave
from skimage.segmentation import find_boundaries
from skimage.color import gray2rgb
from skimage import img_as_ubyte
from skimage.measure import regionprops

# -------------------- User-tunable defaults --------------------
# Adjust these four values as needed for your dataset.
MODEL_TYPE = "cyto3"   # e.g., "cyto3", "cyto2", "nuclei"
DIAMETER = "auto"      # "auto" or a float (e.g., 25.0)
CHANNEL = 0            # 0=grayscale, 1=red, 2=green, 3=blue
USE_GPU = False        # True to use GPU (requires GPU-compatible Cellpose)

IMAGE_EXTS = (".png", ".tif", ".tiff", ".jpg", ".jpeg")


def list_images(folder):
    return sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith(IMAGE_EXTS)
    )


def load_cellpose_model():
    try:
        from cellpose import models
    except Exception as exc:
        print("[ERROR] Cellpose is not available. Install with: pip install cellpose")
        raise exc
    return models.Cellpose(model_type=MODEL_TYPE, gpu=USE_GPU)

def center_label_from_mask(mask):
    """Return the label whose centroid is closest to the image center."""
    labels = np.unique(mask)
    labels = labels[labels != 0]
    if labels.size == 0:
        return None
    H, W = mask.shape
    cx_mid, cy_mid = W / 2.0, H / 2.0
    centroids = {}
    for lab in labels:
        rp = regionprops((mask == lab).astype(np.uint8))
        if rp:
            cy, cx = rp[0].centroid
            centroids[lab] = (cx, cy)
    if not centroids:
        return None
    return min(
        centroids,
        key=lambda k: (centroids[k][0] - cx_mid) ** 2 + (centroids[k][1] - cy_mid) ** 2
    )

def main():
    image_dir = os.getcwd()
    output_dir = os.path.join(image_dir, "preliminary_masks")
    os.makedirs(output_dir, exist_ok=True)

    images = list_images(image_dir)
    if not images:
        print("[WARN] No image files found in the current folder.")
        return 0

    model = load_cellpose_model()
    diam = None if str(DIAMETER).lower() == "auto" else float(DIAMETER)
    channels = [int(CHANNEL), 0]

    for fname in images:
        base, _ = os.path.splitext(fname)
        in_path = os.path.join(image_dir, fname)
        image = imread(in_path)

        # Run Cellpose and keep the full label mask for selection.
        masks, flows, styles, diams = model.eval(
            image, diameter=diam, channels=channels
        )

        # Keep only the object closest to the image center.
        center_label = center_label_from_mask(masks)
        if center_label is None:
            print(f"[WARN] No objects found in {fname}; skipping.")
            continue

        # Renumber the kept object to label 1 for consistency.
        center_mask = np.where(masks == center_label, 1, 0).astype(masks.dtype)

        # Save Cellpose-style *_seg.npy with only the center object.
        seg_path = os.path.join(output_dir, f"{base}_seg.npy")
        np.save(seg_path, {
            "masks": center_mask,
            "model_type": MODEL_TYPE,
            "diameter": DIAMETER,
            "channels": channels,
        })

        outline = find_boundaries(center_mask, mode="outer")
        overlay = gray2rgb(image) if image.ndim == 2 else np.copy(image)
        overlay_uint8 = img_as_ubyte(overlay / (overlay.max() or 1))
        overlay_uint8[outline] = [255, 255, 255]

        overlay_path = os.path.join(output_dir, f"{base}_overlay.png")
        imsave(overlay_path, overlay_uint8)

        print(f"[OK] {fname} -> {os.path.basename(seg_path)}, {os.path.basename(overlay_path)}")

    print(f"[DONE] Wrote masks and overlays to: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
