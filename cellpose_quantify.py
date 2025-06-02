# cellpose_quantify.py
import sys
import os
import csv
import numpy as np
from skimage.io import imread, imsave
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries
from skimage.color import gray2rgb

# Input args
image_dir = sys.argv[1]
mask_dir = sys.argv[2]
output_csv = os.path.join(image_dir, "object_data.csv")

# CSV data list
data = []

# Loop through .npy mask files
for fname in os.listdir(mask_dir):
    if not fname.endswith('_seg.npy'):
        continue

    base = fname.replace('_seg.npy', '')
    mask_path = os.path.join(mask_dir, fname)

    # Look for corresponding image
    image_path_candidates = [
        os.path.join(image_dir, base + ext)
        for ext in ['.png', '.tif', '.jpg']
    ]
    image_path = next((p for p in image_path_candidates if os.path.exists(p)), None)

    if not image_path:
        print(f"No image found for mask {fname}")
        continue

    # Load image and mask
    mask_data = np.load(mask_path, allow_pickle=True).item()
    mask = np.array(mask_data['masks'])  # Ensure it's a proper ndarray
    image = imread(image_path)

    labels = np.unique(mask)
    labels = labels[labels != 0]  # exclude background

    for label in labels:
        binary = (mask == label).astype(np.uint8)
        props = regionprops(binary, intensity_image=image)
        if not props:
            continue
        p = props[0]
        centroid = p.centroid
        area = p.area
        mean = p.mean_intensity
        major = p.major_axis_length
        minor = p.minor_axis_length
        data.append([
            base,
            label,
            f"{int(centroid[1])},{int(centroid[0])}",
            area,
            round(mean[0], 2),
            round(mean[1], 2),
            round(mean[2], 2),
            round(major, 2),
            round(minor, 2)
        ])

    from skimage import img_as_ubyte

    # Generate overlay
    outline = find_boundaries(mask, mode='outer')
    overlay = np.copy(image)

    # Ensure overlay is RGB
    if overlay.ndim == 2 or overlay.shape[2] == 1:
        overlay = gray2rgb(overlay)

    # Normalize overlay to uint8 (0-255)
    overlay_uint8 = img_as_ubyte(overlay / np.max(overlay))

    # Draw white outlines 
    overlay_uint8[outline] = [255, 255, 255]

    overlay_path = os.path.join(mask_dir, f"{base}_overlay.png")
    imsave(overlay_path, overlay_uint8)

# Save CSV
with open(output_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Image', 'Object Number', 'Object Location', 'Area',
                     'Red Avg', 'Green Avg', 'Blue Avg', 'Length', 'Width'])
    writer.writerows(data)

print(f"Saved results to {output_csv}")
