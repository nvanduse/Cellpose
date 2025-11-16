"""
cellpose_quantify_Nppa_v1.py

Measures cytoplasm via two-pass segmentation (cells, nuclei). For each cell:
- Build cytoplasm mask = cell_mask - nuclei_mask
- Report intensity metrics on cytoplasm only
- Replace "Area" with Cytoplasm Area, plus add Cell Area, Nucleus Area(in cell), Cytoplasm Fraction

Outputs:
  output/
    masks_cell/<image>_cell_seg.npy
    masks_nuc/<image>_nuc_seg.npy
    overlays/<image>_overlay.png
    object_data.csv        (per-cell table; cytoplasm-based intensities)
"""

import os, sys, csv, math
import numpy as np
from skimage.io import imread, imsave
from skimage.segmentation import find_boundaries
from skimage.color import gray2rgb
from skimage import img_as_ubyte
from skimage.measure import regionprops

# -------------------- arg parsing --------------------
def parse_args():
    import argparse
    p = argparse.ArgumentParser(description="Cells + nuclei segmentation; cytoplasm-only intensity metrics.")
    p.add_argument('image_dir')
    p.add_argument('output_dir')
    p.add_argument('--cell_channel', default=None)
    p.add_argument('--cell_diameter', default='auto')
    p.add_argument('--cell_model', default='cyto3')

    p.add_argument('--nuc_channel', default=None)
    p.add_argument('--nuc_diameter', default='auto')
    p.add_argument('--nuc_model', default='nuclei')

    p.add_argument('--lwr_filter', type=int, default=0)
    p.add_argument('--lwr_cutoff', type=float, default=0.0)
    p.add_argument('--save_overlays', type=int, default=1)
    return p.parse_args()

def channel_to_index(ch):
    if ch is None or ch == "": return None
    c = str(ch).strip().lower()
    m = {'r':0,'g':1,'b':2}
    if c in m: return m[c]
    try:
        i = int(c)
        return i
    except: return None

def ensure_dir(p):
    os.makedirs(p, exist_ok=True); return p

def load_image(path):
    img = imread(path)
    return img

def select_channel(img, idx):
    if img.ndim == 2 or idx is None: return img
    if idx < 0 or idx >= img.shape[2]:
        raise ValueError(f"channel {idx} out of range for {img.shape}")
    return img[..., idx]

def run_cellpose_single(img2d, model_type="cyto3", diameter='auto'):
    from cellpose import models
    model = models.Cellpose(model_type=model_type, gpu=False)
    diam = None if str(diameter).lower() == 'auto' else float(diameter)
    masks, flows, styles, diams = model.eval(img2d, diameter=diam, channels=[0,0])
    return masks

def p99(v):
    if v.size == 0: return np.nan
    return float(np.percentile(v, 99))

def compute_rgb_means_p99(full_img, coords):
    """Return (mean_r, mean_g, mean_b, p99_r, p99_g, p99_b). If grayscale, repeats the single value."""
    if coords.size == 0:
        return (np.nan,)*6
    if full_img.ndim == 2:
        pix = full_img[coords[:,0], coords[:,1]].astype(float)
        m = float(np.mean(pix)); q = p99(pix)
        return (m,m,m,q,q,q)
    else:
        C = min(3, full_img.shape[2])
        ch = full_img[coords[:,0], coords[:,1], :C].astype(float)
        means = np.mean(ch, axis=0)
        q = np.percentile(ch, 99, axis=0)
        # pad to length 3 if needed
        if C < 3:
            means = np.pad(means, (0,3-C), constant_values=np.nan)
            q     = np.pad(q,     (0,3-C), constant_values=np.nan)
        return (float(means[0]), float(means[1]), float(means[2]),
                float(q[0]), float(q[1]), float(q[2]))

def safe_lw(props):
    maj = float(props.major_axis_length)
    minr = float(props.minor_axis_length) if props.minor_axis_length != 0 else np.nan
    lw = maj/minr if (not np.isnan(minr) and minr != 0) else np.nan
    return maj, minr, lw

# -------------------- main --------------------
def main():
    args = parse_args()

    cell_ch = channel_to_index(args.cell_channel)
    nuc_ch  = channel_to_index(args.nuc_channel)

    outdir = ensure_dir(args.output_dir)
    cell_dir = ensure_dir(os.path.join(outdir, "masks_cell"))
    nuc_dir  = ensure_dir(os.path.join(outdir, "masks_nuc"))
    ovl_dir  = ensure_dir(os.path.join(outdir, "overlays"))

    prompt_log = (f"Prompt_Log: cell_ch={args.cell_channel}; cell_diam={args.cell_diameter}; "
                  f"cell_model={args.cell_model}; nuc_ch={args.nuc_channel}; "
                  f"nuc_diam={args.nuc_diameter}; nuc_model={args.nuc_model}; "
                  f"L/W_filter={args.lwr_filter}; cutoff={args.lwr_cutoff}")

    images = sorted([f for f in os.listdir(args.image_dir)
                     if f.lower().endswith((".png",".tif",".tiff",".jpg",".jpeg",".bmp"))])

    rows = []
    for fname in images:
        path = os.path.join(args.image_dir, fname)
        base = os.path.splitext(fname)[0]
        print(f"[INFO] {fname}")

        img = load_image(path)
        cell_img2d = select_channel(img, cell_ch)
        nuc_img2d  = select_channel(img, nuc_ch)

        # --- segment ---
        cell_masks = run_cellpose_single(cell_img2d, model_type=args.cell_model, diameter=args.cell_diameter)
        nuc_masks  = run_cellpose_single(nuc_img2d,  model_type=args.nuc_model,  diameter=args.nuc_diameter)

        # save npy bundles
        np.save(os.path.join(cell_dir, f"{base}_cell_seg.npy"), {"masks":cell_masks})
        np.save(os.path.join(nuc_dir,  f"{base}_nuc_seg.npy"),  {"masks":nuc_masks})

        # label lists
        cell_labels = np.unique(cell_masks); cell_labels = cell_labels[cell_labels!=0]

        # optional L/W filter cache
        cell_prop_cache = {}

        for lbl in cell_labels:
            bin_cell = (cell_masks==lbl).astype(np.uint8)
            props = regionprops(bin_cell, intensity_image=cell_img2d)
            if not props: continue
            p = props[0]
            cell_prop_cache[lbl] = p

            maj, minr, lw = safe_lw(p)
            keep = True
            if args.lwr_filter and args.lwr_cutoff>0 and not np.isnan(lw):
                keep = (lw > args.lwr_cutoff)
            if not keep:
                continue

            # areas
            cell_area = int(p.area)

            # cytoplasm = cell minus ANY nucleus pixels
            nuc_any = (nuc_masks>0)
            cyto_bin = (bin_cell==1) & (~nuc_any)
            cyto_area = int(np.count_nonzero(cyto_bin))

            # nucleus area within this cell (for bookkeeping)
            nuc_in_cell = int(np.count_nonzero((bin_cell==1) & nuc_any))

            # cytoplasm coords
            cyto_coords = np.column_stack(np.nonzero(cyto_bin))

            # intensities (computed on cytoplasm only)
            r_avg, g_avg, b_avg, r_p99, g_p99, b_p99 = compute_rgb_means_p99(img, cyto_coords)

            cyto_frac = (cyto_area/float(cell_area)) if cell_area>0 else np.nan
            cy = int(p.centroid[0]); cx = int(p.centroid[1])

            # NOTE: "Area" column below is CYTOPLASM area (as requested)
            rows.append([
                base, int(lbl), f"{cx},{cy}",
                int(cyto_area),
                round(r_avg if not np.isnan(r_avg) else 0.0, 2),
                round(g_avg if not np.isnan(g_avg) else 0.0, 2),
                round(b_avg if not np.isnan(b_avg) else 0.0, 2),
                round(r_p99 if not np.isnan(r_p99) else 0.0, 2),
                round(g_p99 if not np.isnan(g_p99) else 0.0, 2),
                round(b_p99 if not np.isnan(b_p99) else 0.0, 2),
                round(maj,2), round(minr if not np.isnan(minr) else 0.0,2),
                round(lw if not np.isnan(lw) else 0.0,3),

                # new bookkeeping columns
                int(cell_area),
                int(nuc_in_cell),
                round(cyto_frac if not np.isnan(cyto_frac) else 0.0,4),

                prompt_log
            ])

        # overlays (white = cell outlines; red = nuclei outlines)
        if args.save_overlays:
            ov = img
            if ov.ndim==2 or (ov.ndim==3 and ov.shape[2]==1):
                ov = gray2rgb(ov)
            mx = float(np.max(ov)) if ov.size else 1.0
            ov8 = img_as_ubyte(ov/mx) if mx>0 else np.zeros_like(ov, dtype=np.uint8)

            co = find_boundaries(cell_masks, mode="outer")
            no = find_boundaries(nuc_masks,  mode="outer")
            ov8[co] = [255,255,255]
            ov8[no] = [255, 64, 64]
            imsave(os.path.join(ovl_dir, f"{base}_overlay.png"), ov8)

    # write CSV
    out_csv = os.path.join(outdir, "object_data.csv")
    header = [
        'Image','Object Number','Object Location','Area',           # Area = CYTOPLASM area
        'Red Avg','Green Avg','Blue Avg','Red p99','Green p99','Blue p99',
        'Length','Width','L/W Ratio',
        'Cell Area (total)','Nucleus Area (in cell)','Cytoplasm Fraction',
        'Prompt_Log'
    ]
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

    print(f"[DONE] Per-cell (cytoplasm) table -> {out_csv}")

if __name__ == "__main__":
    main()
