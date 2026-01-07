"""
cellpose_segmentation_for_Nppa_confocal40x_1024x1024px.py

Segmentation-only pass for per-channel TIFFs using the Cellpose CLI.
Generates Cellpose masks for cells and nuclei so they can be edited in the GUI.

Usage:
python cellpose_segmentation_for_Nppa_confocal40x_1024x1024px.py <image_dir> <output_dir>
# Prompts for cell/nucleus channels (with tags) and optional diameter inputs.
# Use --cell_channel/--nuc_channel to skip channel prompts.
# Use --cell_diameter/--nuc_diameter to skip diameter prompts.
# If overlays are present, you can select the c1-4 overlay as a segmentation input.
# Optional overlay styling: --cell_outline_color R,G,B --nuc_outline_color R,G,B --outline_thickness N

Expected inputs (per sample):
  *_c1.tif, *_c2.tif, *_c3.tif, *_c4.tif  (or custom tags via --channel-tags)
  *_c1-4.tif (optional overlay, also usable as a segmentation input)

Outputs:
  output/
    masks_cell/<base>_cell_seg.npy
    masks_cell/<base>_cell.tif
    masks_nuc/<base>_nuc_seg.npy
    masks_nuc/<base>_nuc.tif
    overlays/<base>_overlay.png (2-color composite: cell channel + nuclei channel)
    segmentation_log.txt
"""

import os, sys, glob, warnings, shutil, subprocess, tempfile, datetime, shlex
import numpy as np
from skimage.io import imread, imsave
from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation, square
from skimage import img_as_ubyte

# -------------------- arg parsing --------------------
def parse_args():
    import argparse
    p = argparse.ArgumentParser(description="Segment cells + nuclei for manual GUI edits.")
    p.add_argument('image_dir', help="Directory with per-channel TIFFs and optional overlay TIFFs.")
    p.add_argument('output_dir')
    p.add_argument('--cell_channel', default=None,
                   help="Channel name for cells (far_red, red, green, blue, overlay) or leave blank to be prompted.")
    p.add_argument('--cell_diameter', default=None,
                   help="Cell diameter in pixels (number) or 'auto'. Leave blank to be prompted.")
    p.add_argument('--cell_model', default=None,
                   help="Cellpose model for cells (cyto3, cyto, nuclei). Leave blank to be prompted.")

    p.add_argument('--nuc_channel', default=None,
                   help="Channel name for nuclei (far_red, red, green, blue, overlay) or leave blank to be prompted.")
    p.add_argument('--nuc_diameter', default=None,
                   help="Nucleus diameter in pixels (number) or 'auto'. Leave blank to be prompted.")
    p.add_argument('--nuc_model', default=None,
                   help="Cellpose model for nuclei (nuclei, cyto3, cyto). Leave blank to be prompted.")

    p.add_argument('--channel-tags', nargs=4, default=["c1", "c2", "c3", "c4"],
                   help="Filename tags for far-red/red/green/blue channels (default: c1 c2 c3 c4).")
    p.add_argument('--save_overlays', type=int, default=1)
    p.add_argument('--cell_outline_color', default="255,255,255",
                   help="Cell outline color as R,G,B (default: 255,255,255).")
    p.add_argument('--nuc_outline_color', default="255,255,0",
                   help="Nucleus outline color as R,G,B (default: 255,255,0).")
    p.add_argument('--outline_thickness', type=int, default=1,
                   help="Outline thickness in pixels (default: 1).")
    return p.parse_args()

def ensure_dir(p):
    os.makedirs(p, exist_ok=True); return p

def load_image(path):
    return imread(path)

def run_cellpose_cli(image_dir: str, model_type: str, diameter: str) -> None:
    if shutil.which("cellpose") is None:
        raise RuntimeError("cellpose CLI not found on PATH. Install Cellpose or update PATH.")
    diam = "0" if str(diameter).lower() == "auto" else str(diameter)
    cmd = [
        "cellpose",
        "--dir", image_dir,
        "--pretrained_model", str(model_type),
        "--chan", "0",
        "--diameter", diam,
        "--verbose",
    ]
    subprocess.run(cmd, check=True)
    return cmd

def find_channel_file(base: str, image_dir: str, tag: str):
    patterns = [
        os.path.join(image_dir, f"{base}*_{tag}.tif"),
        os.path.join(image_dir, f"{base}*_{tag}.tiff"),
    ]
    matches = []
    for pat in patterns:
        matches.extend(sorted(glob.glob(pat)))
    return matches[0] if matches else None

def find_overlay_file(base: str, image_dir: str, overlay_tag: str):
    patterns = [
        os.path.join(image_dir, f"{base}*_{overlay_tag}.tif"),
        os.path.join(image_dir, f"{base}*_{overlay_tag}.tiff"),
    ]
    matches = []
    for pat in patterns:
        matches.extend(sorted(glob.glob(pat)))
    return matches[0] if matches else None

def load_mask(mask_path: str) -> np.ndarray:
    arr = np.load(mask_path, allow_pickle=True)
    try:
        obj = arr.item()
        if isinstance(obj, dict) and "masks" in obj:
            return np.asarray(obj["masks"])
    except Exception:
        pass
    return np.asarray(arr)

def collect_bases(image_dir: str, channel_tags, overlay_tag: str):
    bases = set()
    for tag in channel_tags:
        for ext in ("tif", "tiff"):
            for path in glob.glob(os.path.join(image_dir, f"*_{tag}.{ext}")):
                stem = os.path.splitext(os.path.basename(path))[0]
                suffix = f"_{tag}"
                if stem.endswith(suffix):
                    bases.add(stem[:-len(suffix)])
    for ext in ("tif", "tiff"):
        for path in glob.glob(os.path.join(image_dir, f"*_{overlay_tag}.{ext}")):
            stem = os.path.splitext(os.path.basename(path))[0]
            suffix = f"_{overlay_tag}"
            if stem.endswith(suffix):
                bases.add(stem[:-len(suffix)])
    return sorted(bases)

def prompt_channel(label: str, choices):
    valid = {c[0] for c in choices}
    while True:
        print(f"Select {label} channel:")
        for idx, (_, display) in enumerate(choices, start=1):
            print(f"  [{idx}] {display}")
        resp = input("Enter a channel number: ").strip().lower()
        if resp.isdigit():
            idx = int(resp)
            if 1 <= idx <= len(choices):
                return choices[idx - 1][0]
        if resp in valid:
            return resp
        print("Invalid channel. Try again.")

def prompt_diameter(label: str):
    while True:
        resp = input(f"Specify {label} diameter? [1] Auto [2] Manual: ").strip()
        if resp == "1":
            return "auto", "auto"
        if resp == "2":
            val = input(f"Enter expected {label} diameter (pixels): ").strip()
            if val:
                return val, "manual"
        print("Invalid choice. Try again.")

def prompt_model(label: str, choices):
    while True:
        print(f"Select {label} model:")
        for idx, name in enumerate(choices, start=1):
            print(f"  [{idx}] {name}")
        resp = input("Enter a model number: ").strip()
        if resp.isdigit():
            idx = int(resp)
            if 1 <= idx <= len(choices):
                return choices[idx - 1]
        print("Invalid choice. Try again.")

def normalize_to_uint8(img: np.ndarray) -> np.ndarray:
    arr = np.asarray(img)
    if arr.ndim == 3:
        arr = arr.mean(axis=2)
    arr = arr.astype(np.float32)
    mx = float(np.max(arr)) if arr.size else 0.0
    if mx <= 0:
        return np.zeros(arr.shape, dtype=np.uint8)
    return img_as_ubyte(arr / mx)

def parse_rgb(value: str, default):
    try:
        parts = [int(p) for p in value.split(",")]
        if len(parts) != 3:
            return default
        return tuple(max(0, min(255, p)) for p in parts)
    except Exception:
        return default

def thicken_outline(mask: np.ndarray, thickness: int) -> np.ndarray:
    if thickness <= 1:
        return mask
    return binary_dilation(mask, square(1 + 2 * (thickness - 1)))

# -------------------- main --------------------
def main():
    args = parse_args()

    outdir = ensure_dir(args.output_dir)
    cell_dir = ensure_dir(os.path.join(outdir, "masks_cell"))
    nuc_dir  = ensure_dir(os.path.join(outdir, "masks_nuc"))
    ovl_dir  = ensure_dir(os.path.join(outdir, "overlays"))

    channel_order = ["far_red", "red", "green", "blue"]
    channel_tag_map = dict(zip(channel_order, args.channel_tags))
    overlay_tag = "c1-4"

    bases = collect_bases(args.image_dir, args.channel_tags, overlay_tag)
    if not bases:
        print("[WARN] No input TIFFs found. Check --channel-tags and image_dir.")
        return
    overlay_available = any(find_overlay_file(base, args.image_dir, overlay_tag) for base in bases)
    choices = [
        ("far_red", f"{channel_tag_map['far_red']}_far_red"),
        ("red", f"{channel_tag_map['red']}_red"),
        ("green", f"{channel_tag_map['green']}_green"),
        ("blue", f"{channel_tag_map['blue']}_blue"),
    ]
    if overlay_available:
        choices.append(("overlay", f"{overlay_tag}_all_grayscale"))

    valid_channels = {c[0] for c in choices}
    cell_channel = (args.cell_channel or "").strip().lower()
    nuc_channel = (args.nuc_channel or "").strip().lower()
    if cell_channel not in valid_channels:
        cell_channel = prompt_channel("cell", choices)
    if nuc_channel not in valid_channels:
        nuc_channel = prompt_channel("nucleus", choices)

    cell_diameter = args.cell_diameter
    nuc_diameter = args.nuc_diameter
    if cell_diameter is None:
        cell_diameter, cell_diam_mode = prompt_diameter("cell")
    else:
        cell_diam_mode = "auto" if str(cell_diameter).lower() in {"auto", "0"} else "manual"
    if nuc_diameter is None:
        nuc_diameter, nuc_diam_mode = prompt_diameter("nucleus")
    else:
        nuc_diam_mode = "auto" if str(nuc_diameter).lower() in {"auto", "0"} else "manual"

    if args.cell_model is None:
        cell_model = prompt_model("cell", ["cyto3", "cyto", "nuclei"])
    else:
        cell_model = args.cell_model
    if args.nuc_model is None:
        nuc_model = prompt_model("nucleus", ["nuclei", "cyto3", "cyto"])
    else:
        nuc_model = args.nuc_model

    log_path = os.path.join(outdir, "segmentation_log.txt")
    log_lines = []
    log_lines.append(f"===== Cellpose segmentation log: {datetime.datetime.now().isoformat()} =====")
    log_lines.append(f"Command: {shlex.join(sys.argv)}")
    log_lines.append(f"Image dir: {args.image_dir}")
    log_lines.append(f"Output dir: {outdir}")
    log_lines.append(f"Channel tags: far_red={channel_tag_map['far_red']}, red={channel_tag_map['red']}, green={channel_tag_map['green']}, blue={channel_tag_map['blue']}")
    log_lines.append(f"Overlay tag: {overlay_tag} (available={overlay_available})")
    log_lines.append(f"Cell channel: {cell_channel}")
    log_lines.append(f"Nuc channel: {nuc_channel}")
    log_lines.append(f"Cell model: {cell_model}")
    log_lines.append(f"Nuc model: {nuc_model}")
    log_lines.append(f"Cell diameter: {cell_diameter} (mode={cell_diam_mode})")
    log_lines.append(f"Nuc diameter: {nuc_diameter} (mode={nuc_diam_mode})")
    log_lines.append(f"Bases detected: {len(bases)}")
    log_lines.append("Mask/image pairing: <base>_cell_seg.npy with <base>_cell.tif and <base>_nuc_seg.npy with <base>_nuc.tif")
    cell_outline_color = parse_rgb(args.cell_outline_color, (255, 255, 255))
    nuc_outline_color = parse_rgb(args.nuc_outline_color, (255, 255, 0))
    outline_thickness = max(1, int(args.outline_thickness))
    log_lines.append(f"Cell outline color: {cell_outline_color}")
    log_lines.append(f"Nuc outline color: {nuc_outline_color}")
    log_lines.append(f"Outline thickness: {outline_thickness}")

    def stage_channel_images(selection: str, suffix: str, temp_dir: str):
        staged = {}
        for base in bases:
            if selection == "overlay":
                src = find_overlay_file(base, args.image_dir, overlay_tag)
            else:
                tag = channel_tag_map[selection]
                src = find_channel_file(base, args.image_dir, tag)
            if not src:
                warnings.warn(f"Missing channel image for {base}; skipping.")
                continue
            ext = os.path.splitext(src)[1]
            dst = os.path.join(temp_dir, f"{base}__{suffix}{ext}")
            shutil.copy2(src, dst)
            staged[base] = {"src": src, "staged": dst, "ext": ext}
        return staged

    def move_masks(staged_map, suffix: str, out_dir: str):
        for base, info in staged_map.items():
            src_mask = os.path.join(temp_dir, f"{base}__{suffix}_seg.npy")
            if not os.path.exists(src_mask):
                warnings.warn(f"Missing mask for {base} ({suffix}); skipping.")
                continue
            dst_mask = os.path.join(out_dir, f"{base}_{suffix}_seg.npy")
            shutil.move(src_mask, dst_mask)
            dst_img = os.path.join(out_dir, f"{base}_{suffix}{info['ext']}")
            shutil.copy2(info["src"], dst_img)

    with tempfile.TemporaryDirectory() as temp_dir:
        print("[INFO] Staging cell channel images...")
        staged_cells = stage_channel_images(cell_channel, "cell", temp_dir)
        if staged_cells:
            print("[INFO] Running Cellpose CLI for cells...")
            cell_cmd = run_cellpose_cli(temp_dir, cell_model, cell_diameter)
            move_masks(staged_cells, "cell", cell_dir)
        else:
            print("[WARN] No cell channel images staged.")
            cell_cmd = None

        for fname in os.listdir(temp_dir):
            if fname.endswith((".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp")):
                os.remove(os.path.join(temp_dir, fname))
            if fname.endswith("_seg.npy"):
                os.remove(os.path.join(temp_dir, fname))

        print("[INFO] Staging nucleus channel images...")
        staged_nucs = stage_channel_images(nuc_channel, "nuc", temp_dir)
        if staged_nucs:
            print("[INFO] Running Cellpose CLI for nuclei...")
            nuc_cmd = run_cellpose_cli(temp_dir, nuc_model, nuc_diameter)
            move_masks(staged_nucs, "nuc", nuc_dir)
        else:
            print("[WARN] No nucleus channel images staged.")
            nuc_cmd = None

        if cell_cmd:
            log_lines.append(f"Cellpose cell command: {shlex.join(cell_cmd)}")
        if nuc_cmd:
            log_lines.append(f"Cellpose nuc command: {shlex.join(nuc_cmd)}")
        log_lines.append(f"Cell staged: {len(staged_cells)}")
        log_lines.append(f"Nuc staged: {len(staged_nucs)}")

        # overlays (2-color composite)
        if args.save_overlays:
            for base in bases:
                cell_info = staged_cells.get(base)
                nuc_info = staged_nucs.get(base)
                if not (cell_info and nuc_info):
                    continue
                cell_mask_path = os.path.join(cell_dir, f"{base}_cell_seg.npy")
                nuc_mask_path = os.path.join(nuc_dir, f"{base}_nuc_seg.npy")
                if not (os.path.exists(cell_mask_path) and os.path.exists(nuc_mask_path)):
                    continue
                cell_masks = load_mask(cell_mask_path)
                nuc_masks = load_mask(nuc_mask_path)
                cell_img = normalize_to_uint8(load_image(cell_info["src"]))
                nuc_img = normalize_to_uint8(load_image(nuc_info["src"]))

                rgb = np.zeros((cell_img.shape[0], cell_img.shape[1], 3), dtype=np.uint8)
                if cell_channel == "green":
                    rgb[..., 1] = cell_img
                else:
                    rgb[..., 0] = cell_img
                rgb[..., 2] = nuc_img
                cell_outline = thicken_outline(find_boundaries(cell_masks, mode="outer"), outline_thickness)
                nuc_outline = thicken_outline(find_boundaries(nuc_masks, mode="outer"), outline_thickness)
                rgb[cell_outline] = cell_outline_color
                rgb[nuc_outline] = nuc_outline_color
                imsave(os.path.join(ovl_dir, f"{base}_overlay.png"), rgb)

    next_steps = [
        "===== Next Steps =====",
        "Open the Cellpose GUI (enter 'cellpose' in the command line), and inspect the masks.",
        "Just drag and drop an image into the GUI and the masks should automatically load if the matching *seg.npy mask file is in the same folder as the image.",
        "Manually correct or remove bad masks if necessary (see mask editing controls below).",
        "After correcting the masks, then run the quantification script.",
        "",
        "CELLPOSE GUI SHORTCUTS:",
        "",
        "Select mask: Left-click on a mask.",
        "Delete mask: Ctrl + left-click on the mask.",
        "Start drawing new mask: Right-click to start a mask.",
        "Finish drawing mask: Right-click again or return cursor to the red start circle (depending on if 'single stroke' box is checked).",
        "Merge masks: Alt + left-click (merges the last two masks clicked).",
        "",
        "Undo last mask/stroke: Ctrl+Z.",
        "Clear all masks in the image: Ctrl+0.",
        "Toggle mask visibility: X (turn masks on/off).",
        "Toggle outlines: Z (show/hide outlines).",
        "Change brush size while drawing: , and . to increase/decrease brush size.",
        "",
        "Save current masks to _seg.npy: Ctrl+S.",
        "Save masks as PNG: Ctrl+N.",
        "Save ROIs to ImageJ ROI format: Ctrl+R.",
    ]
    with open(log_path, "w", encoding="ascii") as fh:
        fh.write("\n".join(log_lines) + "\n\n")
        fh.write("\n".join(next_steps) + "\n")
    print(f"[DONE] Segmentation masks saved to {outdir}")

if __name__ == "__main__":
    main()
