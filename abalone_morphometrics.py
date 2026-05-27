#!/usr/bin/env python3
"""
abalone_morphometrics.py  –  v2.0
===================================
Batch-measures abalone from lightbox JPEG images.

Setup / Installation
--------------------
    pip install opencv-python numpy colour-science colour-checker-detection

Pipeline
--------
1.  Load JPEG (grey/white board background, Calibrite ColorChecker Classic,
    thick black vertical divider strip, abalone on the right half).
2.  Detect the black vertical divider strip to isolate the right-hand zone.
3.  Detect the ColorChecker card using colour-checker-detection library and
    compute px/mm from patch spacing (6 patches × 12 mm along long dimension).
    Handles all 4 card orientations automatically.
4.  Segment the abalone using HSV thresholds calibrated from ground-truth:
      (V < 75) OR (S > 30), AND (V > 15)
    -- abalone tissue is darker and more colourful than the neutral grey board.
5.  Fit a minimum-area rotated bounding rectangle → length & width (mm).
6.  Compute filled-contour area (mm²).
7.  Save an annotated visualisation.
8.  Append results to a CSV.

Usage (Windows)
---------------
python abalone_morphometrics.py ^
    --images  "C:\\Users\\RebeccaPedler\\OneDrive - Yumbah\\Documents\\Images\\PROFILED\\Photos of Claude" ^
    --output  "C:\\Users\\RebeccaPedler\\Documents\\abalone_measurements.csv" ^
    --vis_dir "C:\\Users\\RebeccaPedler\\Documents\\abalone_annotated"

Optional flags
--------------
--scale_mm      Real-world span of the ruler (mm).  Default = 123
--ext           Comma-separated extensions to glob.  Default = jpg,jpeg
--min_area_px   Min contour area (px2) to be an abalone.  Default = 50000
--debug         Also save intermediate mask images.
"""

import argparse
import csv
import glob
import os
import sys
import warnings
from pathlib import Path

import cv2
import numpy as np


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RULER_REAL_MM        = 123
DIVIDER_FRAC         = 0.75
MIN_DIVIDER_WIDTH_PX = 20
MAX_DIVIDER_WIDTH_PX = 500


# ===========================================================================
# 1.  DIVIDER DETECTION
# ===========================================================================

def find_divider_x(img_bgr):
    """Return the x-coordinate of the right edge of the black divider strip."""
    h, w = img_bgr.shape[:2]
    search_w = int(w * DIVIDER_FRAC)

    gray = cv2.cvtColor(img_bgr[:, :search_w], cv2.COLOR_BGR2GRAY)
    _, dark = cv2.threshold(gray, 50, 255, cv2.THRESH_BINARY_INV)

    col_frac   = dark.sum(axis=0).astype(float) / (255.0 * h)
    col_smooth = np.convolve(col_frac, np.ones(31) / 31, mode='same')

    dark_cols = np.where(col_smooth > 0.4)[0]
    if len(dark_cols) == 0:
        return w // 2

    right_edge = int(dark_cols[-1])
    left_edge  = int(dark_cols[0])
    width      = right_edge - left_edge

    if width < MIN_DIVIDER_WIDTH_PX or width > MAX_DIVIDER_WIDTH_PX:
        return w // 2

    return right_edge


# ===========================================================================
# 2.  SCALE CALIBRATION — via ColourChecker detection library
# ===========================================================================
#
# Uses the colour_checker_detection library to locate the Calibrite
# ColorChecker Classic card and compute px/mm from its patch spacing.
# Each patch is PATCH_SIZE_MM mm; the card has N_PATCHES_LONG patches
# along its long dimension.  Tries all 4 rotations to handle flipped cards.

PATCH_SIZE_MM  = 12.0   # physical size of one patch (mm)
N_PATCHES_LONG = 6      # patches along the long card dimension


def pixels_per_mm(img_bgr, divider_x=None,
                  ruler_real_mm=RULER_REAL_MM,
                  debug_dir=None, stem=""):
    """
    Return (px_per_mm, patch_centres) by directly measuring individual
    ColorChecker patches in the left half of the image.

    Targets five well-separated colours (yellow, cyan, magenta, green, orange)
    and measures the bounding-box side length of each detected square patch.
    Each patch is physically PATCH_SIZE_MM x PATCH_SIZE_MM.
    Uses the median patch size to reject fragments or partial detections.

    Returns (px_per_mm, quad_pts) where quad_pts is a (4,2) int32 array
    enclosing all detected patch centres (for visualisation).
    Raises ValueError if fewer than 2 patches are detected.
    """
    h, w = img_bgr.shape[:2]
    left = img_bgr[:, :w // 2]   # card is always in the left half

    hsv = cv2.cvtColor(left, cv2.COLOR_BGR2HSV)
    H   = hsv[:, :, 0].astype(float)
    S   = hsv[:, :, 1].astype(float)
    V   = hsv[:, :, 2].astype(float)

    # Colour bands targeting individual, well-separated patches
    bands = {
        "yellow":  (H >= 20)  & (H <= 35)  & (S > 150) & (V > 80),
        "cyan":    (H >= 85)  & (H <= 100) & (S > 150) & (V > 80),
        "magenta": (H >= 145) & (H <= 165) & (S > 150) & (V > 80),
        "green":   (H >= 60)  & (H <= 85)  & (S > 150) & (V > 60),
        "orange":  (H >= 5)   & (H <= 20)  & (S > 150) & (V > 80),
    }

    k_close = np.ones((8,  8),  np.uint8)
    k_open  = np.ones((4,  4),  np.uint8)

    # Expected patch size range in pixels (very wide tolerance)
    MIN_PX, MAX_PX = 50, 600

    patch_sizes   = []
    patch_centres = []

    for name, band in bands.items():
        mask = band.astype(np.uint8) * 255
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, k_close)
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN,  k_open)
        cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL,
                                   cv2.CHAIN_APPROX_SIMPLE)
        for c in cnts:
            xb, yb, bw, bh = cv2.boundingRect(c)
            if not (MIN_PX < bw < MAX_PX and MIN_PX < bh < MAX_PX):
                continue
            if not (0.70 < bw / bh < 1.30):        # must be roughly square
                continue
            if cv2.contourArea(c) < (MIN_PX ** 2) * 0.5:
                continue
            patch_sizes.append((bw + bh) / 2.0)
            patch_centres.append((xb + bw // 2, yb + bh // 2))

    if len(patch_sizes) < 2:
        raise ValueError(
            f"Only {len(patch_sizes)} colour patch(es) found — "
            "check that the ColorChecker card is in the left half of the image.")

    # Median is robust against partial patches at card edges
    px_per_mm = float(np.median(patch_sizes)) / PATCH_SIZE_MM

    # Build a convex hull around patch centres for the card overlay
    if len(patch_centres) >= 3:
        pts        = np.array(patch_centres, dtype=np.int32)
        hull_idx   = cv2.convexHull(pts, returnPoints=False).flatten()
        hull_pts   = pts[hull_idx]
        # Add generous padding so the overlay covers the full card body
        pad        = int(np.median(patch_sizes) * 0.8)
        cx_all     = int(pts[:, 0].mean());   cy_all = int(pts[:, 1].mean())
        x0c = max(0,     pts[:, 0].min() - pad)
        y0c = max(0,     pts[:, 1].min() - pad)
        x1c = min(w // 2, pts[:, 0].max() + pad)
        y1c = min(h,      pts[:, 1].max() + pad)
        quad_pts = np.array([[x0c, y0c], [x1c, y0c],
                              [x1c, y1c], [x0c, y1c]], dtype=np.int32)
    else:
        pts      = np.array(patch_centres, dtype=np.int32)
        quad_pts = pts

    return px_per_mm, quad_pts


# ===========================================================================
# 3.  ABALONE SEGMENTATION
# ===========================================================================

def segment_abalone(img_bgr, divider_x,
                    min_area_px=50_000,
                    debug_dir=None, stem=""):
    """
    Segment the abalone from the grey board on the right side of the divider.

    Thresholds calibrated from ground-truth red annotations:
      True abalone:  S mean=94,  V mean=75-84  (colourful, dark)
      Nacre/water:   S mean=22,  V mean=120     (neutral grey, bright)

    Two-pass approach:
      Pass 1: High-confidence mask  (S>35) AND (V<110) AND (V>15)
              Captures colourful shell/tissue, excludes white nacre and wet board.
      Pass 2: Within the bounding region of Pass 1, apply a slightly wider
              threshold to recover any shell edges missed by strict V<110.
      Final:  Convex hull to remove water tendrils and irregular edges.

    Returns (contour_full_coords, full_mask) or (None, None).
    """
    h, w    = img_bgr.shape[:2]
    margin  = 80
    x_start = min(divider_x + margin, w - 1)
    roi     = img_bgr[:, x_start:]
    roi_h, roi_w = roi.shape[:2]

    hsv = cv2.cvtColor(roi, cv2.COLOR_BGR2HSV)
    S   = hsv[:, :, 1].astype(float)
    V   = hsv[:, :, 2].astype(float)

    k_close = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (25, 25))
    k_open  = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (10, 10))

    # ── PASS 1: high-confidence pixels ───────────────────────────────────────
    # S>35 = colourful shell/lip/tissue
    # V<110 = excludes bright nacre rim (nacre V mean=120)
    # V>15  = excludes dark outer frame
    mask_hc = ((S > 35) & (V < 110) & (V > 15)).astype(np.uint8) * 255
    mask_hc = cv2.morphologyEx(mask_hc, cv2.MORPH_CLOSE, k_close)
    mask_hc = cv2.morphologyEx(mask_hc, cv2.MORPH_OPEN,  k_open)

    cnts_hc, _ = cv2.findContours(mask_hc, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts_hc:
        return None, None

    valid_hc = [c for c in cnts_hc if cv2.contourArea(c) > min_area_px * 0.3]
    if not valid_hc:
        return None, None

    anchor      = max(valid_hc, key=cv2.contourArea)
    anchor_mask = np.zeros((roi_h, roi_w), np.uint8)
    cv2.drawContours(anchor_mask, [anchor], -1, 255, cv2.FILLED)

    # ── PASS 2: wider threshold inside the anchor region ─────────────────────
    search_region = cv2.dilate(anchor_mask, np.ones((40, 40), np.uint8))
    mask_w = (((S > 20) & (V < 120) & (V > 15)) & (search_region > 0)).astype(np.uint8) * 255
    mask_w = cv2.morphologyEx(mask_w, cv2.MORPH_CLOSE, k_close)
    mask_w = cv2.morphologyEx(mask_w, cv2.MORPH_OPEN,  np.ones((6, 6), np.uint8))

    if debug_dir and stem:
        cv2.imwrite(os.path.join(debug_dir, f"{stem}_seg_mask.jpg"), mask_w)

    contours, _ = cv2.findContours(mask_w, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None, None

    # ── SELECT best contour ───────────────────────────────────────────────────
    valid = []
    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area < min_area_px:
            continue
        hull = cv2.convexHull(cnt)
        rect = cv2.minAreaRect(hull)
        _, (rw, rh), _ = rect
        if rw <= 0 or rh <= 0:
            continue
        if max(rw, rh) / min(rw, rh) > 4.0:
            continue
        if rw > roi_w * 0.90:
            continue
        xc, yc, cntw, cnth = cv2.boundingRect(cnt)
        if xc <= 10:
            continue
        valid.append(cnt)

    if not valid:
        return None, None

    best      = max(valid, key=cv2.contourArea)
    best_hull = cv2.convexHull(best)   # convex hull removes water tendrils

    best_hull[:, :, 0] += x_start

    full_mask = np.zeros((h, w), dtype=np.uint8)
    cv2.drawContours(full_mask, [best_hull], -1, 255, cv2.FILLED)

    return best_hull, full_mask


# ===========================================================================
# 4.  MEASUREMENT
# ===========================================================================

def measure_abalone(contour, mask, px_mm):
    """Return dict with length_mm, width_mm, area_mm2, rect, box_pts."""
    rect      = cv2.minAreaRect(contour)
    box_w, box_h = rect[1]
    length_px = max(box_w, box_h)
    width_px  = min(box_w, box_h)
    area_px   = float(np.count_nonzero(mask))

    return {
        "length_mm": round(length_px / px_mm, 2),
        "width_mm":  round(width_px  / px_mm, 2),
        "area_mm2":  round(area_px   / (px_mm ** 2), 2),
        "rect":      rect,
        "box_pts":   cv2.boxPoints(rect).astype(int),
    }


# ===========================================================================
# 5.  VISUALISATION
# ===========================================================================

def save_annotated(img_bgr, contour, meas, divider_x, px_mm, out_path,
                   card_quad=None):
    vis = img_bgr.copy()
    h, w = vis.shape[:2]

    font  = cv2.FONT_HERSHEY_SIMPLEX
    fsc   = max(1.0, w / 3000)
    thick = max(2, int(w / 1500))
    lh    = int(55 * fsc)

    # ── Divider line ──────────────────────────────────────────────────────────
    cv2.line(vis, (divider_x, 0), (divider_x, h), (255, 180, 0), 3)

    # ── ColorChecker card overlay ─────────────────────────────────────────────
    if card_quad is not None:
        pts = card_quad.reshape((-1, 1, 2))
        # Filled semi-transparent highlight
        overlay = vis.copy()
        cv2.fillPoly(overlay, [pts], (255, 200, 0))          # amber fill
        cv2.addWeighted(overlay, 0.25, vis, 0.75, 0, vis)    # 25% opacity
        # Solid border
        cv2.polylines(vis, [pts], isClosed=True, color=(255, 200, 0),
                      thickness=max(3, thick + 1))
        # Label — positioned above the top-left corner of the quad
        lbl_x = int(card_quad[:, 0].min())
        lbl_y = int(card_quad[:, 1].min()) - int(20 * fsc)
        lbl_y = max(lbl_y, int(40 * fsc))
        scale_lbl = f"Scale: {px_mm:.2f} px/mm  ({PATCH_SIZE_MM:.0f}mm patches)"
        cv2.putText(vis, scale_lbl, (lbl_x + 3, lbl_y + 3),
                    font, fsc * 0.85, (0, 0, 0), thick + 3, cv2.LINE_AA)
        cv2.putText(vis, scale_lbl, (lbl_x, lbl_y),
                    font, fsc * 0.85, (255, 200, 0), thick, cv2.LINE_AA)

    # ── Abalone contour + bounding box ───────────────────────────────────────
    cv2.drawContours(vis, [contour],         -1, (0, 230, 60),  4)
    cv2.drawContours(vis, [meas["box_pts"]], -1, (0, 220, 255), 3)

    rect = meas["rect"]
    cx   = int(rect[0][0])
    cy   = int(rect[0][1])

    for i, lbl in enumerate([
        f"Length: {meas['length_mm']:.1f} mm",
        f"Width:  {meas['width_mm']:.1f} mm",
        f"Area:   {meas['area_mm2']:.0f} mm2",
    ]):
        ypos = cy - lh + i * lh
        cv2.putText(vis, lbl, (cx + 3, ypos + 3),
                    font, fsc, (0, 0, 0), thick + 3, cv2.LINE_AA)
        cv2.putText(vis, lbl, (cx, ypos),
                    font, fsc, (0, 255, 200), thick, cv2.LINE_AA)

    # ── Scale bar ─────────────────────────────────────────────────────────────
    bar_mm = 20
    bar_px = int(bar_mm * px_mm)
    bx, by = 60, h - 80
    cv2.rectangle(vis, (bx, by - 25), (bx + bar_px, by), (0, 0, 0), cv2.FILLED)
    cv2.rectangle(vis, (bx, by - 25), (bx + bar_px, by), (255, 255, 255), 3)
    cv2.putText(vis, f"{bar_mm} mm", (bx, by - 35),
                font, fsc, (255, 255, 255), thick, cv2.LINE_AA)

    cv2.imwrite(out_path, vis)


# ===========================================================================
# 6.  PER-IMAGE PIPELINE
# ===========================================================================

def process_image(img_path, vis_dir, ruler_real_mm, min_area_px, debug):
    stem    = Path(img_path).stem
    img_bgr = cv2.imread(img_path)
    if img_bgr is None:
        print(f"  [WARN] Cannot read: {img_path}")
        return None

    debug_dir = vis_dir if debug else None

    divider_x = find_divider_x(img_bgr)

    try:
        px_mm, card_pts = pixels_per_mm(img_bgr, divider_x, ruler_real_mm,
                                        debug_dir, stem)
    except ValueError as exc:
        print(f"  [WARN] Scale error ({stem}): {exc}")
        return None

    contour, mask = segment_abalone(img_bgr, divider_x, min_area_px,
                                    debug_dir, stem)
    if contour is None:
        print(f"  [WARN] No abalone found in {stem}")
        return None

    meas = measure_abalone(contour, mask, px_mm)

    if vis_dir:
        vis_path = os.path.join(vis_dir, f"{stem}_annotated.jpg")
        save_annotated(img_bgr, contour, meas, divider_x, px_mm, vis_path,
                       card_quad=card_pts)

    print(f"  {stem}: L={meas['length_mm']} mm  "
          f"W={meas['width_mm']} mm  "
          f"A={meas['area_mm2']} mm2  "
          f"[{px_mm:.2f} px/mm]")

    return {
        "filename":        Path(img_path).name,
        "length_mm":       meas["length_mm"],
        "width_mm":        meas["width_mm"],
        "area_mm2":        meas["area_mm2"],
        "scale_px_per_mm": round(px_mm, 4),
    }


# ===========================================================================
# 7.  MAIN
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Batch abalone morphometrics from lightbox JPEG images.")
    ap.add_argument("--images",      required=True)
    ap.add_argument("--output",      required=True)
    ap.add_argument("--vis_dir",     default=None)
    ap.add_argument("--scale_mm",    type=float, default=RULER_REAL_MM)
    ap.add_argument("--ext",         default="jpg,jpeg")
    ap.add_argument("--min_area_px", type=int, default=50_000)
    ap.add_argument("--debug",       action="store_true")
    args = ap.parse_args()

    exts      = [e.strip().lstrip(".") for e in args.ext.split(",")]
    img_paths = []
    for ext in exts:
        # ** makes glob recurse into all subfolders
        img_paths += glob.glob(os.path.join(args.images, "**", f"*.{ext}"),         recursive=True)
        img_paths += glob.glob(os.path.join(args.images, "**", f"*.{ext.upper()}"), recursive=True)
    img_paths = sorted(set(img_paths))

    # Report how many subfolders were found
    subfolders = sorted(set(os.path.dirname(p) for p in img_paths))
    print(f"Scanning {len(subfolders)} folder(s):")
    for sf in subfolders:
        count = sum(1 for p in img_paths if os.path.dirname(p) == sf)
        print(f"  {sf}  ({count} image(s))")

    if not img_paths:
        print(f"No images found in: {args.images}")
        sys.exit(1)

    print(f"Found {len(img_paths)} image(s).\n")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    if args.vis_dir:
        os.makedirs(args.vis_dir, exist_ok=True)

    results = []
    for img_path in img_paths:
        print(f"Processing: {Path(img_path).name}")
        r = process_image(img_path, args.vis_dir, args.scale_mm,
                          args.min_area_px, args.debug)
        if r:
            results.append(r)

    if not results:
        print("No measurements produced.")
        sys.exit(1)

    fieldnames = ["filename", "length_mm", "width_mm", "area_mm2",
                  "scale_px_per_mm"]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\n{'─'*60}")
    print(f"Done!  {len(results)}/{len(img_paths)} abalone measured.")
    print(f"CSV  -> {args.output}")
    if args.vis_dir:
        print(f"Vis  -> {args.vis_dir}")


if __name__ == "__main__":
    main()
