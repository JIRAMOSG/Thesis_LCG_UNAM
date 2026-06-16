"""
CELL MOTILITY ANALYSIS 

Author: Johana Itzel Ramos-Galguera

Description: fluorescent cell motility quantification from Zeiss CZI timelapse (confocal microscopy) data
exported as a multi-channel ImageJ TIFF hyperstack (the imput for this analysis was generated using FIJI and thr macro file "Confocalprocessing_brightfield_fluorescense.java".

#################################################################################################

Filtters:
1. TISSUE MASK  — detections outside the zebrafish body are rejected.
   A brightfield-derived fish-body mask (BF < 180, largest component,
   dilated 25 px) is applied to every frame, eliminating background noise
   and debris that would otherwise appear as false positives (~34% of raw
   detections in the example dataset).
2. ROBUST TRACKING  — max linking distance reduced to 40 µm and a 4×
   area-consistency check prevents jumps between two separate nearby cells.
   Candidates are globally sorted by distance before assignment, removing
   order-dependency in the greedy linker (fixes the spurious long-range
   links visible in f1).

Input format:
Your TIF must be a 2-channel ImageJ hyperstack where:
  • Channel 1 (frames 0, 2, 4 … in the TIFF) = fluorescent marker (sparse, dim)
  • Channel 2 (frames 1, 3, 5 … in the TIFF) = brightfield reference (bright)

  # Single file (all parameters auto-read from TIFF metadata)
  # default parameters:
  python cell_motility_analysis.py file.tif 
      --fl-channel 1     # 1-based channel number of the fluorescent marker
      --um-per-px 0.566  # pixel size in microns
      --dt 20            # frame interval in minutes
      --min-area 6       # minimum cell area in µm²
      --max-area 1600    # maximum cell area in µm²
      --max-link 40      # max distance a cell can travel between frames (µm)
      --max-area-ratio 4 # max fold-change in cell area between frames
      --min-track 5      # minimum number of timepoints to keep a track

outputs:
  *_tracks.csv               — one row per track: speed, displacement, directness
  *_msd.csv                  — MSD at each lag for every track
  *_detections.csv           — raw per-frame centroids and areas
  *_cell_detection.png       — two-timepoint detection overlay (red cells on BF)
  *_thesis_figure.png / .pdf — 6-panel publication figure (300 dpi)
  *_motility_video.mp4       — split-panel tracking video

Requirements:
  pip install Pillow numpy matplotlib opencv-python


#################################################################################################

ANALYSIS 

#################################################################################################
1. Read data
   Read pixel size (µm/px) and frame interval (min) from TIFF ImageDescription tag.
   Fall back to command-line values if not present.

2. ROI
   Brightfield tissue pixels are darker than the white background.
   We find the bounding box to avoid processing empty space at image edges.

3. Body mask
   From the first brightfield frame:
   a. Threshold at BF < 180  -> binary tissue map
   b. Morphological close 30 px  -> fill gaps within the fish body
   c. Keep only the LARGEST connected component  -> discard small debris patches
   d. Dilate 25 px  -> extend mask to include cells at the tissue boundary
   This mask is applied to every frame after thresholding (step 4e).

4. Cell detection
   a. Gaussian blur σ=2 px  -> removes single-pixel shot noise
   b. Morphological opening 51 px -> estimates slowly-varying background (erases features smaller than ~29 µm)
   c. Background subtraction  -> reveals cells above local baseline
   d. Otsu threshold  -> automatic signal/noise boundary
   e. Morphological closing 5 px  -> fills small holes in cell bodies
   f. Apply fish body mask   -> discard detections outside tissue  [NEW]
   g. Area filter 20–5000 px² -> removes noise and large artefacts
   Centroids and contours returned for each valid component.

5. Tracking
   For each pair of consecutive frames t-1 -> t:
   a. Build a list of ALL (distance, track_idx, detection_idx) candidates where distance ≤ max_link AND area_ratio ≤ max_area_ratio
   b. Sort candidates by distance (shortest first)
   c. Assign greedily, skipping already-used tracks or detections
   This globally-optimal greedy approach prevents order-dependent errors.

6. Motility data  (per track)
   • Mean speed (µm/min)   = total_path / duration
   • Net displacement (µm) = Euclidean distance from first to last point
   • Total path (µm)       = sum of all step distances
   • Directness (0–1)      = net_displacement / total_path
   • MSD at lag τ          = mean( Δx(τ)² + Δy(τ)² ) over all valid pairs

7. Imaging output
   • Composite image: BF (gray) + fluorescent cells (red, detected pixels only)
   • Track overlay: plasma colormap by mean speed
   • Video: left panel = per-frame detection, right panel = track history
"""

import os, sys, csv, argparse, warnings
import numpy as np
import cv2
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from PIL import Image
import subprocess

warnings.filterwarnings('ignore')
# Default parameters 
P = dict(
    fl_channel      = 1,       # 1-based channel index for fluorescent marker
    n_channels      = 2,       # total channels in the hyperstack
    um_per_px       = None,    # µm/pixel — auto from TIFF, else set here
    dt_min          = None,    # frame interval (minutes) — auto from TIFF
    bg_kernel       = 51,      # rolling-ball kernel size (pixels)
    gauss_sigma     = 2.0,     # Gaussian pre-blur sigma (pixels)
    min_area_um2    = 6.0,     # minimum cell area (µm²)
    max_area_um2    = 1600.0,  # maximum cell area (µm²)
    max_link_um     = 40.0,    # maximum linking distance between frames (µm) [reduced from 80]
    max_area_ratio  = 4.0,     # maximum fold-change in cell area per step [NEW]
    min_track       = 5,       # minimum timepoints per track
    max_lag         = 10,      # maximum lag for MSD curve
    fl_boost        = 3.0,     # brightness multiplier for red overlay
    video_scale     = 0.4,     # downscale factor for video frames
    video_fps       = 3,       # output video frames per second
    # Tissue mask parameters [NEW]
    bf_tissue_thr   = 180,     # BF pixel < this = tissue (was BF < 200 for ROI)
    tissue_close_px = 30,      # morphological close radius for fish mask
    tissue_dilate_px= 25,      # dilation radius to include boundary cells
)


# 
# 1. Load data
# 
def read_tif_metadata(path):
    """Return (um_per_px, dt_min, n_channels, n_timepoints) from ImageJ TIFF tags."""
    with Image.open(path) as img:
        tags    = img.tag_v2 if hasattr(img, 'tag_v2') else {}
        desc    = tags.get(270, '')
        px_res  = tags.get(282, None)
        n_frames = getattr(img, 'n_frames', 1)
    meta = {}
    for line in desc.replace('\r', '\n').split('\n'):
        if '=' in line:
            k, _, v = line.partition('=')
            meta[k.strip()] = v.strip()
    n_ch = int(meta.get('channels', 2))
    n_tp = int(meta.get('frames', n_frames // max(n_ch, 1)))
    fint = float(meta.get('finterval', 0))
    dt   = fint / 60.0 if fint > 0 else None
    um   = None
    if px_res is not None:
        r = float(px_res[0]) if isinstance(px_res, tuple) else float(px_res)
        if r > 0:
            um = 1.0 / r
    return um, dt, n_ch, n_tp


# 
# 2. ROI 
# 
def get_tissue_roi(path, bf_channel, n_ch, n_tp):
    """Find bounding box of zebrafish tissue using brightfield (dark on bright BG)."""
    ch_idx = bf_channel - 1
    with Image.open(path) as img:
        img.seek(ch_idx)
        bf = np.array(img, dtype=np.uint8)
    tissue = bf < 200
    cols   = np.where(tissue.sum(axis=0) > 50)[0]
    rows   = np.where(tissue.sum(axis=1) > 50)[0]
    if len(cols) == 0 or len(rows) == 0:
        return (0, 0, bf.shape[1], bf.shape[0])
    h, w = bf.shape
    return (max(int(cols[0])-10, 0),
            max(int(rows[0])-10, 0),
            min(int(cols[-1])+10, w),
            min(int(rows[-1])+10, h))


# 
# 3. Body mask 
# 
def build_fish_mask(bf_roi, tissue_thr=180, close_px=30, dilate_px=25):
    """
    Build a binary mask of the fish body from a brightfield ROI.

    Steps:
      1. BF < tissue_thr  -> raw tissue pixels (white=tissue, black=background)
      2. Morphological close (close_px)  -> fill internal gaps (e.g. gut, swim bladder)
      3. Keep only the LARGEST connected component  -> remove debris patches
      4. Dilate (dilate_px)  -> extend to include boundary cells

    Returns: uint8 mask (255 = inside fish, 0 = outside)

    Why BF < 180 (not 200)?
    BF < 200 includes 15,000+ disconnected fragments in typical images.
    BF < 180 gives a clean single-body segmentation with far fewer fragments,
    making the largest-component selection reliable.
    """
    tissue_raw = (bf_roi < tissue_thr).astype(np.uint8) * 255
    kern_close = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (close_px, close_px))
    closed     = cv2.morphologyEx(tissue_raw, cv2.MORPH_CLOSE, kern_close)
    n_lbl, labels, stats, _ = cv2.connectedComponentsWithStats(closed)
    if n_lbl < 2:
        return closed  # fallback: no components found
    largest_idx = 1 + stats[1:, cv2.CC_STAT_AREA].argmax()
    fish_only   = (labels == largest_idx).astype(np.uint8) * 255
    kern_dil    = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (dilate_px, dilate_px))
    fish_mask   = cv2.dilate(fish_only, kern_dil)
    return fish_mask


# 
# 4. Cell Detection
# 
def detect_cells_in_roi(fl_roi, min_area_px, max_area_px,
                        fish_mask=None, bg_kernel=51, gauss_sigma=2.0):
    """
    Detect fluorescent cells in a single ROI frame, optionally filtered by fish mask.
    Returns: list of (cx_roi, cy_roi, area_px2),  binary mask (uint8)
    """
    g    = cv2.GaussianBlur(fl_roi, (5, 5), gauss_sigma)
    kern = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (bg_kernel, bg_kernel))
    bg   = cv2.morphologyEx(g, cv2.MORPH_OPEN, kern)
    sub  = cv2.subtract(g, bg)
    _, thr = cv2.threshold(sub, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    cl   = cv2.morphologyEx(thr, cv2.MORPH_CLOSE,
           cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))

    # Apply fish body mask before finding contours
    if fish_mask is not None:
        cl = cv2.bitwise_and(cl, fish_mask)

    cnts, _ = cv2.findContours(cl, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cells = []
    for c in cnts:
        a = cv2.contourArea(c)
        if min_area_px <= a <= max_area_px:
            M = cv2.moments(c)
            if M['m00'] > 0:
                cells.append((M['m10']/M['m00'], M['m01']/M['m00'], a))
    return cells, cl


def detect_all(path, fl_ch, n_ch, n_tp, roi, um_per_px,
               min_area_um2, max_area_um2, bg_kernel, gauss_sigma,
               tissue_thr=180, close_px=30, dilate_px=25):
    """
    Loop over all timepoints, detect cells within the fish body mask.
    Returns list-of-lists of (cx, cy, area) in full-image pixel coordinates.
    """
    min_px = min_area_um2 / um_per_px**2
    max_px = max_area_um2 / um_per_px**2
    ch_idx = fl_ch - 1
    bf_ch  = 2 if fl_ch == 1 else 1
    rx0, ry0, rx1, ry1 = roi

    # Build fish mask from first BF frame
    print('  Building fish body mask from brightfield...')
    with Image.open(path) as img:
        img.seek(bf_ch - 1)  # first BF frame
        bf0 = np.array(img, dtype=np.uint8)[ry0:ry1, rx0:rx1]
    fish_mask = build_fish_mask(bf0, tissue_thr, close_px, dilate_px)
    pct_fish  = 100 * (fish_mask > 0).sum() / fish_mask.size
    print(f'  Fish mask covers {pct_fish:.1f}% of ROI')

    all_det = []
    with Image.open(path) as img:
        for t in range(n_tp):
            img.seek(t * n_ch + ch_idx)
            fl     = np.array(img, dtype=np.uint8)
            fl_roi = fl[ry0:ry1, rx0:rx1]
            cells, _ = detect_cells_in_roi(
                fl_roi, min_px, max_px, fish_mask, bg_kernel, gauss_sigma)
            # Convert ROI-local -> full-image coords
            all_det.append([(cx+rx0, cy+ry0, a) for cx, cy, a in cells])
            if t % 10 == 0:
                print(f'  t={t:3d}/{n_tp}: {len(cells)} cells (tissue-filtered)')
    return all_det, fish_mask


# 
# 5. Tracking  
# 
def link_tracks(all_det, max_link_um, um_per_px, min_track, max_area_ratio=4.0):
    """
    Globally-sorted greedy nearest-neighbour tracking.

    Improvements over f1:
    • max_link reduced to 40 µm (was 80 µm)  — prevents long spurious jumps
    • Area-consistency filter: reject links where cell area changes by >4×
      (prevents linking to debris, adjacent clusters, or cell doublets)
    • Globally-sorted candidates: all (distance, track, detection) pairs are
      sorted by distance before assignment, so the globally shortest link is
      always assigned first regardless of iteration order.

    Returns list of tracks (each track = list of (t, cx_px, cy_px, area_px)).
    Only tracks with >= min_track timepoints are returned.
    """
    max_px = max_link_um / um_per_px
    # Seed tracks from frame 0
    tracks = [[(0, cx, cy, a)] for cx, cy, a in all_det[0]]

    for t in range(1, len(all_det)):
        dets      = list(all_det[t])
        used_trk  = [False] * len(tracks)
        used_det  = [False] * len(dets)

        # Build ALL valid candidate links
        candidates = []
        for ti, tr in enumerate(tracks):
            if tr[-1][0] != t - 1:
                continue  # track not active last frame
            lx, ly, la = tr[-1][1], tr[-1][2], tr[-1][3]
            for di, (cx, cy, a) in enumerate(dets):
                d = np.hypot(cx - lx, cy - ly)
                if d > max_px:
                    continue
                # Area consistency check
                ratio = a / la if la > 0 else 999.0
                if ratio > max_area_ratio or ratio < 1.0 / max_area_ratio:
                    continue
                candidates.append((d, ti, di))

        # Assign shortest-distance first (globally optimal greedy)
        candidates.sort()
        for d, ti, di in candidates:
            if used_trk[ti] or used_det[di]:
                continue
            cx, cy, a = dets[di]
            tracks[ti].append((t, cx, cy, a))
            used_trk[ti] = True
            used_det[di] = True

        # Start new tracks for unmatched detections
        for di, (cx, cy, a) in enumerate(dets):
            if not used_det[di]:
                tracks.append([(t, cx, cy, a)])

    return [tr for tr in tracks if len(tr) >= min_track]


# 
# 6. Motility data
# 
def compute_metrics(tracks, um_per_px, dt_min, max_lag):
    metrics = []
    for tr in tracks:
        tpts = np.array([p[0] for p in tr])
        xs   = np.array([p[1] for p in tr]) * um_per_px
        ys   = np.array([p[2] for p in tr]) * um_per_px
        dx, dy   = np.diff(xs), np.diff(ys)
        steps    = np.sqrt(dx**2 + dy**2)
        net      = float(np.sqrt((xs[-1]-xs[0])**2 + (ys[-1]-ys[0])**2))
        path     = float(steps.sum())
        dirn     = net / path if path > 0 else 0.0
        speed    = float((steps / dt_min).mean())
        duration = float((tpts[-1] - tpts[0]) * dt_min)
        msds = []
        for lag in range(1, min(max_lag+1, len(tr))):
            msds.append(float(np.mean(
                [(xs[j+lag]-xs[j])**2 + (ys[j+lag]-ys[j])**2
                 for j in range(len(tr)-lag)])))
        metrics.append({
            'n': len(tr), 'dur': duration, 'speed': speed,
            'net': net, 'path': path, 'dir': dirn,
            'xs': xs.tolist(), 'ys': ys.tolist(),
            'tpts': tpts.tolist(), 'msds': msds
        })
    return metrics


def save_csvs(metrics, all_det, um_per_px, dt_min, out_prefix):
    with open(out_prefix+'_tracks.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['track_id','n_timepoints','duration_min','mean_speed_um_min',
                    'net_displacement_um','total_path_um','directness'])
        for i, m in enumerate(metrics):
            w.writerow([i, m['n'], round(m['dur'],2), round(m['speed'],4),
                        round(m['net'],2), round(m['path'],2), round(m['dir'],4)])

    with open(out_prefix+'_msd.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['track_id'] + [f'MSD_lag{l}_um2' for l in range(1, 11)])
        for i, m in enumerate(metrics):
            row = [i] + [round(m['msds'][l-1], 3) if l-1 < len(m['msds']) else ''
                         for l in range(1, 11)]
            w.writerow(row)

    with open(out_prefix+'_detections.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['frame','time_min','x_um','y_um','area_um2'])
        for t, dets in enumerate(all_det):
            for cx, cy, a in dets:
                w.writerow([t, round(t*dt_min,1), round(cx*um_per_px,2),
                             round(cy*um_per_px,2), round(a*um_per_px**2,2)])
    print(f'  CSVs saved to {out_prefix}_*.csv')


# 
# 7. Imaging output (BF gray + FL red, cells only)
# 
def make_composite_red(fl_roi, bf_roi, fl_p99, mask=None, fish_mask=None, boost=3.0):
    """
    RGB composite: BF as gray base, fluorescent cells overlaid in red.
    Only pixels that pass Otsu (and optionally the fish mask) are coloured.
    """
    if mask is None:
        g   = cv2.GaussianBlur(fl_roi, (5,5), 2.0)
        bg  = cv2.morphologyEx(g, cv2.MORPH_OPEN,
              cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (51,51)))
        sub = cv2.subtract(g, bg)
        _, mask = cv2.threshold(sub, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE,
               cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)))
        if fish_mask is not None:
            mask = cv2.bitwise_and(mask, fish_mask)
    fl_f  = np.clip(fl_roi.astype(np.float32)/fl_p99, 0, 1) * (mask > 0)
    bf_f  = bf_roi.astype(np.float32) / 255.0
    comp  = np.stack([bf_f, bf_f, bf_f], axis=-1)   # RGB
    comp[:,:,0] = np.clip(comp[:,:,0] + fl_f * boost, 0, 1)   # R channel
    return (comp * 255).astype(np.uint8)


# 
# 8. overlay
# 
def save_detection_image(path, all_det, roi, um_per_px, dt_min,
                         fl_ch, n_ch, fl_p99, fish_mask, out_path,
                         timepoints=(0, 24)):
    """
    Save detection overlay as PNG and PDF.
    Two rows: one per timepoint. Red=cells, cyan=outlines, tissue-filtered.
    """
    ch_idx = fl_ch - 1
    bf_ch  = 2 if fl_ch == 1 else 1
    fig, axes = plt.subplots(len(timepoints), 1,
                              figsize=(22, 5*len(timepoints)), facecolor='black')
    if len(timepoints) == 1:
        axes = [axes]

    with Image.open(path) as img:
        for ax, t in zip(axes, timepoints):
            img.seek(t * n_ch + ch_idx); fl = np.array(img, dtype=np.uint8)
            img.seek(t * n_ch + (bf_ch-1)); bf = np.array(img, dtype=np.uint8)
            fl_roi = fl[roi[1]:roi[3], roi[0]:roi[2]]
            bf_roi = bf[roi[1]:roi[3], roi[0]:roi[2]]

            g   = cv2.GaussianBlur(fl_roi,(5,5),2.0)
            bg  = cv2.morphologyEx(g,cv2.MORPH_OPEN,
                  cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(51,51)))
            sub = cv2.subtract(g,bg)
            _,mask=cv2.threshold(sub,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
            mask=cv2.morphologyEx(mask,cv2.MORPH_CLOSE,
                 cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)))
            if fish_mask is not None:
                mask = cv2.bitwise_and(mask, fish_mask)

            comp_rgb = make_composite_red(fl_roi, bf_roi, fl_p99, mask)
            comp_bgr = cv2.cvtColor(comp_rgb, cv2.COLOR_RGB2BGR)

            cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            valid   = [c for c in cnts if 20 <= cv2.contourArea(c) <= 5000]
            cv2.drawContours(comp_bgr, valid, -1, (0,255,255), 2)

            h, w = comp_bgr.shape[:2]
            sb = int(500/um_per_px)
            cv2.line(comp_bgr,(w-sb-25,h-20),(w-25,h-20),(255,255,255),3)
            cv2.putText(comp_bgr,'500 µm',(w-sb-25,h-6),
                        cv2.FONT_HERSHEY_SIMPLEX,0.6,(255,255,255),1)
            cv2.rectangle(comp_bgr,(0,0),(270,54),(0,0,0),-1)
            cv2.putText(comp_bgr,f't = {t*dt_min:.0f} min',(8,24),
                        cv2.FONT_HERSHEY_SIMPLEX,0.75,(0,200,255),2)
            cv2.putText(comp_bgr,f'{len(valid)} cells (tissue-filtered)',(8,48),
                        cv2.FONT_HERSHEY_SIMPLEX,0.5,(0,255,200),1)

            ax.imshow(cv2.cvtColor(comp_bgr,cv2.COLOR_BGR2RGB),aspect='auto')
            ax.axis('off')
            ax.set_title(f't={t*dt_min:.0f} min  |  red = cells (tissue-only)  |  '
                         f'cyan = detected outline', color='white', fontsize=9, pad=3)

    fig.suptitle('Cell detection — tissue-filtered (background detections removed)',
                 color='white', fontsize=12, y=1.01)
    plt.tight_layout(pad=0.3)
    # Save both PNG and PDF
    base = out_path.rsplit('.', 1)[0]
    for ext in ('png', 'pdf'):
        fig.savefig(f'{base}.{ext}', dpi=150, bbox_inches='tight', facecolor='black')
    plt.close(fig)
    print(f'  Saved: {base}.png / .pdf')


# 
# 9. videos  (three outputs: split-panel, detection-only, tracks-only)
# 
def _stitch_video(frm_dir, out_path, fps):
    """Stitch PNG frames in frm_dir into an mp4 with ffmpeg."""
    cmd = ['ffmpeg', '-y', '-framerate', str(fps), '-i', f'{frm_dir}/frame_%03d.png',
           '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
           '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20', out_path]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode == 0:
        print(f'  Saved: {out_path}')
        import shutil; shutil.rmtree(frm_dir, ignore_errors=True)
    else:
        print(f'  ffmpeg error: {r.stderr.decode()[-300:]}')
        print(f'  Frame PNGs kept at: {frm_dir}')


def render_video(path, all_det, metrics, roi, um_per_px, dt_min,
                 fl_ch, n_ch, n_tp, fl_p99, fish_mask, out_path,
                 scale=0.4, fps=3):
    """
    Render three videos from the same set of frames:
      • <base>_motility_video.mp4    — split panel (detection LEFT | tracks RIGHT)
      • <base>_cell_detection_video.mp4 — detection only
      • <base>_cell_tracks_video.mp4    — tracks only

    All three are produced in a single pass over the frames to avoid re-reading
    the large TIFF multiple times.
    """
    rx0, ry0, rx1, ry1 = roi
    ch_idx = fl_ch - 1
    bf_ch  = 2 if fl_ch == 1 else 1
    pw = int((rx1-rx0)*scale); ph = int((ry1-ry0)*scale)

    speeds  = np.array([m['speed'] for m in metrics])
    sp_min  = speeds.min(); sp_max = np.percentile(speeds, 95)
    cmap_cv = cv2.applyColorMap(np.arange(256,dtype=np.uint8).reshape(1,256),
                                 cv2.COLORMAP_PLASMA)[0,:]
    tdata = []
    for m in metrics:
        sn  = np.clip((m['speed']-sp_min)/(sp_max-sp_min+1e-9), 0, 1)
        ci  = int(sn*255)
        col = (int(cmap_cv[ci][0]), int(cmap_cv[ci][1]), int(cmap_cv[ci][2]))
        tdata.append({'tpts': m['tpts'],
                      'xs_px': [v/um_per_px for v in m['xs']],
                      'ys_px': [v/um_per_px for v in m['ys']],
                      'col': col})

    # Static dim BF background for track panel
    with Image.open(path) as img:
        img.seek(bf_ch - 1)
        bf0 = np.array(img, dtype=np.uint8)[ry0:ry1, rx0:rx1]
    track_bg = cv2.resize(
        (bf0.astype(np.float32)*0.35).clip(0,255).astype(np.uint8), (pw, ph))
    track_bg = cv2.cvtColor(track_bg, cv2.COLOR_GRAY2BGR)

    # Temporary frame directories
    base     = out_path.rsplit('_motility', 1)[0]   # strip suffix if present
    spl_dir  = base + '_tmp_split'
    det_dir  = base + '_tmp_det'
    trk_dir  = base + '_tmp_trk'
    for d in (spl_dir, det_dir, trk_dir):
        os.makedirs(d, exist_ok=True)

    sb = int(500/um_per_px*scale)
    print(f'  Rendering {n_tp} frames (3 video outputs)...')

    with Image.open(path) as img:
        for t in range(n_tp):
            img.seek(t*n_ch + ch_idx); fl = np.array(img, dtype=np.uint8)
            img.seek(t*n_ch + (bf_ch-1)); bf = np.array(img, dtype=np.uint8)
            fl_roi = fl[ry0:ry1, rx0:rx1]; bf_roi = bf[ry0:ry1, rx0:rx1]

            # Detection panel 
            g   = cv2.GaussianBlur(fl_roi,(5,5),2.0)
            bg  = cv2.morphologyEx(g,cv2.MORPH_OPEN,
                  cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(51,51)))
            sub = cv2.subtract(g,bg)
            _,mask = cv2.threshold(sub,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
            mask   = cv2.morphologyEx(mask,cv2.MORPH_CLOSE,
                     cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)))
            if fish_mask is not None:
                mask = cv2.bitwise_and(mask, fish_mask)
            fl_f = np.clip(fl_roi.astype(np.float32)/fl_p99,0,1)*(mask>0)
            bf_f = bf_roi.astype(np.float32)/255.0
            comp = np.stack([bf_f,bf_f,bf_f],axis=-1)
            comp[:,:,2] = np.clip(comp[:,:,2]+fl_f*3.0,0,1)  # BGR: index 2=red
            det = (comp*255).astype(np.uint8)
            cnts,_=cv2.findContours(mask,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            for c in [c for c in cnts if 20<=cv2.contourArea(c)<=5000]:
                cv2.drawContours(det,[(c*scale).astype(np.int32)],-1,(255,255,0),1)
            det = cv2.resize(det,(pw,ph))
            # Labels
            cv2.rectangle(det,(0,0),(185,30),(0,0,0),-1)
            cv2.putText(det,'CELL DETECTION',(5,20),cv2.FONT_HERSHEY_SIMPLEX,0.5,(0,210,255),1)
            cv2.rectangle(det,(0,ph-32),(155,ph),(0,0,0),-1)
            cv2.putText(det,f't={t*dt_min:.0f} min',(5,ph-10),cv2.FONT_HERSHEY_SIMPLEX,0.6,(255,215,0),2)
            cv2.line(det,(pw-sb-22,ph-18),(pw-22,ph-18),(255,255,255),2)
            cv2.putText(det,'500µm',(pw-sb-22,ph-4),cv2.FONT_HERSHEY_SIMPLEX,0.36,(255,255,255),1)

            #  Tracks panel 
            trk = track_bg.copy()
            for td in tdata:
                pts = [(i,ti) for i,ti in enumerate(td['tpts']) if ti<=t]
                for k in range(1,len(pts)):
                    i1,_=pts[k-1]; i2,t2=pts[k]
                    p1=(int((td['xs_px'][i1]-rx0)*scale),int((td['ys_px'][i1]-ry0)*scale))
                    p2=(int((td['xs_px'][i2]-rx0)*scale),int((td['ys_px'][i2]-ry0)*scale))
                    alpha = max(0.15, 1.0-(t-t2)/14.0)
                    cv2.line(trk,p1,p2,tuple(int(c*alpha) for c in td['col']),2,cv2.LINE_AA)
                last=[i for i,ti in enumerate(td['tpts']) if ti<=t]
                if last:
                    li=last[-1]
                    cv2.circle(trk,(int((td['xs_px'][li]-rx0)*scale),
                                    int((td['ys_px'][li]-ry0)*scale)),4,td['col'],-1)
            cv2.rectangle(trk,(0,0),(155,30),(0,0,0),-1)
            cv2.putText(trk,'CELL TRACKS',(5,20),cv2.FONT_HERSHEY_SIMPLEX,0.5,(180,180,255),1)
            cv2.rectangle(trk,(0,ph-32),(155,ph),(0,0,0),-1)
            cv2.putText(trk,f't={t*dt_min:.0f} min',(5,ph-10),cv2.FONT_HERSHEY_SIMPLEX,0.6,(255,215,0),2)
            cv2.line(trk,(pw-sb-22,ph-18),(pw-22,ph-18),(255,255,255),2)
            cv2.putText(trk,'500µm',(pw-sb-22,ph-4),cv2.FONT_HERSHEY_SIMPLEX,0.36,(255,255,255),1)

            # Enforce even dimensions
            dw, dh = det.shape[1]//2*2, det.shape[0]//2*2
            det = det[:dh,:dw]; trk = trk[:dh,:dw]

            # Split panel 
            div  = np.full((dh,4,3),128,dtype=np.uint8)
            spl  = np.concatenate([det, div, trk], axis=1)

            cv2.imwrite(f'{det_dir}/frame_{t:03d}.png', det)
            cv2.imwrite(f'{trk_dir}/frame_{t:03d}.png', trk)
            cv2.imwrite(f'{spl_dir}/frame_{t:03d}.png', spl)

    # Stitch all three
    _stitch_video(spl_dir, out_path,               fps)   # split panel
    _stitch_video(det_dir, base+'_cell_detection_video.mp4', fps)
    _stitch_video(trk_dir, base+'_cell_tracks_video.mp4',    fps)


# 
# 10. .pdf summary fugure per fish
# 
def save_thesis_figure(path, all_det, metrics, roi, um_per_px, dt_min,
                       fl_ch, n_ch, n_tp, fl_p99, fish_mask, out_prefix, min_track):
    rx0,ry0,rx1,ry1 = roi; ch_idx = fl_ch - 1
    bf_ch = 2 if fl_ch == 1 else 1
    speeds = np.array([m['speed'] for m in metrics])
    nets   = np.array([m['net']   for m in metrics])
    dirns  = np.array([m['dir']   for m in metrics])
    counts = [len(d) for d in all_det]; times = np.arange(n_tp)*dt_min

    # Max-projection (frame-by-frame to save RAM)
    fl_max = None; bf_t0 = None
    with Image.open(path) as img:
        for t in range(n_tp):
            img.seek(t*n_ch + ch_idx)
            fl = np.array(img, dtype=np.uint8)[ry0:ry1, rx0:rx1]
            fl_max = fl if fl_max is None else np.maximum(fl_max, fl)
            if t == 0:
                img.seek(bf_ch - 1)
                bf_t0 = np.array(img, dtype=np.uint8)[ry0:ry1, rx0:rx1]

    # Cell mask on max-projection
    g   = cv2.GaussianBlur(fl_max,(5,5),2.0)
    bg  = cv2.morphologyEx(g,cv2.MORPH_OPEN,
          cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(51,51)))
    sub = cv2.subtract(g,bg)
    _,mask = cv2.threshold(sub,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    if fish_mask is not None:
        mask = cv2.bitwise_and(mask, fish_mask)

    comp_rgb = make_composite_red(fl_max, bf_t0, fl_p99, mask)

    # Draw tracks (plasma by speed)
    norm     = Normalize(vmin=speeds.min(), vmax=np.percentile(speeds,95))
    cmap_plt = cm.plasma
    for m in metrics:
        xs_r = np.array(m['xs'])/um_per_px - rx0
        ys_r = np.array(m['ys'])/um_per_px - ry0
        c    = cmap_plt(norm(m['speed']))
        rc   = (int(c[0]*255), int(c[1]*255), int(c[2]*255))
        for k in range(1,len(xs_r)):
            cv2.line(comp_rgb,(int(xs_r[k-1]),int(ys_r[k-1])),
                     (int(xs_r[k]),int(ys_r[k])),rc,2,cv2.LINE_AA)
        cv2.circle(comp_rgb,(int(xs_r[0]),int(ys_r[0])),5,rc,-1)
        cv2.drawMarker(comp_rgb,(int(xs_r[-1]),int(ys_r[-1])),rc,
                       cv2.MARKER_TRIANGLE_UP,10,2)

    plt.rcParams.update({'font.size':9,'axes.linewidth':0.8,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig = plt.figure(figsize=(18/2.54,22/2.54), facecolor='white')
    gs  = gridspec.GridSpec(3,3,figure=fig,left=0.10,right=0.97,top=0.93,
                            bottom=0.09,wspace=0.45,hspace=0.55)

    ax_a = fig.add_subplot(gs[0,:])
    ax_a.imshow(comp_rgb, aspect='auto', interpolation='bilinear')
    sm  = plt.cm.ScalarMappable(cmap='plasma', norm=norm); sm.set_array([])
    cb  = plt.colorbar(sm, ax=ax_a, fraction=0.012, pad=0.01)
    cb.set_label('Speed (µm/min)', fontsize=7.5); cb.ax.tick_params(labelsize=7)
    rw=rx1-rx0; rh=ry1-ry0; sb=int(500/um_per_px)
    ax_a.plot([rw-sb-40,rw-40],[rh-40,rh-40],'w-',lw=2.5)
    ax_a.text(rw-sb//2-40,rh-60,'500 µm',color='w',ha='center',fontsize=7.5,fontweight='bold')
    ax_a.set_xticks([]); ax_a.set_yticks([])
    ax_a.set_title('Max-projection: BF (gray) + FL cells (red, tissue-only) + tracks (plasma=speed)',
                   fontsize=8, pad=4)
    ax_a.text(-0.03,1.06,'A',transform=ax_a.transAxes,fontsize=12,fontweight='bold',va='top')

    for ax,data,xlabel,title,letter,col in [
        (fig.add_subplot(gs[1,0]),speeds,'Speed (µm/min)','Speed','B','#4C72B0'),
        (fig.add_subplot(gs[1,1]),nets,'Net displacement (µm)','Net displacement','C','#55A868'),
        (fig.add_subplot(gs[1,2]),dirns,'Directness index','Directness','D','#8172B2')]:
        ax.hist(data,bins=12,color=col,edgecolor='white',linewidth=0.4)
        ax.axvline(data.mean(),color='#C44E52',lw=1.3,ls='--',
                   label=f'Mean={data.mean():.2f}')
        ax.axvline(np.median(data),color='#DD8452',lw=1.3,ls=':',
                   label=f'Med={np.median(data):.2f}')
        ax.set_xlabel(xlabel); ax.set_ylabel('Tracks'); ax.set_title(title,fontsize=8.5)
        ax.legend(fontsize=6.5,frameon=False)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.text(-0.22,1.08,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',va='top')

    ax_e = fig.add_subplot(gs[2,0:2])
    lag_t = np.arange(1,11)*dt_min
    vm = [m['msds'] for m in metrics if len(m['msds'])>=8]
    if vm:
        mat = np.array([v[:8] for v in vm])
        for r in mat: ax_e.plot(lag_t[:8],r,color='#4C72B0',alpha=0.15,lw=0.7)
        mn = mat.mean(0); se = mat.std(0)/np.sqrt(len(mat))
        ax_e.fill_between(lag_t[:8],mn-se,mn+se,color='#C44E52',alpha=0.25)
        ax_e.plot(lag_t[:8],mn,'o-',color='#C44E52',lw=1.8,ms=4,
                  label=f'Mean±SEM (n={len(mat)})')
    ax_e.set_xlabel('Lag time (min)'); ax_e.set_ylabel('MSD (µm²)')
    ax_e.set_title('Mean squared displacement',fontsize=8.5)
    ax_e.legend(fontsize=7,frameon=False)
    ax_e.text(-0.12,1.08,'E',transform=ax_e.transAxes,fontsize=12,fontweight='bold',va='top')

    ax_f = fig.add_subplot(gs[2,2])
    ax_f.plot(times,counts,'o-',color='#C44E52',lw=1.4,ms=3.5,
              markerfacecolor='white',markeredgewidth=1.2)
    ax_f.set_xlabel('Time (min)'); ax_f.set_ylabel('Cells per frame')
    ax_f.set_title('Cells detected over time',fontsize=8.5)
    ax_f.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_f.text(-0.22,1.08,'F',transform=ax_f.transAxes,fontsize=12,fontweight='bold',va='top')

    cap = (f"Figure X. Fluorescent cell motility in zebrafish. "
           "(A) Maximum-intensity projection: brightfield (gray), fluorescent cells "
           "(red, tissue-filtered), tracks coloured by mean speed (plasma). "
           "●=start ▲=end. Scale bar 500 µm. "
           f"(B–D) Speed, net displacement and directness (n={len(metrics)} tracks "
           f"≥{min_track} timepoints). Dashed=mean, dotted=median. "
           "(E) MSD vs lag; blue=individual, red=mean±SEM. "
           "(F) Cells detected per frame. "
           "Background detections removed by BF tissue mask. "
           "Tracking: max_link=40 µm, area_ratio≤4, globally-sorted greedy.")
    fig.text(0.01,0.002,cap,fontsize=6.0,color='#333',ha='left',va='bottom',wrap=True)

    for ext in ('png','pdf'):
        fig.savefig(out_prefix+f'_thesis_figure.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Figure saved: {out_prefix}_thesis_figure.png/pdf')



def parse_args():
    p = argparse.ArgumentParser(
        description='Cell motility analysis (f2) — tissue-filtered + robust tracking')
    p.add_argument('tif', help='Path to ImageJ TIFF hyperstack')
    p.add_argument('--fl-channel',    type=int,   default=None, dest='fl_ch')
    p.add_argument('--um-per-px',     type=float, default=None)
    p.add_argument('--dt',            type=float, default=None)
    p.add_argument('--min-area',      type=float, default=None)
    p.add_argument('--max-area',      type=float, default=None)
    p.add_argument('--max-link',      type=float, default=None)
    p.add_argument('--max-area-ratio',type=float, default=None)
    p.add_argument('--min-track',     type=int,   default=None)
    return p.parse_args()


def main():
    args     = parse_args()
    tif_path = os.path.abspath(args.tif)
    if not os.path.isfile(tif_path):
        sys.exit(f'File not found: {tif_path}')

    out_dir    = os.path.dirname(tif_path)
    base       = os.path.splitext(os.path.basename(tif_path))[0]
    out_prefix = os.path.join(out_dir, base)

    print(f'\n{"═"*65}')
    print(f'  Cell Motility Analysis  [f2 — tissue-filtered + robust tracking]')
    print(f'  {tif_path}')
    print(f'{"═"*65}')

    # Read metadata
    um_meta, dt_meta, n_ch_meta, n_tp = read_tif_metadata(tif_path)
    um     = args.um_per_px      or um_meta   or 0.5657
    dt     = args.dt             or dt_meta   or 20.0
    fl_ch  = args.fl_ch          or P['fl_channel']
    n_ch   = n_ch_meta           or P['n_channels']
    mina   = args.min_area       or P['min_area_um2']
    maxa   = args.max_area       or P['max_area_um2']
    maxl   = args.max_link       or P['max_link_um']
    maxar  = args.max_area_ratio or P['max_area_ratio']
    mint   = args.min_track      or P['min_track']

    print(f'  Pixel size: {um:.4f} µm | Interval: {dt:.1f} min | '
          f'Channels: {n_ch} | Timepoints: {n_tp}')
    print(f'  FL channel: {fl_ch} | max_link: {maxl} µm | max_area_ratio: {maxar}')

    # ROI
    bf_ch = 2 if fl_ch == 1 else 1
    roi   = get_tissue_roi(tif_path, bf_ch, n_ch, n_tp)
    print(f'  Tissue ROI: x={roi[0]}-{roi[2]}, y={roi[1]}-{roi[3]}')

    # FL brightness reference
    with Image.open(tif_path) as img:
        img.seek(fl_ch - 1)
        fl0 = np.array(img, dtype=np.uint8)[roi[1]:roi[3], roi[0]:roi[2]]
    fl_p99 = float(np.percentile(fl0[fl0>5], 99)) if fl0.max() > 5 else 50.0

    # Detect with tissue mask
    print('\nDetecting cells (with tissue filter)...')
    all_det, fish_mask = detect_all(
        tif_path, fl_ch, n_ch, n_tp, roi, um, mina, maxa,
        P['bg_kernel'], P['gauss_sigma'],
        P['bf_tissue_thr'], P['tissue_close_px'], P['tissue_dilate_px'])

    # Track with improved linker
    print('\nTracking (globally-sorted greedy, area consistency)...')
    long_tracks = link_tracks(all_det, maxl, um, mint, maxar)
    print(f'  Long tracks (≥{mint} tp): {len(long_tracks)}')

    # Metrics
    metrics = compute_metrics(long_tracks, um, dt, P['max_lag'])
    sp = np.array([m['speed'] for m in metrics])
    print(f'\n  Speed:       {sp.mean():.3f} ± {sp.std():.3f} µm/min')
    print(f'  Net disp:    {np.array([m["net"] for m in metrics]).mean():.1f} µm')
    print(f'  Directness:  {np.array([m["dir"] for m in metrics]).mean():.3f}')

    # Save
    print('\nSaving outputs...')
    save_csvs(metrics, all_det, um, dt, out_prefix)
    save_detection_image(
        tif_path, all_det, roi, um, dt, fl_ch, n_ch, fl_p99, fish_mask,
        out_prefix+'_cell_detection.png')
    save_thesis_figure(
        tif_path, all_det, metrics, roi, um, dt, fl_ch, n_ch, n_tp,
        fl_p99, fish_mask, out_prefix, mint)
    render_video(
        tif_path, all_det, metrics, roi, um, dt, fl_ch, n_ch, n_tp,
        fl_p99, fish_mask, out_prefix+'_motility_video.mp4',   # split panel
        scale=P['video_scale'], fps=P['video_fps'])
    # cell_detection_video.mp4 and cell_tracks_video.mp4 are auto-saved by render_video

    print(f'\nAll outputs saved to: {out_dir}')
    print('Done.')


if __name__ == '__main__':
    main()
