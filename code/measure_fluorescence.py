import tifffile
import numpy as np
import pandas as pd
from skimage import io, filters, morphology, measure

# Parameters
tif_path = "time_lapse_stack.tif"  # multi-page tiff with frames (or z-stacks per time -> adapt)
output_csv = "mean_intensity_vs_frame.csv"

# Load stack (frames, y, x) — if 3D per frame, you'll need to reshape/maximum-project
stack = tifffile.imread(tif_path)  # shape: (n_frames, height, width) or (n_time, z, y, x)

# If you have z-stacks per time, convert to max projection per time:
if stack.ndim == 4:
    # stack shape (n_time, z, y, x)
    stack = stack.max(axis=1)  # now (n_time, y, x)

n_frames = stack.shape[0]
means = []

for i in range(n_frames):
    frame = stack[i].astype(np.float32)

    # Simple background removal + threshold for a rough mask (replace with segmentation)
    bg = filters.gaussian(frame, sigma=10)
    frame_corr = frame - bg
    thresh = filters.threshold_otsu(frame_corr.clip(min=0))
    mask = frame_corr > thresh
    mask = morphology.remove_small_objects(mask, min_size=64)
    mask = morphology.binary_closing(mask, morphology.disk(3))

    # Measure mean intensity within mask
    if mask.sum() == 0:
        mean_intensity = np.nan
    else:
        mean_intensity = frame[mask].mean()

    means.append(mean_intensity)

# Save results
df = pd.DataFrame({"frame": np.arange(n_frames), "mean_intensity": means})
df.to_csv(output_csv, index=False)
print("Saved", output_csv)
