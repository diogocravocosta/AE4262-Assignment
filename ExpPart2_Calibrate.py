import os
import sys
import numpy as np
import struct
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d
import matplotlib.pyplot as plt

# ── path setup ────────────────────────────────────────────────────────────────
ROOT = r"C:\Users\Yasmine\OneDrive - Delft University of Technology\Desktop\Q3 Projects\Combustion\Assignment Documents\Combustion Assignment Data"
sys.path.append(ROOT)

cal_path = os.path.join(ROOT, "Calibration", "CH4", "B16")


#── read function ──────────────────────────────────────────────────────
def readB16(dirPath, filename=None):

    # Parse optional input
    if filename is not None:
        filename = os.path.join(dirPath, filename)
    else:
        filename = dirPath

    # Open the file
    try:
        f = open(filename, "rb")
    except OSError:
        raise FileNotFoundError(f"Could not open file: {filename}")

    with f:
        # Check that it is a PCO file
        filetype = f.read(4)
        if filetype != b"PCO-":
            try:
                ft = filetype.decode("ascii", errors="replace")
            except Exception:
                ft = str(filetype)
            raise ValueError(f"Wrong filetype: {ft}")

        # Get image dimensions:
        fileSize, headLength, imgWidth, imgHeight = struct.unpack("<4i", f.read(16))

        # Look into the extended header, throw error if color image
        extHeader = struct.unpack("<i", f.read(4))[0]
        if extHeader == -1:
            colorMode = struct.unpack("<i", f.read(4))[0]
            if colorMode != 0:
                raise ValueError(
                    "Color image detected. Only b/w images have been tested with this function"
                )

        # Get the image (match MATLAB column-major fread + transpose)
        f.seek(headLength, os.SEEK_SET)

        # Read exactly width*height pixels as uint16
        n = imgWidth * imgHeight
        data = np.fromfile(f, dtype=np.dtype("<u2"), count=n)
        if data.size != n:
            raise EOFError("File ended before all pixels were read.")

        # MATLAB: fread(...,[imgWidth,imgHeight],'uint16')'  -> (imgHeight x imgWidth)
        image = data.reshape((imgWidth, imgHeight), order="F").T
        print(image)
        print(image.shape)
    return image

# ── calibration function ──────────────────────────────────────────────────────
def get_scale(cal_image):
    results = {}
    fig, axes = plt.subplots(1, 2, figsize=(12, 8))

    for i, (axis, label) in enumerate([(0, "horizontal"), (1, "vertical")]):
        profile        = cal_image.mean(axis=axis).astype(float)
        profile_inv    = profile.max() - profile
        profile_smooth = uniform_filter1d(profile_inv, size=5)
        peaks, _       = find_peaks(profile_smooth, distance=5,
                                    prominence=profile_smooth.max() * 0.03)
        spacing_px     = np.median(np.diff(peaks))
        mm_per_px      = 1.0 / spacing_px   # major lines = 1 mm apart


        print(f"  {label}: spacing = {spacing_px:.1f} px → {mm_per_px:.5f} mm/px")
        results[label] = {"mm_per_px": mm_per_px, "peaks": peaks}

        axes[i].imshow(cal_image, cmap="gray")
        H, W = cal_image.shape
        axes[i].set_xlim(0, W)
        axes[i].set_ylim(H, 0)  # note: reversed because image y goes top→bottom

        for p in peaks:
            if axis == 0:
                axes[i].axhline(p, color="red", linewidth=0.8, alpha=0.7)
            else:
                axes[i].axvline(p, color="cyan", linewidth=0.8, alpha=0.7)
        axes[i].set_title(f"{label} lines — {spacing_px:.1f} px → {mm_per_px:.4f} mm/px")
        axes[i].axis("off")

    plt.tight_layout()
    plt.show()

    # average both axes for final scale
    mm_per_px = np.mean([results["horizontal"]["mm_per_px"],
                         results["vertical"]["mm_per_px"]])
    print(f"\n  Final scale (mean): {mm_per_px:.5f} mm/px")
    return mm_per_px, results

# ── run ───────────────────────────────────────────────────────────────────────
b16_files = [f for f in os.listdir(cal_path) if f.endswith(".b16")]
print(f"Found {len(b16_files)} calibration file(s)\n")

scales = []
for fname in b16_files:
    print(f"Processing: {fname}")
    img = readB16(cal_path, fname)
    mm_per_px, results = get_scale(img)   # ← unpack, no axis argument
    scales.append(mm_per_px)

final_scale = np.mean(scales)
print(f"\nFinal scale: {final_scale:.5f} mm/px  ({1 / final_scale:.2f} px/mm)")

np.save("scale_mm_per_px.npy", final_scale)

