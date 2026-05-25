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
def get_scale(cal_image, out_filename):
    results = {}
    
    # 1. Create the figure with exact pixel dimensions
    H, W = cal_image.shape
    dpi = 100
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    
    # Fill the entire canvas, no borders
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    
    # Render the base image with pixel-perfect interpolation
    ax.imshow(cal_image, cmap="gray", interpolation="nearest")

    # 2. Process both axes and overlay lines
    # axis=0 -> averages rows -> array of size W (X-coordinates) -> Vertical lines
    # axis=1 -> averages cols -> array of size H (Y-coordinates) -> Horizontal lines
    for axis, label in [(0, "vertical"), (1, "horizontal")]:
        profile        = cal_image.mean(axis=axis).astype(float)
        profile_inv    = profile.max() - profile
        profile_smooth = uniform_filter1d(profile_inv, size=5)
        peaks, _       = find_peaks(profile_smooth, distance=5,
                                    prominence=profile_smooth.max() * 0.03)
        spacing_px     = np.median(np.diff(peaks))
        mm_per_px      = 1.0 / spacing_px   # major lines = 1 mm apart

        print(f"  {label}: spacing = {spacing_px:.1f} px → {mm_per_px:.5f} mm/px")
        results[label] = {"mm_per_px": mm_per_px, "peaks": peaks}

        # Overlay the lines directly onto the pixel-perfect axes
        # Using a thin line (linewidth=0.5) so it doesn't obscure the pixels beneath
        for p in peaks:
            if axis == 0:
                ax.axvline(p, color="cyan", linewidth=0.5, alpha=0.8)  # X-coords
            else:
                ax.axhline(p, color="red", linewidth=0.5, alpha=0.8)   # Y-coords

    # 3. Save the image natively and cleanly close the figure
    plt.savefig(out_filename, dpi=dpi, pad_inches=0)
    plt.close(fig) 

    # 4. Average both axes for final scale
    mm_per_px_mean = np.mean([results["horizontal"]["mm_per_px"],
                              results["vertical"]["mm_per_px"]])
    print(f"  Final scale (mean): {mm_per_px_mean:.5f} mm/px\n")
    
    return mm_per_px_mean, results

# ── run ───────────────────────────────────────────────────────────────────────
b16_files = [f for f in os.listdir(cal_path) if f.endswith(".b16")]
print(f"Found {len(b16_files)} calibration file(s)\n")

scales = []
for fname in b16_files:
    print(f"Processing: {fname}")
    img = readB16(cal_path, fname)
    
    # Generate an output filename so each image saves uniquely
    out_file = os.path.join(cal_path, fname.replace(".b16", "_pixel_perfect.png"))
    
    # Pass the out_file to the function
    mm_per_px, results = get_scale(img, out_file)
    scales.append(mm_per_px)

final_scale = np.mean(scales)
print(f"=========================================")
print(f"Total Final scale: {final_scale:.5f} mm/px  ({1 / final_scale:.2f} px/mm)")
print(f"=========================================")

np.save("scale_mm_per_px.npy", final_scale)