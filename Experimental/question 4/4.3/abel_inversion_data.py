"""
Abel inversion of a PCO B16 image file — data-only, no plotting.

Usage as a module:
    from abel_inversion_data import process_b16
    img_bg, recon, vmax = process_b16("path/to/file.b16")

Usage from the command line (saves arrays to a .npz file):
    python abel_inversion_data.py <path_to_b16_file>

Requirements:
    pip install numpy PyAbel
"""

import os
import sys
import numpy as np
import abel
from readB16 import readB16


def process_b16(filepath):
    """
    Read a B16 file, subtract background, and apply Hansen-Law Abel inversion.

    Returns
    -------
    img_bg : ndarray
        Background-subtracted original image (float64, clipped ≥ 0).
    recon : ndarray
        Abel-inverted radial emissivity (same shape as img_bg).
    vmax : float
        99th-percentile of |recon|, suitable as the color scale upper limit.
    """
    img = readB16(filepath).astype(np.float64)

    bg = np.median(img[:50, :])
    img_bg = np.clip(img - bg, 0, None)

    # Centre the image on the flame axis before inverting.
    # PyAbel expects the symmetry axis at the image centre (col width//2).
    # Use centre-of-mass: the cone axis is dark so argmax lands on a cone wall.
    col_sums = img_bg.sum(axis=0)
    cols = np.arange(img_bg.shape[1])
    center_col = int(round(np.average(cols, weights=col_sums)))
    half_w = min(center_col, img_bg.shape[1] - center_col - 1)
    img_centered = img_bg[:, center_col - half_w : center_col + half_w + 1]

    result = abel.Transform(
        img_centered,
        direction="inverse",
        method="hansenlaw",
        symmetry_axis=0,
        verbose=False,
    )
    recon = result.transform
    vmax = float(np.percentile(np.abs(recon), 99))

    return img_bg, recon, vmax


def main():
    if len(sys.argv) < 2:
        print(f"Usage: python {os.path.basename(__file__)} <path_to_b16_file>")
        sys.exit(1)

    filepath = sys.argv[1]
    if not os.path.isfile(filepath):
        print(f"Error: file not found: {filepath}")
        sys.exit(1)

    img_bg, recon, vmax = process_b16(filepath)

    base = os.path.splitext(os.path.basename(filepath))[0]
    out_path = f"{base}_abel.npz"
    np.savez(out_path, img_bg=img_bg, recon=recon, vmax=vmax)
    print(f"Saved arrays to {out_path}  (keys: img_bg, recon, vmax)")


if __name__ == "__main__":
    main()
