"""
Abel inversion of a PCO B16 image file.

Usage:
    python abel_b16.py <path_to_b16_file>

Requirements:
    pip install numpy matplotlib PyAbel
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import abel
from readB16 import readB16


def main():
    if len(sys.argv) < 2:
        print(f"Usage: python {os.path.basename(__file__)} <path_to_b16_file>")
        sys.exit(1)

    filepath = sys.argv[1]
    if not os.path.isfile(filepath):
        print(f"Error: file not found: {filepath}")
        sys.exit(1)

    # --- Read ---
    print(f"Reading {filepath} ...")
    img = readB16(filepath).astype(np.float64)
    print(f"  Image size: {img.shape[1]} x {img.shape[0]}, "
          f"range: [{img.min():.0f}, {img.max():.0f}]")

    # --- Background subtraction ---
    bg = np.median(img[:50, :])
    img_bg = np.clip(img - bg, 0, None)
    print(f"  Background estimate: {bg:.1f}")

    # --- Centre on flame axis (PyAbel expects symmetry axis at image centre) ---
    # Use centre-of-mass: for a cone flame the axis is dark so argmax would
    # land on a cone wall instead of the true centreline.
    col_sums = img_bg.sum(axis=0)
    cols = np.arange(img_bg.shape[1])
    center_col = int(round(np.average(cols, weights=col_sums)))
    half_w = min(center_col, img_bg.shape[1] - center_col - 1)
    img_centered = img_bg[:, center_col - half_w : center_col + half_w + 1]

    # --- Abel inversion ---
    print("Running Abel inversion (Hansen-Law) ...")
    result = abel.Transform(
        img_centered,
        direction="inverse",
        method="hansenlaw",
        symmetry_axis=0,
        verbose=False,
    )
    recon = result.transform
    print(f"  Result range: [{recon.min():.2f}, {recon.max():.2f}]")

    # --- Plot ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))

    im1 = ax1.imshow(img_centered, cmap="hot", aspect="auto")
    ax1.set_title("Original (bg subtracted)")
    plt.colorbar(im1, ax=ax1, shrink=0.5)

    vmax = np.percentile(np.abs(recon), 99)
    im2 = ax2.imshow(recon, cmap="hot", aspect="auto", vmin=0, vmax=vmax)
    ax2.set_title("Abel Inversion (Hansen-Law)")
    plt.colorbar(im2, ax=ax2, shrink=0.5)

    base = os.path.splitext(os.path.basename(filepath))[0]
    fig.suptitle(f"{base} — Abel Inversion", fontsize=14, fontweight="bold")
    plt.tight_layout()

    out_path = f"{base}_abel.png"
    plt.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")
    plt.show()


if __name__ == "__main__":
    main()