"""
Display OH* Raw, OH* Abel-Inverted, and OH-PLIF for each equivalence ratio.

Usage:
    python display_phi_trio.py

Configure the two base paths below before running.
readB16.py must be in the same directory.

Requirements:
    pip install numpy matplotlib PyAbel
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import abel
from readB16 import readB16

# ─── CONFIGURE THESE ────────────────────────────────────────────────
BASE_OHSTAR = r"C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 4\OHstar"
BASE_PLIF   = r"C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\LIF"
PHI_VALUES  = np.arange(0.5, 1.4, 0.1)
# ────────────────────────────────────────────────────────────────────


def abel_invert(img):
    """Center on flame axis, then inverse Abel transform (Hansen-Law)."""
    col_sums = img.sum(axis=0)
    center_col = int(np.round(np.average(np.arange(img.shape[1]), weights=col_sums)))
    half_w = min(center_col, img.shape[1] - center_col - 1)
    img_centered = img[:, center_col - half_w : center_col + half_w + 1]
    return abel.Transform(
        img_centered, direction="inverse", method="hansenlaw",
        symmetry_axis=0, verbose=False,
    ).transform


OUT_DIR = os.path.join(os.path.dirname(__file__), "Phi_trio")


def plot_image(img, title, save_path):
    """Save a single image with white theme and per-image color scale."""
    fig, ax = plt.subplots(figsize=(5, 9))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    vmin = 0
    vmax = np.percentile(img[img > 0], 99) if np.any(img > 0) else 1
    im = ax.imshow(img, cmap="hot", aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_title(title, color="black", fontsize=12)
    ax.axis("off")
    cb = fig.colorbar(im, ax=ax, shrink=0.6)
    cb.ax.yaxis.set_tick_params(color="black")
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="black")
    plt.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    for phi in PHI_VALUES:
        tag = f"H2_PHI{phi:.1f}"
        oh_path   = os.path.join(BASE_OHSTAR, tag, "B16", "B0001.b16")
        plif_path = os.path.join(BASE_PLIF,   tag, "B16", "B0001.b16")

        oh  = readB16(oh_path).astype(np.float64)
        lif = readB16(plif_path).astype(np.float64)
        inv = abel_invert(oh)

        slug = f"phi{phi:.1f}"
        plot_image(oh,  f"OH* Raw  \u03c6={phi:.1f}",           os.path.join(OUT_DIR, f"{slug}_OH_raw.png"))
        plot_image(inv, f"OH* Abel Inverted  \u03c6={phi:.1f}", os.path.join(OUT_DIR, f"{slug}_OH_abel.png"))
        plot_image(lif, f"OH-PLIF  \u03c6={phi:.1f}",           os.path.join(OUT_DIR, f"{slug}_PLIF.png"))

        print(f"  \u03c6 = {phi:.1f} \u2713  \u2192 saved to Phi_trio/")


if __name__ == "__main__":
    main()
