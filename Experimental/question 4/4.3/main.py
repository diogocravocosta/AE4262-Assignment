import os
import numpy as np
import matplotlib.pyplot as plt
from readB16 import readB16
from New_abel_imversion import abel_invert_b16

# --- Settings ---
SAVE_FIGURES = True
SAVE_DIR     = os.path.join(os.path.dirname(__file__), "Trio figures")

BASE_OHSTAR = r"C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 4\OHstar"
BASE_PLIF   = r"C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\LIF"

phi_values = np.arange(0.5, 1.4, 0.1)

def get_ohstar_paths():
    return [fr"{BASE_OHSTAR}\H2_PHI{phi:.1f}\B16\B0001.b16" for phi in phi_values]

def get_plif_paths():
    return [fr"{BASE_PLIF}\H2_PHI{phi:.1f}\B16\B0001.b16" for phi in phi_values]

# --- Load images ---
ohstar_paths = get_ohstar_paths()
plif_paths   = get_plif_paths()

ohstar_raw = [readB16(p).astype(np.float64) for p in ohstar_paths]
recon_all  = [abel_invert_b16(p)            for p in ohstar_paths]
plif_raw   = [readB16(p).astype(np.float64) for p in plif_paths]

# --- Shared color limits per image type ---
def clim(images):
    return np.min([img.min() for img in images]), np.max([img.max() for img in images])

ohstar_clim = clim(ohstar_raw)
recon_clim  = clim(recon_all)
plif_clim   = clim(plif_raw)

if SAVE_FIGURES:
    os.makedirs(SAVE_DIR, exist_ok=True)

# --- One figure per phi, 3 subplots ---
for i, phi in enumerate(phi_values):
    fig, axes = plt.subplots(1, 3, figsize=(15, 8))
    fig.patch.set_facecolor('white')
    fig.suptitle(f'φ = {phi:.1f}', color='black', fontsize=14)

    data = [
        (ohstar_raw[i], ohstar_clim, f'OHstar Raw  [φ={phi:.1f}]'),
        (recon_all[i],  recon_clim,  f'OHstar After Inversion  [φ={phi:.1f}]'),
        (plif_raw[i],   plif_clim,   f'PLIF  [φ={phi:.1f}]'),
    ]

    for ax, (img, (vmin, vmax), title) in zip(axes, data):
        ax.set_facecolor('white')
        im = ax.imshow(img, cmap='hot', aspect='auto', vmin=vmin, vmax=vmax)
        ax.set_title(title, color='black')
        ax.axis('off')
        cb = fig.colorbar(im, ax=ax, shrink=0.6)
        cb.ax.yaxis.set_tick_params(color='black')
        plt.setp(cb.ax.yaxis.get_ticklabels(), color='black')

    plt.tight_layout()

    if SAVE_FIGURES:
        fig.savefig(os.path.join(SAVE_DIR, f"phi_{phi:.1f}.png"), dpi=150, bbox_inches='tight')

plt.show()
