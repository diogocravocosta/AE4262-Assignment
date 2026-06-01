import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# ── inputs ────────────────────────────────────────────────
U0      = 0.5    # inlet (average) velocity [m/s]
R       = 0.05   # pipe radius [m]
nu      = 1e-6   # kinematic viscosity [m²/s]  (water ~ 1e-6)
L_ratio = 40     # pipe length in diameters
iso_u   = 0.30   # isoline velocity [m/s]
# ──────────────────────────────────────────────────────────

D  = 2 * R
Re = U0 * D / nu
Le = 0.06 * Re * D   # hydrodynamic entrance length (Langhaar estimate)
L  = L_ratio * D

print(f"Re              = {Re:.0f}  {'(laminar ✓)' if Re < 2300 else '(WARNING: turbulent, model invalid)'}")
print(f"Entrance length = {Le:.4f} m  ({Le/D:.1f} diameters)")
print(f"Pipe length     = {L:.4f} m  ({L_ratio} diameters)")
print(f"u_max (CL, FD)  = {2*U0:.3f} m/s")

NX, NR = 500, 200
x = np.linspace(0, L, NX)
r = np.linspace(-R, R, NR)
X, Rv = np.meshgrid(x, r)


def velocity(x, r, R, U0, Le):
    """
    Approximate velocity field for laminar pipe flow including entrance region.

    Parameters
    ----------
    x  : axial position [m]
    r  : radial position [m]  (signed, 0 = centreline)
    R  : pipe radius [m]
    U0 : mean inlet velocity [m/s]
    Le : hydrodynamic entrance length [m]

    Returns
    -------
    u  : axial velocity [m/s]
    """
    xi   = np.clip(x / Le, 0.0, 1.0)   # normalised entrance coordinate
    rn   = np.abs(r) / R                # normalised radius
    uFD  = 2.0 * U0 * (1.0 - rn**2)    # fully-developed Hagen-Poiseuille
    uIn  = U0 * np.ones_like(rn)        # uniform plug profile at inlet
    blend = xi**0.5                     # smooth 0→1 transition
    return uIn * (1.0 - blend) + uFD * blend


U = velocity(X, Rv, R, U0, Le)

# ── plot ──────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(13, 7),
                         gridspec_kw={'height_ratios': [3, 1]})

# --- top: 2-D colour field ---
ax = axes[0]
cf = ax.contourf(X * 100, Rv * 100, U, levels=80, cmap='viridis')
cbar = plt.colorbar(cf, ax=ax, pad=0.01)
cbar.set_label('velocity  u  [m/s]', fontsize=10)

# isoline
if 0 < iso_u < U.max():
    cs = ax.contour(X * 100, Rv * 100, U,
                    levels=[iso_u], colors='white', linewidths=2.0)
    ax.clabel(cs, fmt=f'u = {iso_u:.2f} m/s', fontsize=9, colors='white')
else:
    print(f"Note: iso_u = {iso_u} m/s is outside the velocity range "
          f"[0, {U.max():.3f}] — no isoline drawn.")

# entrance-length marker
if Le < L:
    ax.axvline(Le * 100, color='white', lw=1.2, ls='--', alpha=0.75)
    ax.text(Le * 100 + L * 0.005 * 100, 0,
            'fully developed →', color='white', fontsize=8,
            va='center', ha='left')

ax.set_xlabel('axial position  x  [cm]', fontsize=10)
ax.set_ylabel('radial position  r  [cm]', fontsize=10)
ax.set_title(
    f'Laminar pipe flow velocity field\n'
    f'Re = {Re:.0f}   U₀ = {U0} m/s   R = {R*100:.1f} cm   '
    f'ν = {nu:.2e} m²/s',
    fontsize=11)
ax.set_aspect('auto')

# --- bottom: axial profiles at selected x stations ---
ax2 = axes[1]
stations = [0.0, 0.1, 0.25, 0.5, 1.0]  # fractions of pipe length
r_plot = np.linspace(-R, R, 300)
colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(stations)))

for frac, col in zip(stations, colors):
    xi_pos = frac * L
    u_prof = velocity(xi_pos, r_plot, R, U0, Le)
    lbl = f'x = {frac:.0%}·L'
    ax2.plot(r_plot * 100, u_prof, color=col, lw=1.8, label=lbl)

if 0 < iso_u < 2 * U0:
    ax2.axhline(iso_u, color='gray', ls=':', lw=1.2, label=f'isoline ({iso_u:.2f} m/s)')

ax2.set_xlabel('r  [cm]', fontsize=10)
ax2.set_ylabel('u  [m/s]', fontsize=10)
ax2.set_title('Velocity profiles at selected axial stations', fontsize=10)
ax2.legend(fontsize=8, ncol=len(stations) + 1, loc='lower center')
ax2.set_xlim(-R * 100, R * 100)
ax2.set_ylim(0, 2.1 * U0)

plt.tight_layout()
plt.savefig('pipe_flow.png', dpi=150, bbox_inches='tight')
print("Figure saved to pipe_flow.png")
plt.show()