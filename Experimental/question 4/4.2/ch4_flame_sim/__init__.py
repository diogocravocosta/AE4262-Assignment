"""
ch4_flame_sim/__init__.py
=========================
Cantera 1-D freely-propagating CH₄/air flame simulation
with RADIS emission spectrum synthesis.

PURPOSE
-------
Runs 1-D freely-propagating premixed CH₄/air flames at φ = 0.8, 1.0, 1.2
using the GRI-Mech 3.0 mechanism (53 species, 325 reactions – bundled with
Cantera). For each φ it:

  1. Solves the flame and extracts spatial profiles of T and all species.
  2. Plots temperature + key species mole fractions vs. distance from the
     reaction zone.
  3. Synthesises a line-by-line emission spectrum (300–1000 nm) for the
     burnt-gas composition using RADIS, combining contributions from the
     major radiating species present in methane combustion:
       OH, H2O, CO, CO2, CH (if present), NO

DEPENDENCIES
------------
    pip install cantera matplotlib numpy radis

MECHANISM
---------
GRI-Mech 3.0 ('gri30.yaml') – always bundled with Cantera.
Contains all relevant CH₄ combustion species including CH, CH2, CH2O,
C2H2, C2H4, NO, N2O, HCN etc.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ── output directory ──────────────────────────────────────────────────────────
RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── simulation settings ───────────────────────────────────────────────────────

#: Equivalence ratios
PHI_VALUES: List[float] = [0.8, 1.0, 1.2]

#: GRI-Mech 3.0 – comprehensive CH4 mechanism, always bundled with Cantera
MECHANISM: str = "gri30.yaml"

#: Fuel / oxidiser
FUEL: str = "CH4"
OXIDISER: str = "O2:1, N2:3.76"

#: Inlet conditions
T_INLET: float = 300.0    # [K]
P_INLET: float = 101325.0 # [Pa]  1 atm

#: Initial domain width [m]
DOMAIN_WIDTH: float = 0.05  # 5 cm – wider than H2 flames

#: Spatial coordinate: 'x_shifted' | 'x_raw' | 'c'
COORD: str = "x_shifted"

#: Species to track in profiles – key CH4 combustion species
SPECIES_TO_PLOT: List[str] = [
    "CH4", "O2", "H2O", "CO2", "CO",
    "OH",  "H",  "O",   "H2",  "HO2",
    "CH2O","CH", "C2H2","NO",
]

#: Species RADIS can synthesise spectra for (subset with HITRAN/HITEMP data)
#: These will be attempted; any that RADIS can't handle are skipped gracefully.
SPECTRO_SPECIES: List[str] = ["OH", "H2O", "CO", "CO2", "NO"]

#: Wavelength range for emission spectrum [nm]
WL_MIN: float = 300.0
WL_MAX: float = 750.0
WL_STEP: float = 0.05   # nm  (fine enough to resolve main bands)

#: Scale factors for plotting (species with very small mole fractions)
SCALE: Dict[str, int] = {
    "CH4": 1, "O2": 1, "H2O": 1, "CO2": 1, "CO": 1,
    "OH": 50, "H": 100, "O": 100, "H2": 1,
    "HO2": 500, "CH2O": 200, "CH": 1000, "C2H2": 200, "NO": 500,
}

#: Figure DPI
FIG_DPI: int = 150

# ── colours for species ───────────────────────────────────────────────────────
COLOURS: Dict[str, str] = {
    "CH4":  "#1f77b4",
    "O2":   "#9467bd",
    "H2O":  "#aec7e8",
    "CO2":  "#ffbb78",
    "CO":   "#d62728",
    "OH":   "#2ca02c",
    "H":    "#17becf",
    "O":    "#e377c2",
    "H2":   "#ff7f0e",
    "HO2":  "#8c564b",
    "CH2O": "#bcbd22",
    "CH":   "#7f7f7f",
    "C2H2": "#e8a838",
    "NO":   "#c49c94",
}


# ── data container ────────────────────────────────────────────────────────────

@dataclass
class SimResult:
    """Container for one φ simulation result."""
    phi: float
    grid: np.ndarray            # raw Cantera grid [m]
    T: np.ndarray               # temperature [K]
    X: Dict[str, np.ndarray]    # mole fractions keyed by species name
    hrr: np.ndarray             # heat release rate [W/m³]
    sl: float                   # laminar flame speed [m/s]
    x_flame: float              # axial position of max-HRR [m]
    coord_values: np.ndarray = field(default_factory=lambda: np.array([]))
    coord_label: str = ""

    def to_dataframe(self):
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("Install pandas to use to_dataframe().")
        data = {
            "x_m": self.grid,
            "x_shifted_m": self.grid - self.x_flame,
            "T_K": self.T,
            "HRR_W_m3": self.hrr,
        }
        for sp, xsp in self.X.items():
            data[f"X_{sp}"] = xsp
        return pd.DataFrame(data)

    def save_csv(self, path: Optional[Path] = None) -> Path:
        if path is None:
            path = RESULTS_DIR / f"phi_{self.phi:.1f}.csv"
        self.to_dataframe().to_csv(path, index=False)
        return path


# ── core simulation ───────────────────────────────────────────────────────────

def _load_gas():
    """Load GRI-Mech 3.0 gas object."""
    import cantera as ct
    return ct.Solution(MECHANISM)


def simulate_phi(phi: float, loglevel: int = 0) -> SimResult:
    """
    Run a 1-D freely-propagating CH₄/air flame at equivalence ratio *phi*.

    Parameters
    ----------
    phi : float
        Equivalence ratio (0.8, 1.0, 1.2 for the methane range).
    loglevel : int
        Cantera solver verbosity (0 = silent).

    Returns
    -------
    SimResult
    """
    import cantera as ct

    gas = _load_gas()
    gas.set_equivalence_ratio(phi, FUEL, OXIDISER)
    gas.TP = T_INLET, P_INLET

    flame = ct.FreeFlame(gas, width=DOMAIN_WIDTH)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
    flame.set_max_grid_points(flame.flame, 1000)
    flame.solve(loglevel=loglevel, auto=True)

    grid = flame.grid.copy()
    T    = flame.T.copy()
    hrr  = flame.heat_release_rate.copy()
    sl   = float(flame.velocity[0])

    X: Dict[str, np.ndarray] = {}
    for sp in SPECIES_TO_PLOT:
        try:
            idx = gas.species_index(sp)
            X[sp] = flame.X[idx].copy()
        except ValueError:
            X[sp] = np.zeros_like(grid)

    i_flame = int(np.argmax(hrr))
    x_flame = grid[i_flame]

    result = SimResult(
        phi=phi, grid=grid, T=T, X=X,
        hrr=hrr, sl=sl, x_flame=x_flame,
    )
    result.coord_values, result.coord_label = _make_coord(result)

    # ── burnt gas summary ─────────────────────────────────────────────────────
    i_post  = int(0.8 * len(grid))
    T_burnt = float(np.mean(T[i_post:]))
    print(f"    T_burnt   = {T_burnt:.0f} K")
    for sp in ["OH", "CH", "C2H2"]:
        xsp = X.get(sp, np.zeros_like(grid))
        print(f"    X({sp:<4s})  peak = {xsp.max():.4e}   burnt gas = {float(np.mean(xsp[i_post:])):.4e}")

    return result


def _make_coord(res: SimResult) -> Tuple[np.ndarray, str]:
    if COORD == "x_shifted":
        return (res.grid - res.x_flame) * 1e3, "Distance from reaction zone [mm]"
    elif COORD == "x_raw":
        return res.grid * 1e3, "Axial position x [mm]"
    elif COORD == "c":
        x_h2o   = res.X.get("H2O", np.zeros_like(res.grid))
        x_h2o_u = x_h2o[0]
        x_h2o_b = x_h2o.max()
        if (x_h2o_b - x_h2o_u) < 1e-10:
            warnings.warn("H2O range too small; falling back to x_shifted.", UserWarning)
            return (res.grid - res.x_flame) * 1e3, "x − x_flame [mm]"
        c = (x_h2o - x_h2o_u) / (x_h2o_b - x_h2o_u)
        return c, r"Progress variable $c$ [-]"
    else:
        raise ValueError(f"Unknown COORD: '{COORD}'")


# ── flame profile plots ───────────────────────────────────────────────────────

def plot_single_phi(res: SimResult, save: bool = True):
    """Full-domain species + temperature profile for one φ."""
    import matplotlib.pyplot as plt

    fig, ax_sp = plt.subplots(figsize=(10, 5))
    ax_T = ax_sp.twinx()
    xc = res.coord_values

    for sp, xsp in res.X.items():
        f = SCALE.get(sp, 1)
        if np.max(xsp) * f < 1e-8:
            continue
        label = f"X({sp}) ×{f}" if f > 1 else f"X({sp})"
        ax_sp.plot(xc, xsp * f,
                   color=COLOURS.get(sp, "grey"),
                   linewidth=2.0 if sp in ("OH", "CO", "CH4") else 1.3,
                   label=label)

    ax_T.plot(xc, res.T, color="black", linewidth=2.0,
              linestyle="--", label="T [K]")

    ax_sp.set_xlabel(res.coord_label, fontsize=11)
    ax_sp.set_ylabel("Mole fraction [-] (scaled)", fontsize=11)
    ax_sp.set_ylim(bottom=0)
    ax_T.set_ylabel("Temperature [K]", fontsize=11)
    ax_T.set_ylim(bottom=0)
    ax_sp.set_title(
        f"CH₄/air freely-propagating flame  |  φ = {res.phi:.1f}  |"
        f"  SL = {res.sl*100:.1f} cm/s",
        fontsize=11,
    )
    l1, b1 = ax_sp.get_legend_handles_labels()
    l2, b2 = ax_T.get_legend_handles_labels()
    ax_sp.legend(l1 + l2, b1 + b2, loc="upper left",
                 fontsize=7, ncol=3, framealpha=0.85)
    ax_sp.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"flame_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_reaction_zone(res: SimResult, window_mm: float = 5.0, save: bool = True):
    """Zoomed ±window_mm plot centred on the reaction zone."""
    import matplotlib.pyplot as plt

    fig, ax_sp = plt.subplots(figsize=(9, 5))
    ax_T = ax_sp.twinx()

    xc   = res.coord_values
    mask = (xc >= -window_mm) & (xc <= window_mm)

    # In the zoomed view only plot the radical / intermediate species
    zoom_species = ["CH4", "O2", "CO", "OH", "H", "O", "HO2", "CH2O", "CH", "NO"]

    for sp in zoom_species:
        xsp = res.X.get(sp, np.zeros_like(xc))
        f   = SCALE.get(sp, 1)
        if np.max(xsp[mask]) * f < 1e-8:
            continue
        label = f"X({sp}) ×{f}" if f > 1 else f"X({sp})"
        ax_sp.plot(xc[mask], xsp[mask] * f,
                   color=COLOURS.get(sp, "grey"),
                   linewidth=2.0 if sp == "OH" else 1.4,
                   label=label)

    ax_T.plot(xc[mask], res.T[mask], color="black",
              linewidth=2.0, linestyle="--", label="T [K]")
    ax_sp.axvline(0.0, color="grey", linestyle=":", linewidth=1.0)

    ax_sp.set_xlabel(res.coord_label, fontsize=11)
    ax_sp.set_ylabel("Mole fraction [-] (scaled)", fontsize=11)
    ax_sp.set_ylim(bottom=0)
    ax_T.set_ylabel("Temperature [K]", fontsize=11)
    ax_T.set_ylim(bottom=0)
    ax_sp.set_title(
        f"Reaction zone  |  φ = {res.phi:.1f}  |  ±{window_mm} mm",
        fontsize=11)

    l1, b1 = ax_sp.get_legend_handles_labels()
    l2, b2 = ax_T.get_legend_handles_labels()
    ax_sp.legend(l1 + l2, b1 + b2, loc="upper left",
                 fontsize=8, ncol=2, framealpha=0.85)
    ax_sp.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"zoom_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_phi_sweep(results: List[SimResult], save: bool = True):
    """Overlay of OH and CO mole fractions + HRR for all φ."""
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
    ax_oh, ax_co, ax_hrr = axes

    cmap     = cm.plasma
    phi_vals = [r.phi for r in results]
    norm     = plt.Normalize(vmin=min(phi_vals), vmax=max(phi_vals))

    for res in results:
        colour = cmap(norm(res.phi))
        xc = res.coord_values
        ax_oh.plot(xc, res.X.get("OH",  np.zeros_like(xc)) * SCALE.get("OH",  1),
                   color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_co.plot(xc, res.X.get("CO",  np.zeros_like(xc)),
                   color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_hrr.plot(xc, res.hrr / 1e9,
                    color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")

    for ax, ylabel, title in [
        (ax_oh,  f"X(OH) ×{SCALE['OH']}  [-]",      "OH radical"),
        (ax_co,  "X(CO)  [-]",                       "CO intermediate"),
        (ax_hrr, "Heat release rate  [GW/m³]",       "Heat release rate"),
    ]:
        ax.axvline(0.0, color="grey", linestyle=":", linewidth=1.0)
        ax.set_xlabel(results[0].coord_label, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "CH₄/air premixed flames — φ sweep  |  1-D Cantera FreeFlame (GRI-Mech 3.0)",
        fontsize=11,
    )
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / "phi_sweep.png"
        fig.savefig(fname, dpi=FIG_DPI)
        print(f"  Saved → {fname}")

    return fig


# ── spectroscopy ──────────────────────────────────────────────────────────────

def _burnt_gas_composition(res: SimResult) -> Dict[str, float]:
    """
    Extract the mean post-flame (burnt gas) mole fractions for spectroscopy.
    Averages over the last 20% of the domain (fully burnt side).
    Returns a dict {species: mole_fraction} for species with X > 1e-8.
    """
    n_pts  = len(res.grid)
    i_post = int(0.8 * n_pts)   # start of post-flame averaging window
    composition = {}
    for sp, xsp in res.X.items():
        x_mean = float(np.mean(xsp[i_post:]))
        if x_mean > 1e-8:
            composition[sp] = x_mean
    return composition


def compute_emission_spectrum(res: SimResult) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    """
    Compute a synthetic line-by-line emission spectrum of the burnt gas
    using RADIS for the radiating species present.

    Compatible with RADIS 0.17+.
    Key API differences from older versions:
      - 'species'  instead of 'molecule'
      - pressure in bar  (not atm)
      - no 'quantity' argument; emission is retrieved via s.get('emisscoeff')
        after the spectrum is computed in default (transmittance+radiance) mode

    Parameters
    ----------
    res : SimResult

    Returns
    -------
    wl : np.ndarray
        Wavelength axis [nm].
    total : np.ndarray
        Total integrated emission intensity [a.u.].
    per_species : dict
        Individual species contributions {species_name: intensity_array}.
    """
    try:
        from radis import calc_spectrum
    except ImportError:
        raise ImportError(
            "Install RADIS to compute emission spectra:\n"
            "    pip install radis"
        )

    composition = _burnt_gas_composition(res)
    T_burnt     = float(np.mean(res.T[int(0.8 * len(res.T)):]))  # mean post-flame T [K]
    P_bar       = P_INLET / 1e5                                   # Pa → bar

    wl = np.arange(WL_MIN, WL_MAX, WL_STEP)
    total       = np.zeros_like(wl)
    per_species = {}

    print(f"    Burnt gas T = {T_burnt:.0f} K,  P = {P_bar:.3f} bar")

    for sp in SPECTRO_SPECIES:
        x_sp = composition.get(sp, 0.0)
        if x_sp < 1e-8:
            print(f"    {sp}: not present in burnt gas – skipping")
            continue

        try:
            s = calc_spectrum(
                wmin         = WL_MIN,
                wmax         = WL_MAX,
                wunit        = "nm",
                species      = sp,
                isotope      = "1",
                pressure     = P_bar,
                Tgas         = T_burnt,
                mole_fraction= x_sp,
                path_length  = 1.0,          # cm
                databank     = "hitemp",
                cutoff       = 0,            # keep ALL lines, no strength threshold
                verbose      = 0,
            )

            # Retrieve radiance (emission) from the Spectrum object.
            # 'radiance_noslit' is the emitted spectral radiance [W/cm²/sr/nm]
            # and is always computed in equilibrium mode.
            wl_rad, I_rad = s.get("radiance_noslit", wunit="nm")
            I_interp = np.interp(wl, wl_rad, I_rad, left=0.0, right=0.0)
            per_species[sp] = I_interp
            total += I_interp
            print(f"    {sp}: peak = {I_interp.max():.3e}  (X={x_sp:.4f})")

        except Exception as exc:
            print(f"    {sp}: RADIS failed ({exc}) – skipping")

    return wl, total, per_species


def plot_emission_spectrum(res: SimResult, save: bool = True):
    """
    Plot the synthetic emission spectrum of the burnt gas at one φ.
    Shows total spectrum plus individual species contributions stacked below.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    print(f"  Computing emission spectrum for φ = {res.phi:.1f} …")
    wl, total, per_species = compute_emission_spectrum(res)

    n_sp   = len(per_species)
    n_rows = max(2, 1 + n_sp)  # at least 2 rows so sharex works
    fig, axes = plt.subplots(n_rows, 1,
                             figsize=(12, 2.5 * n_rows),
                             sharex=True)
    if n_rows == 1:
        axes = [axes]

    # ── top panel: total spectrum ─────────────────────────────────────────────
    ax0 = axes[0]
    total_norm = total / total.max() if total.max() > 0 else total
    ax0.fill_between(wl, total_norm, alpha=0.15, color="white")
    ax0.plot(wl, total_norm, color="white", linewidth=0.8)
    ax0.set_xlim(WL_MIN, WL_MAX)

    # Colour the background by approximate visible wavelength
    _shade_visible(ax0, wl)

    ax0.set_ylabel("Norm. intensity [-]", fontsize=10)
    ax0.set_title(
        f"Synthetic emission spectrum  |  CH₄/air  |  φ = {res.phi:.1f}  |"
        f"  T_burnt ≈ {float(np.mean(res.T[int(0.8*len(res.T)):])):.0f} K",
        fontsize=11,
    )
    ax0.set_ylim(0, 1.15)
    ax0.grid(True, alpha=0.25, color="grey")

    # ── per-species panels ────────────────────────────────────────────────────
    sp_colours = {"OH": "#2ca02c", "H2O": "#1f77b4", "CO": "#d62728",
                  "CO2": "#ff7f0e", "NO": "#c49c94", "CH": "#7f7f7f"}

    for i, (sp, I) in enumerate(per_species.items(), start=1):
        ax = axes[i]
        col = sp_colours.get(sp, "grey")
        ax.fill_between(wl, I / I.max() if I.max() > 0 else I,
                        alpha=0.35, color=col)
        ax.plot(wl, I / I.max() if I.max() > 0 else I,
                color=col, linewidth=0.9, label=sp)
        ax.set_ylabel(f"X({sp})\nnorm. [-]", fontsize=9)
        ax.set_ylim(0, 1.15)
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(True, alpha=0.25)
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    axes[-1].set_xlabel("Wavelength [nm]", fontsize=11)

    fig.suptitle(
        f"RADIS line-by-line emission  |  HITRAN database  |  P = {P_INLET/101325:.1f} atm",
        fontsize=10, y=1.001,
    )
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"spectrum_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight",
                    facecolor="black")
        print(f"  Saved → {fname}")

    return fig


def _shade_visible(ax, wl: np.ndarray):
    """
    Shade the background of the spectrum plot by wavelength colour.
    UV (< 380 nm) is dark purple, visible (380–780 nm) is rainbow,
    near-IR (> 780 nm) is dark red.
    """
    from matplotlib.colors import ListedColormap

    all_wl = np.linspace(wl[0], wl[-1], 800)
    colors = [_wl_to_rgb(w) for w in all_wl]
    cmap   = ListedColormap(colors)

    ax.imshow(
        np.linspace(0, 1, 800).reshape(1, -1),
        aspect="auto",
        extent=[wl[0], wl[-1], -0.02, 0.0],
        cmap=cmap,
        transform=ax.transData,
        zorder=3,
    )
    ax.set_facecolor("black")


def _wl_to_rgb(wl: float) -> Tuple[float, float, float]:
    """Approximate wavelength (nm) → (R, G, B) in [0,1]."""
    if 380 <= wl < 440:
        r, g, b = (440 - wl) / 60, 0.0, 1.0
    elif 440 <= wl < 490:
        r, g, b = 0.0, (wl - 440) / 50, 1.0
    elif 490 <= wl < 510:
        r, g, b = 0.0, 1.0, (510 - wl) / 20
    elif 510 <= wl < 580:
        r, g, b = (wl - 510) / 70, 1.0, 0.0
    elif 580 <= wl < 645:
        r, g, b = 1.0, (645 - wl) / 65, 0.0
    elif 645 <= wl <= 780:
        r, g, b = 1.0, 0.0, 0.0
    else:
        r, g, b = 0.0, 0.0, 0.0
    # Intensity falloff at UV and deep red edges
    if wl < 420:
        factor = 0.3 + 0.7 * (wl - 380) / 40
    elif wl > 700:
        factor = 0.3 + 0.7 * (780 - wl) / 80
    else:
        factor = 1.0
    return (r * factor, g * factor, b * factor)


# ── public API ────────────────────────────────────────────────────────────────

def run_all_phi(
    phi_values: Optional[List[float]] = None,
    loglevel: int = 0,
    save_csv: bool = False,
) -> List[SimResult]:
    """
    Run simulations for all equivalence ratios, produce flame profile plots
    and synthetic emission spectra.
    """
    if phi_values is None:
        phi_values = PHI_VALUES

    results: List[SimResult] = []
    print(f"Running CH₄/air FreeFlame simulations for φ = {phi_values}")
    print(f"Mechanism : {MECHANISM}")
    print(f"Coordinate: {COORD}")
    print("-" * 60)

    for phi in phi_values:
        print(f"\n  φ = {phi:.1f} … ", end="", flush=True)
        try:
            res = simulate_phi(phi, loglevel=loglevel)
            print(f"SL = {res.sl * 100:.2f} cm/s")
            plot_single_phi(res, save=True)
            plot_reaction_zone(res, save=True)
            plot_emission_spectrum(res, save=True)
            if save_csv:
                res.save_csv()
            results.append(res)
        except Exception as exc:
            print(f"FAILED: {exc}")
            warnings.warn(f"Simulation at φ={phi:.1f} failed: {exc}", RuntimeWarning)

    if results:
        print("\nGenerating φ-sweep overlay …")
        plot_phi_sweep(results, save=True)

    print("\nDone.  Output written to:", RESULTS_DIR.resolve())
    return results


__all__ = [
    "PHI_VALUES", "MECHANISM", "COORD", "SPECIES_TO_PLOT",
    "SimResult", "simulate_phi",
    "plot_single_phi", "plot_reaction_zone",
    "plot_emission_spectrum", "plot_phi_sweep",
    "run_all_phi",
]