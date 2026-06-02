"""
No_flame_sim/__init__.py
========================
AE4262 Combustion – Assignment 2026, Section 4.3
Cantera 1-D freely-propagating H₂/air flame simulation
-------------------------------------------------------

PURPOSE
-------
This module runs a set of 1-D freely-propagating premixed H₂/air flame
simulations over the experimental equivalence-ratio range φ = 0.5 … 1.3
(step 0.1, matching the lab setpoints in Table 1 of the Laboratory Exercise).

For each φ it extracts and plots – on a single figure per φ – the spatial
profiles of:
  • Temperature  T  [K]
  • Mole fractions of key species: OH, H₂O, H₂, O₂, H, O, HO₂, NO
  • Heat release rate  HRR  [W/m³]

NO (nitric oxide) is extracted alongside the combustion species to allow
direct comparison of its spatial profile with temperature and radicals.
Because NO formation via the thermal (Zeldovich) mechanism is kinetically
slow (high activation energy), its profile peaks well downstream of the
reaction zone — a clear illustration of why equilibrium overpredicts NO.

─────────────────────────────────────────────────────────────────────────────
MECHANISM
─────────────────────────────────────────────────────────────────────────────
GRI-Mech 3.0 ('gri30.yaml') is used by default because it includes full
thermal-NO chemistry (Zeldovich + prompt + N₂O pathways), making it the
natural choice for NOx estimation in H₂/air flames.

If you prefer the leaner Burke et al. H₂/O₂ mechanism ('h2o2.yaml'),
set MECHANISM = 'h2o2.yaml' — but note that this mechanism does NOT contain
nitrogen species, so NO will be zero everywhere.

─────────────────────────────────────────────────────────────────────────────
SPATIAL COORDINATE
─────────────────────────────────────────────────────────────────────────────
  'x_shifted'  Distance from peak-HRR (reaction zone), [mm].  DEFAULT.
  'x_raw'      Raw Cantera grid [mm].
  'c'          Progress variable c = X_H2O / X_H2O,b ∈ [0,1].

─────────────────────────────────────────────────────────────────────────────
MODULE STRUCTURE
─────────────────────────────────────────────────────────────────────────────
  __init__.py      ← this file; public API and all simulation logic
  __main__.py      ← entry point:  python -m No_flame_sim
  results/         ← output figures and CSV files are written here
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ── output directory ─────────────────────────────────────────────────────────
RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── simulation settings ───────────────────────────────────────────────────────

#: Equivalence ratios to simulate
PHI_VALUES: List[float] = [0.7, 1.0, 1.3]

#: Mechanism file.
#: GRI-Mech 3.0 is required for NO chemistry.
#: Swap to 'h2o2.yaml' for a leaner H₂ mechanism (NO will be zero).
MECHANISM: str = "gri30.yaml"

#: Fuel and oxidiser strings for set_equivalence_ratio()
FUEL: str = "H2"
OXIDISER: str = "O2:1, N2:3.76"

#: Baseline conditions (1 atm, 300 K unburnt)
T_INLET: float = 300.0    # [K]
P_INLET: float = 100000.0  # [Pa]

#: Initial domain width for the flame solver [m]
DOMAIN_WIDTH: float = 0.04  # 4 cm

#: Spatial coordinate for plotting
COORD: str = "x_shifted"

#: Species to include in the per-φ profile plots
SPECIES_TO_PLOT: List[str] = ["OH", "H2O", "H2", "O2", "H", "O", "HO2", "NO"]

#: Scale factors so every species is visible on the same axis
#: OH, H, O are O(0.01) → ×50 brings them to O(0.5)
#: HO2 is O(5e-4) → ×500
#: NO is O(1e-6 … 1e-5) → ×50 000 brings it to O(0.05)
SCALE: Dict[str, int] = {
    "OH":  50,
    "H2O": 1,
    "H2":  1,
    "O2":  1,
    "H":   50,
    "O":   50,
    "HO2": 500,
    "NO":  50_000,
}

#: Colours per species
COLOURS: Dict[str, str] = {
    "OH":  "#2ca02c",
    "H2O": "#1f77b4",
    "H2":  "#ff7f0e",
    "O2":  "#9467bd",
    "H":   "#17becf",
    "O":   "#e377c2",
    "HO2": "#8c564b",
    "NO":  "#d62728",  # red — stands out as the emissions species
}

#: Figure DPI
FIG_DPI: int = 150


# ── data container ────────────────────────────────────────────────────────────

@dataclass
class SimResult:
    """Container for a single φ simulation result."""

    phi: float
    grid: np.ndarray           # raw Cantera grid [m]
    T: np.ndarray              # temperature [K]
    X: Dict[str, np.ndarray]  # mole fractions, keyed by species name
    hrr: np.ndarray            # heat release rate [W/m³]
    sl: float                  # laminar flame speed [m/s]
    x_flame: float             # axial position of max-HRR [m]

    coord_values: np.ndarray = field(default_factory=lambda: np.array([]))
    coord_label: str = ""

    def to_dataframe(self):
        """Return a pandas DataFrame with all profiles (requires pandas)."""
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
        """Save profiles to CSV."""
        if path is None:
            path = RESULTS_DIR / f"phi_{self.phi:.1f}.csv"
        self.to_dataframe().to_csv(path, index=False)
        return path

    @property
    def no_peak_ppm(self) -> float:
        """Peak NO mole fraction in ppm."""
        x_no = self.X.get("NO", np.zeros_like(self.T))
        return float(x_no.max()) * 1e6

    @property
    def no_peak_position_mm(self) -> float:
        """Position of peak NO relative to reaction zone [mm]."""
        x_no = self.X.get("NO", np.zeros_like(self.T))
        i_peak = int(np.argmax(x_no))
        return float((self.grid[i_peak] - self.x_flame) * 1e3)


# ── core simulation function ──────────────────────────────────────────────────

def _load_gas() -> "cantera.Solution":  # type: ignore[name-defined]
    """Load the gas object, falling back to gri30.yaml if the primary is absent."""
    import cantera as ct
    try:
        gas = ct.Solution(MECHANISM)
        print(f"Loaded mechanism: {MECHANISM}  ({gas.n_species} species, {gas.n_reactions} reactions)")
    except Exception:
        fallback = "gri30.yaml"
        warnings.warn(
            f"Mechanism '{MECHANISM}' not found; falling back to '{fallback}'.",
            RuntimeWarning,
            stacklevel=2,
        )
        gas = ct.Solution(fallback)
        print(f"Loaded fallback: {fallback}  ({gas.n_species} species, {gas.n_reactions} reactions)")
    # Warn if NO is absent (happens with h2o2.yaml)
    try:
        gas.species_index("NO")
    except ValueError:
        warnings.warn(
            "Species 'NO' not found in mechanism. NO profiles will be zero. "
            "Switch to gri30.yaml for NOx predictions.",
            UserWarning,
            stacklevel=2,
        )
    return gas


def simulate_phi(phi: float, loglevel: int = 0) -> SimResult:
    """
    Run a 1-D freely-propagating H₂/air flame at equivalence ratio *phi*.

    Parameters
    ----------
    phi : float
        Equivalence ratio.
    loglevel : int
        Cantera solver verbosity (0 = silent, 1 = progress).

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
    T = flame.T.copy()
    hrr = flame.heat_release_rate.copy()
    sl = float(flame.velocity[0])

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
        phi=phi, grid=grid, T=T, X=X, hrr=hrr, sl=sl, x_flame=x_flame,
    )
    result.coord_values, result.coord_label = _make_coord(result)
    return result


def _make_coord(res: SimResult) -> Tuple[np.ndarray, str]:
    """Build the plotting coordinate from a SimResult."""
    if COORD == "x_shifted":
        return (res.grid - res.x_flame) * 1e3, "Distance from reaction zone [mm]"
    elif COORD == "x_raw":
        return res.grid * 1e3, "Axial position  x  [mm]"
    elif COORD == "c":
        x_h2o = res.X.get("H2O", np.zeros_like(res.grid))
        x_h2o_u = x_h2o[0]
        x_h2o_b = x_h2o.max()
        if (x_h2o_b - x_h2o_u) < 1e-10:
            warnings.warn("H2O range too small; falling back to x_shifted.", UserWarning)
            return (res.grid - res.x_flame) * 1e3, "x − x_flame  [mm]"
        c = (x_h2o - x_h2o_u) / (x_h2o_b - x_h2o_u)
        return c, r"Progress variable  $c = X_{\mathrm{H_2O}} / X_{\mathrm{H_2O,b}}$  [-]"
    else:
        raise ValueError(f"Unknown COORD: '{COORD}'. Choose 'x_shifted', 'x_raw', or 'c'.")


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_single_phi(res: SimResult, save: bool = True) -> "matplotlib.figure.Figure":
    """
    Full-domain flame profile:
      Left axis  – species mole fractions (scaled for visibility)
      Right axis – temperature [K]
    NO is highlighted in red with its scale factor shown in the legend.
    """
    import matplotlib.pyplot as plt

    fig, ax_sp = plt.subplots(figsize=(9, 5))
    ax_T = ax_sp.twinx()
    xc = res.coord_values

    for sp, xsp in res.X.items():
        f = SCALE.get(sp, 1)
        if np.max(xsp) * f < 1e-8:
            continue
        label = f"X({sp}) ×{f}" if f > 1 else f"X({sp})"
        lw = 2.5 if sp in ("OH", "NO") else 1.5
        ls = "--" if sp == "NO" else "-"
        ax_sp.plot(xc, xsp * f, color=COLOURS.get(sp, "grey"),
                   linewidth=lw, linestyle=ls, label=label)

    ax_T.plot(xc, res.T, color="black", linewidth=2.0, linestyle=":", label="T [K]")
    ax_sp.axvline(0.0, color="grey", linestyle=":", linewidth=1.0, alpha=0.6)

    ax_sp.set_xlabel(res.coord_label, fontsize=11)
    ax_sp.set_ylabel("Mole fraction  [-]  (scaled)", fontsize=11)
    ax_sp.set_ylim(bottom=0)
    ax_T.set_ylabel("Temperature  [K]", fontsize=11)
    ax_T.set_ylim(bottom=0)
    ax_sp.set_title(
        f"H₂/air freely-propagating flame  |  φ = {res.phi:.1f}  |  "
        f"SL = {res.sl*100:.1f} cm/s  |  NO_peak = {res.no_peak_ppm:.2f} ppm",
        fontsize=10,
    )

    l1, b1 = ax_sp.get_legend_handles_labels()
    l2, b2 = ax_T.get_legend_handles_labels()
    ax_sp.legend(l1 + l2, b1 + b2, loc="upper left", fontsize=8,
                 ncol=2, framealpha=0.85)
    ax_sp.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"flame_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_reaction_zone(res: SimResult, window_mm: float = 3.0, save: bool = True):
    """
    Zoomed plot centred on the reaction zone (x = 0 ± window_mm).
    Shows radicals (OH, H, O, HO2), reactants (H2, O2), and NO.
    Temperature on right axis.
    """
    import matplotlib.pyplot as plt

    fig, ax_sp = plt.subplots(figsize=(8, 5))
    ax_T = ax_sp.twinx()

    xc = res.coord_values
    mask = (xc >= -window_mm) & (xc <= window_mm)

    # Exclude H2O from the zoomed view — its large scale crushes radical detail
    species_zoom = [sp for sp in SPECIES_TO_PLOT if sp != "H2O"]

    for sp in species_zoom:
        xsp = res.X.get(sp, np.zeros_like(xc))
        f = SCALE.get(sp, 1)
        if np.max(xsp[mask]) * f < 1e-8:
            continue
        label = f"X({sp}) ×{f}" if f > 1 else f"X({sp})"
        lw = 2.0 if sp in ("OH", "NO") else 1.4
        ls = "--" if sp == "NO" else "-"
        ax_sp.plot(xc[mask], xsp[mask] * f, color=COLOURS.get(sp, "grey"),
                   linewidth=lw, linestyle=ls, label=label)

    ax_T.plot(xc[mask], res.T[mask], color="black", linewidth=2.0,
              linestyle=":", label="T [K]")
    ax_sp.axvline(0.0, color="grey", linestyle=":", linewidth=1.0, alpha=0.6)

    ax_sp.set_xlabel(res.coord_label, fontsize=11)
    ax_sp.set_ylabel("Mole fraction [-] (scaled)", fontsize=11)
    ax_sp.set_ylim(bottom=0)
    ax_T.set_ylabel("Temperature [K]", fontsize=11)
    ax_T.set_ylim(bottom=0)
    ax_sp.set_title(
        f"Reaction zone detail  |  φ = {res.phi:.1f}  |  ±{window_mm} mm  |  "
        f"NO_peak @ {res.no_peak_position_mm:.1f} mm",
        fontsize=10)

    l1, b1 = ax_sp.get_legend_handles_labels()
    l2, b2 = ax_T.get_legend_handles_labels()
    ax_sp.legend(l1 + l2, b1 + b2, loc="upper left", fontsize=8,
                 ncol=2, framealpha=0.85)
    ax_sp.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"zoom_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_no_profile(res: SimResult, save: bool = True):
    """
    Dedicated NO profile plot with T on the right axis.
    Shows NO in ppm (not scaled) so absolute concentrations are clear.
    Annotates the NO peak position relative to the reaction zone.
    """
    import matplotlib.pyplot as plt

    fig, ax_no = plt.subplots(figsize=(8, 4))
    ax_T = ax_no.twinx()

    xc = res.coord_values
    x_no = res.X.get("NO", np.zeros_like(xc))

    ax_no.plot(xc, x_no * 1e6, color="#d62728", linewidth=2.0, label="NO [ppm]")
    ax_T.plot(xc, res.T, color="black", linewidth=1.5, linestyle=":", label="T [K]")

    ax_no.axvline(0.0, color="grey", linestyle=":", linewidth=1.0, alpha=0.6,
                  label="Reaction zone (peak HRR)")

    # Annotate peak NO
    i_peak = int(np.argmax(x_no))
    ax_no.annotate(
        f"peak NO = {x_no[i_peak]*1e6:.2f} ppm\n@ {xc[i_peak]:.1f} mm",
        xy=(xc[i_peak], x_no[i_peak] * 1e6),
        xytext=(xc[i_peak] + 1.5, x_no[i_peak] * 1e6 * 0.85),
        arrowprops=dict(arrowstyle="->", color="#d62728"),
        fontsize=8, color="#d62728",
    )

    ax_no.set_xlabel(res.coord_label, fontsize=11)
    ax_no.set_ylabel("NO mole fraction  [ppm]", fontsize=11)
    ax_no.set_ylim(bottom=0)
    ax_T.set_ylabel("Temperature  [K]", fontsize=11)
    ax_T.set_ylim(bottom=0)
    ax_no.set_title(
        f"NO spatial profile  |  φ = {res.phi:.1f}  |  SL = {res.sl*100:.1f} cm/s",
        fontsize=11)

    l1, b1 = ax_no.get_legend_handles_labels()
    l2, b2 = ax_T.get_legend_handles_labels()
    ax_no.legend(l1 + l2, b1 + b2, loc="upper left", fontsize=9, framealpha=0.85)
    ax_no.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"NO_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_phi_sweep(results: List[SimResult], save: bool = True):
    """
    Overlay: OH mole fraction, NO mole fraction [ppm], and HRR for all φ.
    Three-panel figure using a plasma colormap to distinguish φ values.
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
    ax_oh, ax_no, ax_hrr = axes

    cmap = cm.plasma
    phi_vals = [r.phi for r in results]
    norm = plt.Normalize(vmin=min(phi_vals), vmax=max(phi_vals))

    for res in results:
        colour = cmap(norm(res.phi))
        xc = res.coord_values
        oh = res.X.get("OH", np.zeros_like(xc))
        no = res.X.get("NO", np.zeros_like(xc))

        ax_oh.plot(xc, oh, color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_no.plot(xc, no * 1e6, color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_hrr.plot(xc, res.hrr / 1e9, color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")

    for ax, ylabel, title in [
        (ax_oh,  "OH mole fraction [-]",    "OH radical"),
        (ax_no,  "NO mole fraction [ppm]",  "NO emissions"),
        (ax_hrr, "Heat release rate [GW/m³]", "Heat release rate"),
    ]:
        ax.axvline(0.0, color="grey", linestyle=":", linewidth=1.0, alpha=0.6)
        ax.set_xlabel(results[0].coord_label, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=2)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "H₂/air premixed flames  –  φ sweep  |  1-D Cantera FreeFlame (GRI-Mech 3.0)",
        fontsize=11,
    )
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / "phi_sweep_OH_NO_HRR.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_phi_sweep_zoom(results: List[SimResult], window_mm: float = 3.0, save: bool = True):
    """Zoomed reaction-zone overlay: OH, NO, and HRR for all φ."""
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
    ax_oh, ax_no, ax_hrr = axes

    cmap = cm.plasma
    phi_vals = [r.phi for r in results]
    norm = plt.Normalize(vmin=min(phi_vals), vmax=max(phi_vals))

    for res in results:
        colour = cmap(norm(res.phi))
        xc = res.coord_values
        mask = (xc >= -window_mm) & (xc <= window_mm)
        oh = res.X.get("OH", np.zeros_like(xc))
        no = res.X.get("NO", np.zeros_like(xc))

        ax_oh.plot(xc[mask], oh[mask], color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_no.plot(xc[mask], no[mask] * 1e6, color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_hrr.plot(xc[mask], res.hrr[mask] / 1e9, color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")

    for ax, ylabel, title in [
        (ax_oh,  "OH mole fraction [-]",      "OH radical"),
        (ax_no,  "NO mole fraction [ppm]",    "NO emissions"),
        (ax_hrr, "Heat release rate [GW/m³]", "Heat release rate"),
    ]:
        ax.axvline(0.0, color="grey", linestyle=":", linewidth=1.0, alpha=0.6)
        ax.set_xlabel(results[0].coord_label, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=2)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"H₂/air flames φ sweep  –  reaction zone detail ±{window_mm} mm  |  GRI-Mech 3.0",
        fontsize=11,
    )
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / "phi_sweep_zoom_OH_NO_HRR.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_no_sweep(results: List[SimResult], save: bool = True):
    """
    Summary bar chart: peak NO [ppm] vs φ.
    Quick visual showing how NO emissions vary across equivalence ratios.
    """
    import matplotlib.pyplot as plt

    phi_vals = [r.phi for r in results]
    no_peaks = [r.no_peak_ppm for r in results]

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(phi_vals, no_peaks, width=0.07, color="#d62728",
                  edgecolor="black", linewidth=0.7, alpha=0.85)

    for bar, val in zip(bars, no_peaks):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                f"{val:.2f}", ha="center", va="bottom", fontsize=8)

    ax.set_xlabel("Equivalence ratio  φ  [-]", fontsize=11)
    ax.set_ylabel("Peak NO  [ppm]", fontsize=11)
    ax.set_title("Peak NO vs φ  |  H₂/air 1-D flame  |  GRI-Mech 3.0", fontsize=11)
    ax.set_xticks(phi_vals)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / "NO_sweep_summary.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


# ── public API ────────────────────────────────────────────────────────────────

def run_all_phi(
    phi_values: Optional[List[float]] = None,
    loglevel: int = 0,
    save_csv: bool = False,
) -> List[SimResult]:
    """
    Run simulations for all equivalence ratios and produce all plots.

    Parameters
    ----------
    phi_values : list of float, optional
        Equivalence ratios to simulate.  Defaults to PHI_VALUES.
    loglevel : int
        Cantera solver verbosity (0 = silent).
    save_csv : bool
        If True, save each result as a CSV file in results/.

    Returns
    -------
    list of SimResult
    """
    if phi_values is None:
        phi_values = PHI_VALUES

    results: List[SimResult] = []
    print(f"Running H₂/air FreeFlame simulations for φ = {phi_values}")
    print(f"Mechanism : {MECHANISM}")
    print(f"Coordinate: {COORD}")
    print("-" * 60)

    for phi in phi_values:
        print(f"  φ = {phi:.1f} … ", end="", flush=True)
        try:
            res = simulate_phi(phi, loglevel=loglevel)
            print(f"SL = {res.sl*100:.2f} cm/s  |  NO_peak = {res.no_peak_ppm:.2f} ppm")
            plot_single_phi(res, save=True)
            plot_reaction_zone(res, save=True)
            plot_no_profile(res, save=True)
            if save_csv:
                res.save_csv()
            results.append(res)
        except Exception as exc:
            print(f"FAILED: {exc}")
            warnings.warn(f"Simulation at φ={phi:.1f} failed: {exc}", RuntimeWarning)

    if results:
        print("\nGenerating φ-sweep overlays …")
        plot_phi_sweep(results, save=True)
        plot_phi_sweep_zoom(results, save=True)
        plot_no_sweep(results, save=True)

    print("\nDone.  Output written to:", RESULTS_DIR.resolve())
    return results


# ── convenience re-exports ────────────────────────────────────────────────────
__all__ = [
    "PHI_VALUES",
    "MECHANISM",
    "COORD",
    "SPECIES_TO_PLOT",
    "SCALE",
    "SimResult",
    "simulate_phi",
    "plot_single_phi",
    "plot_reaction_zone",
    "plot_no_profile",
    "plot_phi_sweep",
    "plot_phi_sweep_zoom",
    "plot_no_sweep",
    "run_all_phi",
]