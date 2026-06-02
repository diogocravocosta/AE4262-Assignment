"""
h2_flame_sim/__init__.py
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
  • Mole fractions of key species: OH, H₂O, H₂, O₂, H, O, HO₂

The comparison is the core physical argument requested in section 4.3 of the
annotated assignment guide:

    "A 1D flame simulation comparing the spatial profiles of heat release
     rate and OH mole fraction confirms that heat release (and therefore
     OH* emission) peaks sharply at the flame front, while the OH
     concentration profile has a broader distribution extending downstream."

─────────────────────────────────────────────────────────────────────────────
CHOICE OF SPATIAL COORDINATE
─────────────────────────────────────────────────────────────────────────────
A freely-propagating flame in Cantera is solved on a 1-D domain.  The natural
output coordinate is the **physical axial distance x [m]**, stored in
`flame.grid`.  However, this coordinate is mechanism- and φ-dependent
(the domain width and grid spacing adapt automatically during refinement),
so direct comparison of profiles across different φ values on the same x-axis
can be misleading because the flame sits at different absolute positions.

Two better choices exist:

  1. **x − x_flame  (distance from flame front) [m]**
     Shift every profile so that x = 0 corresponds to the location of
     maximum heat release rate (the reaction zone).  This is the most
     physically transparent axis for the OH* vs OH-PLIF comparison
     because it makes the peak HRR always coincide at the origin,
     clearly showing that OH extends into the post-flame side while
     HRR does not.  ← **DEFAULT CHOICE IN THIS MODULE** (coord='x_shifted')

  2. **Progress variable c  [-]**
     The premixed analogue of mixture fraction Z (which is uniform in a
     premixed flame and carries no spatial information).
     Defined as c = (X_H2O - X_H2O,u) / (X_H2O,b - X_H2O,u), using H2O
     as the tracking species (recommended in Lecture 9 for H2/air flames).
     Ranges from 0 (unburnt reactants) to 1 (fully burnt products).

     ┌───────────────────────────────────────────────────────────────────┐
     │  NOTE: 'x_shifted' is clearest for the §4.3 spatial argument;    │
     │  'c' is useful for comparing flame structure across φ in a        │
     │  normalised 0→1 frame independent of physical thickness.          │
     └───────────────────────────────────────────────────────────────────┘

  3. **Raw grid x [m]** – simplest; use coord='x_raw' for quick checks.

Set the module-level constant `COORD` to 'x_shifted' (default), 'x_raw',
or 'c' to switch behaviour globally.

─────────────────────────────────────────────────────────────────────────────
MECHANISM
─────────────────────────────────────────────────────────────────────────────
The annotated guide explicitly states:
    "a suitable H₂/O₂ mechanism such as the Li et al. or Burke et al.
     mechanisms"

This module uses **Burke et al. (2012)** (9 species, 21 reactions), widely
available as 'h2o2.yaml' in Cantera's bundled data.  If you have the Li et
al. mechanism installed separately, change MECHANISM below.

Fallback: if 'h2o2.yaml' is not found, the module tries 'gri30.yaml', which
is always bundled with Cantera but contains many CH₄ species irrelevant to
H₂ combustion (heavier and slower to solve).

─────────────────────────────────────────────────────────────────────────────
WHAT NEEDS TO BE DONE (step-by-step checklist)
─────────────────────────────────────────────────────────────────────────────
□  1. Install dependencies
       pip install cantera matplotlib numpy

□  2. Verify mechanism file availability
       python -c "import cantera as ct; ct.Solution('h2o2.yaml')"
       If this fails, install the extended Cantera data package or point
       MECHANISM to the path of your downloaded .yaml file.

□  3. Run the module
       python -m h2_flame_sim          # runs __main__.py (see below)
       – or –
       from h2_flame_sim import run_all_phi
       run_all_phi()

□  4. Interpret the plots
       For each φ a figure is produced with two y-axes:
         • Left  axis : Temperature [K]  and mole fractions (dimensionless)
         • Right axis : Heat release rate [W/m³] (dashed, red)
       The OH mole fraction and HRR are the critical profiles:
         • HRR peaks sharply at x = 0 → this is the reaction zone where
           OH* chemiluminescence originates.
         • OH X_OH is still significant well downstream (post-flame side)
           → this is what OH-PLIF measures (ground-state OH).
       This spatial separation directly explains why Abel-inverted OH* images
       are thinner / sharper than OH-PLIF images of the same flame.

□  5. Compile multi-φ summary plot
       Call `plot_phi_sweep()` for an overlay plot of OH mole fraction and
       HRR across all φ on a single figure (useful for the assignment
       discussion of how the post-flame OH tail changes with equivalence
       ratio).

□  6. (Optional) Save data to CSV
       Each SimResult object exposes `.to_dataframe()` for downstream
       analysis in Excel / MATLAB.

─────────────────────────────────────────────────────────────────────────────
PHYSICAL BACKGROUND (to annotate your report)
─────────────────────────────────────────────────────────────────────────────
OH* (electronically excited OH) is produced mainly by the reaction
    CH + O₂  →  OH* + CO          (dominant in hydrocarbons)
    H + O + M →  OH* + M          (relevant in H₂ flames)
Its radiative lifetime is ~700 ns, so it emits at the point of formation.
Formation rate ∝ heat release rate, so OH* emission is a marker of the
reaction zone (sharply peaked).

Ground-state OH is produced in the reaction zone but also persists in the
hot post-flame gas at super-equilibrium concentrations that relax slowly.
It therefore has a broader spatial distribution (reaction zone + post-flame).

OH-PLIF probes ground-state OH → broader signal.
OH* chemiluminescence probes excited OH → narrow, reaction-zone signal.

─────────────────────────────────────────────────────────────────────────────
MODULE STRUCTURE
─────────────────────────────────────────────────────────────────────────────
  __init__.py      ← this file; public API and all simulation logic
  __main__.py      ← entry point:  python -m h2_flame_sim
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

#: Equivalence ratios to simulate (matches Table 1, H₂ premixed, Lab Exercise)
PHI_VALUES: List[float] = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]

#: Mechanism file.  Burke et al. H₂/O₂ mechanism (bundled with Cantera ≥ 3.0)
#: Change to a full path string if using a custom mechanism:
#:   MECHANISM = "/path/to/li_h2.yaml"
MECHANISM: str = "h2o2.yaml"

#: Fallback mechanism (always bundled, but includes many CH₄ species)
MECHANISM_FALLBACK: str = "gri30.yaml"

#: Fuel and oxidiser strings for set_equivalence_ratio()
FUEL: str = "H2"
OXIDISER: str = "O2:1, N2:3.76"

#: Baseline conditions (1 atm, 300 K unburnt – matches assignment baseline)
T_INLET: float = 300.0   # [K]
P_INLET: float = 100000.0  # [Pa]

#: Initial domain width for the flame solver [m]
#: Cantera refines the grid automatically; this is just the starting guess.
DOMAIN_WIDTH: float = 0.04  # 4 cm is sufficient for H₂ flames at 1 atm

#: Spatial coordinate to use for plotting.
#:   'x_shifted' – physical distance measured from the peak-HRR location [m]
#:                 (RECOMMENDED – most physically meaningful for this analysis)
#:   'x_raw'     – raw Cantera grid [m]  (quick debug)
#:   'c'         – progress variable c = (X_H2O - X_H2O,u)/(X_H2O,b - X_H2O,u)
#:                 ranges 0 (reactants) → 1 (products); premixed analogue of Z
COORD: str = "x_shifted"

#: Species to include in the per-φ profile plots
SPECIES_TO_PLOT: List[str] = ["OH", "H2O", "H2", "O2", "H", "O", "HO2"]

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

    # derived coordinate (set after construction)
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


# ── core simulation function ──────────────────────────────────────────────────

def _load_gas() -> "cantera.Solution":  # type: ignore[name-defined]
    """Load the gas object, falling back to GRI-Mech if h2o2.yaml is absent."""
    import cantera as ct  # local import so the module can be imported without cantera
    try:
        gas = ct.Solution(MECHANISM)
    except Exception:
        warnings.warn(
            f"Mechanism '{MECHANISM}' not found; falling back to '{MECHANISM_FALLBACK}'. "
            "Install the Cantera data package for the Burke H₂/O₂ mechanism.",
            RuntimeWarning,
            stacklevel=2,
        )
        gas = ct.Solution(MECHANISM_FALLBACK)
    return gas


def simulate_phi(phi: float, loglevel: int = 0) -> SimResult:
    """
    Run a 1-D freely-propagating H₂/air flame at equivalence ratio *phi*.

    Parameters
    ----------
    phi : float
        Equivalence ratio (0.5 … 1.3 for the H₂ experimental range).
    loglevel : int
        Cantera solver verbosity (0 = silent, 1 = progress, 2 = verbose).

    Returns
    -------
    SimResult
        Object containing grid, temperature, species mole fractions,
        heat release rate, and laminar flame speed.

    Notes
    -----
    The solver sequence follows the tutorial given in the 21 May lecture:
      1. Create gas, set composition and state.
      2. Construct FreeFlame on an initial coarse grid.
      3. Call flame.solve(auto=True) to iteratively refine grid and
         tighten tolerances until convergence.
      4. Extract profiles from flame.grid, flame.T, flame.X, flame.heat_release_rate.
    """
    import cantera as ct

    gas = _load_gas()
    gas.set_equivalence_ratio(phi, FUEL, OXIDISER)
    gas.TP = T_INLET, P_INLET

    # Build initial grid (Cantera will refine this automatically)
    flame = ct.FreeFlame(gas, width=DOMAIN_WIDTH)

    # Convergence settings – start loose, tighten during auto-solve
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
    flame.set_max_grid_points(flame.flame, 1000)

    # Solve with automatic grid refinement (recommended in lecture tutorial)
    flame.solve(loglevel=loglevel, auto=True)

    # ── extract profiles ──────────────────────────────────────────────────
    grid = flame.grid.copy()                      # [m]
    T = flame.T.copy()                            # [K]
    hrr = flame.heat_release_rate.copy()          # [W/m³]
    sl = float(flame.velocity[0])                 # [m/s] inlet velocity at convergence

    # Species mole fractions
    X: Dict[str, np.ndarray] = {}
    for sp in SPECIES_TO_PLOT:
        try:
            idx = gas.species_index(sp)
            X[sp] = flame.X[idx].copy()
        except ValueError:
            # species not in this mechanism – fill with zeros
            X[sp] = np.zeros_like(grid)

    # Location of peak HRR (= reaction zone / flame front)
    i_flame = int(np.argmax(hrr))
    x_flame = grid[i_flame]

    result = SimResult(
        phi=phi,
        grid=grid,
        T=T,
        X=X,
        hrr=hrr,
        sl=sl,
        x_flame=x_flame,
    )

    # Attach coordinate values based on module-level COORD setting
    result.coord_values, result.coord_label = _make_coord(result)

    return result


def _make_coord(res: SimResult) -> Tuple[np.ndarray, str]:
    """
    Build the plotting coordinate from a SimResult according to COORD setting.

    Returns
    -------
    values : np.ndarray
        The coordinate values along the flame profile.
    label : str
        Axis label string for the plot.
    """
    if COORD == "x_shifted":
        return (res.grid - res.x_flame) * 1e3, "Distance from reaction zone [mm]"
    elif COORD == "x_raw":
        return res.grid * 1e3, "Axial position  x  [mm]"
    elif COORD == "c":
        # Progress variable using H2O as tracking species (Lecture 9 recommendation)
        # c = (X_H2O - X_H2O,u) / (X_H2O,b - X_H2O,u)
        # X_H2O,u = 0 (no water in fresh H2/air mixture)
        # X_H2O,b = equilibrium value = max of the H2O profile
        x_h2o = res.X.get("H2O", np.zeros_like(res.grid))
        x_h2o_u = x_h2o[0]               # unburnt value (inlet)
        x_h2o_b = x_h2o.max()            # burnt value (equilibrium)
        if (x_h2o_b - x_h2o_u) < 1e-10:
            warnings.warn("H2O range too small to compute progress variable; falling back to x_shifted.", UserWarning)
            return (res.grid - res.x_flame) * 1e3, "x − x_flame  [mm]"
        c = (x_h2o - x_h2o_u) / (x_h2o_b - x_h2o_u)
        return c, r"Progress variable  $c = X_{\mathrm{H_2O}} / X_{\mathrm{H_2O,b}}$  [-] "
    else:
        raise ValueError(f"Unknown COORD value: '{COORD}'. Choose 'x_shifted', 'x_raw', or 'c'.")


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_single_phi(res: SimResult, save: bool = True) -> "matplotlib.figure.Figure":  # type: ignore[name-defined]
    """
    One plot per φ: species mole fractions on left axis, temperature on right axis.
    """
    import matplotlib.pyplot as plt

    fig, ax_sp = plt.subplots(figsize=(9, 5))
    ax_T = ax_sp.twinx()

    xc = res.coord_values

    colours = {
        "OH":  "#2ca02c",
        "H2O": "#1f77b4",
        "H2":  "#ff7f0e",
        "O2":  "#9467bd",
        "H":   "#17becf",
        "O":   "#e377c2",
        "HO2": "#8c564b",
    }
    # Scale small species so every line is clearly visible.
    # H2O, H2, O2 are O(0.1–0.3): no scaling needed.
    # OH, H, O are O(0.01): ×10 brings them into view.
    # HO2 is O(5e-4): ×100.
    scale = {"OH": 50, "H2O": 1, "H2": 1, "O2": 1, "H": 50, "O": 50, "HO2": 500}

    for sp, xsp in res.X.items():
        if np.max(xsp) < 1e-8:
            continue
        f = scale.get(sp, 1)
        label = f"X({sp}) ×{f}" if f > 1 else f"X({sp})"
        ax_sp.plot(xc, xsp * f, color=colours.get(sp, "grey"),
                   linewidth=2.5 if sp == "OH" else 1.2,
                   label=label)

    ax_T.plot(xc, res.T, color="black", linewidth=2.0, linestyle="--", label="T [K]")

    ax_sp.set_xlabel(res.coord_label, fontsize=11)
    ax_sp.set_ylabel("Mole fraction  [-]", fontsize=11)
    ax_sp.set_ylim(bottom=0)
    ax_T.set_ylabel("Temperature  [K]", fontsize=11)
    ax_T.set_ylim(bottom=0)
    ax_sp.set_title(
        f"H₂/air freely-propagating flame  |  φ = {res.phi:.1f}  |  SL = {res.sl*100:.1f} cm/s",
        fontsize=11,
    )

    lines1, labs1 = ax_sp.get_legend_handles_labels()
    lines2, labs2 = ax_T.get_legend_handles_labels()
    ax_sp.legend(lines1 + lines2, labs1 + labs2,
                 loc="upper left", fontsize=8, ncol=2, framealpha=0.85)
    ax_sp.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"flame_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved → {fname}")

    return fig


def plot_phi_sweep(results: List[SimResult], save: bool = True) -> "matplotlib.figure.Figure":  # type: ignore[name-defined]
    """
    Overlay plot: OH mole fraction and HRR for all φ on a single figure.

    This 'sweep' plot is useful for the assignment discussion:
      – How does the post-flame OH tail change with φ?
      – At which φ is the reaction zone (peak HRR) sharpest?
      – Confirm that the OH profile always extends further downstream
        than the HRR profile, regardless of φ.
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    ax_oh, ax_hrr = axes

    cmap = cm.plasma
    phi_vals = [r.phi for r in results]
    norm = plt.Normalize(vmin=min(phi_vals), vmax=max(phi_vals))

    for res in results:
        colour = cmap(norm(res.phi))
        xc = res.coord_values

        ax_oh.plot(xc, res.X.get("OH", np.zeros_like(xc)),
                   color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_hrr.plot(xc, res.hrr / 1e9,
                    color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")

    for ax, ylabel, title in [
        (ax_oh,  "Mole fraction X(OH)  [-]",   "OH mole fraction"),
        (ax_hrr, "Heat release rate  [GW/m³]", "Heat release rate"),
    ]:
        ax.axvline(0.0, color="grey", linestyle=":", linewidth=1.0)
        ax.set_xlabel(results[0].coord_label, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    #fig.colorbar(sm, ax=axes, label="Equivalence ratio  φ  [-]", shrink=0.8)

    fig.suptitle(
        "H₂/air premixed flames φ sweep  |  1-D Cantera FreeFlame\n",
        fontsize=11,
    )
    #"Left: OH-PLIF analogue  ·  Right: OH* chemiluminescence analogue"
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / "phi_sweep_OH_HRR.png"
        fig.savefig(fname, dpi=FIG_DPI)
        print(f"  Saved → {fname}")

    return fig


def plot_phi_sweep_zoom(results: List[SimResult], window_mm: float = 3.0, save: bool = True):
    """Zoomed overlay: OH and HRR for all φ, centred on reaction zone."""
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    ax_oh, ax_hrr = axes

    cmap = cm.plasma
    phi_vals = [r.phi for r in results]
    norm = plt.Normalize(vmin=min(phi_vals), vmax=max(phi_vals))

    for res in results:
        colour = cmap(norm(res.phi))
        xc = res.coord_values
        mask = (xc >= -window_mm) & (xc <= window_mm)
        ax_oh.plot(xc[mask], res.X.get("OH", np.zeros_like(xc))[mask],
                   color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")
        ax_hrr.plot(xc[mask], res.hrr[mask] / 1e9,
                    color=colour, linewidth=1.5, label=f"φ={res.phi:.1f}")

    for ax, ylabel, title in [
        (ax_oh,  "Mole fraction X(OH)  [-]",  "OH mole fraction"),
        (ax_hrr, "Heat release rate  [GW/m³]", "Heat release rate"),
    ]:
        ax.axvline(0.0, color="grey", linestyle=":", linewidth=1.0)
        ax.set_xlabel(results[0].coord_label, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"H₂/air premixed flames φ sweep  |  reaction zone detail  ±{window_mm} mm",
        fontsize=11,
    )
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / "phi_sweep_zoom_OH_HRR.png"
        fig.savefig(fname, dpi=FIG_DPI)
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
        Equivalence ratios to simulate.  Defaults to PHI_VALUES (0.5 … 1.3).
    loglevel : int
        Cantera solver verbosity (0 = silent).
    save_csv : bool
        If True, save each result as a CSV file in results/.

    Returns
    -------
    list of SimResult
        One SimResult per φ, in the same order as phi_values.
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
            print(f"SL = {res.sl*100:.2f} cm/s")
            plot_single_phi(res, save=True)
            plot_reaction_zone(res, save=True)
            if save_csv:
                res.save_csv()
            results.append(res)
        except Exception as exc:
            print(f"FAILED: {exc}")
            warnings.warn(f"Simulation at φ={phi:.1f} failed: {exc}", RuntimeWarning)

    if results:
        print("\nGenerating φ-sweep overlay …")
        plot_phi_sweep(results, save=True)
        plot_phi_sweep_zoom(results, save=True)  

    print("\nDone.  Output written to:", RESULTS_DIR.resolve())
    return results


# ── convenience re-exports ────────────────────────────────────────────────────
__all__ = [
    "PHI_VALUES",
    "MECHANISM",
    "COORD",
    "SPECIES_TO_PLOT",
    "SimResult",
    "simulate_phi",
    "plot_single_phi",
    "plot_phi_sweep",
    "run_all_phi",
]


def plot_reaction_zone(res: "SimResult", window_mm: float = 3.0, save: bool = True):
    """
    Zoomed plot centred on the reaction zone (x = 0 ± window_mm).
    Radicals (OH, H, O, HO2) and reactants (H2, O2) only — H2O excluded
    so the radical detail is not crushed by the product scale.
    Temperature on right axis.
    """
    import matplotlib.pyplot as plt

    fig, ax_sp = plt.subplots(figsize=(8, 5))
    ax_T = ax_sp.twinx()

    xc = res.coord_values
    mask = (xc >= -window_mm) & (xc <= window_mm)

    colours = {"OH": "#2ca02c", "H2": "#ff7f0e", "O2": "#9467bd",
               "H": "#17becf", "O": "#e377c2", "HO2": "#8c564b"}
    scale_z = {"OH": 50, "H2": 1, "O2": 1, "H": 50, "O": 50, "HO2": 500}

    for sp, col in colours.items():
        xsp = res.X.get(sp, np.zeros_like(xc))
        if np.max(xsp[mask]) < 1e-8:
            continue
        f = scale_z.get(sp, 1)
        label = f"X({sp}) \u00d7{f}" if f > 1 else f"X({sp})"
        ax_sp.plot(xc[mask], xsp[mask] * f, color=col,
                   linewidth=2.0 if sp == "OH" else 1.4, label=label)

    ax_T.plot(xc[mask], res.T[mask], color="black", linewidth=2.0,
              linestyle="--", label="T [K]")
    ax_sp.axvline(0.0, color="grey", linestyle=":", linewidth=1.0)

    ax_sp.set_xlabel(res.coord_label, fontsize=11)
    ax_sp.set_ylabel("Mole fraction [-] (scaled)", fontsize=11)
    ax_sp.set_ylim(bottom=0)
    ax_T.set_ylabel("Temperature [K]", fontsize=11)
    ax_T.set_ylim(bottom=0)
    ax_sp.set_title(
        f"Reaction zone detail  |  \u03c6 = {res.phi:.1f}  |  \u00b1{window_mm} mm",
        fontsize=11)

    l1, b1 = ax_sp.get_legend_handles_labels()
    l2, b2 = ax_T.get_legend_handles_labels()
    ax_sp.legend(l1 + l2, b1 + b2, loc="upper left", fontsize=8, ncol=2, framealpha=0.85)
    ax_sp.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = RESULTS_DIR / f"zoom_phi_{res.phi:.1f}.png"
        fig.savefig(fname, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved \u2192 {fname}")

    return fig