# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

AE4262 combustion engineering assignment (Delft University of Technology). Combines numerical thermochemical simulations (Python/Cantera) with experimental data analysis (Python image processing + MATLAB flame imaging).

## Running Scripts

No build system — scripts run independently per question:

```powershell
# Numerical thermochemical calculations
cd "Numerical Part"
python part2.py        # adiabatic flame temperature & NO emissions (outputs Part2a–d PDFs)
python part3a.py       # two-step NO kinetics (mass fractions)
python part3appm.py    # same in ppm units

# Experimental Q1: stoichiometry & flow rate tables
cd "Experimental/question 1"
python Chemical_reactions.py   # stoichiometric coefficients
python phitable.py             # phi vs. flow rate tables (CSV + PNG)

# Experimental Q4: spectroscopy & Abel inversion
cd "Experimental/question 4"
python abel_inversion.py       # processes B16 images, applies Abel inversion

# Experimental Q2: flame front analysis (requires MATLAB)
cd "Experimental/question 2/UNIFORM_ATEMP"
# run main_q2.m in MATLAB  ← current pipeline (auto-detection, uniform params)
# root-level main_q2.m     ← legacy pipeline (manual per-image tuning)
```

## Key Dependencies

```
cantera         # chemical kinetics & equilibrium (GRI-Mech 3.0, gri30.yaml)
numpy, scipy    # numerical computation
matplotlib      # plotting (outputs as PDF/PNG)
CoolProp        # fluid properties (density, viscosity, MW)
PyAbel          # Abel inversion (Hansen-Law method) for radial reconstruction
```

No `requirements.txt` exists — install these manually via pip/conda.

## Architecture

**Numerical Part** (`Numerical Part/`): Cantera equilibrium calculations parametrized over equivalence ratio φ, pressure, and inlet temperature. Uses `gri30.yaml` mechanism. Pre-computed thermodynamic data cached in `cantera_cp_table.csv` at the repo root.

**Experimental Q1** (`Experimental/question 1/`): `Chemical_reactions.py` provides stoichiometry logic consumed by `phitable.py` and the flame-speed scripts. Outputs CSV lookup tables mapping φ → mass flow rates and Reynolds numbers for CH₄ and H₂.

**Experimental Q2** (`Experimental/question 2/`): MATLAB pipeline for Bunsen flame speed measurement via OH-PLIF imaging. Two implementations exist:

**Current pipeline** (`UNIFORM_ATEMP/`): Refactored, parameter-uniform approach with automatic boundary detection.

Pipeline flow in `UNIFORM_ATEMP/main_q2.m`:
1. **Calibration** — `calibratePixelPerMM` returns `px_per_mm_x/y` for each fuel (separate CH₄ and H₂ calibration images).
2. **Image locations & velocities** — `OHPLIF_locations` and `Get_experimental_data` (same as legacy).
3. **Flame front detection** — Both fuels call `flameFrontDetect` via `flamefront_sf_H2` / `flamefront_sf_CH4`. Three shared parameters replace all manual per-image tuning:
   - `threshold` (0.25) — local normalised intensity threshold for edge crossing
   - `drop_frac` (0.80) — gradient drops to this fraction of peak → working point
   - `trunc_frac` (0.10) — fraction trimmed from each fit end
4. **Flame speed**: Sₗ = U_mean · sin(α/2), stored as `SF_from_angle`.
5. **Plot**: `plotSfAngle` (shared with legacy).

`flameFrontDetect.m` is the core algorithm: intensity-weighted auto-centering, 5%-peak row-average threshold for automatic `row_min`/`row_max`, gradient-guided edge detection (robust to OH double-peak profiles in rich H₂), linear fit below the cone tip.

**Legacy pipeline** (root of `Experimental/question 2/`): Uses manually-tuned per-image arrays (`baselines`, `centers`, `tops`, `thresholds`) and separate detection functions (`flameFrontAnglePlot` for H₂, `flameFrontStats` for CH₄). Kept for reference.

**Critical**: `OHPLIF_locations.m` and `main_q2.m` (in both pipeline roots) contain absolute paths hardcoded to `C:\Users\franc\...`. These must be updated when running on a different machine.

**Experimental Q4** (`Experimental/question 4/`): Python pipeline.
- `readB16.py` — parses the PCO B16 binary image format; call `read_b16(filepath)` to get a NumPy array.
- `abel_inversion.py` — loads B16 images from `OHstar/` (H₂, φ = 0.5–1.2) and `Spectrometer/` (CH₄, φ = 0.8–1.2), applies PyAbel Hansen-Law inversion to reconstruct radial emission profiles, and plots results.

**Data locations:**
- `Experimental/question 4/OHstar/` — 10 subdirectories of H₂ OH* images
- `Experimental/question 4/Spectrometer/` — 3 CH₄ spectrometer images
- `Experimental/question 2/calibration/`, `LIF/` — MATLAB input images

## Cantera Usage Pattern

All numerical scripts follow this pattern:

```python
import cantera as ct
gas = ct.Solution('gri30.yaml')
gas.TP = T, P
gas.set_equivalence_ratio(phi, fuel, oxidizer)
gas.equilibrate('HP')   # or 'TP' depending on the problem
```

Fuels studied: CH₄ (`CH4`) and H₂ (`H2`). Oxidizer is air.
