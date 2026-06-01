# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Pipeline

All scripts run in MATLAB. The preferred entry point is the unified pipeline:

```matlab
% Preferred — gradient-guided, no manual tuning
cd UNIFORM_ATEMP
run main_q2.m

% Legacy — manual per-image parameter arrays
run main_q2.m   % from question 2 root
```

Before running on a new machine, update the absolute paths hardcoded in `OHPLIF_locations.m` and both `main_q2.m` files (`C:\Users\franc\...`).

## Architecture

There are two parallel implementations sharing the same support files:

**Root (legacy):** Each fuel uses a separate detection function (`flameFrontAnglePlot.m` for H2, `flameFrontStats.m` for CH4) driven by manually-tuned per-image arrays (`baselines`, `centers`, `tops`, `thresholds`).

**`UNIFORM_ATEMP/` (preferred):** Both fuels route through a single `flameFrontDetect.m` with three shared parameters and automatic boundary detection. The batch wrappers (`flamefront_sf_H2.m`, `flamefront_sf_CH4.m`) are identical in both directories — the `UNIFORM_ATEMP` versions simply call `flameFrontDetect` instead of the fuel-specific functions.

### Shared support files (root directory, used by both pipelines)

| File | Role |
|---|---|
| `readB16.m` | Parses PCO binary camera format; returns `uint16` image array |
| `calibratePixelPerMM.m` | Grid-paper calibration → `px_per_mm_x/y` (~29.7 px/mm) |
| `OHPLIF_locations.m` | Hardcoded cell arrays of B16 file paths for all 12 images |
| `Get_experimental_data.m` | Hardcoded bulk velocity data [phi, v_mm/s] for H2 (9) and CH4 (3) |
| `plotSfAngle.m` | Plots experimental vs. literature vs. Cantera results for both fuels |
| `literature_data_{H2,CH4}.m` | Hardcoded literature flame speed reference data |
| `cantera_data_{H2,CH4}.m` | Hardcoded Cantera numerical results |

### Core detection: `UNIFORM_ATEMP/flameFrontDetect.m`

1. Reads and Gaussian-smooths the image (15×15 kernel, σ=3)
2. Auto-detects `row_min`/`row_max` from the row-averaged intensity profile (no manual input)
3. Per row, scans outward from the flame axis; locates the steepest gradient rise then finds where it drops to `drop_frac × peak` → **working point**; re-normalises locally; takes first threshold crossing as the edge
4. Detects flame tip as the narrowest width row (incrementally relaxes threshold if needed)
5. Fits lines to left/right edges below the tip; returns opening angle `alpha_deg`, `Dc_mm`, and per-row edge positions

The gradient-guided working point is the key innovation — it handles the OH double-peak profiles that appear in rich H2 flames and would fool a plain threshold scan.

**Parameters (set in `UNIFORM_ATEMP/main_q2.m`):**

| Parameter | Value | Effect |
|---|---|---|
| `threshold` | 0.25 | Normalised intensity crossing (0–1) |
| `drop_frac` | 0.90 | Gradient must drop to 90 % of peak to set working point |
| `trunc_frac` | 0.10 | Trim 10 % from each end of the edge set before fitting |

### Flame speed formula

```
Sf = v_bulk * sin(alpha_deg / 2)
```

where `alpha_deg` is the full cone opening angle and `v_bulk` is the bulk mean velocity from `Get_experimental_data`.

## Image Data Layout

```
LIF/
  H2_PHI{0.5,0.6,...,1.3}/B16/B0001.b16   ← 9 H2 images
  CH4_PHI{0.8,1.0,1.2}/B16/B0001.b16      ← 3 CH4 images
calibration/
  Calibration/H2/B16/B0001.b16
  Calibration/CH4/B16/B0001.b16
```

## `verbose` Diagnostic Plots

Pass `verbose = true` to `flameFrontDetect` to get four diagnostic figures per image:
1. Row-averaged intensity with auto-detected flame boundaries
2. Raw intensity profiles (right half) with detected edges overlaid
3. Gradient profiles (right half) with working points
4. Full flame image with both edge fits and angle annotation
