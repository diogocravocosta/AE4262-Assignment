# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Pipeline

```matlab
% From this directory in MATLAB:
run main_q2.m
```

Before running on a new machine: update the two absolute paths in `main_q2.m` (calibration `.b16` files) and in `OHPLIF_locations.m` (all LIF image paths). Both are hardcoded to `C:\Users\franc\...`.

## File Roles

| File | Role |
|---|---|
| `main_q2.m` | Entry point — calibration, detection parameters, batch runs, final plot |
| `flameFrontDetect.m` | Core algorithm: gradient-guided edge detection for a single image |
| `flamefront_sf_H2.m` | Batch wrapper: loops over 9 H2 images, calls `flameFrontDetect` |
| `flamefront_sf_CH4.m` | Batch wrapper: loops over 3 CH4 images, calls `flameFrontDetect` |
| `calibratePixelPerMM.m` | Returns `px_per_mm_x/y` from a grid-paper calibration image |
| `OHPLIF_locations.m` | Returns cell arrays of `.b16` paths for all 12 flame images |
| `Get_experimental_data.m` | Returns `[phi, v_bulk_m/s]` arrays for H2 (9 rows) and CH4 (3 rows) |
| `plotSfAngle.m` | Final plot: experimental vs. literature vs. Cantera flame speeds |
| `readB16.m` | Parses PCO binary format; returns `uint16` image array |

Reference data files (`literature_data_H2/CH4.m`, `cantera_data_H2/CH4.m`) contain hardcoded arrays consumed only by `plotSfAngle`.

## Detection Algorithm (`flameFrontDetect.m`)

**Key insight:** Uses the smoothed-intensity *gradient* rather than raw intensity to locate the flame edge. This handles OH double-peak profiles that appear in rich H2 flames and would fool a plain threshold scan.

**Per-row edge detection (`edge_from_gradient`):**
1. Gradient of smoothed row profile → `peak_i` (steepest rise)
2. After `peak_i`: find `x_drop` where gradient ≤ `drop_frac × peak` (robust to plateaus)
3. Re-normalise raw intensity locally around `x_drop`
4. First threshold crossing → flame edge pixel

**Auto-boundary detection (no manual input needed):**
- `row_min`: first row where row-averaged intensity exceeds 35% of max
- `row_max`: row with steepest negative gradient in the row-average profile below `row_min`

**Tip detection:** scans backward from bottom; finds first row where width ≤ 5% of base width (threshold relaxes incrementally if not met).

**Output struct:** `alpha_deg`, `Dc_mm`, `tip_row/col`, `rows`, `left_edge`, `right_edge`, `c_left/c_right` (pixel-space linear fits), `row_min/row_max`.

## Tunable Parameters (set in `main_q2.m`)

| Parameter | Current value | Effect |
|---|---|---|
| `threshold` | 0.25 | Local normalised intensity threshold for edge crossing (0–1) |
| `drop_frac` | 0.90 | Gradient must drop to 90% of its peak to set the working point |
| `trunc_frac` | 0.10 | Trims 10% from each end of below-tip rows before linear fitting |

## Flame Speed Formula

```
Sf = v_bulk [mm/s] × sin(alpha_deg / 2)
```

`v_bulk` comes from `Get_experimental_data`; `alpha_deg` is the full cone opening angle from `flameFrontDetect`.

## Diagnostics

Pass `verbose = true` (default in `main_q2.m`) to `flameFrontDetect` to get four figures per image:
1. Row-average intensity with detected `row_min`/`row_max`
2. Right-half intensity profiles with detected edge points
3. Right-half gradient profiles with working points
4. Full image overlay with edge fits and angle annotation

Set `verbose = false` in `main_q2.m` to suppress all per-image diagnostic plots and only see the final `plotSfAngle` figure.
