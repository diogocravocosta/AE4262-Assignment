# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this directory contains

Sub-task 4.3 of the AE4262 Combustion Lab (TU Delft). Two parallel implementations for Abel inversion of OH\* chemiluminescence B16 images — one in Python, one in MATLAB — plus a CH₄ spectra overlay generator.

## Running the scripts

```powershell
# Python: Abel-invert a single B16 file (outputs <name>_abel.png)
python abel_inversion.py <path_to_file.b16>

# Python: Overlay CH4 spectra at phi = 0.8 / 1.0 / 1.2 (outputs CH4_spectra_overlay.png)
python CH4_spectra_overlay.py
```

```matlab
% MATLAB: Read a B16 file, display it, and show Abel-inverted result
process_OHstar_Abel('path\to\file.b16')

% MATLAB: Test the Fourier-based Abel inversion algorithm standalone
abel_inversion()          % no args → runs generate_test_data internally
```

## Architecture

**Python stack**

- `readB16.py` — PCO B16 parser; API is `readB16(filepath)` → NumPy `uint16` array. Note: the root-level `readB16.py` at `Experimental/question 4/readB16.py` has a slightly different internal layout but the same signature.
- `abel_inversion.py` — CLI wrapper: background-subtracts (median of top 50 rows), then calls `PyAbel` Hansen-Law inversion (`abel.Transform`). This is a **single-file tool**, unlike the parent directory's batch processor.
- `CH4_spectra_overlay.py` — composites three CH₄ spectrometer PNGs using painter's algorithm with per-φ tints. **Contains hardcoded absolute paths** (`C:\Users\franc\...`) — update the `paths` dict before running on another machine.

**MATLAB stack**

- `readB16.m` — MATLAB port of the B16 parser (Carl Hall, 2016). Same semantics as the Python version.
- `process_OHstar_Abel.m` — top-level function: reads B16, finds symmetry axis by column-sum argmax, averages left/right halves, calls `abel_inversion()` row by row, mirrors result. This calls the local `abel_inversion.m` (not PyAbel).
- `abel_inversion.m` — Fourier-based Abel inversion (Pretzler 1991, [Z. Naturforsch. 46a, 639]) implemented by C. Killer (MATLAB File Exchange #43639). Controlled by `upf` (number of cosine expansions, default 10) — lower values act as a low-pass filter.
- `compute_expansion.m`, `solve_lsq.m`, `generate_test_data.m` — helpers for `abel_inversion.m`; not called independently.

## Key differences between the two implementations

| | Python (`abel_inversion.py`) | MATLAB (`process_OHstar_Abel.m`) |
|---|---|---|
| Algorithm | Hansen-Law (PyAbel) | Pretzler 1991 Fourier cosine series |
| Symmetry axis | `symmetry_axis=0` (hardcoded left edge) | argmax of column sums |
| Background | Median of top 50 rows | None |
| Output | PNG saved to working directory | Two MATLAB figures |

## Hardcoded paths to update

- `CH4_spectra_overlay.py` lines 4–8: `paths` dict points to `C:\Users\franc\...\Spectrometer\CH4_PHI*.png`
