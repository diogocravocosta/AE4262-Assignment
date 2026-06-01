function results = flamefront_sf(filenames, px2mm_x, px2mm_y, threshold, drop_frac, trunc_tip, trunc_base, velocities, verbose_diag, verbose_overlay, fit_lines_only)
% FLAMEFRONT_SF  Batch flame speed using unified gradient-guided detection.
%
%   Works for any fuel — pass the appropriate calibration and velocity data.
%
%   Inputs:
%       filenames       - cell array of .b16 file paths
%       px2mm_x         - mm/px in x  (1 / px_per_mm_x from calibratePixelPerMM)
%       px2mm_y         - mm/px in y
%       threshold       - local normalised threshold (default 0.25)
%       drop_frac       - gradient drop fraction for working-point detection (default 0.80)
%       trunc_tip       - fraction trimmed from the tip end of below-tip rows (default 0.10)
%       trunc_base      - fraction trimmed from the base end of below-tip rows (default 0.10)
%       velocities      - mean bulk velocity per image [mm/s]
%       verbose_diag    - row-avg, intensity, gradient figures (default true)
%       verbose_overlay - flame image with front overlay and fit (default true)
%       fit_lines_only  - overlay shows only fitted lines (default false)
%
%   Output fields per image:
%       filename, flame_left, flame_right, angle_deg, alpha_half_deg, SF_from_angle

    if nargin < 9  || isempty(verbose_diag),    verbose_diag    = true;  end
    if nargin < 10 || isempty(verbose_overlay), verbose_overlay = true;  end
    if nargin < 11 || isempty(fit_lines_only),  fit_lines_only  = false; end

    n       = length(filenames);
    results = struct('filename', cell(1,n), 'flame_left', cell(1,n), ...
                     'flame_right', cell(1,n), 'angle_deg', cell(1,n), ...
                     'alpha_half_deg', cell(1,n), 'SF_from_angle', cell(1,n));

    for i = 1:n
        stats = flameFrontDetect(filenames{i}, ...
            1/px2mm_x, 1/px2mm_y, threshold, drop_frac, trunc_tip, trunc_base, verbose_diag, verbose_overlay, fit_lines_only);

        results(i).filename       = filenames{i};
        results(i).flame_left     = stats.left_edge;
        results(i).flame_right    = stats.right_edge;
        results(i).V_bulk         = velocities(i);
        results(i).angle_deg      = stats.alpha_deg;
        results(i).alpha_half_deg = stats.alpha_deg / 2;
        results(i).SF_from_angle  = velocities(i) * sin(deg2rad(stats.alpha_deg) / 2);
    end
end
