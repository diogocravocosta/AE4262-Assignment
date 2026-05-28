function [angle_deg, line_left, line_right, r_left_trunc, r_right_trunc, z_left_trunc, z_right_trunc] = flameFrontAngle(flame_left, flame_right, flame_center, px2mm_x, px2mm_y, trunc_frac)
% FLAMEFRONTANGLE  Compute the angle between left and right flame front lines.
%
%   [angle_deg, line_left, line_right, r_left_trunc, r_right_trunc, z_left_trunc, z_right_trunc] = ...
%       flameFrontAngle(flame_left, flame_right, flame_center, px2mm_x, px2mm_y, trunc_frac)
%
%   Inputs:
%       flame_left   - Nx1 array of left front x-positions (pixels, NaN where undetected)
%       flame_right  - Nx1 array of right front x-positions (pixels, NaN where undetected)
%       flame_center - scalar, mean flame center x-position (pixels)
%       px2mm_x      - pixel-to-mm conversion factor in x (mm/px)
%       px2mm_y      - pixel-to-mm conversion factor in y (mm/px)
%       trunc_frac   - fraction to keep from each end (e.g. 0.2 keeps only the 20%-80% range)
%
%   Outputs:
%       angle_deg      - angle between the two fitted lines [degrees]
%       line_left      - [slope, intercept] of the left front fit (r = slope*z + intercept)
%       line_right     - [slope, intercept] of the right front fit
%       r_left_trunc   - left radii used for fit [mm]
%       r_right_trunc  - right radii used for fit [mm]
%       z_left_trunc   - z coordinates of valid left points after truncation [mm]
%       z_right_trunc  - z coordinates of valid right points after truncation [mm]

    % --- Validation (same as flameFrontArea) ---
    valid = ~isnan(flame_left) | ~isnan(flame_right);
    row_min = find(valid, 1, 'first');
    row_max = find(valid, 1, 'last');
    rows = (row_min:row_max)';
    left  = flame_left(row_min:row_max);
    right = flame_right(row_min:row_max);

    % Radii in mm
    r_left  = abs(left  - flame_center) * px2mm_x;
    r_right = abs(right - flame_center) * px2mm_x;

    % Axial coordinate in mm
    z_mm = (rows - row_min) * px2mm_y;

    % --- Truncate each side independently based on its own valid points ---
    [z_left_trunc,  r_left_trunc]  = truncate_valid(z_mm, r_left,  trunc_frac);
    [z_right_trunc, r_right_trunc] = truncate_valid(z_mm, r_right, trunc_frac);

    % --- Linear fits ---
    if length(r_left_trunc) < 2 || length(r_right_trunc) < 2
        error('Not enough valid points after truncation for a linear fit.');
    end

    p_left  = polyfit(z_left_trunc,  r_left_trunc,  1);
    p_right = polyfit(z_right_trunc, r_right_trunc, 1);

    line_left  = p_left;
    line_right = p_right;

    % --- Angle between the two lines ---
    theta_left  = atan(abs(p_left(1)));
    theta_right = atan(abs(p_right(1)));
    angle_deg   = rad2deg(theta_left + theta_right);
end


function [z_trunc, r_trunc] = truncate_valid(z_mm, r, trunc_frac)
% TRUNCATE_VALID  Keep only the middle portion of valid (non-NaN) points.
%   trunc_frac = 0.2 means keep points between the 20th and 80th percentile
%   of the valid z range.

    valid_mask = ~isnan(r);
    z_valid = z_mm(valid_mask);
    r_valid = r(valid_mask);

    N = length(z_valid);
    idx_start = max(1,  floor(trunc_frac * N) + 1);
    idx_end   = min(N, floor((1 - trunc_frac) * N));

    z_trunc = z_valid(idx_start:idx_end);
    r_trunc = r_valid(idx_start:idx_end);
end