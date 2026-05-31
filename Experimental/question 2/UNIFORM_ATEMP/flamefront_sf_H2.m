function results = flamefront_sf_H2(filenames, px2mm_x, px2mm_y, threshold, drop_frac, trunc_frac, velocities, verbose)
% FLAMEFRONT_SF_H2  Batch flame speed for H2 using unified gradient-guided detection.
%
%   Inputs:
%       filenames  - cell array of .b16 file paths
%       px2mm_x    - mm/px in x  (1 / px_per_mm_x from calibratePixelPerMM)
%       px2mm_y    - mm/px in y
%       threshold  - local normalised threshold (default 0.25)
%       drop_frac  - gradient drop fraction for working-point detection (default 0.80)
%       trunc_frac - fraction trimmed from each fit end (default 0.10)
%       velocities - mean bulk velocity per image [mm/s]
%       verbose    - passed to flameFrontDetect
%
%   Output fields per image:
%       filename, flame_left, flame_right, angle_deg, SF_from_angle

    if nargin < 8
        verbose = true;
    end

    n       = length(filenames);
    results = struct('filename', cell(1,n), 'flame_left', cell(1,n), ...
                     'flame_right', cell(1,n), 'angle_deg', cell(1,n), ...
                     'SF_from_angle', cell(1,n));

    for i = 1:n
        stats = flameFrontDetect(filenames{i}, ...
            1/px2mm_x, 1/px2mm_y, threshold, drop_frac, trunc_frac, verbose);

        results(i).filename      = filenames{i};
        results(i).flame_left    = stats.left_edge;
        results(i).flame_right   = stats.right_edge;
        results(i).angle_deg     = stats.alpha_deg;
        results(i).SF_from_angle = velocities(i) * sin(deg2rad(stats.alpha_deg) / 2);
    end
end
