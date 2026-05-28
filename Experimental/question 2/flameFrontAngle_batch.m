function results = flameFrontAngle_batch(filenames, px2mm_x, px2mm_y, limit, trunc_frac, verbose)
% FLAMEFRONTANGLE_BATCH  Batch flame front angle detection and fitting.
%
%   results = flameFrontAngle_batch(filenames, baseline, centers, tops, ...
%       thresholds, px2mm_x, px2mm_y, trunc_frac, verbose)
%
%   Inputs:
%       filenames   - cell array of B16 file paths
%       baseline    - Nx1 array of row_min per image
%       centers     - Nx1 array of center_col per image
%       tops        - Nx1 array of row_max per image
%       thresholds  - Nx1 array of intensity thresholds per image
%       px2mm_x     - mm/px in x (scalar, same for all)
%       px2mm_y     - mm/px in y (scalar, same for all)
%       trunc_frac  - fraction to remove from each end (e.g. 0.2 keeps 20%-80%)
%       verbose     - logical, if true plot results (default: true)
%
%   Outputs:
%       results - struct array with fields:
%           filename     - file path
%           flame_left   - Nx1 left front x-positions (pixels)
%           flame_center - mean flame center (pixels)
%           flame_right  - Nx1 right front x-positions (pixels)
%           angle_deg    - opening angle [degrees]
%           line_left    - [slope, intercept] of left fit in mm
%           line_right   - [slope, intercept] of right fit in mm

    if nargin < 9
        verbose = true;
    end

    n = length(filenames);
    results = struct('filename', {}, 'flame_left', {}, 'flame_center', {}, ...
                     'flame_right', {}, 'angle_deg', {}, 'line_left', {}, 'line_right', {});

    for i = 1:n

        stats = flameFrontStats(filenames{i}, 1/px2mm_x, 1/px2mm_y, limit, trunc_frac, verbose);
    
        results(i).filename     = filenames{i};
        results(i).flame_left   = stats.left_edge ;
        results(i).flame_right  = stats.right_edge;
        results(i).angle_deg = stats.alpha_deg;
    end
end