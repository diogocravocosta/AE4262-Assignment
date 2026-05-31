function results = flamefront_sf_CH4(filenames, px2mm_x, px2mm_y, threshold, trunc_frac, velocities, verbose)
% FLAMEFRONT_BATCH_AREA  Detect flame fronts, compute surface areas and flame speed.
%
%   results = flamefront_batch_area(filenames, baseline, centers, tops, thresholds, ...
%                                   px2mm_x, px2mm_y, tube_dia_mm, velocities, verbose)
%
%   Inputs:
%       filenames    - cell array of file paths
%       baseline     - array of row_min values per image
%       centers      - array of center_col values per image
%       tops         - array of row_max values per image
%       thresholds   - array of threshold values per image
%       px2mm_x      - pixel-to-mm factor in x (scalar or per image)
%       px2mm_y      - pixel-to-mm factor in y (scalar or per image)
%       tube_dia_mm  - inner tube diameter [mm] (scalar)
%       velocities   - array of mean flow velocities per image [mm/s]
%       verbose      - logical, passed to flameFront (default: true)
%
%   Output fields per image:
%       filename, flame_left, flame_center, flame_right,
%       area, r_avg, z_mm, Sf

    if nargin < 10
        verbose = true;
    end

    % Run the existing batch detection



    n = length(filenames);
    results = struct('filename', {});




    for i = 1:n

        stats = flameFrontStats(filenames{i}, 1/px2mm_x, 1/px2mm_y, threshold, trunc_frac, verbose);
    
        results(i).filename     = filenames{i};
        results(i).flame_left   = stats.left_edge ;
        results(i).flame_right  = stats.right_edge;
        results(i).angle_deg = stats.alpha_deg;
     
        results(i).SF_from_angle = velocities(i) * sin(deg2rad(results(i).angle_deg)/2);

    end
end