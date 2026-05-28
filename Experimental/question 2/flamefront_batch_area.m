function results = flamefront_batch_area(filenames, baseline, centers, tops, thresholds, px2mm_x, px2mm_y, tube_dia_mm, velocities, verbose)
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
    results = flamefront_batch(filenames, baseline, centers, tops, thresholds, false);

    results_angle = flameFrontAngle_batch(filenames, baseline, centers, tops, thresholds, px2mm_x, px2mm_y, 0.2, verbose);


    n = length(results);

    % Allow scalar px2mm values
    if isscalar(px2mm_x), px2mm_x = repmat(px2mm_x, n, 1); end
    if isscalar(px2mm_y), px2mm_y = repmat(px2mm_y, n, 1); end

    tube_area = pi * (tube_dia_mm / 2)^2;  % [mm^2]

    for i = 1:n
        %fprintf('\n=== Area & flame speed: image %d / %d ===\n', i, n);

        [area, r_avg, z_mm, area_integration, area_cap, flame_lengh] = flameFrontArea(results(i).flame_left, ...
                                              results(i).flame_right, ...
                                              results(i).flame_center, ...
                                              px2mm_x(i), px2mm_y(i));

        V = velocities(i) * tube_area;  % [mm^3/s]

        
        results(i).r_avg = r_avg;    % [mm]
        results(i).z_mm  = z_mm;     % [mm]
        results(i).vel = velocities(i); % [mm/s]
        results(i).area_int = area_integration;
        results(i).area_cap = area_cap;
        results(i).area  = area;     % [mm^2]
        results(i).Sf_areas    = V / area; % [mm/s]
        results(i).SF_angle = velocities(i) * sin(deg2rad(results_angle(i).angle_deg)/2);
        results(i).SF_perfect_tecnic = (4 * V * sin(deg2rad(results_angle(i).angle_deg)/2)) / (pi  * (2*r_avg (end))^2 );

        %fprintf('Velocity = %.2f mm/s | Area = %.2f mm^2 | Sf = %.4f mm/s\n', ...
        %       velocities(i), area, results(i).Sf);
    end
end