function [area, r_avg, z_mm, area_integration, area_cap, avg_line_length] = flameFrontArea(flame_left, flame_right, flame_center, px2mm_x, px2mm_y)
% FLAMEFRONTAREA  Compute flame front surface area assuming axisymmetry.
%
%   Outputs:
%       ...
%       avg_line_length  - arc length of the average flame front profile [mm]

    % Find rows where at least one side has a detection
    valid = ~isnan(flame_left) | ~isnan(flame_right);
    row_min = find(valid, 1, 'first');
    row_max = find(valid, 1, 'last');
    rows = (row_min:row_max)';
    left  = flame_left(row_min:row_max);
    right = flame_right(row_min:row_max);

    % Radii in mm
    r_left  = abs(left  - flame_center) * px2mm_x;
    r_right = abs(right - flame_center) * px2mm_x;

    % Average radius (ignoring NaNs where only one side was detected)
    r_avg = mean([r_left, r_right], 2, 'omitnan');

    % Axial coordinate in mm
    z_mm = (rows - row_min) * px2mm_y;

    % Remove any interior rows where neither side was detected
    interior_valid = ~isnan(r_avg);
    r_avg = r_avg(interior_valid);
    z_mm  = z_mm(interior_valid);

    % dr/dz via central differences
    dr_dz = gradient(r_avg, z_mm);

    % Surface area of revolution
    integrand = 2 * pi * r_avg .* sqrt(1 + dr_dz.^2);
    area_integration = trapz(z_mm, integrand);

    % Circular cap at the open top
    area_cap = pi * r_avg(1)^2;

    % Total area
    area = area_integration + area_cap;

    % Arc length of the average flame front line: L = integral( sqrt(dr^2 + dz^2) )
    dr = diff(r_avg);
    dz = diff(z_mm);
    avg_line_length = sum(sqrt(dr.^2 + dz.^2));
end