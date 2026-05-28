function [angle_deg, line_left, line_right, flame_left, flame_right, flame_center] = flameFrontAnglePlot(filename, center_col, limit, row_min, row_max, px2mm_x, px2mm_y, trunc_frac, verbose)
% FLAMEFRONTANGLEPLOT  Detect flame fronts, fit lines, and compute opening angle.
%
%   [angle_deg, line_left, line_right, flame_left, flame_right, flame_center] = ...
%       flameFrontAnglePlot(filename, center_col, limit, row_min, row_max, ...
%       px2mm_x, px2mm_y, trunc_frac, verbose)
%
%   Inputs:
%       filename     - path to B16 image file
%       center_col   - center column of flame (pixels)
%       limit        - intensity threshold for flame detection
%       row_min      - first row to scan
%       row_max      - last row to scan
%       px2mm_x      - mm/px in x
%       px2mm_y      - mm/px in y
%       trunc_frac   - fraction to remove from each end (e.g. 0.2 keeps 20%-80%)
%       verbose      - logical, if true plot results (default: true)
%
%   Outputs:
%       angle_deg    - opening angle between left and right fitted lines [degrees]
%       line_left    - [slope, intercept] of left fit in mm (r = slope*z + intercept)
%       line_right   - [slope, intercept] of right fit in mm
%       flame_left   - Nx1 left front x-positions (pixels, NaN where undetected)
%       flame_right  - Nx1 right front x-positions (pixels, NaN where undetected)
%       flame_center - scalar, mean flame center x-position (pixels)

    if nargin < 9
        verbose = true;
    end

    % --- Read image ---
    image = readB16(filename);
    image_d = double(image);
    nRows = size(image, 1);

    % --- Detect flame fronts (same as flameFront) ---
    flame_left  = NaN(nRows, 1);
    flame_right = NaN(nRows, 1);


    if verbose
    parts = strsplit(filename, {'\', '/'});
    short_name = parts{end-2};
    figure('Name', ['Intensity profiles - ', short_name]);
    hold on;
    xlabel('Pixel index from center');
    ylabel('Intensity');
    title({'Right-half intensity profile for each row', short_name}, 'FontSize', 8);
    yline(limit, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');

    % Set number of highlighted lines (0 = all faint blue, >0 = that many coloured)
    n_highlight = 5;
    if n_highlight > 0
        highlight_rows = round(linspace(row_min, row_max, n_highlight));
        cmap = lines(n_highlight);
    else
        highlight_rows = [];
    end
    end

    for row = row_min:row_max
        line = image_d(row, :);

        % Right side
        right_half = line(center_col:end);
        idx = find(right_half > limit, 1, 'first');
        if ~isempty(idx)
            flame_right(row) = center_col + idx - 1;
        end

        if verbose
            hi = find(highlight_rows == row, 1);
            if ~isempty(hi)
                plot(1:length(right_half), right_half, 'Color', cmap(hi,:), ...
                    'LineWidth', 1.2, 'DisplayName', sprintf('Row %d', row));
            else
                plot(1:length(right_half), right_half, 'Color', [0 0 1 0.1], ...
                    'HandleVisibility', 'off');
            end
            if ~isempty(idx)
                plot(idx, right_half(idx), 'ro', 'MarkerSize', 3, 'HandleVisibility', 'off');
            end
        end

        % Left side
        left_half = fliplr(line(1:center_col));
        idx = find(left_half > limit, 1, 'first');
        if ~isempty(idx)
            flame_left(row) = center_col - idx + 1;
        end
    end

    flame_center = mean([flame_left; flame_right], 'omitnan');

    % --- Validation ---
    valid = ~isnan(flame_left) | ~isnan(flame_right);
    first_valid = find(valid, 1, 'first');
    last_valid  = find(valid, 1, 'last');
    rows = (first_valid:last_valid)';
    left  = flame_left(first_valid:last_valid);
    right = flame_right(first_valid:last_valid);

    % --- Truncate each side independently on its valid points ---
    [rows_left_trunc,  left_trunc]  = truncate_valid(rows, left,  trunc_frac);
    [rows_right_trunc, right_trunc] = truncate_valid(rows, right, trunc_frac);

    if length(left_trunc) < 2 || length(right_trunc) < 2
        error('Not enough valid points after truncation for a linear fit.');
    end

    % --- Scale to mm ---
    r_left_mm  = abs(left_trunc  - flame_center) * px2mm_x;
    r_right_mm = abs(right_trunc - flame_center) * px2mm_x;
    z_left_mm  = (rows_left_trunc  - first_valid) * px2mm_y;
    z_right_mm = (rows_right_trunc - first_valid) * px2mm_y;

    % --- Linear fits in mm space ---
    p_left  = polyfit(z_left_mm,  r_left_mm,  1);
    p_right = polyfit(z_right_mm, r_right_mm, 1);

    line_left  = p_left;
    line_right = p_right;

    % --- Angle ---
    theta_left  = atan(abs(p_left(1)));
    theta_right = atan(abs(p_right(1)));
    angle_deg   = rad2deg(theta_left + theta_right);

    % --- Plotting ---
    if verbose
            z_plot_L = linspace(z_left_mm(1), z_left_mm(end), 100)';
            r_fit_L  = polyval(p_left, z_plot_L);
            row_plot_L = z_plot_L / px2mm_y + first_valid;
            x_plot_L   = flame_center - r_fit_L / px2mm_x;
    
            z_plot_R = linspace(z_right_mm(1), z_right_mm(end), 100)';
            r_fit_R  = polyval(p_right, z_plot_R);
            row_plot_R = z_plot_R / px2mm_y + first_valid;
            x_plot_R   = flame_center + r_fit_R / px2mm_x;
    
            parts = strsplit(filename, {'\', '/'});
            short_name = parts{end-2};
    
            figure('Name', ['Flame Front Angle - ', short_name]);
            imagesc(image); colormap(gray); colorbar; hold on;
    
            plot(flame_left,  1:nRows, 'r', 'LineWidth', 1, 'DisplayName', 'Left front');
            plot(flame_right, 1:nRows, 'b', 'LineWidth', 1, 'DisplayName', 'Right front');
    
            plot(x_plot_L, row_plot_L, 'g-', 'LineWidth', 2, 'DisplayName', ...
                sprintf('Left fit (%.2f°)', rad2deg(theta_left)));
            plot(x_plot_R, row_plot_R, 'm-', 'LineWidth', 2, 'DisplayName', ...
                sprintf('Right fit (%.2f°)', rad2deg(theta_right)));
    
            xline(flame_center, 'c--', 'LineWidth', 1, 'DisplayName', 'Center');
    
            yline(rows_left_trunc(1),   'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
            yline(rows_left_trunc(end), 'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
            yline(rows_right_trunc(1),   'm--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
            yline(rows_right_trunc(end), 'm--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
            axis image;
            title({short_name, sprintf('Angle: %.1f°', angle_deg)}, 'FontSize', 10);
            legend('Location', 'best');
            hold off;
    
            fprintf('--- Flame Front Angle ---\n');
            fprintf('File: %s\n', filename);
            fprintf('Truncation: %.0f%% - %.0f%% of valid points\n', trunc_frac*100, (1-trunc_frac)*100);
            fprintf('Left  fit: r = %.4f * z + %.4f mm  (angle: %.2f°)\n', p_left(1), p_left(2), rad2deg(theta_left));
            fprintf('Right fit: r = %.4f * z + %.4f mm  (angle: %.2f°)\n', p_right(1), p_right(2), rad2deg(theta_right));
            fprintf('Opening angle: %.2f°\n', angle_deg);
            fprintf('Left  points used: %d\n', length(r_left_mm));
            fprintf('Right points used: %d\n', length(r_right_mm));
     end
end


function [rows_trunc, vals_trunc] = truncate_valid(rows, vals, trunc_frac)
% Keep only the middle portion of valid (non-NaN) points.
% trunc_frac = 0.2 keeps the 20%-80% range.
    valid_mask = ~isnan(vals);
    rows_valid = rows(valid_mask);
    vals_valid = vals(valid_mask);

    N = length(vals_valid);
    idx_start = max(1, floor(trunc_frac * N) + 1);
    idx_end   = min(N, floor((1 - trunc_frac) * N));

    rows_trunc = rows_valid(idx_start:idx_end);
    vals_trunc = vals_valid(idx_start:idx_end);
end