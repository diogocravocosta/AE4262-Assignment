function stats = flameFrontStats(filepath, px_per_mm_x, px_per_mm_y, threshold, trunc_frac,verbose)
% FLAMEFRONTSTATS  Detect flame front and return all quantities needed
%                  for both flame-speed calculation methods.
%
%   stats = flameFrontStats(filepath)
%   stats = flameFrontStats(filepath, px_per_mm_x, px_per_mm_y, threshold, trunc_frac)
%
%   Inputs
%   ------
%   filepath     - path to .b16 image file
%   px_per_mm_x  - calibration: pixels per mm in x (horizontal)  (default 29.706)
%   px_per_mm_y  - calibration: pixels per mm in y (vertical)    (default 29.695)
%   threshold    - normalised intensity threshold 0-1             (default 0.25)
%   trunc_frac   - fraction trimmed from each end of
%                  the below-tip region before fitting            (default 0.10)
%
%   Output fields (struct)
%   ----------------------
%   alpha_half_deg   half-angle alpha/2  [deg]  (in physical mm space)
%   alpha_deg        full cone angle alpha = 2*(alpha/2)  [deg]
%   Dc_mm            flame-cone base width  [mm]
%                    (mean inner-front separation over bottom 10% of rows)
%   tip_row          tip row index  [px]
%   tip_col          tip column index  [px]
%   rows             row indices where both fronts were found  [Nx1]
%   left_edge        left  front x-positions [px, Nx1]
%   right_edge       right front x-positions [px, Nx1]
%   c_left           [slope intercept] of left  edge fit (pixel space)
%   c_right          [slope intercept] of right edge fit (pixel space)

    % --- defaults ---
    if nargin < 2 || isempty(px_per_mm_x), px_per_mm_x = 29.706; end
    if nargin < 3 || isempty(px_per_mm_y), px_per_mm_y = 29.695; end
    if nargin < 4 || isempty(threshold),   threshold   = 0.25;   end
    if nargin < 5 || isempty(trunc_frac),  trunc_frac  = 0.10;   end
    if nargin < 6 || isempty(verbose),     verbose     = true;   end

    % --- load & normalise ---
    image = readB16(filepath);
    img_d = double(image);
    img_n = (img_d - min(img_d(:))) / (max(img_d(:)) - min(img_d(:)));

    % gaussian smooth (sigma = 3 px)
    h     = fspecial('gaussian', [15 15], 3);
    img_s = imfilter(img_n, h, 'replicate');

    [nRows, nCols] = size(img_s);
    % --- Calculate intensity-weighted center column ---
    col_indices = 1:nCols;
    % Sum the brightness of each column across all rows
    col_sums = sum(img_s, 1); 
    % Calculate the weighted average (center of mass)
    center_col = round(sum(col_indices .* col_sums) / sum(col_sums));

    % --- walk outward from centre for every row ---
    left_all  = NaN(nRows, 1);
    right_all = NaN(nRows, 1);

    for r = 51 : nRows - 50
        row = img_s(r, :);
        if max(row) < 0.10, continue; end

        % left: walk from centre toward col 1
        for c = center_col : -1 : 2
            if row(c) > threshold
                left_all(r) = c;
                break;
            end
        end

        % right: walk from centre toward last col
        for c = center_col : nCols - 1
            if row(c) > threshold
                right_all(r) = c;
                break;
            end
        end
    end

    % keep only rows where BOTH fronts were found
    valid      = ~isnan(left_all) & ~isnan(right_all) & (right_all > left_all);
    rows       = find(valid);
    left_edge  = left_all(rows);
    right_edge = right_all(rows);
    widths_px  = right_edge - left_edge;

    % --- tip: narrowest width ---
    [~, tip_idx] = min(widths_px);
    tip_row = rows(tip_idx);
    tip_col = round((left_edge(tip_idx) + right_edge(tip_idx)) / 2);

    % --- work only with rows BELOW the tip ---
    below    = rows > tip_row;
    rows_b   = rows(below);
    left_b   = left_edge(below);
    right_b  = right_edge(below);
    widths_b = widths_px(below);

    n  = numel(rows_b);
    i0 = max(1, floor(trunc_frac * n) + 1);
    i1 = min(n, floor((1 - trunc_frac) * n));

    % --- linear fits in pixel space ---
    c_left  = polyfit(rows_b(i0:i1), left_b(i0:i1),  1);  % [slope intercept]
    c_right = polyfit(rows_b(i0:i1), right_b(i0:i1), 1);

    % Convert pixel-space slope to physical slope before taking arctan:
    %   slope_px = dx_px / dy_px
    %   slope_mm = (dx_px / px_per_mm_x) / (dy_px / px_per_mm_y)
    %            = slope_px * (px_per_mm_y / px_per_mm_x)
    slope_correction = px_per_mm_y / px_per_mm_x;

    a_left  = rad2deg(atan(abs(c_left(1))  * slope_correction));
    a_right = rad2deg(atan(abs(c_right(1)) * slope_correction));
    alpha_half = (a_left + a_right) / 2;

    % --- cone base width: mean of bottom 10% of below-tip rows ---
    %     width is horizontal -> use px_per_mm_x
    base_start = floor(0.90 * n) + 1;
    Dc_mm = mean(widths_b(base_start:end)) / px_per_mm_x;

% --- pack output ---
    alpha_deg = 2 * alpha_half;
    stats.alpha_half_deg = alpha_half;
    stats.alpha_deg      = alpha_deg;
    stats.Dc_mm          = Dc_mm;
    stats.tip_row        = tip_row;
    stats.tip_col        = tip_col;
    stats.rows           = rows;
    stats.left_edge      = left_edge;
    stats.right_edge     = right_edge;
    stats.c_left         = c_left;
    stats.c_right        = c_right;

    % =========================================================================
    % --- PLOTTING & VERBOSE OUTPUT (Conditional) ---
    % =========================================================================
    if verbose
        parts = strsplit(filepath, {'\', '/'});
        short_name = parts{end-2};
        
        % --- Figure 1: Intensity Profiles ---
        figure('Name', ['Intensity profiles - ', short_name]);
        hold on;
        xlabel('Pixel index from center');
        ylabel('Intensity');
        title({'Right-half intensity profile for each row', short_name}, 'FontSize', 8);
        yline(threshold, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');
        
        % Highlight specific rows for clean visualization
        n_highlight = 5;
        highlight_rows = round(linspace(rows_b(i0), rows_b(i1), n_highlight));
        cmap = lines(n_highlight);
        
        for r = 51 : nRows - 50
            row_data = img_s(r, :);
            if max(row_data) < 0.10, continue; end
            right_half = row_data(center_col:end);
            
            hi = find(highlight_rows == r, 1);
            if ~isempty(hi)
                plot(1:length(right_half), right_half, 'Color', cmap(hi,:), ...
                    'LineWidth', 1.2, 'DisplayName', sprintf('Row %d', r));
            else
                plot(1:length(right_half), right_half, 'Color', [0 0 1 0.1], ...
                    'HandleVisibility', 'off');
            end
            
            % Mark detected right-edge pixel crossing
            if ~isnan(right_all(r))
                idx = right_all(r) - center_col + 1;
                if idx > 0 && idx <= length(right_half)
                    plot(idx, right_half(idx), 'ro', 'MarkerSize', 3, 'HandleVisibility', 'off');
                end
            end
        end
        hold off;
    
        % --- Figure 2: Flame Front Angle Image Overlay ---
        rows_fit = (rows_b(i0):rows_b(i1))';
        
        % Reconstruct fits in pixel space
        x_plot_L = polyval(c_left,  rows_fit);
        x_plot_R = polyval(c_right, rows_fit);
    
        figure('Name', ['Flame Front Angle - ', short_name]);
        imagesc(image); colormap(hot); colorbar; hold on;
    
        % Plot raw tracked edges across ALL valid rows
        plot(left_edge,  rows, 'r', 'LineWidth', 1, 'DisplayName', 'Left front');
        plot(right_edge, rows, 'b', 'LineWidth', 1, 'DisplayName', 'Right front');
    
        % Calculate individual angles for the legend
        a_left  = rad2deg(atan(abs(c_left(1))  * slope_correction));
        a_right = rad2deg(atan(abs(c_right(1)) * slope_correction));
    
        % Plot the trimmed linear regression lines
        plot(x_plot_L, rows_fit, 'g-', 'LineWidth', 2, 'DisplayName', ...
            sprintf('Left fit (%.2f°)', a_left));
        plot(x_plot_R, rows_fit, 'm-', 'LineWidth', 2, 'DisplayName', ...
            sprintf('Right fit (%.2f°)', a_right));
    
        % Reference markers
        xline(center_col, 'c--', 'LineWidth', 1, 'DisplayName', 'Center');
        yline(rows_b(i0), 'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        yline(rows_b(i1), 'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
        axis image;
        title({short_name, sprintf('Angle: %.1f°', alpha_deg)}, 'FontSize', 10);
        legend('Location', 'best');
        hold off;
    
        % --- Terminal Summary Output ---
        fprintf('--- Flame Front Angle (Code 1 Processing Engine) ---\n');
        fprintf('File: %s\n', filepath);
        fprintf('Truncation: %.0f%% - %.0f%% of below-tip regions\n', trunc_frac*100, (1-trunc_frac)*100);
        fprintf('Left  fit: x_px = %.4f * row + %.4f\n', c_left(1), c_left(2));
        fprintf('Right fit: x_px = %.4f * row + %.4f\n', c_right(1), c_right(2));
        fprintf('Opening angle (Full Cone): %.2f°\n', alpha_deg);
        fprintf('Points used in fit: %d rows\n', length(rows_fit));
    end
end