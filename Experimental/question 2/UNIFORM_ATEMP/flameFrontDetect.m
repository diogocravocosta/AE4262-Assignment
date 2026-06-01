function stats = flameFrontDetect(filepath, px_per_mm_x, px_per_mm_y, threshold, drop_frac, trunc_tip, trunc_base, verbose_diag, verbose_overlay, fit_lines_only)
% FLAMEFRONTDETECT  Unified gradient-guided flame front detection for CH4 and H2.
%
%   Handles flames with double OH-brightness peaks (rich H2) by using the
%   smoothed-intensity gradient rather than the raw intensity peak:
%
%   Per row (within auto-detected flame boundaries), scanning outward from axis:
%     1. Row-average of smoothed image → locate flame top and bottom rows
%     2. Gaussian-smoothed gradient → steepest rising point (peak)
%     3. After the peak: find where gradient drops to drop_frac of its max
%        (x_drop) — robust to OH plateaus and double-peak profiles
%     4. Re-normalise RAW intensity locally in window around x_drop,
%        find first crossing of threshold → precise flame edge
%
%   Inputs
%   ------
%   filepath        - path to .b16 image file
%   px_per_mm_x     - pixels per mm in x  (default 29.706)
%   px_per_mm_y     - pixels per mm in y  (default 29.695)
%   threshold       - local normalised threshold 0-1  (default 0.25)
%   drop_frac       - gradient drop fraction to locate working point (default 0.80)
%   trunc_frac      - fraction trimmed from each end of below-tip rows (default 0.10)
%   verbose_diag    - if true: plot row-avg, intensity, and gradient figures (Figs 1-3)
%   verbose_overlay - if true: plot flame image with front overlay and fit (Fig 4)
%
%   Output struct fields  (identical to flameFrontStats)
%   ----------------------------------------------------
%   alpha_half_deg, alpha_deg, Dc_mm
%   tip_row, tip_col
%   rows, left_edge, right_edge
%   c_left, c_right   ([slope intercept] pixel-space fits)
%   row_min, row_max  (detected flame boundaries)

    if nargin < 2 || isempty(px_per_mm_x),    px_per_mm_x    = 29.706; end
    if nargin < 3 || isempty(px_per_mm_y),    px_per_mm_y    = 29.695; end
    if nargin < 4 || isempty(threshold),       threshold      = 0.25;   end
    if nargin < 5 || isempty(drop_frac),       drop_frac      = 0.80;   end
    if nargin < 6  || isempty(trunc_tip),       trunc_tip       = 0.10;  end
    if nargin < 7  || isempty(trunc_base),     trunc_base      = 0.10;  end
    if nargin < 8  || isempty(verbose_diag),   verbose_diag    = true;  end
    if nargin < 9  || isempty(verbose_overlay), verbose_overlay = true; end
    if nargin < 10 || isempty(fit_lines_only), fit_lines_only  = false; end

    % --- load & normalise ---
    image      = readB16(filepath);
    img_d      = double(image);
    img_s      = (img_d - min(img_d(:))) / (max(img_d(:)) - min(img_d(:)));
    h_filt     = fspecial('gaussian', [15 15], 3);
    img_smooth = imfilter(img_s, h_filt, 'replicate');

    [nRows, nCols] = size(img_s);

    % --- intensity-weighted center column ---
    col_idx    = 1:nCols;
    col_sums   = sum(img_smooth, 1);
    center_col = round(sum(col_idx .* col_sums) / sum(col_sums));

    % --- flame row boundaries from row-average intensity ---
    row_avg            = mean(img_smooth, 2);
    row_thresh_frac_up = 0.35;
    row_thresh_up      = row_thresh_frac_up * max(row_avg);
    flame_rows_up      = find(row_avg > row_thresh_up);
    
    if isempty(flame_rows_up)
        error('flameFrontDetect: no flame rows detected for top boundary in %s', filepath);
    end
    row_min = max(2, flame_rows_up(1));
    
    % row_max: steepest negative gradient of row-average, detected independently
    % for each half; the one further up (smaller row index) is used as the limit.
    search_region  = row_min : nRows;

    row_avg_left   = mean(img_smooth(:, 1:center_col),    2);
    row_avg_right  = mean(img_smooth(:, center_col:end),  2);

    [~, idx_left]  = max(-gradient(row_avg_left(search_region)));
    [~, idx_right] = max(-gradient(row_avg_right(search_region)));

    row_max_left   = search_region(idx_left);
    row_max_right  = search_region(idx_right);
    row_max        = min(nRows-1, min(row_max_left, row_max_right));

    % --- per-row detection (within flame boundaries only) ---
    left_all  = NaN(nRows, 1);
    right_all = NaN(nRows, 1);

    for r = row_min : row_max
        % right side: scan outward from center toward right edge
        right_raw    = img_s(r, center_col:end);
        right_smooth = img_smooth(r, center_col:end);
        k = edge_from_gradient(right_raw, right_smooth, threshold, drop_frac);
        if ~isnan(k)
            right_all(r) = center_col + k - 1;
        end

        % left side: flip so index 1 = center, higher index = further left
        left_raw    = fliplr(img_s(r, 1:center_col));
        left_smooth = fliplr(img_smooth(r, 1:center_col));
        k = edge_from_gradient(left_raw, left_smooth, threshold, drop_frac);
        if ~isnan(k)
            left_all(r) = center_col - k + 1;
        end
    end

    % keep only rows where BOTH fronts were found
    valid      = ~isnan(left_all) & ~isnan(right_all) & (right_all > left_all);
    rows       = find(valid);
    left_edge  = left_all(rows);
    right_edge = right_all(rows);
    widths_px  = right_edge - left_edge;

    if isempty(rows)
        error('flameFrontDetect: no valid rows found in %s', filepath);
    end

   % --- tip: 5% of original baseline width scanning backward ---
    % Original distance is measured at the widest/lowest valid rows of the flame
    original_distance = widths_px(end); 
    
    tip_found = false;
    pct_threshold = 0.05; % Start at 5%
    
    while ~tip_found && pct_threshold <= 1.0
        % Scan in reverse: from end (bottom) to beginning (top)
        for idx = length(widths_px) : -1 : 1
            if widths_px(idx) <= pct_threshold * original_distance
                tip_idx = idx;
                tip_found = true;
                break; % Found the first crossing moving backward!
            end
        end
        
        if ~tip_found
            % If no row is narrow enough, incrementally ease the threshold
            pct_threshold = pct_threshold + 0.01; 
        end
    end
    
    % Fallback if something goes completely wrong
    if ~tip_found
        [~, tip_idx] = min(widths_px); % Fallback to absolute minimum
        warning('flameFrontDetect: Could not find 5%% threshold. Using absolute minimum.');
    elseif pct_threshold > 0.05
        warning('flameFrontDetect: 5%% tip criteria not met. Dynamically increased threshold to %.0f%%.', pct_threshold * 100);
    end

    tip_row = rows(tip_idx);
    tip_col = round((left_edge(tip_idx) + right_edge(tip_idx)) / 2);

    % --- work only with rows BELOW the tip ---
    below    = rows > tip_row;
    rows_b   = rows(below);
    left_b   = left_edge(below);
    right_b  = right_edge(below);
    widths_b = widths_px(below);

    nb  = numel(rows_b);
    i0 = max(1, floor(trunc_tip  * nb) + 1);
    i1 = min(nb, floor((1 - trunc_base) * nb));


    if nb < 5
        warning('flameFrontDetect: only %d rows below cone tip — lower row_thresh_frac or check image in %s', nb, filepath);

    end


    % --- linear fits in pixel space ---
    c_left  = polyfit(rows_b(i0:i1), left_b(i0:i1),  1);
    c_right = polyfit(rows_b(i0:i1), right_b(i0:i1), 1);

    % convert pixel-space slopes to physical angles
    slope_corr = px_per_mm_y / px_per_mm_x;
    a_left  = rad2deg(atan(abs(c_left(1))  * slope_corr));
    a_right = rad2deg(atan(abs(c_right(1)) * slope_corr));
    alpha_half = (a_left + a_right) / 2;
    alpha_deg  = 2 * alpha_half;

    % cone base width (bottom 10% of below-tip rows)
    base_start = floor(0.90 * nb) + 1;
    Dc_mm = mean(widths_b(base_start:end)) / px_per_mm_x;

    % --- pack output ---
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
    stats.row_min        = row_min;
    stats.row_max        = row_max;

    % =========================================================
    % --- verbose plots ---
    % =========================================================
    if verbose_diag || verbose_overlay
        parts      = strsplit(filepath, {'\', '/'});
        short_name = parts{end-2};

        % highlight rows evenly spaced by count within the fitted region
        n_highlight    = 5;
        fit_rows       = rows_b(i0:i1);
        hi_idx         = round(linspace(1, length(fit_rows), n_highlight));
        highlight_rows = fit_rows(hi_idx);
        cmap           = lines(n_highlight);
    end

    if verbose_diag
        % Figure 1: row-average profile with detected flame boundaries
        figure('Name', ['Row average - ', short_name]);
        plot(1:nRows, row_avg,       'b',  'LineWidth', 1.2, 'DisplayName', 'Full width'); hold on;
        plot(1:nRows, row_avg_left,  '--', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.2, 'DisplayName', 'Left half');
        plot(1:nRows, row_avg_right, '--', 'Color', [0.8 0.4 0.0], 'LineWidth', 1.2, 'DisplayName', 'Right half');
        yline(row_thresh_up, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('%.0f%% threshold (row\\_min)', row_thresh_frac_up*100));
        xline(row_min,       'b-',  'LineWidth', 1.5, 'DisplayName', sprintf('Flame top (Row %d)', row_min));
        xline(row_max_left,  'Color', [0.2 0.7 0.2], 'LineWidth', 1.5, 'DisplayName', sprintf('row\\_max left (Row %d)',  row_max_left));
        xline(row_max_right, 'Color', [0.8 0.4 0.0], 'LineWidth', 1.5, 'DisplayName', sprintf('row\\_max right (Row %d)', row_max_right));
        xline(row_max,       'm-',  'LineWidth', 2.0, 'DisplayName', sprintf('row\\_max used (Row %d)', row_max));
        xlabel('Row index'); ylabel('Mean intensity (smoothed)');
        title({'Row-average intensity — flame boundaries', short_name}, 'FontSize', 8);
        legend('Location', 'best'); grid on; hold off;

        % Figure 2: right-half intensity profiles (unfiltered, detection dot)
        h_intens = figure('Name', ['Intensity profiles - ', short_name]);
        hold on;
        xlabel('Pixel index from center');
        ylabel('Intensity (normalised, unfiltered)');
        title({'Right-half intensity profile for each row', short_name}, 'FontSize', 8);
        yline(threshold, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');

        % Figure 3: right-half gradient profiles (smoothed)
        h_grad = figure('Name', ['Intensity gradient - ', short_name]);
        hold on;
        xlabel('Pixel index from center');
        ylabel('d(Intensity)/dx  (smoothed)');
        title({'Right-half intensity gradient for each row', short_name}, 'FontSize', 8);
        yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        for r = row_min : row_max
            right_raw  = img_s(r, center_col:end);
            right_sm   = img_smooth(r, center_col:end);
            grad_right = gradient(double(right_sm));
            hi         = find(highlight_rows == r, 1);

            figure(h_intens);
            if ~isempty(hi)
                plot(1:length(right_raw), right_raw, 'Color', cmap(hi,:), ...
                    'LineWidth', 1.2, 'DisplayName', sprintf('Row %d', r));
            % else
            %     plot(1:length(right_raw), right_raw, 'Color', [0 0 1 0.1], ...
            %         'HandleVisibility', 'off');
            end
            if ~isnan(right_all(r))
                idx = right_all(r) - center_col + 1;
                if idx >= 1 && idx <= length(right_raw)
                    if ~isempty(hi)
                        plot(idx, right_raw(idx), 's', 'Color', cmap(hi,:), ...
                            'MarkerFaceColor', cmap(hi,:), 'MarkerSize', 7, 'HandleVisibility', 'off');
                        [~, xd] = edge_from_gradient(right_raw, right_sm, threshold, drop_frac);
                        if ~isnan(xd) && xd >= 1 && xd <= length(right_raw)
                            plot(xd, right_raw(xd), '^', 'Color', cmap(hi,:), ...
                                'MarkerFaceColor', cmap(hi,:), 'MarkerSize', 7, 'HandleVisibility', 'off');
                        end
                    % else
                    %     plot(idx, right_raw(idx), 'ro', 'MarkerSize', 3, 'HandleVisibility', 'off');
                    end
                end
            end

            figure(h_grad);
            if ~isempty(hi)
                plot(1:length(grad_right), grad_right, 'Color', cmap(hi,:), ...
                    'LineWidth', 1.2, 'DisplayName', sprintf('Row %d', r));
            % else
            %     plot(1:length(grad_right), grad_right, 'Color', [0 0 1 0.1], ...
            %         'HandleVisibility', 'off');
            end
            if ~isnan(right_all(r))
                idx = right_all(r) - center_col + 1;
                if idx >= 1 && idx <= length(grad_right)
                    if ~isempty(hi)
                        plot(idx, grad_right(idx), 's', 'Color', cmap(hi,:), ...
                            'MarkerFaceColor', cmap(hi,:), 'MarkerSize', 7, 'HandleVisibility', 'off');
                    % else
                    %     plot(idx, grad_right(idx), 'ro', 'MarkerSize', 3, 'HandleVisibility', 'off');
                    end
                end
            end
        end
        figure(h_intens); hold off;
        figure(h_grad);   hold off;

        fprintf('--- flameFrontDetect (gradient-guided) ---\n');
        fprintf('File:              %s\n', filepath);
        fprintf('Flame rows:        %d  to  %d\n', row_min, row_max);
        fprintf('Gradient drop:     %.0f%% of gradient peak\n', drop_frac*100);
        fprintf('Local threshold:   %.2f (normalised)\n', threshold);
        fprintf('Truncation tip:    %.0f%%\n', trunc_tip*100);
        fprintf('Truncation base:   %.0f%%\n', trunc_base*100);
        fprintf('Left  fit:  x = %.4f*row + %.4f  (%.2f°)\n', c_left(1),  c_left(2),  a_left);
        fprintf('Right fit:  x = %.4f*row + %.4f  (%.2f°)\n', c_right(1), c_right(2), a_right);
        fprintf('Opening angle:     %.2f°\n', alpha_deg);
        fprintf('Rows used in fit:  %d\n', numel(fit_rows));
    end

    if verbose_overlay
        % Figure 4: flame front overlay on image
        rows_fit = (rows_b(i0):rows_b(i1))';
        x_plot_L = polyval(c_left,  rows_fit);
        x_plot_R = polyval(c_right, rows_fit);

        figure('Name', ['Flame Front Angle - ', short_name]);
        imagesc(image); colormap(hot); colorbar; hold on;
        if ~fit_lines_only
            plot(left_b(i0:i1),  rows_b(i0:i1), 'g', 'LineWidth', 1, 'DisplayName', 'Left front');
            plot(right_b(i0:i1), rows_b(i0:i1), 'b', 'LineWidth', 1, 'DisplayName', 'Right front');
        end
        plot(x_plot_L, rows_fit, 'g-', 'LineWidth', 2, 'DisplayName', ...
            sprintf('Left fit (%.2f°)', a_left));
        plot(x_plot_R, rows_fit, 'm-', 'LineWidth', 2, 'DisplayName', ...
            sprintf('Right fit (%.2f°)', a_right));
        if ~fit_lines_only
            xline(center_col,  'c--', 'LineWidth', 1,   'DisplayName', 'Center');
            yline(row_min,     'w--', 'LineWidth', 1,   'DisplayName', 'Flame top');
            yline(row_max,     'w:',  'LineWidth', 1,   'DisplayName', 'Flame bottom');
            yline(rows_b(i0),  'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
            yline(rows_b(i1),  'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
        axis image;
        title({short_name, sprintf('Angle: %.1f°', alpha_deg)}, 'FontSize', 10);
        legend('Location', 'best');
        hold off;

    end
end


% =========================================================
function [idx, x_drop] = edge_from_gradient(raw_outward, smooth_outward, threshold, drop_frac)
% EDGE_FROM_GRADIENT  Detect flame edge in a 1D outward-scanning profile.
%
%   1. Gradient of smoothed profile → peak_i (steepest rising point)
%   2. After peak_i: x_drop where gradient <= drop_frac * gradient_peak
%      (robust to OH plateaus and double-peak profiles)
%   3. Re-normalise raw data in window around x_drop locally
%   4. First crossing of threshold scanning outward → flame edge index
%
%   Returns index into raw_outward (NaN if not found), and x_drop index.

    idx    = NaN;
    x_drop = NaN;
    grad   = gradient(double(smooth_outward));

    g_global_max = max(grad);
    if g_global_max <= 0, return; end
    [~, locs] = findpeaks(grad, 'MinPeakHeight', 0.10 * g_global_max);
    if isempty(locs), return; end
    peak_i = locs(1);
    g_max  = grad(peak_i);

    % locate x_drop: first point after peak where gradient falls to drop_frac
    post  = grad(peak_i:end);
    d_rel = find(post <= drop_frac * g_max, 1, 'first');
    if isempty(d_rel)
        x_drop = peak_i;
    else
        x_drop = peak_i + d_rel - 2;
        x_drop = max(peak_i, x_drop);
    end

    % window of raw data around x_drop for local re-normalisation
    win_start = max(1, x_drop - 20);
    win_end   = min(length(raw_outward), x_drop + 10);
    seg       = raw_outward(win_start:win_end);

    s_min = min(seg);
    s_max = max(seg);
    if s_max <= s_min, return; end

    seg_norm = (seg - s_min) / (s_max - s_min);

    % first crossing of threshold scanning outward (win_start → win_end)
    for c = 1 : length(seg_norm)
        if seg_norm(c) >= threshold
            idx = win_start + c - 1;
            return;
        end
    end
end
