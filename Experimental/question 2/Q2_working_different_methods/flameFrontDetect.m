function stats = flameFrontDetect(filepath, px_per_mm_x, px_per_mm_y, threshold, drop_frac, trunc_frac, verbose)
% FLAMEFRONTDETECT  Unified gradient-guided flame front detection for CH4 and H2.
%
%   Handles flames with double OH-brightness peaks (rich H2) by using the
%   smoothed-intensity gradient rather than the raw intensity peak:
%
%   Per row, scanning outward from the flame axis:
%     1. Gaussian-smoothed gradient → locate steepest rising point (peak)
%     2. After the peak, find where gradient drops to drop_frac of its max
%        → this is the working point (x_drop), robust to OH plateaus / double peaks
%     3. In a window around x_drop, re-normalise the RAW unfiltered intensity
%        locally and find the first crossing of threshold → precise flame edge
%
%   Inputs
%   ------
%   filepath     - path to .b16 image file
%   px_per_mm_x  - pixels per mm in x  (default 29.706)
%   px_per_mm_y  - pixels per mm in y  (default 29.695)
%   threshold    - local normalised threshold 0-1  (default 0.25)
%   drop_frac    - gradient drop fraction used to locate working point (default 0.80)
%   trunc_frac   - fraction trimmed from each end of below-tip rows (default 0.10)
%   verbose      - if true: plot intensity, gradient, and angle figures (default true)
%
%   Output struct fields  (identical to flameFrontStats)
%   ----------------------------------------------------
%   alpha_half_deg, alpha_deg, Dc_mm
%   tip_row, tip_col
%   rows, left_edge, right_edge
%   c_left, c_right   ([slope intercept] pixel-space fits)

    if nargin < 2 || isempty(px_per_mm_x), px_per_mm_x = 29.706; end
    if nargin < 3 || isempty(px_per_mm_y), px_per_mm_y = 29.695; end
    if nargin < 4 || isempty(threshold),   threshold   = 0.25;   end
    if nargin < 5 || isempty(drop_frac),   drop_frac   = 0.80;   end
    if nargin < 6 || isempty(trunc_frac),  trunc_frac  = 0.10;   end
    if nargin < 7 || isempty(verbose),     verbose     = true;   end

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

    % --- per-row detection ---
    left_all  = NaN(nRows, 1);
    right_all = NaN(nRows, 1);

    for r = 51 : nRows - 50
        if max(img_s(r, :)) < 0.10, continue; end

        % right side: scan from RIGHT EDGE inward toward center (outside-in)
        % This ensures the outer flame front is detected first, even when
        % a brighter inner OH zone exists (rich H2 double-peak case).
        right_half    = img_s(r, center_col:end);
        right_half_sm = img_smooth(r, center_col:end);
        L_right       = length(right_half);
        k = edge_from_gradient(fliplr(right_half), fliplr(right_half_sm), threshold, drop_frac);
        if ~isnan(k)
            right_all(r) = center_col + (L_right - k);
        end

        % left side: scan from LEFT EDGE inward toward center (outside-in)
        left_half    = img_s(r, 1:center_col);
        left_half_sm = img_smooth(r, 1:center_col);
        k = edge_from_gradient(left_half, left_half_sm, threshold, drop_frac);
        if ~isnan(k)
            left_all(r) = k;
        end
    end

    % keep only rows where BOTH fronts detected and right > left
    valid      = ~isnan(left_all) & ~isnan(right_all) & (right_all > left_all);
    rows       = find(valid);
    left_edge  = left_all(rows);
    right_edge = right_all(rows);
    widths_px  = right_edge - left_edge;

    if isempty(rows)
        error('flameFrontDetect: no valid rows found in %s', filepath);
    end

    % --- tip: row with narrowest flame width ---
    [~, tip_idx] = min(widths_px);
    tip_row = rows(tip_idx);
    tip_col = round((left_edge(tip_idx) + right_edge(tip_idx)) / 2);

    % --- fit only rows BELOW the tip ---
    below    = rows > tip_row;
    rows_b   = rows(below);
    left_b   = left_edge(below);
    right_b  = right_edge(below);
    widths_b = widths_px(below);

    nb = numel(rows_b);
    i0 = max(1, floor(trunc_frac * nb) + 1);
    i1 = min(nb, floor((1 - trunc_frac) * nb));

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

    % =========================================================
    % --- verbose plots ---
    % =========================================================
    if verbose
        parts      = strsplit(filepath, {'\', '/'});
        short_name = parts{end-2};

        % highlight rows evenly spaced in the fitted region
        n_highlight    = 5;
        highlight_rows = round(linspace(rows_b(i0), rows_b(i1), n_highlight));
        cmap           = lines(n_highlight);

        % Figure 1: intensity profiles (smoothed for display, raw dot at edge)
        h_intens = figure('Name', ['Intensity profiles - ', short_name]);
        hold on;
        xlabel('Pixel index from center');
        ylabel('Intensity (normalised, unfiltered)');
        title({'Right-half intensity profile for each row', short_name}, 'FontSize', 8);
        yline(threshold, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');

        % Figure 2: gradient profiles
        h_grad = figure('Name', ['Intensity gradient - ', short_name]);
        hold on;
        xlabel('Pixel index from center');
        ylabel('d(Intensity)/dx  (smoothed)');
        title({'Right-half intensity gradient for each row', short_name}, 'FontSize', 8);
        yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        for r = 51 : nRows - 50
            if max(img_s(r, :)) < 0.05, continue; end
            right_raw  = img_s(r, center_col:end);
            right_sm   = img_smooth(r, center_col:end);
            grad_right = gradient(double(right_sm));

            hi = find(highlight_rows == r, 1);

            figure(h_intens);
            if ~isempty(hi)
                plot(1:length(right_raw), right_raw, 'Color', cmap(hi,:), ...
                    'LineWidth', 1.2, 'DisplayName', sprintf('Row %d', r));
            else
                plot(1:length(right_raw), right_raw, 'Color', [0 0 1 0.1], ...
                    'HandleVisibility', 'off');
            end
            if ~isnan(right_all(r))
                idx = right_all(r) - center_col + 1;
                if idx >= 1 && idx <= length(right_raw)
                    plot(idx, right_raw(idx), 'ro', 'MarkerSize', 3, 'HandleVisibility', 'off');
                end
            end

            figure(h_grad);
            if ~isempty(hi)
                plot(1:length(grad_right), grad_right, 'Color', cmap(hi,:), ...
                    'LineWidth', 1.2, 'DisplayName', sprintf('Row %d', r));
            else
                plot(1:length(grad_right), grad_right, 'Color', [0 0 1 0.1], ...
                    'HandleVisibility', 'off');
            end
            if ~isnan(right_all(r))
                idx = right_all(r) - center_col + 1;
                if idx >= 1 && idx <= length(grad_right)
                    plot(idx, grad_right(idx), 'ro', 'MarkerSize', 3, 'HandleVisibility', 'off');
                end
            end
        end
        figure(h_intens); hold off;
        figure(h_grad);   hold off;

        % Figure 3: flame front overlay
        rows_fit = (rows_b(i0):rows_b(i1))';
        x_plot_L = polyval(c_left,  rows_fit);
        x_plot_R = polyval(c_right, rows_fit);

        figure('Name', ['Flame Front Angle - ', short_name]);
        imagesc(image); colormap(hot); colorbar; hold on;
        plot(left_edge,  rows, 'r',  'LineWidth', 1,   'DisplayName', 'Left front');
        plot(right_edge, rows, 'b',  'LineWidth', 1,   'DisplayName', 'Right front');
        plot(x_plot_L, rows_fit, 'g-', 'LineWidth', 2, 'DisplayName', ...
            sprintf('Left fit (%.2f°)', a_left));
        plot(x_plot_R, rows_fit, 'm-', 'LineWidth', 2, 'DisplayName', ...
            sprintf('Right fit (%.2f°)', a_right));
        xline(center_col,  'c--', 'LineWidth', 1,   'DisplayName', 'Center');
        yline(rows_b(i0),  'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        yline(rows_b(i1),  'g--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        axis image;
        title({short_name, sprintf('Angle: %.1f°', alpha_deg)}, 'FontSize', 10);
        legend('Location', 'best');
        hold off;

        fprintf('--- flameFrontDetect (gradient-guided) ---\n');
        fprintf('File:              %s\n', filepath);
        fprintf('Gradient drop:     %.0f%% of gradient peak\n', drop_frac*100);
        fprintf('Local threshold:   %.2f (normalised)\n', threshold);
        fprintf('Truncation:        %.0f%% - %.0f%%\n', trunc_frac*100, (1-trunc_frac)*100);
        fprintf('Left  fit:  x = %.4f * row + %.4f  (%.2f°)\n', c_left(1),  c_left(2),  a_left);
        fprintf('Right fit:  x = %.4f * row + %.4f  (%.2f°)\n', c_right(1), c_right(2), a_right);
        fprintf('Opening angle:     %.2f°\n', alpha_deg);
        fprintf('Rows used in fit:  %d\n', length(rows_fit));
    end
end


% =========================================================
function idx = edge_from_gradient(raw_outward, smooth_outward, threshold, drop_frac)
% EDGE_FROM_GRADIENT  Detect flame edge in a 1D outward-scanning profile.
%
%   1. Gradient of smoothed profile → peak_i (steepest rising point)
%   2. After peak_i: find x_drop where gradient <= drop_frac * gradient_peak
%      (robust to OH plateaus and double-peak profiles)
%   3. Re-normalise raw data in window [x_drop-20, x_drop+10] locally
%   4. First crossing of threshold scanning outward → flame edge index
%
%   Returns index into raw_outward (NaN if not found).

    idx  = NaN;
    grad = gradient(double(smooth_outward));

    [g_max, peak_i] = max(grad);
    if g_max <= 0, return; end   % no rising edge in this row

    % locate x_drop: first point after peak where gradient falls to drop_frac
    post  = grad(peak_i:end);
    d_rel = find(post <= drop_frac * g_max, 1, 'first');
    if isempty(d_rel)
        x_drop = peak_i;           % gradient never drops — use peak position
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
