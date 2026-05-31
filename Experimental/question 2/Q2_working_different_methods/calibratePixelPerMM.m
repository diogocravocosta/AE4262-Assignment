function [px_per_mm_x, px_per_mm_y] = calibratePixelPerMM(Path, verbose)
% CALIBRATEPIXELPERMM  Compute pixel-per-mm scale from millimetric paper image.
%   Uses peak detection on row/column-averaged intensity profiles to find
%   the 1 mm grid line spacing.
%
%   [px_per_mm_x, px_per_mm_y] = calibratePixelPerMM(dirPath, filename, verbose)
%
%   Inputs:
%       dirPath  – folder containing the .b16 calibration file
%       filename – name of the .b16 file (e.g. 'B0001.b16')
%       verbose  – (optional, default false) if true, print results and show plots
%
%   Outputs:
%       px_per_mm_x – pixels per millimetre in the X (horizontal) direction
%       px_per_mm_y – pixels per millimetre in the Y (vertical) direction

    if nargin < 2
        verbose = false;
    end

    %% 1. Read the image
    image = readB16(Path);

    %% 2. Select a clean region (avoid bright object at bottom)
    rowRange = 150:1100;
    colRange = 150:900;

    %% 3a. X-direction: detect vertical lines
    profile_x = mean(image(rowRange, colRange), 1);
    profile_x_inv = -profile_x;
    profile_x_inv = profile_x_inv - min(profile_x_inv);

    [~, locs_x] = findpeaks(profile_x_inv, 'MinPeakDistance', 15, 'MinPeakProminence', 1);
    
    locs_x =locs_x(2:end-1);
    spacings_x = diff(locs_x);
    px_per_mm_x = mean(spacings_x);

    %% 3b. Y-direction: detect horizontal lines
    profile_y = mean(image(rowRange, colRange), 2);
    profile_y_inv = -profile_y;
    profile_y_inv = profile_y_inv - min(profile_y_inv);

    [~, locs_y] = findpeaks(profile_y_inv, 'MinPeakDistance', 15, 'MinPeakProminence', 1);

    locs_y =locs_y(2:end-1);
    spacings_y = diff(locs_y);
    px_per_mm_y = mean(spacings_y);

    %% 4. Print and plot if verbose
    if verbose
        fprintf('========== PEAK DETECTION ==========\n');
        fprintf('X-direction: %d lines, spacings: ', length(locs_x));
        fprintf('%.0f ', spacings_x);
        fprintf('\nX mean spacing: %.2f px/mm\n\n', px_per_mm_x);

        fprintf('Y-direction: %d lines, spacings: ', length(locs_y));
        fprintf('%.0f ', spacings_y);
        fprintf('\nY mean spacing: %.2f px/mm\n\n', px_per_mm_y);

        % Profiles with detected peaks
        fig1 = figure;
        set(fig1, 'Position', [100 100 1400 500]);

        subplot(1,2,1);
        plot(profile_x, 'b'); hold on;
        plot(locs_x, profile_x(locs_x), 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        xlabel('X (pixels)'); ylabel('Mean intensity');
        title(sprintf('X-profile: mean spacing = %.1f px', px_per_mm_x));
        legend('Profile', 'Detected lines');
        grid on;

        subplot(1,2,2);
        plot(profile_y, 'b'); hold on;
        plot(locs_y, profile_y(locs_y), 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        xlabel('Y (pixels)'); ylabel('Mean intensity');
        title(sprintf('Y-profile: mean spacing = %.1f px', px_per_mm_y));
        legend('Profile', 'Detected lines');
        grid on;

        % Detected lines overlaid on image
        fig2 = figure;
        imagesc(image); colormap(gray); colorbar;
        hold on;
        for i = 1:length(locs_x)
            xline(locs_x(i) + colRange(1) - 1, 'r', 'LineWidth', 0.5);
        end
        for i = 1:length(locs_y)
            yline(locs_y(i) + rowRange(1) - 1, 'g', 'LineWidth', 0.5);
        end
        title('Detected grid lines overlaid on image');
    end
end