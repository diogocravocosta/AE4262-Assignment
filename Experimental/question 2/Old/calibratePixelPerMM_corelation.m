function [px_per_mm_x, px_per_mm_y] = calibratePixelPerMM(dirPath, filename, verbose)
% CALIBRATEPIXELPERMM  Determine pixel-per-mm from millimetric paper image.
%
%   [px_per_mm_x, px_per_mm_y] = calibratePixelPerMM(dirPath, filename, verbose)
%
%   Inputs:
%       dirPath  - folder containing the .b16 file
%       filename - name of the .b16 file (e.g. 'B0001.b16')
%       verbose  - logical; if true, prints results and shows plots
%
%   Outputs:
%       px_per_mm_x - pixels per mm in the X-direction
%       px_per_mm_y - pixels per mm in the Y-direction

if nargin < 3
    verbose = false;
end

%% 1. Read the image
image = readB16(dirPath, filename);

%% 2. Select a clean region
rowRange = 100:1100;
colRange = 100:900;

%% 3a. X-direction autocorrelation
profile_x = mean(image(rowRange, colRange), 1);
profile_x = profile_x - mean(profile_x);

[acor_x, lags_x] = xcorr(profile_x, 'coeff');

half_x = acor_x(lags_x > 10 & lags_x < 200);
lags_x_pos = lags_x(lags_x > 10 & lags_x < 200);

[~, loc_x] = findpeaks(half_x, 'MinPeakDistance', 15, 'NPeaks', 1);
px_per_mm_x = lags_x_pos(loc_x);

%% 3b. Y-direction autocorrelation
profile_y = mean(image(rowRange, colRange), 2)';
profile_y = profile_y - mean(profile_y);

[acor_y, lags_y] = xcorr(profile_y, 'coeff');

half_y = acor_y(lags_y > 10 & lags_y < 200);
lags_y_pos = lags_y(lags_y > 10 & lags_y < 200);

[~, loc_y] = findpeaks(half_y, 'MinPeakDistance', 15, 'NPeaks', 1);
px_per_mm_y = lags_y_pos(loc_y);

%% 4. Verbose output
if verbose
    fprintf('========== CALIBRATION (Autocorrelation) ==========\n');
    fprintf('X-direction: 1 mm = %d pixels\n', px_per_mm_x);
    fprintf('Y-direction: 1 mm = %d pixels\n', px_per_mm_y);

    fig1 = figure;
    set(fig1, 'Position', [100 100 1400 500]);

    subplot(1,2,1);
    plot(lags_x_pos, half_x, 'b'); hold on;
    xline(px_per_mm_x, 'r--', 'LineWidth', 2);
    xlabel('Lag (pixels)'); ylabel('Autocorrelation');
    title(sprintf('X autocorrelation: 1 mm = %d px', px_per_mm_x));
    grid on;

    subplot(1,2,2);
    plot(lags_y_pos, half_y, 'b'); hold on;
    xline(px_per_mm_y, 'r--', 'LineWidth', 2);
    xlabel('Lag (pixels)'); ylabel('Autocorrelation');
    title(sprintf('Y autocorrelation: 1 mm = %d px', px_per_mm_y));
    grid on;

    % Overlay on image
    fig2 = figure;
    imagesc(image); colormap(gray); colorbar;
    title(sprintf('Calibration image — %.0f x %.0f px/mm', px_per_mm_x, px_per_mm_y));
end

end