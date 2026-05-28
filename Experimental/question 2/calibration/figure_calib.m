%% Calibration: Pixel spacing from millimetric paper
clear; clc; close all;

%% 1. Read the image
dirPath = 'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\calibration\Calibration\CH4\B16'
%dirPath = 'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\calibration\Calibration\H2\B16'

filename = 'B0001.b16';
image = readB16(dirPath, filename);

%% 2. Select a clean region (avoid bright object at bottom)
rowRange = 200:1100;
colRange = 200:900;

%% ============================================================
%% METHOD 1: Peak detection
%% ============================================================

%% 3a. X-direction: detect vertical lines
profile_x = mean(image(rowRange,colRange ), 1);
profile_x_inv = -profile_x;
profile_x_inv = profile_x_inv - min(profile_x_inv);

[pks_x, locs_x] = findpeaks(profile_x_inv, 'MinPeakDistance', 15, 'MinPeakProminence', 1);

spacings_x = diff(locs_x);
mean_spacing_x = mean(spacings_x);

%% 3b. Y-direction: detect horizontal lines
profile_y = mean(image(rowRange, colRange), 2);
profile_y_inv = -profile_y;
profile_y_inv = profile_y_inv - min(profile_y_inv);

[pks_y, locs_y] = findpeaks(profile_y_inv, 'MinPeakDistance', 15, 'MinPeakProminence', 1);

spacings_y = diff(locs_y);
mean_spacing_y = mean(spacings_y);

%% 4. Print peak detection results
fprintf('========== PEAK DETECTION ==========\n');
fprintf('X-direction: %d lines, spacings: ', length(locs_x));
fprintf('%.0f ', spacings_x);
fprintf('\nX mean spacing: %.2f pixels\n\n', mean_spacing_x);

fprintf('Y-direction: %d lines, spacings: ', length(locs_y));
fprintf('%.0f ', spacings_y);
fprintf('\nY mean spacing: %.2f pixels\n\n', mean_spacing_y);

%% 5. Plot profiles with detected peaks
fig1 = figure;
set(fig1, 'Position', [100 100 1400 500]);

subplot(1,2,1);
plot(profile_x, 'b'); hold on;
plot(locs_x, profile_x(locs_x), 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('X (pixels)'); ylabel('Mean intensity');
title(sprintf('X-profile: mean spacing = %.1f px', mean_spacing_x));
legend('Profile', 'Detected lines');
grid on;

subplot(1,2,2);
plot(profile_y, 'b'); hold on;
plot(locs_y, profile_y(locs_y), 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Y (pixels)'); ylabel('Mean intensity');
title(sprintf('Y-profile: mean spacing = %.1f px', mean_spacing_y));
legend('Profile', 'Detected lines');
grid on;

%% 6. Overlay detected lines on image
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

% %% ============================================================
% %% METHOD 2: Autocorrelation (more robust)
% %% ============================================================
% 
% %% 7a. X-direction autocorrelation
% profile_x_clean = mean(image(rowRange, colRange), 1);
% profile_x_clean = profile_x_clean - mean(profile_x_clean);
% 
% [acor_x, lags_x_ac] = xcorr(profile_x_clean, 'coeff');
% 
% half_x = acor_x(lags_x_ac > 10 & lags_x_ac < 200);
% lags_x_pos = lags_x_ac(lags_x_ac > 10 & lags_x_ac < 200);
% 
% [~, loc_ac_x] = findpeaks(half_x, 'MinPeakDistance', 15, 'NPeaks', 1);
% spacing_ac_x = lags_x_pos(loc_ac_x);
% 
% %% 7b. Y-direction autocorrelation
% profile_y_clean = mean(image(rowRange, colRange), 2)';
% profile_y_clean = profile_y_clean - mean(profile_y_clean);
% 
% [acor_y, lags_y_ac] = xcorr(profile_y_clean, 'coeff');
% 
% half_y = acor_y(lags_y_ac > 10 & lags_y_ac < 200);
% lags_y_pos = lags_y_ac(lags_y_ac > 10 & lags_y_ac < 200);
% 
% [~, loc_ac_y] = findpeaks(half_y, 'MinPeakDistance', 15, 'NPeaks', 1);
% spacing_ac_y = lags_y_pos(loc_ac_y);
% 
% %% 8. Print autocorrelation results
% fprintf('========== AUTOCORRELATION ==========\n');
% fprintf('X-direction: 1mm = %d pixels\n', spacing_ac_x);
% fprintf('Y-direction: 1mm = %d pixels\n', spacing_ac_y);
% 
% %% 9. Plot autocorrelation
% fig3 = figure;
% set(fig3, 'Position', [100 100 1400 500]);
% 
% subplot(1,2,1);
% plot(lags_x_pos, half_x, 'b'); hold on;
% xline(spacing_ac_x, 'r--', 'LineWidth', 2);
% xlabel('Lag (pixels)'); ylabel('Autocorrelation');
% title(sprintf('X autocorrelation: 1mm = %d px', spacing_ac_x));
% grid on;
% 
% subplot(1,2,2);
% plot(lags_y_pos, half_y, 'b'); hold on;
% xline(spacing_ac_y, 'r--', 'LineWidth', 2);
% xlabel('Lag (pixels)'); ylabel('Autocorrelation');
% title(sprintf('Y autocorrelation: 1mm = %d px', spacing_ac_y));
% grid on;