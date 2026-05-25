%% Calibration: Pixel spacing from millimetric paper
% This script reads the B0001.b16 image and calculates the pixel spacing
% between grid lines in both x and y directions.

clear; clc; close all;

%% 1. Read the image
dirPath = 'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\calibration\Calibration\CH4\B16\';
filename = 'B0001.b16';
image = readB16(dirPath, filename);

%% 2. Select a clean region (avoid the bright object at bottom)
rowRange = 100:1000;   % rows to use (avoid bottom bright region)
colRange = 100:1000;   % columns to use

%% 3. X-direction: detect vertical lines
% Average intensity along rows to get a 1D profile vs x
profile_x = mean(image(rowRange, :), 1);

% The grid lines appear darker, so invert the profile to find peaks
profile_x_inv = -profile_x;
profile_x_inv = profile_x_inv - min(profile_x_inv); % shift to positive

% Find peaks (tune MinPeakDistance and MinPeakProminence as needed)
[pks_x, locs_x] = findpeaks(profile_x_inv, 'MinPeakDistance', 30, 'MinPeakProminence', 5);

% Pixel spacings between consecutive vertical lines
spacings_x = diff(locs_x);
mean_spacing_x = mean(spacings_x);

%% 4. Y-direction: detect horizontal lines
% Average intensity along columns to get a 1D profile vs y
profile_y = mean(image(:, colRange), 2);

% Invert
profile_y_inv = -profile_y;
profile_y_inv = profile_y_inv - min(profile_y_inv);

% Find peaks
[pks_y, locs_y] = findpeaks(profile_y_inv, 'MinPeakDistance', 30, 'MinPeakProminence', 5);

% Pixel spacings between consecutive horizontal lines
spacings_y = diff(locs_y);
mean_spacing_y = mean(spacings_y);

%% 5. Display results
fprintf('=== X-direction (vertical lines) ===\n');
fprintf('Number of lines detected: %d\n', length(locs_x));
fprintf('Individual spacings (pixels): ');
fprintf('%.1f  ', spacings_x);
fprintf('\nMean spacing: %.2f pixels\n\n', mean_spacing_x);

fprintf('=== Y-direction (horizontal lines) ===\n');
fprintf('Number of lines detected: %d\n', length(locs_y));
fprintf('Individual spacings (pixels): ');
fprintf('%.1f  ', spacings_y);
fprintf('\nMean spacing: %.2f pixels\n\n', mean_spacing_y);

%% 6. Plot profiles with detected peaks
figure('Position', [100 100 1400 500]);

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

%% 7. Overlay detected lines on image
figure;
imagesc(image); colormap(gray); colorbar;
hold on;
for i = 1:length(locs_x)
    xline(locs_x(i), 'r', 'LineWidth', 0.5);
end
for i = 1:length(locs_y)
    yline(locs_y(i), 'g', 'LineWidth', 0.5);
end
title('Detected grid lines overlaid on image');
legend('Vertical lines (red)', 'Horizontal lines (green)');