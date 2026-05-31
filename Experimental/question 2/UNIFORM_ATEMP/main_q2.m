clc ; clear; close all

% CH4 calibration
[px_mm_x_CH4, px_mm_y_CH4] = calibratePixelPerMM( ...
    'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\calibration\Calibration\CH4\B16\B0001.b16', ...
    false);

% H2 calibration
[px_mm_x_H2, px_mm_y_H2] = calibratePixelPerMM( ...
    'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\calibration\Calibration\H2\B16\B0001.b16', ...
    false);

[H2_image_locs, CH4_image_loc] = OHPLIF_locations();

[data_H2, data_CH4] = Get_experimental_data();

% --- Detection parameters (shared by both fuels) ---
threshold  = 0.25;   % local normalised intensity threshold
drop_frac  = 0.90;   % gradient drops to this fraction of its peak → working point
trunc_frac = 0.10;   % fraction trimmed from each end of below-tip rows for fitting

% Old manual H2 tuning arrays kept here for reference (no longer used):
%   baselines_H2  = [400,500,625,720,750,780,800,820,820];
%   centers_H2    = [477,477,477,477,477,477,477,477,477];
%   tops_H2       = [925,925,925,925,915,915,915,915,915];
%   thersholds_H2 = [1100,1100,1100,1100,1100,600,360,300,300];

results_H2 = flamefront_sf_H2(H2_image_locs, ...
    1/px_mm_x_H2, 1/px_mm_y_H2, threshold, drop_frac, trunc_frac, ...
    data_H2(:,2)'*1000, ...
    true);

results_CH4 = flamefront_sf_CH4(CH4_image_loc, ...
    1/px_mm_x_CH4, 1/px_mm_y_CH4, threshold, drop_frac, trunc_frac, ...
    data_CH4(:,2)'*1000, ...
    true);

plotSfAngle(data_CH4(:,1)', data_H2(:,1)', results_H2, results_CH4)
