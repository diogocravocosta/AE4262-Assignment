clc ; clear; close all

% CH4 calibration
[px_mm_x_CH4, px_mm_y_CH4] = calibratePixelPerMM( ...
    'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\calibration\Calibration\CH4\B16\B0001.b16', ...
    false);

% H2 calibration
[px_mm_x_H2, px_mm_y_H2] = calibratePixelPerMM( ...
    'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\calibration\Calibration\H2\B16\B0001.b16' ...
    ,false);


[H2_image_locs, CH4_image_loc] = OHPLIF_locations();

[data_H2,data_CH4] = Get_experimental_data();

%sizes_H2 = showB16(H2_image_locs);

%sizes_ch4 = showB16(CH4_image_loc);

%flameFront(H2_image_locs{1}, 480, 500,200,900, true)



n = length(H2_image_locs);
baselines_H2 = [400,500,625,720,750,780,800,820,820];
centers_H2 = [477,477,477,477,477,477,477,477,477];
tops_H2 = [925,925,925,925,915,915,915,915,915]; % Define tops for H2 images
thersholds_H2 = [1100,1100,1100,1100,1100,600,360,300,300];

%flamefront_batch(H2_image_locs,baselines_H2,centers_H2,tops_H2,thersholds_H2,true)


n = length(CH4_image_loc);
baselines_CH4 = [580,820,670];
centers_CH4 = [500,500,500];
tops_CH4 = [1320,1320,1315]; % Define tops for images
thersholds_CH4 = [350,550,390];

%flamefront_batch(CH4_image_loc,baselines_CH4,centers_CH4,tops_CH4,thersholds_CH4,true)



results_H2 = flamefront_sf_H2(H2_image_locs, baselines_H2, centers_H2, tops_H2, thersholds_H2, ...
    1/px_mm_x_H2, 1/px_mm_y_H2, ...
     data_H2(:,2)'*1000, ...
    true);


results_CH4 = flamefront_sf_CH4(CH4_image_loc, ...
    1/px_mm_x_CH4, 1/px_mm_y_CH4, 0.25, 0.1, ...
    data_CH4(:,2)'*1000, ...
    true);


plotSfAngle(data_CH4(:,1)', data_H2(:,1)', results_H2, results_CH4)



