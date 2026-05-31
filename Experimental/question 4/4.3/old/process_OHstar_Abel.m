function process_OHstar_Abel(filepath)
% PROCESS_OHSTAR_ABEL  Read an OH* chemiluminescence .b16 file, display it,
%   apply Abel inversion, and show the reconstructed radial emissivity
%   with a 'hot' colormap.
%
%   process_OHstar_Abel(filepath)
%
%   INPUT
%       filepath  -  full path to an OH* .b16 image file
%
%   REQUIRES (in the same folder / on the MATLAB path)
%       readB16.m                   -  provided with the data package
%       abel.m  (Carsten Killer)    -  MATLAB File Exchange #43639
%
%   AE4262 Combustion Lab - TU Delft

    %% ---- 1. Read the B16 file ----
    img = double(readB16(filepath));

    %% ---- 2. Display the raw OH* image (showB16 style) ----
    figure('Name','Raw OH* Chemiluminescence','Color','k');
    imagesc(img);
    axis image; axis off;
    colormap(hot);
    cb = colorbar;
    cb.Color = 'w';
    ylabel(cb, 'Intensity [counts]', 'Color', 'w');
    title('OH* Chemiluminescence (line-of-sight)', 'Color', 'w');
    set(gca, 'Color', 'k');

    %% ---- 3. Abel inversion ----
    [nRows, nCols] = size(img);

    % Find the axis of symmetry (column with max total intensity)
    [~, cIdx] = max(sum(img, 1));

    % Build symmetric half-image by averaging left & right halves
    halfW  = min(cIdx - 1, nCols - cIdx);
    leftH  = img(:, cIdx - (1:halfW));      % flipped left half
    rightH = img(:, cIdx + (1:halfW));
    halfImg = 0.5 * (leftH + rightH);

    % Apply Abel inversion row by row using abel() from Carsten Killer
    abelHalf = zeros(nRows, halfW);
    for row = 1:nRows
        profile = halfImg(row, :);
        abelHalf(row, :) = abel_inversion(profile);
    end

    % Clamp negative values (noise artefacts)
    abelHalf(abelHalf < 0) = 0;

    % Mirror to reconstruct full image
    abel_img = [fliplr(abelHalf), abelHalf];

    %% ---- 4. Display the Abel-inverted image ----
    figure('Name','Abel-Inverted OH* Emissivity','Color','k');
    imagesc(abel_img);
    axis image; axis off;
    colormap(hot);
    cb2 = colorbar;
    cb2.Color = 'w';
    ylabel(cb2, 'Local emissivity [a.u.]', 'Color', 'w');
    title('Abel-Inverted OH* (radial emissivity)', 'Color', 'w');
    set(gca, 'Color', 'k');

end