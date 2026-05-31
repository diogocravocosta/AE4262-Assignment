function display_batch_images(images, label, phi_values)
% DISPLAY_BATCH_IMAGES  Display a batch of images with a shared hot color scale.
%
%   display_batch_images(images, label, phi_values)
%
%   INPUT
%       images     - N×1 cell array of 2-D double arrays
%       label      - string prefix shown in every figure title (e.g. 'OHstar Raw')
%       phi_values - (optional) 1×N array of phi values for titles

    if nargin < 3
        phi_values = [];
    end

    % Shared color limits across all images in this batch
    gmin = inf;  gmax = -inf;
    for i = 1:numel(images)
        gmin = min(gmin, min(images{i}(:)));
        gmax = max(gmax, max(images{i}(:)));
    end

    for i = 1:numel(images)
        if ~isempty(phi_values)
            fig_title = sprintf('%s  [\\phi = %.1f]', label, phi_values(i));
        else
            fig_title = sprintf('%s  [%d]', label, i);
        end

        figure('Color', 'k');
        imagesc(images{i}, [gmin gmax]);
        axis image; axis off;
        colormap(hot);
        cb = colorbar;
        cb.Color = 'w';
        title(fig_title, 'Color', 'w', 'Interpreter', 'tex');
        set(gca, 'Color', 'k');
    end
end
