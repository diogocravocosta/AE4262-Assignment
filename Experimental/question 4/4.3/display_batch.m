function display_batch(paths)
% DISPLAY_BATCH  Display raw B16 images with hot colormap.
%
%   display_batch(paths)
%
%   INPUT
%       paths - N×1 cell array of .b16 file paths (from get_paths)

    for i = 1:numel(paths)
        img = double(readB16(paths{i}));

        figure('Name', paths{i}, 'Color', 'k');
        imagesc(img);
        axis image; axis off;
        colormap(hot);
        cb = colorbar;
        cb.Color = 'w';
        ylabel(cb, 'Intensity [counts]', 'Color', 'w');
        title(paths{i}, 'Color', 'w', 'Interpreter', 'none');
        set(gca, 'Color', 'k');
    end
end
