function [flame_left, flame_center, flame_right] = flameFront(filename, center_col, limit, row_min, row_max, verbose)
if nargin < 6
    verbose = true;
end

parts = strsplit(filename, {'\', '/'});
short_name = parts{end-2};

image = readB16(filename);
image_d = double(image);
bg = mean([image_d(1:50,1:50), image_d(1:50,end-50:end)], 'all');
thresh = limit;
nRows = size(image, 1);
flame_left = NaN(nRows, 1);
flame_right = NaN(nRows, 1);
if verbose
    figure('Name', ['Intensity profiles - ', short_name]);
    hold on;
    xlabel('Pixel index from center');
    ylabel('Intensity');
    title({'Right-half intensity profile for each row', short_name}, 'FontSize', 8);
    yline(thresh, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');
end
for row = row_min:row_max
    line = image_d(row,:);
    % Right side
    right_half = line(center_col:end);
    idx = find(right_half > thresh, 1, 'first');
    if ~isempty(idx)
        flame_right(row) = center_col + idx - 1;
    end
    % Plot the right-half profile for this row
    if verbose
        plot(1:length(right_half), right_half, 'Color', [0 0 1 0.1], 'HandleVisibility', 'off');
        if ~isempty(idx)
            plot(idx, right_half(idx), 'ro', 'MarkerSize', 3, 'HandleVisibility', 'off');
        end
    end
    % Left side
    left_half = fliplr(line(1:center_col));
    idx = find(left_half > thresh, 1, 'first');
    if ~isempty(idx)
        flame_left(row) = center_col - idx + 1;
    end
end
% Mean center from all detected front points
flame_center = mean([flame_left; flame_right], 'omitnan');
if verbose
    legend('Threshold');
    hold off;
    fprintf('Flame front detection for: %s\n', filename);
    fprintf('Center column: %d\n', center_col);
    fprintf('Row range: %d to %d\n', row_min, row_max);
    fprintf('Estimated background: %.1f\n', bg);
    fprintf('Threshold: %.1f\n', thresh);
    fprintf('Image size: %d x %d\n', size(image,1), size(image,2));
    fprintf('Left front detected in %d / %d rows\n', sum(~isnan(flame_left)), row_max - row_min + 1);
    fprintf('Right front detected in %d / %d rows\n', sum(~isnan(flame_right)), row_max - row_min + 1);
    fprintf('Mean flame center: %.1f px\n', flame_center);
    figure('Name', ['Flame front - ', short_name]);
    imagesc(image); colormap(gray); colorbar; hold on;
    plot(flame_left, 1:nRows, 'r', 'LineWidth', 1.5);
    plot(flame_right, 1:nRows, 'r', 'LineWidth', 1.5);
    xline(center_col, 'b--', 'LineWidth', 1);
    xline(flame_center, 'm--', 'LineWidth', 1.5);
    yline(row_min, 'g--', 'LineWidth', 1);
    yline(row_max, 'g--', 'LineWidth', 1);
    axis image;
    title({'Flame front overlaid on image', short_name}, 'FontSize', 8);
    legend('Left front', 'Right front', 'Center col', 'Mean center', 'Row range');
end
end