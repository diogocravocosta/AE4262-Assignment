function display_phi_trio(ohstar_raw, recon_all, plif_raw, phi_values)
% DISPLAY_PHI_TRIO  One figure per phi value with 3 subplots:
%   OHstar Raw | OHstar After Inversion | PLIF
%   All images of the same type share a common color scale.

    % Shared color limits per image type
    ohstar_clim = [min(cellfun(@(x) min(x(:)), ohstar_raw)), ...
                   max(cellfun(@(x) max(x(:)), ohstar_raw))];
    recon_clim  = [min(cellfun(@(x) min(x(:)), recon_all)),  ...
                   max(cellfun(@(x) max(x(:)), recon_all))];
    plif_clim   = [min(cellfun(@(x) min(x(:)), plif_raw)),   ...
                   max(cellfun(@(x) max(x(:)), plif_raw))];

    for i = 1:numel(phi_values)
        phi_str = sprintf('\\phi = %.1f', phi_values(i));

        figure('Color', 'k', 'Name', sprintf('phi=%.1f', phi_values(i)));
        colormap(gcf, hot);

        subplot(1, 3, 1);
        imagesc(ohstar_raw{i}, ohstar_clim);
        axis image; axis off;
        colorbar('Color', 'w');
        title(['OHstar Raw  [' phi_str ']'], 'Color', 'w', 'Interpreter', 'tex');
        set(gca, 'Color', 'k');

        subplot(1, 3, 2);
        imagesc(recon_all{i}, recon_clim);
        axis image; axis off;
        colorbar('Color', 'w');
        title(['OHstar After Inversion  [' phi_str ']'], 'Color', 'w', 'Interpreter', 'tex');
        set(gca, 'Color', 'k');

        subplot(1, 3, 3);
        imagesc(plif_raw{i}, plif_clim);
        axis image; axis off;
        colorbar('Color', 'w');
        title(['PLIF  [' phi_str ']'], 'Color', 'w', 'Interpreter', 'tex');
        set(gca, 'Color', 'k');
    end
end
