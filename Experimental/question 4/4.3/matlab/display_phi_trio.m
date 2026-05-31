%% display_phi_trio.m
%  Display OH* Raw, OH* Abel-Inverted, and OH-PLIF for each phi.
%  Uses readB16.m for file reading and Python (abel_invert.py) for the
%  Abel inversion.
%
%  BEFORE FIRST RUN — execute once in the Command Window:
%    pyenv('Version', 'C:\Users\franc\AppData\Local\Programs\Python\Python311\python.exe')
%
%  Place readB16.m and abel_invert.py in the same folder as this script.

clc; clear; close all;

%% ─── CONFIGURE THESE ───────────────────────────────────────────────
BASE_OHSTAR = 'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 4\OHstar';
BASE_PLIF   = 'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\LIF';
phi_values  = 0.5:0.1:1.3;
%% ───────────────────────────────────────────────────────────────────

% Make sure Python can find abel_invert.py in this folder
script_dir = fileparts(mfilename('fullpath'));
if count(py.sys.path, script_dir) == 0
    py.sys.path().insert(int32(0), script_dir);
end
mod = py.importlib.import_module('abel_invert');
py.importlib.reload(mod);

%% Loop over phi values
for k = 1:numel(phi_values)
    phi = phi_values(k);
    tag = sprintf('H2_PHI%.1f', phi);

    oh_path   = fullfile(BASE_OHSTAR, tag, 'B16', 'B0001.b16');
    plif_path = fullfile(BASE_PLIF,   tag, 'B16', 'B0001.b16');

    % Read images
    oh_raw  = double(readB16(oh_path));
    plif    = double(readB16(plif_path));

    % Abel inversion via Python
    py_result = mod.abel_invert(py.numpy.array(oh_raw));
    oh_inv    = double(py_result);

    % --- OH* Raw ---
    figure('Color', 'k', 'Name', sprintf('OH* Raw  phi=%.1f', phi));
    imagesc(oh_raw);
    axis image; axis off; colormap(hot);
    caxis([0, prctile(oh_raw(oh_raw > 0), 99)]);
    cb = colorbar; cb.Color = 'w';
    title(sprintf('OH* Raw  \\phi = %.1f', phi), 'Color', 'w', 'Interpreter', 'tex');
    set(gca, 'Color', 'k');

    % --- OH* Abel Inverted ---
    figure('Color', 'k', 'Name', sprintf('OH* Abel  phi=%.1f', phi));
    imagesc(oh_inv);
    axis image; axis off; colormap(hot);
    pos_vals = oh_inv(oh_inv > 0);
    if ~isempty(pos_vals)
        caxis([0, prctile(pos_vals, 99)]);
    end
    cb = colorbar; cb.Color = 'w';
    title(sprintf('OH* Abel Inverted  \\phi = %.1f', phi), 'Color', 'w', 'Interpreter', 'tex');
    set(gca, 'Color', 'k');

    % --- OH-PLIF ---
    figure('Color', 'k', 'Name', sprintf('PLIF  phi=%.1f', phi));
    imagesc(plif);
    axis image; axis off; colormap(hot);
    caxis([0, prctile(plif(plif > 0), 99)]);
    cb = colorbar; cb.Color = 'w';
    title(sprintf('OH-PLIF  \\phi = %.1f', phi), 'Color', 'w', 'Interpreter', 'tex');
    set(gca, 'Color', 'k');

    fprintf('  phi = %.1f  done\n', phi);
end
