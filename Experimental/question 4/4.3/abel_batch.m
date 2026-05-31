function recon_all = abel_batch(paths)
% ABEL_BATCH  Run Abel inversion on each path using New_abel_imversion.py.
%
%   recon_all = abel_batch(paths)
%
%   INPUT
%       paths     - N×1 cell array of .b16 file paths (from get_paths)
%
%   OUTPUT
%       recon_all - N×1 cell array; recon_all{i} is the Abel-inverted
%                   array (double) for paths{i}

    mod = py.importlib.import_module('New_abel_imversion');
    py.importlib.reload(mod);

    n = numel(paths);
    recon_all = cell(n, 1);

    for i = 1:n
        fprintf('Processing %d/%d: %s\n', i, n, paths{i});
        result = mod.abel_invert_b16(paths{i});
        recon_all{i} = double(result);
    end
end
