clc; clear; close all;

% Run once per MATLAB session before any py.* call:
%   pyenv('Version', 'C:\Users\franc\AppData\Local\Programs\Python\Python311\python.exe')

% Add this directory to Python's path so abel_inversion_data.py is importable
script_dir = fileparts(mfilename('fullpath'));
if count(py.sys.path, script_dir) == 0
    py.sys.path().insert(int32(0), script_dir);
end

phi_values = 0.5:0.1:1.3;

%% 1. Get paths
ohstar_paths    = get_paths();
[plif_paths, ~] = OHPLIF_locations();

%% 2. Load raw OHstar
ohstar_raw = cell(numel(ohstar_paths), 1);
for i = 1:numel(ohstar_paths)
    ohstar_raw{i} = double(readB16(ohstar_paths{i}));
end

%% 3. Abel inversion on OHstar only
recon_all = abel_batch(ohstar_paths);

%% 4. Load raw PLIF
plif_raw = cell(numel(plif_paths), 1);
for i = 1:numel(plif_paths)
    plif_raw{i} = double(readB16(plif_paths{i}));
end

%% 5. Display — one figure per phi, ordered 0.5 → 1.3
display_phi_trio(ohstar_raw, recon_all, plif_raw, phi_values);
