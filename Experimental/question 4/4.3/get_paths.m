function paths = get_paths()
% GET_PATHS  Return cell array of B16 image paths for H2 OH* data (phi 0.5–1.3).
    base = 'C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 4\OHstar';
    phi_values = 0.5:0.1:1.3;
    paths = cell(numel(phi_values), 1);
    for i = 1:numel(phi_values)
        phi_str = sprintf('%.1f', phi_values(i));
        paths{i} = fullfile(base, ['H2_PHI' phi_str], 'B16', 'B0001.b16');
    end
end
