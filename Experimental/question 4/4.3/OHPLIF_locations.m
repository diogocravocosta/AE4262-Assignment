function [H2_image_locs, CH4_image_loc] = OHPLIF_locations()
%returns location of all figures

phi_values = 0.5:0.1:1.3;
H2_image_locs = arrayfun(@(phi) ['C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\LIF\H2_PHI' num2str(phi,'%.1f') '\B16\B0001.b16'], phi_values, 'UniformOutput', false);

phi_values = 0.8:0.2:1.2;
CH4_image_loc = arrayfun(@(phi) ['C:\Users\franc\Documents\GitHub\AE4262-Assignment\Experimental\question 2\LIF\CH4_PHI' num2str(phi,'%.1f') '\B16\B0001.b16'], phi_values, 'UniformOutput', false);

end