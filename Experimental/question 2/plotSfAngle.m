function plotSfAngle(phi_CH4, phi_H2, results_H2, results_CH4)
% PLOTSFANGLE  Plot flame speed angle Sf vs equivalence ratio for CH4 and H2.
%
%   plotSfAngle(phi_CH4, phi_H2, results)
%
%   Inputs:
%       phi_CH4  - array of equivalence ratios for methane cases
%       phi_H2   - array of equivalence ratios for hydrogen cases
%       results  - struct array where results(i).Sf_angle is the flame angle

    Sf_H2 = [results_H2.SF_angle];
    Sf_CH4 = [results_CH4.SF_angle];



% --- CH4 plot ---
    figure('Name', 'Sf vs phi - CH4');
    hold on; % Allows plotting multiple datasets on the same figure
    
    % Dataset 1: Experiments
    plot(phi_CH4, Sf_CH4, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.2 0.4 0.8]);
    
    % Dataset 2: Literature
    [phi_CH4_lit,Sf_CH4_lit] = literature_data_CH4();
    plot(phi_CH4_lit, Sf_CH4_lit, 's--', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
    
    % Labels and formatting
    xlabel('\phi  (equivalence ratio)');
    ylabel('S_f  (°)');
    title('Flame angle S_f vs \phi  —  CH_4');
    
    % Add legend to distinguish the data
    legend('Experiments', 'Literature', 'Location', 'best');
    
    grid on;
    hold off;

% --- H2 plot ---
    figure('Name', 'Sf vs phi - H2');
    hold on; % Allows plotting multiple datasets on the same figure
    
    % Dataset 1: Experiments
    plot(phi_H2, Sf_H2, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
    
    % Dataset 2: Literature (Adjust variable names 'phi_H2_lit' and 'Sf_H2_lit' as needed)
    [phi_H2_lit,Sf_H2_lit] = literature_data_H2();
    plot(phi_H2_lit, Sf_H2_lit, 'd--', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.1 0.6 0.3]);
    
    % Labels and formatting
    xlabel('\phi  (equivalence ratio)');
    ylabel('S_f  (°)');
    title('Flame angle S_f vs \phi  —  H_2');
    
    % Add legend to distinguish the data
    legend('Experiments', 'Literature', 'Location', 'best');
    
    grid on;
    hold off;

end