function plotSfAngle(phi_CH4, phi_H2, results_H2, results_CH4)
% PLOTSFANGLE  Plot flame speed angle Sf vs equivalence ratio for CH4 and H2.
%
%   plotSfAngle(phi_CH4, phi_H2, results_H2, results_CH4)
%
%   Inputs:
%       phi_CH4     - array of equivalence ratios for methane cases
%       phi_H2      - array of equivalence ratios for hydrogen cases
%       results_H2  - struct array for H2 where results_H2(i).SF_angle is the flame angle
%       results_CH4 - struct array for CH4 where results_CH4(i).SF_angle is the flame angle

    Sf_H2 = [results_H2.SF_from_angle];
    Sf_CH4 = [results_CH4.SF_from_angle];

    % ==========================================
    % --- CH4 plot ---
    % ==========================================
    figure('Name', 'Sf vs phi - CH4');
    hold on; % Allows plotting multiple datasets on the same figure
    
    % Dataset 1: Experiments
    plot(phi_CH4, Sf_CH4, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.2 0.4 0.8]);
    
    % Dataset 2: Literature
    [phi_CH4_lit, Sf_CH4_lit] = literature_data_CH4();
    plot(phi_CH4_lit, Sf_CH4_lit, 's--', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
    
    % Dataset 3: Cantera Numerical Results
    [phi_CH4_cantera, Sf_CH4_cantera] = cantera_data_CH4();
    plot(phi_CH4_cantera, Sf_CH4_cantera, '^:', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.5 0.1 0.7]);
    
    % Labels and formatting
    xlabel('\phi  (equivalence ratio)');
    ylabel('S_f  (°)');
    title('Flame angle S_f vs \phi  —  CH_4');
    
    % Add legend to distinguish the data
    legend('Experiments', 'Literature', 'Cantera', 'Location', 'best');
    
    grid on;
    hold off;

    % ==========================================
    % --- H2 plot ---
    % ==========================================
    figure('Name', 'Sf vs phi - H2');
    hold on; % Allows plotting multiple datasets on the same figure
    
    % Dataset 1: Experiments
    plot(phi_H2, Sf_H2, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
    
    % Dataset 2: Literature
    [phi_H2_lit, Sf_H2_lit] = literature_data_H2();
    plot(phi_H2_lit, Sf_H2_lit, 'd--', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.1 0.6 0.3]);
    
    % Dataset 3: Cantera Numerical Results
    [phi_H2_cantera, Sf_H2_cantera] = cantera_data_H2();
    plot(phi_H2_cantera, Sf_H2_cantera, '^:', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.9 0.5 0.1]);
    
    % Labels and formatting
    xlabel('\phi  (equivalence ratio)');
    ylabel('S_f  (°)');
    title('Flame angle S_f vs \phi  —  H_2');
    
    % Add legend to distinguish the data
    legend('Experiments', 'Literature', 'Cantera', 'Location', 'best');
    
    grid on;
    hold off;
end