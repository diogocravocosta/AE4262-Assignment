%% Boundary Layer Flashback Analysis
clear; clc; close all;

%% ── Air properties at 300 K, 1 atm ──────────────────────────────────────
nu_air  = 1.57e-5;
mu_air  = 1.85e-5;
rho_air = 1.17;

%% ── Geometry ─────────────────────────────────────────────────────────────
D_H2  = 4e-3;
D_CH4 = 21e-3;
LD    = 50;

L_H2  = LD * D_H2;
L_CH4 = LD * D_CH4;

dq_H2  = 0.64e-3;
dq_CH4 = 2.5e-3;

%% ── Load experimental and Cantera data ──────────────────────────────────
[data_H2, data_CH4] = Get_experimental_data();

phi_H2  = data_H2(:,1);   V_H2  = data_H2(:,2);   rho_H2  = data_H2(:,3);
phi_CH4 = data_CH4(:,1);  V_CH4 = data_CH4(:,2);  rho_CH4 = data_CH4(:,3);

[phi_can_H2,  sf_H2]  = cantera_data_H2();
[phi_can_CH4, sf_CH4] = cantera_data_CH4();

if ~isequal(phi_H2, phi_can_H2)
    error('H2: phi vectors from Get_experimental_data and cantera_data_H2 do not match.');
end
if ~isequal(phi_CH4, phi_can_CH4)
    error('CH4: phi vectors from Get_experimental_data and cantera_data_CH4 do not match.');
end

%% ── Wall velocity gradients ──────────────────────────────────────────────
gF_H2_P  = 8 * V_H2  / D_H2;
gF_CH4_P = 8 * V_CH4 / D_CH4;

gF_H2_D  = arrayfun(@(V) gF_shah(V, D_H2,  L_H2,  nu_air), V_H2);
gF_CH4_D = arrayfun(@(V) gF_shah(V, D_CH4, L_CH4, nu_air), V_CH4);

gc_H2  = (sf_H2  * 1e-3) / (dq_H2  / 3);
gc_CH4 = (sf_CH4 * 1e-3) / (dq_CH4 / 3);

gc_H2_nd  = (sf_H2  * 1e-3) / dq_H2;
gc_CH4_nd = (sf_CH4 * 1e-3) / dq_CH4;

%% ── Diagnostic: print L+ values ─────────────────────────────────────────
fprintf('\n── Dimensionless tube length  L+ = (L/D) / Re ──────────────────\n');
fprintf('   H2  tube (D = %g mm, L/D = %d):\n', D_H2*1e3, LD);
Re_H2  = V_H2  * D_H2  / nu_air;
Re_CH4 = V_CH4 * D_CH4 / nu_air;
Lp_H2  = (L_H2 /D_H2 ) ./ Re_H2;
Lp_CH4 = (L_CH4/D_CH4) ./ Re_CH4;
fprintf('   Re range: %.0f – %.0f,  L+ range: %.4f – %.4f\n', ...
        min(Re_H2), max(Re_H2), min(Lp_H2), max(Lp_H2));
fprintf('   CH4 tube (D = %g mm, L/D = %d):\n', D_CH4*1e3, LD);
fprintf('   Re range: %.0f – %.0f,  L+ range: %.4f – %.4f\n', ...
        min(Re_CH4), max(Re_CH4), min(Lp_CH4), max(Lp_CH4));
fprintf('   (Fully developed at L+ > 0.06; entry region for smaller values)\n\n');

%% ── Console tables ───────────────────────────────────────────────────────
print_table('H2',  phi_H2,  gF_H2_P,  gF_H2_D,  gc_H2,  gc_H2_nd);
print_table('CH4', phi_CH4, gF_CH4_P, gF_CH4_D, gc_CH4, gc_CH4_nd);

%% ── Critical bulk velocity (flashback limit) ─────────────────────────────
V_lim_H2  = zeros(size(phi_H2));
V_lim_CH4 = zeros(size(phi_CH4));

for i = 1:length(phi_H2)
    target = gc_H2(i);
    V_lim_H2(i) = fzero(@(V) gF_shah(V, D_H2, L_H2, nu_air) - target, [0.01, 50]);
end

for i = 1:length(phi_CH4)
    target = gc_CH4(i);
    V_lim_CH4(i) = fzero(@(V) gF_shah(V, D_CH4, L_CH4, nu_air) - target, [0.01, 50]);
end
%% ── Colour scheme ────────────────────────────────────────────────────────
col_pois    = [0.400 0.620 0.850];
col_dev     = [0.059 0.388 0.651];
col_crit    = [0.950 0.500 0.500];
col_crit_nd = [0.780 0.160 0.160];
col_H2      = [0.216 0.494 0.722];
col_CH4     = [0.894 0.447 0.024];

%% ── Figure 1: H2 flashback ───────────────────────────────────────────────
fig1 = figure('Position', [80 80 800 540], 'Color', 'w');
flashback_panel(fig1, phi_H2, gF_H2_P, gF_H2_D, gc_H2, gc_H2_nd, ...
                col_pois, col_dev, col_crit, col_crit_nd, 1.55, ...
    sprintf('H_2/air  |  D = %g mm,  \\delta_q = %g mm,  L/D = %d', ...
            D_H2*1e3, dq_H2*1e3, LD));

%% ── Figure 2: CH4 flashback ──────────────────────────────────────────────
fig2 = figure('Position', [120 120 800 540], 'Color', 'w');
flashback_panel(fig2, phi_CH4, gF_CH4_P, gF_CH4_D, gc_CH4, gc_CH4_nd, ...
                col_pois, col_dev, col_crit, col_crit_nd, 1.20, ...
    sprintf('CH_4/air  |  D = %g mm,  \\delta_q = %g mm,  L/D = %d', ...
            D_CH4*1e3, dq_CH4*1e3, LD));

%% ── Figure 3: Cantera laminar flame speeds ───────────────────────────────
fig3 = figure('Position', [160 160 700 480], 'Color', 'w');
figure(fig3);
hold on; grid on; box on;

plot(phi_H2,  sf_H2,  '-o', 'Color', col_H2,  'LineWidth', 2.0, ...
     'MarkerSize', 7, 'MarkerFaceColor', col_H2,  ...
     'DisplayName', 'H_2/air  (Cantera)');

plot(phi_CH4, sf_CH4, '-s', 'Color', col_CH4, 'LineWidth', 2.0, ...
     'MarkerSize', 7, 'MarkerFaceColor', col_CH4, ...
     'DisplayName', 'CH_4/air  (Cantera)');

xlabel('\phi  (equivalence ratio)', 'FontSize', 12);
ylabel('Laminar flame speed  S_L  [mm/s]', 'FontSize', 12);
title('Cantera Laminar Flame Speeds', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11);

%% ── Figure 4: V_bulk flashback limit – H2 ───────────────────────────────
fig4 = figure('Position', [200 200 700 460], 'Color', 'w');
hold on; grid on; box on;
plot(phi_H2, V_lim_H2, '-o', 'Color', col_H2, 'LineWidth', 2.0, ...
     'MarkerSize', 7, 'MarkerFaceColor', col_H2, 'DisplayName', 'V_{lim} [Shah]');
plot(phi_H2, V_H2, '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2.0, ...
     'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'V_{bulk} [exp]');
xline(1.55, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.4, 'HandleVisibility', 'off');
text(1.55 + 0.02, max([V_lim_H2; V_H2])*0.93, '\phi_{exp} = 1.55', ...
     'FontSize', 9, 'Color', [0.25 0.25 0.25]);
xlabel('\phi  (equivalence ratio)', 'FontSize', 12);
ylabel('V_{bulk}  [m/s]', 'FontSize', 12);
title('H_2/air  – Flashback limit bulk velocity', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);

%% ── Figure 5: V_bulk flashback limit – CH4 ──────────────────────────────
fig5 = figure('Position', [240 240 700 460], 'Color', 'w');
hold on; grid on; box on;
plot(phi_CH4, V_lim_CH4, '-s', 'Color', col_CH4, 'LineWidth', 2.0, ...
     'MarkerSize', 7, 'MarkerFaceColor', col_CH4, 'DisplayName', 'V_{lim} [Shah]');
plot(phi_CH4, V_CH4, '-s', 'Color', [0.5 0.5 0.5], 'LineWidth', 2.0, ...
     'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'V_{bulk} [exp]');
xline(1.20, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.4, 'HandleVisibility', 'off');
text(1.20 + 0.02, max([V_lim_CH4; V_CH4])*0.93, '\phi_{exp} = 1.20', ...
     'FontSize', 9, 'Color', [0.25 0.25 0.25]);
xlabel('\phi  (equivalence ratio)', 'FontSize', 12);
ylabel('V_{bulk}  [m/s]', 'FontSize', 12);
title('CH_4/air  – Flashback limit bulk velocity', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function gF = gF_shah(V, D, L, nu)
% GF_SHAH  Wall velocity gradient using Shah (1978) apparent friction factor.
    fRe_fd = 16;
    K_inf  = 1.25;
    C      = 0.000212;

    Re     = V * D / nu;
    x_plus = (L / D) / Re;

    term_short = 3.44 / sqrt(x_plus);
    term_long  = fRe_fd + K_inf / (4 * x_plus) - term_short;
    fRe_app    = term_short + term_long / (1 + C * x_plus^(-2));

    gF = fRe_app * V / (2 * D);
end


function print_table(fuel, phi, gF_P, gF_D, gc, gc_nd)
% PRINT_TABLE  Print flashback summary table to console.
    fprintf('── %s ──────────────────────────────────────────────────────────────────\n', fuel);
    fprintf('  phi    Re      L+       gF_Pois  gF_Shah   gc/3   gc   FB_Pois  FB_Shah\n');
    nu = 1.57e-5;
    D  = 4e-3 * strcmp(fuel,'H2') + 21e-3 * strcmp(fuel,'CH4');
    L  = 50 * D;
    for i = 1:length(phi)
        Re   = gF_P(i) * D / 8 * D / nu;
        Lp   = (L/D) / Re;
        fp = 'no'; fd = 'no';
        if gF_P(i) < gc(i), fp = 'YES'; end
        if gF_D(i) < gc(i), fd = 'YES'; end
        fprintf('  %.2f  %6.0f  %.4f  %12.1f  %12.1f  %10.1f  %10.1f   %-5s    %-5s\n', ...
                phi(i), Re, Lp, gF_P(i), gF_D(i), gc(i), gc_nd(i), fp, fd);
    end
    fprintf('\n');
end


function flashback_panel(fig, phi, gF_P, gF_D, gc, gc_nd, col_pois, col_dev, col_crit, col_crit_nd, phi_fb, ttl)
% FLASHBACK_PANEL  Plot wall gradients and critical gradients for one fuel.

    figure(fig);
    hold on; grid on; box on;

    % ── Shaded band between the two g_c lines (red) ───────────────────────
    fill([phi; flipud(phi)], [gc; flipud(gc_nd)], ...
         [0.980 0.820 0.820], 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
         'HandleVisibility', 'off');

    % ── Shaded band between the two g_F lines (blue) ──────────────────────
    gF_upper = max([gF_P, gF_D], [], 2);
    gF_lower = min([gF_P, gF_D], [], 2);
    fill([phi; flipud(phi)], [gF_upper; flipud(gF_lower)], ...
         [0.780 0.880 0.960], 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
         'HandleVisibility', 'off');

    % ── g_F curves (blue) ─────────────────────────────────────────────────
    h1 = plot(phi, gF_P, '-s', 'Color', col_pois, 'LineWidth', 2.0, ...
              'MarkerSize', 7, 'MarkerFaceColor', col_pois, ...
              'DisplayName', 'Fully developed [Poiseuille]');

    h2 = plot(phi, gF_D, '--^', 'Color', col_dev, 'LineWidth', 2.0, ...
              'MarkerSize', 7, 'MarkerFaceColor', col_dev, ...
              'DisplayName', 'Developing [Shah 1978]');

    % ── g_c curves (red) ──────────────────────────────────────────────────
    h3 = plot(phi, gc,    '-s', 'Color', col_crit, 'LineWidth', 2.0, ...
              'MarkerSize', 7, 'MarkerFaceColor', col_crit, ...
              'DisplayName', '\delta_p = \delta_q / 3  [Lewis & von Elbe]');

    h4 = plot(phi, gc_nd, '--^', 'Color', col_crit_nd, 'LineWidth', 2.0, ...
              'MarkerSize', 7, 'MarkerFaceColor', col_crit_nd, ...
              'DisplayName', '\delta_p = \delta_q');

    % ── Fixed experimental flashback vertical line ─────────────────────────
    y_top = max([gF_P; gF_D; gc; gc_nd]) * 1.08;

    xline(phi_fb, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.4, ...
          'HandleVisibility', 'off');
    text(phi_fb + 0.02, y_top * 0.93, ...
         sprintf('\\phi_{exp} = %.2f', phi_fb), ...
         'Color', [0.25 0.25 0.25], 'FontSize', 9, 'HorizontalAlignment', 'left');

    % ── Labels and formatting ─────────────────────────────────────────────
    xlabel('\phi  (equivalence ratio)', 'FontSize', 12);
    ylabel('Velocity gradient',  'FontSize', 12);
    title(ttl, 'FontSize', 12, 'FontWeight', 'bold');
    legend([h1 h2 h3 h4], 'Location', 'northwest', 'FontSize', 9, 'Box', 'on');
    set(gca, 'FontSize', 11);
    ylim([0 y_top]);
end