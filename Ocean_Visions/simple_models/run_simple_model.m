% RUN_SIMPLE_MODEL  Master script: choose the complexity level and the physics,
%                   draw the configuration, run the "steady state -> emissivity
%                   perturbation" experiment, and print where the forcing went.
%
%   >> run_simple_model
%
% Everything adjustable is in the CHOICES block.  The experiment is the one in
% the 26 Oct 2026 notes: spin the column up under a liquid cloud until it is in
% steady state, then change the cloud's emissivity (glaciation, with or without
% the precipitation loss of [V22]; or complete dissipation; or a prescribed flux
% change) and follow the skin, the air, the ice temperatures and the thickness.
% The partition of the forcing -- less upwelling LW, more sensible heat, latent
% heat, conduction out of the ice; and conduction into cooling the ice versus
% making ice -- is printed and drawn.  With COMPARE_LEVELS the same experiment
% is repeated at every complexity level so the hierarchy answers the question
% "what does each added layer change".
%
% Figures go to ./figures.  Run time: ~1 min at one level, ~3 min for all four.

clear; clc;
here = fileparts(mfilename('fullpath'));  if isempty(here), here = pwd; end
addpath(here);
fig_dir = fullfile(here, 'figures');  if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ========================== CHOICES ====================================
LEVEL = 3;                       % 1 | 2 | 3 | 4   (see model_config / seb_column_params)
PHYSICS = { ...
    'lw_emission',   'planck', ...        % 'planck' sigma T^4 | 'linearized' tangent (North et al. 1981)
    'turb_scheme',   'bulk_const', ...    % 'linear' | 'bulk_const' | 'most' (IFS Monin-Obukhov)
    'latent_scheme', 'bulk', ...          % 'bulk' | 'off'
    'atm_profile',   'auto', ...          % 'auto' (level 4: 'lapse', else 'isothermal') | 'isothermal' | 'lapse'
    'lapse_kind',    'dry', ...           % 'dry' | 'moist' | number [K/m], negative = inversion
    'ice_growth',    'stefan', ...        % 'stefan' | 'fixed' (ERA5-like constant thickness)
    'atm_forcing',   'relaxation', ...    % 'relaxation' (to T_bg on tau_adv) | 'heating' (prescribed D, Thorndike 1992)
    'D_adv_W_m2',    100, ...             % the heat convergence for 'heating' [W m-2]
    'h_init_m',      1.5, ...             % initial ice thickness [m]
    'h_snow_m',      0.0, ...             % snow on the ice [m] (levels 2-4)
    'use_brine',     false, ...           % Bitz & Lipscomb (1999) brine corrections
    'WP_g_m2',       30, ...              % cloud water path [g m-2]
    'phase',         'liquid', ...        % 'liquid' | 'ice' | 'none'
    'tau_adv_s',     86400};              % advective relaxation time [s]
PERTURBATION = 'glaciate_precip';  % 'glaciate_precip' (phase -> ice, 95 % of mass precipitates [V22])
                                   % 'glaciate' (phase -> ice, same mass) | 'dissipate' (eps -> 0)
                                   % 'flux' (prescribed change F_EXTRA of absorbed flux at the surface)
F_REMAIN     = 0.05;               % condensate remaining after 'glaciate_precip' [V22 Table 2]
R_E_ICE_UM   = 30;                 % ice effective radius after glaciation
F_EXTRA_W_M2 = -10;                % for PERTURBATION = 'flux'
SPIN_DAYS    = 30;
BRANCH_DAYS  = 30;
COMPARE_LEVELS = true;             % repeat at levels 1..4 with the same physics
DPI = 150;
%% ========================================================================

% ------------------------------------------------------------------------
% 1. configuration and its picture
% ------------------------------------------------------------------------
P = model_config(LEVEL, PHYSICS{:});
[~, png] = draw_model_diagram(P, fig_dir);
fprintf('diagram -> %s\n', png);

% ------------------------------------------------------------------------
% 2. the experiment at the chosen level
% ------------------------------------------------------------------------
[D, Rc, Rp, Rs, tag] = run_perturbation_experiment(P, PERTURBATION, F_REMAIN, R_E_ICE_UM, F_EXTRA_W_M2, SPIN_DAYS, BRANCH_DAYS);
fprintf('\n=== level %d: %s ===\n', LEVEL, tag);
print_column_state(Rs, 'spun-up state (end of spin-up)');
print_partition(D, BRANCH_DAYS);
plot_experiment(Rc, Rp, D, P, tag, fig_dir, DPI);

% ------------------------------------------------------------------------
% 3. the same experiment at every level
% ------------------------------------------------------------------------
if COMPARE_LEVELS
    levels = 1:4;  rows = [];
    fprintf('\n=== the hierarchy: %s at every level (same physics options) ===\n', tag);
    fprintf('  %5s %7s %8s %8s %8s %8s %8s %8s %8s %9s %9s %8s %8s\n', 'level', 'dLWD', 'T_sk,ctl', 'dT_sk', 'f_LWU', 'f_SH', 'f_LH', 'f_cool', 'f_grow', 'dh [cm]', 'bound', 'eff %', 'alpha');
    for lev = levels
        Pl = model_config(lev, PHYSICS{:}, 'quiet', true);
        draw_model_diagram(Pl, fig_dir);
        [Dl, Rcl, ~, ~] = run_perturbation_experiment(Pl, PERTURBATION, F_REMAIN, R_E_ICE_UM, F_EXTRA_W_M2, SPIN_DAYS, BRANCH_DAYS);
        Ef = Dl.E_forcing(end);
        rows(end+1,:) = [lev, Dl.dLWD(1), Rcl.T_sk_K(end), Dl.dT_sk(end), Dl.E_LWU(end)/Ef, Dl.E_SH(end)/Ef, Dl.E_LH(end)/Ef, ...
                         Dl.E_icecool(end)/Ef, Dl.E_growth(end)/Ef, 100*Dl.dh(end), 100*Dl.dh_zapalac(end), 100*Dl.dh(end)/Dl.dh_zapalac(end), Dl.alpha_air]; %#ok<SAGROW>
        fprintf('  %5d %7.1f %8.1f %8.1f %8.2f %8.2f %8.2f %8.2f %8.2f %9.2f %9.1f %8.1f %8.2f\n', rows(end,:));
    end
    fprintf('  f_* = share of the cumulative forcing over %d days; f_cool + f_grow = conduction share; bound = Zapalac (all forcing -> freezing)\n', BRANCH_DAYS);
    plot_level_comparison(rows, tag, fig_dir, DPI);
end
fprintf('\nfigures in %s\n', fig_dir);

% ========================================================================
% local functions
% ========================================================================
function [D, Rc, Rp, Rs, tag] = run_perturbation_experiment(P, kind, f_remain, r_e_ice, F_extra, spin_days, branch_days)
% Spin up, then branch a control and a perturbed run from the same state.
P.num.t_end_days = spin_days;
Rs = seb_column_run(P, 'verbose', false);
Pc = P;  Pc.num.t_end_days = branch_days;
Rc = seb_column_run(Pc, 'init_state', Rs.final_state, 'verbose', false);
Pp = Pc;
switch lower(kind)
    case 'glaciate_precip'
        Pp.cloud.phase = 'ice';  Pp.cloud.WP_g_m2 = f_remain * P.cloud.WP_g_m2;  Pp.cloud.r_e_ice_um = r_e_ice;
        tag = sprintf('glaciate+precipitate (LWP %g -> IWP %g g/m2)', P.cloud.WP_g_m2, Pp.cloud.WP_g_m2);
    case 'glaciate'
        Pp.cloud.phase = 'ice';  Pp.cloud.r_e_ice_um = r_e_ice;
        tag = sprintf('glaciate, same mass (%g g/m2)', P.cloud.WP_g_m2);
    case 'dissipate'
        Pp.cloud.phase = 'none';
        tag = 'cloud dissipates';
    case 'flux'
        Pp.atm.F_extra_down_W_m2 = F_extra;
        tag = sprintf('prescribed %+g W/m2 at the surface', F_extra);
    otherwise
        error('unknown PERTURBATION "%s"', kind);
end
Pp = seb_column_finalize(Pp);             % recompute the cloud and slab emissivities
Rp = seb_column_run(Pp, 'init_state', Rs.final_state, 'verbose', false);
D  = response_partition(Rp, Rc);
end

function print_partition(D, days)
fprintf('  forcing at t = 0+: dLWD = %+.1f W/m2\n', D.dLWD(1));
fprintf('  instantaneous shares (sum = 1):\n  %8s %8s %8s %8s %8s %9s %9s\n', 't', 'f_LWU', 'f_SH', 'f_LH', 'f_G', 'dT_sk', 'dT_air');
for th = [1 6 24 72 240 24*days]
    [~, i] = min(abs(D.t_h - th));
    fprintf('  %6.0f h %8.3f %8.3f %8.3f %8.3f %9.2f %9.2f\n', th, D.f_LWU(i), D.f_SH(i), D.f_LH(i), D.f_G(i), D.dT_sk(i), D.dT_air(i));
end
Ef = D.E_forcing(end);
fprintf('  cumulative after %d days: forcing %.1f MJ/m2 -> LW up %.0f %%, sensible %.0f %%, latent %.0f %%, conduction %.0f %% (cooling the ice %.0f %%, making ice %.0f %%)\n', ...
    days, Ef, 100*D.E_LWU(end)/Ef, 100*D.E_SH(end)/Ef, 100*D.E_LH(end)/Ef, 100*D.E_G(end)/Ef, 100*D.E_icecool(end)/Ef, 100*D.E_growth(end)/Ef);
fprintf('  extra ice %.2f cm vs %.1f cm if all forcing froze water (Zapalac 2023): efficiency %.0f %% | alpha = dT_air/dT_sk = %.2f\n', ...
    100*D.dh(end), 100*D.dh_zapalac(end), 100*D.dh(end)/D.dh_zapalac(end), D.alpha_air);
end

function plot_experiment(Rc, Rp, D, P, tag, fig_dir, dpi)
t = Rc.t_days;  nL = P.atm.n_layers;  N = size(Rc.T_ice_K, 1);
cc = [0.35 0.35 0.35];  cp = [0.80 0.15 0.15];
hf = new_figure_clean('Position', [50 50 1400 900]);
tl = tiledlayout(hf, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, Rc.LWD_sfc_W_m2, '--', 'Color', cc, 'LineWidth', 1.5, 'DisplayName', 'LWD control');
plot(ax, t, Rp.LWD_sfc_W_m2, '-', 'Color', cp, 'LineWidth', 1.8, 'DisplayName', 'LWD perturbed');
plot(ax, t, Rc.LWU_sfc_W_m2, ':', 'Color', cc, 'LineWidth', 1.5, 'DisplayName', 'LWU control');
plot(ax, t, Rp.LWU_sfc_W_m2, ':', 'Color', cp, 'LineWidth', 1.8, 'DisplayName', 'LWU perturbed');
ylabel(ax, 'W m^{-2}'); title(ax, '(a) longwave at the surface'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax, 'on');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, Rc.T_sk_K, '--', 'Color', cc, 'LineWidth', 1.5, 'DisplayName', 'T_{skin} control');
plot(ax, t, Rp.T_sk_K, '-', 'Color', cp, 'LineWidth', 1.8, 'DisplayName', 'T_{skin}');
cl = lines(nL);
for k = 1:nL, plot(ax, t, Rp.T_atm_K(k,:), '-', 'Color', cl(k,:), 'LineWidth', 1.3, 'DisplayName', P.atm.layers(k).name); end
ylabel(ax, 'K'); title(ax, '(b) temperatures (perturbed run)'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax, 'on');
ax = nexttile(tl); hold(ax, 'on');
ci = copper(N+2);
for l = 1:N
    plot(ax, t, Rc.T_ice_K(l,:), '--', 'Color', ci(l+1,:), 'LineWidth', 1.3, 'HandleVisibility', 'off');
    plot(ax, t, Rp.T_ice_K(l,:), '-', 'Color', ci(l+1,:), 'LineWidth', 1.8, 'DisplayName', sprintf('ice layer %d', l));
end
yline(ax, P.ocean.T_f_K, 'k:', 'T_f', 'HandleVisibility', 'off');
ylabel(ax, 'K'); title(ax, '(c) ice temperatures (dashed = control)'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax, 'on');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, 100*D.dh, '-', 'Color', cp, 'LineWidth', 2, 'DisplayName', 'perturbed - control');
plot(ax, t, 100*D.dh_zapalac, 'k:', 'LineWidth', 1.4, 'DisplayName', 'Zapalac bound: all forcing -> freezing');
ylabel(ax, 'cm'); title(ax, '(d) extra ice thickness'); legend(ax, 'Location', 'northwest', 'FontSize', 8); grid(ax, 'on');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, 100*86400*Rc.dhdt_bot_m_s, '--', 'Color', cc, 'LineWidth', 1.5, 'DisplayName', 'control');
plot(ax, t, 100*86400*Rp.dhdt_bot_m_s, '-', 'Color', cp, 'LineWidth', 1.8, 'DisplayName', 'perturbed');
xlabel(ax, 'days after perturbation'); ylabel(ax, 'cm day^{-1}'); title(ax, '(e) basal growth rate'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax, 'on');
ax = nexttile(tl); hold(ax, 'on');
tc = [0.30 0.45 0.69; 0.33 0.66 0.41; 0.51 0.45 0.70];
plot(ax, t, Rc.SH_down_W_m2, '--', 'Color', tc(1,:), 'LineWidth', 1.2, 'HandleVisibility', 'off'); plot(ax, t, Rp.SH_down_W_m2, '-', 'Color', tc(1,:), 'LineWidth', 1.8, 'DisplayName', 'SH (down)');
plot(ax, t, Rc.LH_down_W_m2, '--', 'Color', tc(2,:), 'LineWidth', 1.2, 'HandleVisibility', 'off'); plot(ax, t, Rp.LH_down_W_m2, '-', 'Color', tc(2,:), 'LineWidth', 1.8, 'DisplayName', 'LH (down)');
plot(ax, t, Rc.F_cond_top_up_W_m2, '--', 'Color', tc(3,:), 'LineWidth', 1.2, 'HandleVisibility', 'off'); plot(ax, t, Rp.F_cond_top_up_W_m2, '-', 'Color', tc(3,:), 'LineWidth', 1.8, 'DisplayName', 'conduction out of the ice (top)');
plot(ax, t, Rp.F_cond_bot_up_W_m2, ':', 'Color', tc(3,:), 'LineWidth', 1.6, 'DisplayName', 'conduction at the base');
xlabel(ax, 'days after perturbation'); ylabel(ax, 'W m^{-2}'); title(ax, '(f) surface fluxes (dashed = control)'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax, 'on');
title(tl, sprintf('Level %d: %s, \\DeltaLWD = %+.0f W m^{-2}  [%s, %s, %s]', P.opt.complexity, tag, D.dLWD(1), P.opt.lw_emission, P.opt.turb_scheme, P.opt.atm_profile));
exportgraphics(hf, fullfile(fig_dir, sprintf('level%d_timeseries.png', P.opt.complexity)), 'Resolution', dpi); close(hf);

hf = new_figure_clean('Position', [50 50 1300 500]);
tl = tiledlayout(hf, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, D.t_h, D.f_LWU, '-', 'Color', [0.70 0.09 0.17], 'LineWidth', 2, 'DisplayName', 'upwelling LW');
plot(ax, D.t_h, D.f_SH, '-', 'Color', tc(1,:), 'LineWidth', 2, 'DisplayName', 'sensible heat');
plot(ax, D.t_h, D.f_LH, '-', 'Color', tc(2,:), 'LineWidth', 2, 'DisplayName', 'latent heat');
plot(ax, D.t_h, D.f_G, '-', 'Color', tc(3,:), 'LineWidth', 2, 'DisplayName', 'conduction out of the ice');
set(ax, 'XScale', 'log'); xlim(ax, [D.t_h(1) D.t_h(end)]); ylim(ax, [-0.1 1]); grid(ax, 'on');
xlabel(ax, 'hours after perturbation'); ylabel(ax, 'share of \DeltaLWD'); title(ax, '(a) instantaneous partition'); legend(ax, 'Location', 'east', 'FontSize', 8);
ax = nexttile(tl); hold(ax, 'on');
E = -[D.E_forcing(end), D.E_LWU(end), D.E_SH(end), D.E_LH(end), D.E_icecool(end), D.E_growth(end), D.E_ocean(end)];
labels = {'forcing', 'LW up', 'sensible', 'latent', 'ice cooling', 'ice growth', 'ocean flux'};
bc = [0.3 0.3 0.3; 0.70 0.09 0.17; tc(1,:); tc(2,:); 0.6 0.55 0.8; 0.35 0.25 0.55; 0.6 0.6 0.6];
for i = 1:numel(E)
    bar(ax, i, E(i), 0.7, 'FaceColor', bc(i,:), 'EdgeColor', 'none');
    text(ax, i, E(i) + 0.02*max(abs(E))*sign(E(i)+eps), sprintf('%.1f', E(i)), 'HorizontalAlignment', 'center', 'FontSize', 9);
end
yline(ax, 0, 'k-'); set(ax, 'XTick', 1:numel(E), 'XTickLabel', labels); ylabel(ax, 'MJ m^{-2} removed from the surface'); grid(ax, 'on');
title(ax, sprintf('(b) cumulative: %.1f cm of extra ice vs %.1f cm bound', 100*D.dh(end), 100*D.dh_zapalac(end)));
title(tl, sprintf('Level %d: where does the forcing go?  (%s)', P.opt.complexity, tag));
exportgraphics(hf, fullfile(fig_dir, sprintf('level%d_partition.png', P.opt.complexity)), 'Resolution', dpi); close(hf);
end

function plot_level_comparison(rows, tag, fig_dir, dpi)
hf = new_figure_clean('Position', [50 50 1300 480]);
tl = tiledlayout(hf, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
shares = rows(:, 5:9);          % LWU SH LH cool grow
cols = [0.70 0.09 0.17; 0.30 0.45 0.69; 0.33 0.66 0.41; 0.6 0.55 0.8; 0.35 0.25 0.55];
b = bar(ax, rows(:,1), shares, 'stacked');
for i = 1:5, b(i).FaceColor = cols(i,:); b(i).EdgeColor = 'none'; end
legend(ax, {'LW up', 'sensible', 'latent', 'ice cooling', 'ice growth'}, 'Location', 'eastoutside', 'FontSize', 8);
xlabel(ax, 'complexity level'); ylabel(ax, 'share of the cumulative forcing'); ylim(ax, [0 1.05]); set(ax, 'XTick', rows(:,1)); grid(ax, 'on');
title(ax, '(a) partition of the forcing by level');
ax = nexttile(tl); hold(ax, 'on');
yyaxis(ax, 'left'); bar(ax, rows(:,1), rows(:,10), 0.5, 'FaceColor', [0.35 0.25 0.55], 'EdgeColor', 'none'); ylabel(ax, 'extra ice after the run [cm]');
yyaxis(ax, 'right'); plot(ax, rows(:,1), rows(:,2), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); ylabel(ax, '\DeltaLWD at t = 0^+ [W m^{-2}]');
xlabel(ax, 'complexity level'); set(ax, 'XTick', rows(:,1)); grid(ax, 'on'); title(ax, '(b) extra ice and the forcing by level');
title(tl, sprintf('The hierarchy: %s', tag));
exportgraphics(hf, fullfile(fig_dir, 'levels_comparison.png'), 'Resolution', dpi); close(hf);
end
