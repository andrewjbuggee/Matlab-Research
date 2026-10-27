% RUN_GLACIATION_EXPERIMENT  What happens to the sea ice when the cloud above it glaciates?
%
%   >> run_glaciation_experiment
%
% THE EXPERIMENT.  Spin the column up for 30 days under a liquid-bearing cloud
% (LWP = 30 g m-2, the winter mean over Arctic sea ice [L26]), then branch:
%
%   control                      the liquid cloud stays
%   glaciate, same mass          the condensate turns to ice (IWP = LWP, r_e 30 um):
%                                emissivity falls from 0.98 to ~0.74 -- the pure
%                                phase effect
%   glaciate + precipitate       the LES result of Villanueva et al. (2022, Table 2):
%                                at a 1 %/h freezing rate 95 % of the condensate
%                                precipitates out within 18 h; IWP = 0.05*LWP, the
%                                LW effect of the cloud nearly vanishes.  THIS IS THE
%                                HEADLINE CASE -- what MCT is actually meant to do.
%   cloud dissipates             eps_c = 0 (the "thinned out" end member of [L26])
%   ... + fixed atmosphere       the headline case with the atmosphere held at its
%                                cloudy state: the alpha = 0 limit of [LSM], i.e. the
%                                UPPER bound on the turbulent response
%   ... + 3x C_H                 stronger skin-air coupling, same forcing
%   ... + 20 cm snow             the real Arctic surface rather than ERA5's bare ice
%
% WHAT TO LOOK AT.
%   Fig. 1  time series: the forcing (LWD), the skin and air temperatures, the
%           three ice layers cooling from the top down, the thickness anomaly, the
%           growth rate, and the surface fluxes.
%   Fig. 2  WHERE THE FORCING GOES.  The reduction in LWD is balanced at every
%           instant by less upwelling LW (the skin cools), more downward sensible
%           heat (the air is now warmer than the skin), a little latent heat, and
%           more heat conducted up out of the ice.  The conductive share is then
%           split into COOLING THE ICE (sensible heat) and MAKING ICE (latent heat
%           at the base).  Ray's statement at the OV meeting -- "more energy goes
%           into cooling the ice than making more ice" -- is a number here, and so
%           is Zapalac's upper bound dh = dF*t/(rho*L), which assumes ALL of the
%           forcing goes into freezing.
%   Fig. 3  the same after 30 days as a function of initial ice thickness, with
%           and without snow: Matt Shupe's "growth is fastest for the thinnest ice".
%   Fig. 4  ice temperature profiles at successive times after glaciation: how the
%           surface cooling propagates down to the base -- the question in the
%           22 Oct 2026 notes.
%
% Figures are written to ./figures as PNG (invisible figures, exportgraphics).
% Run time: about a minute.

clear; clc;
here = fileparts(mfilename('fullpath'));  if isempty(here), here = pwd; end
addpath(here);
fig_dir = fullfile(here, 'figures');  if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ------------------------------------------------------------------------
% 0.  Configuration
% ------------------------------------------------------------------------
spin_days   = 30;
branch_days = 30;
LWP0_g_m2   = 30;              % [L26] winter LWP over sea ice
f_remain_V22 = 0.05;           % condensate remaining after glaciation, [V22] Table 2 (0.009/0.170)
r_e_ice_um  = 30;

P0 = seb_column_params('phase', 'liquid', 'WP_g_m2', LWP0_g_m2, 't_end_days', spin_days);
rho_i = P0.ice.rho_kg_m3;

% ------------------------------------------------------------------------
% 1.  Spin-up under the liquid cloud
% ------------------------------------------------------------------------
fprintf('\n=== 1. spin-up, %d days, liquid cloud LWP = %g g/m2 (eps_c = %.3f) ===\n', ...
    spin_days, LWP0_g_m2, P0.cloud.emis);
Rspin = seb_column_run(P0);
S0 = Rspin.final_state;
print_state(Rspin, 'spun-up cloudy state');

% ------------------------------------------------------------------------
% 2.  Branch
% ------------------------------------------------------------------------
fprintf('\n=== 2. branch runs, %d days from the same state ===\n', branch_days);
sc = struct('name', {}, 'P', {}, 'R', {}, 'style', {});
sc(end+1) = mk('control',                 P0, {});
sc(end+1) = mk('glaciate, same mass',     P0, {'phase','ice','WP_g_m2',LWP0_g_m2});
sc(end+1) = mk('glaciate + precipitate',  P0, {'phase','ice','WP_g_m2',f_remain_V22*LWP0_g_m2});
sc(end+1) = mk('cloud dissipates',        P0, {'phase','none'});
sc(end+1) = mk('glac.+precip., fixed atmosphere', P0, {'phase','ice','WP_g_m2',f_remain_V22*LWP0_g_m2,'mode','prescribed'});
sc(end+1) = mk('glac.+precip., 3x C_H',   P0, {'phase','ice','WP_g_m2',f_remain_V22*LWP0_g_m2,'C_H',3*P0.atm.C_H,'C_E',3*P0.atm.C_E});
for i = 1:numel(sc)
    sc(i).P.num.t_end_days = branch_days;
    sc(i).R = seb_column_run(sc(i).P, 'init_state', S0, 'verbose', false);
    fprintf('  %-34s eps_c = %.3f  LWD(0+) = %6.1f  T_sk(30d) = %6.1f  dh(30d) = %+5.2f cm\n', ...
        sc(i).name, sc(i).P.cloud.emis, sc(i).R.LWD_sfc_W_m2(1), sc(i).R.T_sk_K(end), ...
        100*(sc(i).R.h_m(end) - sc(1).R.h_m(end)));
end
% the 3x C_H case must be compared with its own control (the coefficient is a
% property of the column, not of the forcing)
P_ctrl3 = seb_column_params('phase','liquid','WP_g_m2',LWP0_g_m2,'t_end_days',spin_days,'C_H',3*P0.atm.C_H,'C_E',3*P0.atm.C_E);
Rspin3  = seb_column_run(P_ctrl3, 'verbose', false);
P_ctrl3.num.t_end_days = branch_days;
Rctrl3  = seb_column_run(P_ctrl3, 'init_state', Rspin3.final_state, 'verbose', false);
sc(6).R = seb_column_run(sc(6).P, 'init_state', Rspin3.final_state, 'verbose', false);
sc(6).control = Rctrl3;

% ------------------------------------------------------------------------
% 3.  Where does the forcing go?   (headline case vs control)
% ------------------------------------------------------------------------
ctrl = sc(1).R;  pert = sc(3).R;  fixed = sc(5).R;
D = response_budget(pert, ctrl, rho_i);
Dfix = response_budget(fixed, ctrl, rho_i);
D3  = response_budget(sc(6).R, Rctrl3, rho_i);
fprintf('\n=== 3. response partition, "glaciate + precipitate" minus control ===\n');
fprintf('  forcing at t = 0+: dLWD = %+.1f W/m2  (cloud eps %.3f -> %.3f)\n', D.dLWD(1), sc(1).P.cloud.emis, sc(3).P.cloud.emis);
fprintf('  instantaneous shares of dLWD (sum = 1 by the skin balance):\n');
fprintf('  %8s %8s %8s %8s %8s %8s | %8s %8s   %s\n', 't', 'f_LWU', 'f_SH', 'f_LH', 'f_G', 'dT_sk', 'f_SH', 'f_G', '<- fixed atmosphere');
for th = [1 6 24 72 240 720]
    [~, i] = min(abs(D.t_h - th));
    fprintf('  %6.0f h %8.3f %8.3f %8.3f %8.3f %8.2f | %8.3f %8.3f\n', th, D.f_LWU(i), D.f_SH(i), D.f_LH(i), D.f_G(i), D.dT_sk(i), Dfix.f_SH(i), Dfix.f_G(i));
end
fprintf('\n  cumulative energy after %d days [MJ/m2]  (negative = energy removed from the surface):\n', branch_days);
fprintf('  %-46s %8.2f\n', 'forcing  integral(dLWD)',                  D.E_forcing(end));
fprintf('  %-46s %8.2f  (%4.0f %%)\n', '  less upwelling LW (skin cooled)',   D.E_LWU(end),  100*D.E_LWU(end)/D.E_forcing(end));
fprintf('  %-46s %8.2f  (%4.0f %%)\n', '  more downward sensible heat',       D.E_SH(end),   100*D.E_SH(end)/D.E_forcing(end));
fprintf('  %-46s %8.2f  (%4.0f %%)\n', '  latent heat',                       D.E_LH(end),   100*D.E_LH(end)/D.E_forcing(end));
fprintf('  %-46s %8.2f  (%4.0f %%)\n', '  more conduction out of the ice, of which:', D.E_G(end), 100*D.E_G(end)/D.E_forcing(end));
fprintf('  %-46s %8.2f  (%4.0f %%)\n', '      cooling the ice (sensible)',   D.E_icecool(end), 100*D.E_icecool(end)/D.E_forcing(end));
fprintf('  %-46s %8.2f  (%4.0f %%)\n', '      making ice (latent, at the base)', D.E_growth(end), 100*D.E_growth(end)/D.E_forcing(end));
fprintf('  %-46s %8.2f\n', '      ocean flux change',                     D.E_ocean(end));
fprintf('  extra ice after %d days: %.2f cm;  Zapalac (2023) upper bound (all forcing -> freezing): %.2f cm;  efficiency %.0f %%\n', ...
    branch_days, 100*D.dh(end), 100*D.dh_zapalac(end), 100*D.dh(end)/D.dh_zapalac(end));
fprintf('  with 3x C_H: extra ice %.2f cm, SH share at 1 day %.2f, at 30 days %.2f\n', ...
    100*D3.dh(end), interp1(D3.t_h, D3.f_SH, 24), D3.f_SH(end));

% ------------------------------------------------------------------------
% 3b. Two honesty checks on the headline number
%     (i)  vertical resolution: 12 equal ice layers instead of 3
%     (ii) the air-coupling knob: tau_adv = 5 days instead of 1 (the air is
%          freer to cool with the surface -> alpha closer to 1, less SH)
% ------------------------------------------------------------------------
fprintf('\n=== 3b. sensitivity of the headline case ===\n');
checks = {{'N_layers', 12, 'layer_fractions', []}, {'tau_adv_s', 5*86400}};
names  = {'12 equal ice layers', 'tau_adv = 5 days'};
for ic = 1:numel(checks)
    Pa = seb_column_params('phase','liquid','WP_g_m2',LWP0_g_m2,'t_end_days',spin_days, checks{ic}{:});
    Ra = seb_column_run(Pa, 'verbose', false);
    Pa.num.t_end_days = branch_days;
    Rc = seb_column_run(Pa, 'init_state', Ra.final_state, 'verbose', false);
    Pb = Pa;  Pb.cloud.phase = 'ice';  Pb.cloud.WP_g_m2 = f_remain_V22*LWP0_g_m2;
    Pb = seb_column_finalize(Pb);
    Rp = seb_column_run(Pb, 'init_state', Ra.final_state, 'verbose', false);
    Dx = response_budget(Rp, Rc, rho_i);
    alpha = (Rp.T_la_K(end) - Rc.T_la_K(end)) / (Rp.T_sk_K(end) - Rc.T_sk_K(end));
    fprintf('  %-20s dLWD %+.1f | 30-day shares: LWU %.2f SH %.2f LH %.2f cond %.2f (cooling %.2f, growth %.2f) | extra ice %.2f cm | dT_sk(30d) %.1f K, alpha = dT_air/dT_sk = %.2f\n', ...
        names{ic}, Dx.dLWD(1), Dx.E_LWU(end)/Dx.E_forcing(end), Dx.E_SH(end)/Dx.E_forcing(end), Dx.E_LH(end)/Dx.E_forcing(end), ...
        Dx.E_G(end)/Dx.E_forcing(end), Dx.E_icecool(end)/Dx.E_forcing(end), Dx.E_growth(end)/Dx.E_forcing(end), 100*Dx.dh(end), Dx.dT_sk(end), alpha);
end
alpha0 = (pert.T_la_K(end) - ctrl.T_la_K(end)) / (pert.T_sk_K(end) - ctrl.T_sk_K(end));
fprintf('  %-20s (reference: default 3 layers, tau_adv = 1 day) alpha = %.2f\n', 'default', alpha0);

% ------------------------------------------------------------------------
% 4.  Thickness dependence, with and without snow
% ------------------------------------------------------------------------
fprintf('\n=== 4. extra ice after %d days vs initial thickness (glaciate + precipitate) ===\n', branch_days);
h_list = [0.5 1.0 1.5 2.5];
snow_list = [0 0.20];
dh30 = zeros(numel(snow_list), numel(h_list));  eff30 = dh30;  fG30 = dh30;  fcool30 = dh30;
Tsk_ctrl = dh30;  dLWD0 = dh30;
for js = 1:numel(snow_list)
    for jh = 1:numel(h_list)
        Pa = seb_column_params('phase','liquid','WP_g_m2',LWP0_g_m2,'t_end_days',spin_days, ...
                               'h_init_m', h_list(jh), 'h_snow_m', snow_list(js));
        Ra = seb_column_run(Pa, 'verbose', false);
        Pa.num.t_end_days = branch_days;
        Rc = seb_column_run(Pa, 'init_state', Ra.final_state, 'verbose', false);
        Pb = Pa;  Pb.cloud.phase = 'ice';  Pb.cloud.WP_g_m2 = f_remain_V22*LWP0_g_m2;
        Pb = seb_column_finalize(Pb);
        Rp = seb_column_run(Pb, 'init_state', Ra.final_state, 'verbose', false);
        Dx = response_budget(Rp, Rc, rho_i);
        dh30(js,jh) = 100*Dx.dh(end);  eff30(js,jh) = Dx.dh(end)/Dx.dh_zapalac(end);
        fG30(js,jh) = Dx.E_G(end)/Dx.E_forcing(end);  fcool30(js,jh) = Dx.E_icecool(end)/Dx.E_forcing(end);
        Tsk_ctrl(js,jh) = Rc.T_sk_K(end);  dLWD0(js,jh) = Dx.dLWD(1);
        fprintf('  h0 = %.1f m, snow = %4.2f m: control T_sk %.1f K, dLWD %+.1f, extra ice %.2f cm (%.0f %% of bound), conduction share %.2f (cooling %.2f, growth %.2f)\n', ...
            h_list(jh), snow_list(js), Tsk_ctrl(js,jh), dLWD0(js,jh), dh30(js,jh), 100*eff30(js,jh), fG30(js,jh), fcool30(js,jh), fG30(js,jh)-fcool30(js,jh));
    end
end

% ------------------------------------------------------------------------
% 5.  Figures
% ------------------------------------------------------------------------
t = ctrl.t_days;  N = P0.ice.N_layers;
col_ctrl = [0.35 0.35 0.35];  col_pert = [0.80 0.15 0.15];  col_same = [0.90 0.55 0.10];  col_none = [0.20 0.45 0.80];
lay_cols = [0.85 0.33 0.10; 0.47 0.67 0.19; 0.00 0.45 0.74];   % top / middle / bottom ice layer

% ---- Fig. 1: time series --------------------------------------------------
hf = new_figure_clean('Position', [50 50 1400 900]);
tl = tiledlayout(hf, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax,'on');
plot(ax, t, ctrl.LWD_sfc_W_m2, '--', 'Color', col_ctrl, 'LineWidth', 1.5, 'DisplayName', 'LWD control');
plot(ax, t, pert.LWD_sfc_W_m2, '-',  'Color', col_pert, 'LineWidth', 1.8, 'DisplayName', 'LWD glaciate + precipitate');
plot(ax, t, sc(2).R.LWD_sfc_W_m2, '-', 'Color', col_same, 'LineWidth', 1.2, 'DisplayName', 'LWD glaciate, same mass');
plot(ax, t, ctrl.LWU_sfc_W_m2, ':', 'Color', col_ctrl, 'LineWidth', 1.5, 'DisplayName', 'LWU control');
plot(ax, t, pert.LWU_sfc_W_m2, ':', 'Color', col_pert, 'LineWidth', 1.8, 'DisplayName', 'LWU glaciate + precipitate');
ylabel(ax, 'W m^{-2}'); title(ax, '(a) longwave at the surface'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
plot(ax, t, ctrl.T_sk_K, '--', 'Color', col_ctrl, 'LineWidth', 1.5, 'DisplayName', 'T_{skin} control');
plot(ax, t, pert.T_sk_K, '-',  'Color', col_pert, 'LineWidth', 1.8, 'DisplayName', 'T_{skin}');
plot(ax, t, pert.T_la_K, '-',  'Color', [0.1 0.6 0.3], 'LineWidth', 1.4, 'DisplayName', 'T_{sub-cloud air}');
plot(ax, t, pert.T_cl_K, '-',  'Color', [0.5 0.3 0.7], 'LineWidth', 1.4, 'DisplayName', 'T_{cloud layer}');
plot(ax, t, fixed.T_sk_K, '-.', 'Color', col_pert, 'LineWidth', 1.0, 'DisplayName', 'T_{skin}, fixed atmosphere');
ylabel(ax, 'K'); title(ax, '(b) temperatures, glaciate + precipitate'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
names = {'top (ice-air)', 'middle', 'bottom (ice-ocean)'};
for l = 1:N
    nm = sprintf('layer %d', l); if N == 3, nm = names{l}; end
    plot(ax, t, ctrl.T_ice_K(l,:), '--', 'Color', lay_cols(min(l,3),:), 'LineWidth', 1.3, 'HandleVisibility', 'off');
    plot(ax, t, pert.T_ice_K(l,:), '-',  'Color', lay_cols(min(l,3),:), 'LineWidth', 1.8, 'DisplayName', nm);
end
yline(ax, P0.ocean.T_f_K, 'k:', 'T_f', 'LineWidth', 1.0, 'HandleVisibility', 'off');
ylabel(ax, 'K'); title(ax, '(c) ice layer temperatures (dashed = control)'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
plot(ax, t, 100*(sc(2).R.h_m - ctrl.h_m), '-', 'Color', col_same, 'LineWidth', 1.6, 'DisplayName', 'glaciate, same mass');
plot(ax, t, 100*(pert.h_m - ctrl.h_m),    '-', 'Color', col_pert, 'LineWidth', 2.0, 'DisplayName', 'glaciate + precipitate');
plot(ax, t, 100*(sc(4).R.h_m - ctrl.h_m), '-', 'Color', col_none, 'LineWidth', 1.6, 'DisplayName', 'cloud dissipates');
plot(ax, t, 100*(fixed.h_m - ctrl.h_m),   '-.', 'Color', col_pert, 'LineWidth', 1.0, 'DisplayName', 'glac.+precip., fixed atmosphere');
plot(ax, t, 100*D.dh_zapalac, ':', 'Color', 'k', 'LineWidth', 1.4, 'DisplayName', 'Zapalac bound: all forcing -> freezing');
ylabel(ax, 'cm'); title(ax, '(d) extra ice thickness relative to control'); legend(ax, 'Location', 'northwest', 'FontSize', 8); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
plot(ax, t, 100*86400*ctrl.dhdt_bot_m_s, '--', 'Color', col_ctrl, 'LineWidth', 1.5, 'DisplayName', 'control');
plot(ax, t, 100*86400*pert.dhdt_bot_m_s, '-',  'Color', col_pert, 'LineWidth', 1.8, 'DisplayName', 'glaciate + precipitate');
xlabel(ax, 'days after glaciation'); ylabel(ax, 'cm day^{-1}'); title(ax, '(e) basal growth rate'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
plot(ax, t, ctrl.SH_down_W_m2, '--', 'Color', [0.30 0.45 0.69], 'LineWidth', 1.3, 'HandleVisibility', 'off');
plot(ax, t, pert.SH_down_W_m2, '-',  'Color', [0.30 0.45 0.69], 'LineWidth', 1.8, 'DisplayName', 'SH (down)');
plot(ax, t, ctrl.LH_down_W_m2, '--', 'Color', [0.33 0.66 0.41], 'LineWidth', 1.3, 'HandleVisibility', 'off');
plot(ax, t, pert.LH_down_W_m2, '-',  'Color', [0.33 0.66 0.41], 'LineWidth', 1.8, 'DisplayName', 'LH (down)');
plot(ax, t, ctrl.F_cond_top_up_W_m2, '--', 'Color', [0.51 0.45 0.70], 'LineWidth', 1.3, 'HandleVisibility', 'off');
plot(ax, t, pert.F_cond_top_up_W_m2, '-',  'Color', [0.51 0.45 0.70], 'LineWidth', 1.8, 'DisplayName', 'conduction up to skin');
plot(ax, t, pert.F_cond_bot_up_W_m2, ':',  'Color', [0.51 0.45 0.70], 'LineWidth', 1.6, 'DisplayName', 'conduction at the base');
xlabel(ax, 'days after glaciation'); ylabel(ax, 'W m^{-2}'); title(ax, '(f) surface fluxes (dashed = control)'); legend(ax, 'Location', 'east', 'FontSize', 8); grid(ax,'on');
title(tl, sprintf('Cloud glaciation over %.1f m of sea ice (no snow): LWP %g -> IWP %g g m^{-2}, \\Delta LWD = %+.0f W m^{-2}', ...
    S0.h_m, LWP0_g_m2, f_remain_V22*LWP0_g_m2, D.dLWD(1)));
exportgraphics(hf, fullfile(fig_dir, 'fig1_glaciation_timeseries.png'), 'Resolution', 150); close(hf);

% ---- Fig. 2: where the forcing goes --------------------------------------
hf = new_figure_clean('Position', [50 50 1300 520]);
tl = tiledlayout(hf, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax,'on');
tc = [0.70 0.09 0.17; 0.30 0.45 0.69; 0.33 0.66 0.41; 0.51 0.45 0.70];
plot(ax, D.t_h, D.f_LWU, '-', 'Color', tc(1,:), 'LineWidth', 2, 'DisplayName', 'upwelling LW');
plot(ax, D.t_h, D.f_SH,  '-', 'Color', tc(2,:), 'LineWidth', 2, 'DisplayName', 'sensible heat');
plot(ax, D.t_h, D.f_LH,  '-', 'Color', tc(3,:), 'LineWidth', 2, 'DisplayName', 'latent heat');
plot(ax, D.t_h, D.f_G,   '-', 'Color', tc(4,:), 'LineWidth', 2, 'DisplayName', 'conduction out of the ice');
plot(ax, Dfix.t_h, Dfix.f_SH, '--', 'Color', tc(2,:), 'LineWidth', 1.2, 'DisplayName', 'sensible, fixed atmosphere');
plot(ax, Dfix.t_h, Dfix.f_G,  '--', 'Color', tc(4,:), 'LineWidth', 1.2, 'DisplayName', 'conduction, fixed atmosphere');
set(ax, 'XScale', 'log'); xlim(ax, [1 branch_days*24]); ylim(ax, [-0.05 1]); grid(ax, 'on');
xlabel(ax, 'hours after glaciation'); ylabel(ax, 'share of \DeltaLWD');
title(ax, '(a) instantaneous partition of the forcing'); legend(ax, 'Location', 'east', 'FontSize', 8);

ax = nexttile(tl); hold(ax,'on');
E = -[D.E_forcing(end), D.E_LWU(end), D.E_SH(end), D.E_LH(end), D.E_icecool(end), D.E_growth(end), D.E_ocean(end)];
labels = {'forcing', 'LW up', 'sensible', 'latent', 'ice cooling', 'ice growth', 'ocean flux'};
bc = [0.3 0.3 0.3; tc(1,:); tc(2,:); tc(3,:); 0.6 0.55 0.8; 0.35 0.25 0.55; 0.6 0.6 0.6];
for i = 1:numel(E)
    bar(ax, i, E(i), 0.7, 'FaceColor', bc(i,:), 'EdgeColor', 'none');
    text(ax, i, E(i) + 0.4*sign(E(i)), sprintf('%.1f', E(i)), 'HorizontalAlignment', 'center', 'FontSize', 9);
end
yline(ax, 0, 'k-');
set(ax, 'XTick', 1:numel(E), 'XTickLabel', labels); ylabel(ax, 'MJ m^{-2} removed from the surface over 30 days');
title(ax, sprintf('(b) cumulative budget: %.1f cm of extra ice vs %.1f cm if ALL forcing froze water', 100*D.dh(end), 100*D.dh_zapalac(end)));
grid(ax, 'on');
title(tl, 'Where does a -\DeltaLWD forcing go?  (glaciate + precipitate minus control)');
exportgraphics(hf, fullfile(fig_dir, 'fig2_response_partition.png'), 'Resolution', 150); close(hf);

% ---- Fig. 3: thickness dependence --------------------------------------
hf = new_figure_clean('Position', [50 50 1300 480]);
tl = tiledlayout(hf, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
mk_s = {'-o', '-s'}; sn = {'no snow (ERA5-like)', '20 cm snow'};
ax = nexttile(tl); hold(ax,'on');
for js = 1:2, plot(ax, h_list, dh30(js,:), mk_s{js}, 'LineWidth', 1.8, 'MarkerFaceColor', 'auto', 'DisplayName', sn{js}); end
xlabel(ax, 'initial ice thickness [m]'); ylabel(ax, 'cm'); title(ax, sprintf('(a) extra ice after %d days', branch_days)); legend(ax, 'Location', 'northeast'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for js = 1:2, plot(ax, h_list, 100*eff30(js,:), mk_s{js}, 'LineWidth', 1.8, 'DisplayName', sn{js}); end
xlabel(ax, 'initial ice thickness [m]'); ylabel(ax, '%'); title(ax, '(b) growth efficiency: extra ice / Zapalac bound'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for js = 1:2
    plot(ax, h_list, fG30(js,:), mk_s{js}, 'Color', tc(4,:), 'LineWidth', 1.8, 'DisplayName', ['conduction share, ' sn{js}]);
    plot(ax, h_list, fcool30(js,:), mk_s{js}, 'Color', [0.6 0.55 0.8], 'LineWidth', 1.2, 'DisplayName', ['  of which ice cooling, ' sn{js}]);
end
xlabel(ax, 'initial ice thickness [m]'); ylabel(ax, 'share of cumulative forcing'); title(ax, '(c) conductive share over 30 days'); legend(ax, 'Location', 'best', 'FontSize', 8); grid(ax,'on');
title(tl, 'Thickness dependence of the response to glaciation (+ precipitation)');
exportgraphics(hf, fullfile(fig_dir, 'fig3_thickness_dependence.png'), 'Resolution', 150); close(hf);

% ---- Fig. 4: profile evolution ------------------------------------------
hf = new_figure_clean('Position', [50 50 1000 520]);
tl = tiledlayout(hf, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
times_h = [0 1 6 24 72 240 720];
cmap = copper(numel(times_h)+1);
for panel = 1:2
    ax = nexttile(tl); hold(ax,'on');
    if panel == 1, Rx = pert; ttl = '(a) glaciate + precipitate'; else, Rx = ctrl; ttl = '(b) control'; end
    for it = 1:numel(times_h)
        [~, i] = min(abs(Rx.t_days*24 - times_h(it)));
        h = Rx.h_m(i);  zc = P0.ice.centre_fractions(:) * h;   % depth of layer centres [m]
        Tz = [Rx.T_sk_K(i); Rx.T_ice_K(:,i); P0.ocean.T_f_K];
        zz = [0; zc; h];
        plot(ax, Tz, -zz, '-o', 'Color', cmap(it,:), 'LineWidth', 1.6, 'MarkerSize', 4, 'MarkerFaceColor', cmap(it,:), ...
            'DisplayName', sprintf('%g h', times_h(it)));
    end
    xlabel(ax, 'temperature [K]'); ylabel(ax, 'depth below surface [m]'); title(ax, ttl); grid(ax,'on');
    if panel == 1, legend(ax, 'Location', 'southwest', 'FontSize', 8); end
end
title(tl, 'Ice temperature profiles after glaciation: skin, layer centres, base at T_f');
exportgraphics(hf, fullfile(fig_dir, 'fig4_profile_evolution.png'), 'Resolution', 150); close(hf);

fprintf('\nfigures written to %s\n', fig_dir);

% ========================================================================
% helpers
% ========================================================================
function s = mk(name, P, overrides)
% Copy P and apply name/value overrides, recomputing the cloud emissivity.
for i = 1:2:numel(overrides)
    nm = overrides{i}; val = overrides{i+1};
    groups = fieldnames(P); done = false;
    for g = 1:numel(groups)
        if isstruct(P.(groups{g})) && isfield(P.(groups{g}), nm), P.(groups{g}).(nm) = val; done = true; end
    end
    if ~done, error('unknown parameter %s', nm); end
end
P = seb_column_finalize(P);                      % cloud and slab emissivities follow the edited phase / water path
s = struct('name', name, 'P', P, 'R', [], 'style', '-');
end

function D = response_budget(Rp, Rc, rho_i)
% Perturbed minus control, and the partition of the LWD forcing.
% Skin balance: LWD - LWU + SH_dn + LH_dn + F_top_up = 0 at every instant, so
%   dLWD = dLWU - dSH - dLH - dF_top  and the shares below sum to one.
% Ice column: F_top_up = F_ocean_up + rho L dh/dt - dE_sens/dt, so the extra
% conduction out of the ice is extra freezing + extra cooling + ocean change.
dt_s = (Rp.t_days(2) - Rp.t_days(1)) * 86400;
D.t_h   = Rp.t_days * 24;
D.dLWD  = Rp.LWD_sfc_W_m2 - Rc.LWD_sfc_W_m2;
D.dLWU  = Rp.LWU_sfc_W_m2 - Rc.LWU_sfc_W_m2;
D.dSH   = Rp.SH_down_W_m2 - Rc.SH_down_W_m2;
D.dLH   = Rp.LH_down_W_m2 - Rc.LH_down_W_m2;
D.dG    = Rp.F_cond_top_up_W_m2 - Rc.F_cond_top_up_W_m2;
D.dT_sk = Rp.T_sk_K - Rc.T_sk_K;
D.dh    = Rp.h_m - Rc.h_m;
D.f_LWU =  D.dLWU ./ D.dLWD;
D.f_SH  = -D.dSH  ./ D.dLWD;
D.f_LH  = -D.dLH  ./ D.dLWD;
D.f_G   = -D.dG   ./ D.dLWD;
% cumulative energies [MJ/m2]
cs = @(x) cumsum(x) * dt_s * 1e-6;
D.E_forcing = cs(D.dLWD);
D.E_LWU     = cs(D.dLWU);
D.E_SH      = -cs(D.dSH);
D.E_LH      = -cs(D.dLH);
D.E_G       = -cs(D.dG);
D.E_growth  = -rho_i * Rp.L_eff_J_kg * D.dh * 1e-6;                         % latent heat released by extra freezing (negative: removed from surface energy)
D.E_icecool = (Rp.E_ice_sens_J_m2 - Rc.E_ice_sens_J_m2) * 1e-6;             % change in ice sensible energy (negative: ice colder)
D.E_ocean   = -cs(Rp.F_ocean_up_W_m2 - Rc.F_ocean_up_W_m2);
% Zapalac's bound: every W m-2 of forcing freezes water, dh = |E_forcing|/(rho L)
D.dh_zapalac = -D.E_forcing * 1e6 / (rho_i * Rp.L_eff_J_kg);
end

function print_state(R, label)
fprintf('  %s:\n', label);
fprintf('    T_sk %.1f K | T_la %.1f | T_cl %.1f | T_ua %.1f | T_ml - T_f %.3f K | h %.3f m\n', ...
    R.T_sk_K(end), R.T_la_K(end), R.T_cl_K(end), R.T_ua_K(end), R.T_ml_K(end)-R.P.ocean.T_f_K, R.h_m(end));
fprintf('    LWD %.1f | LWU %.1f | SH_dn %+.1f | LH_dn %+.1f | cond up to skin %.1f | cond at base %.1f | ocean %.1f W/m2 | growth %.2f cm/day\n', ...
    R.LWD_sfc_W_m2(end), R.LWU_sfc_W_m2(end), R.SH_down_W_m2(end), R.LH_down_W_m2(end), ...
    R.F_cond_top_up_W_m2(end), R.F_cond_bot_up_W_m2(end), R.F_ocean_up_W_m2(end), 100*86400*R.dhdt_bot_m_s(end));
h_eff = R.h_m(end) + R.P.ice.h_snow_m * R.P.ice.k_W_m_K / R.P.ice.k_snow_W_m_K;   % [S25 eq. 10]
fprintf('    damping coefficients: lambda_LW %.2f, lambda_SH %.2f, lambda_LH %.2f W/m2/K; steady lambda_G = k/h_eff = %.2f (h_eff %.2f m)\n', ...
    R.lambda_LW_W_m2_K(end), R.lambda_SH_W_m2_K(end), R.lambda_LH_W_m2_K(end), R.P.ice.k_W_m_K/h_eff, h_eff);
end
