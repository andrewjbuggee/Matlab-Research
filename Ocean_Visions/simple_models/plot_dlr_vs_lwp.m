% PLOT_DLR_VS_LWP  Why does the LWP at which DLR saturates move from month to month?
%
%   >> plot_dlr_vs_lwp
%
% Jenny's result (OV meeting, 22 Oct 2026): at Utqiagvik the downwelling
% longwave saturates in LWP, as Gillespie had shown, but BOTH the asymptote
% and the LWP where it is reached change with month -- near 10 g m-2 in the cold
% months, near 30 g m-2 in September-October.
%
% In the grey-slab model the surface DLR is
%
%     DLR = eps_la sigma T_la^4 + (1 - eps_la) [ eps_c sigma T_c^4 + (1 - eps_c) D_above ]
%
% with eps_c = 1 - exp(-k_abs LWP).  Two separate things follow:
%
%  1. The LWP at which DLR reaches a given fraction of its asymptote depends
%     ONLY on k_abs:   LWP_95 = -ln(0.05)/k_abs = 3/k_abs.   With Stephens' 0.130
%     m2 g-1 that is 23 g m-2; it is independent of the cloud temperature, the
%     water vapour below, or the season -- UNLESS k_abs changes.  k_abs scales as
%     1/r_e (geometric optics, see cloud_lw_emissivity), so larger droplets push
%     the saturation to larger LWP: r_e = 15 um gives LWP_95 = 35 g m-2.  A
%     seasonal change in droplet size (or in the liquid/ice split of the
%     retrieved water path) is the physical candidate for the month-to-month
%     shift.
%  2. The AMPLITUDE of the increase from clear to black,
%        (1 - eps_la) (sigma T_c^4 - D_above),
%     does change with season: a moister, more emissive sub-cloud layer (larger
%     eps_la) and a cloud closer to the clear-sky effective temperature both
%     shrink it.  A smaller amplitude makes the curve LOOK saturated earlier
%     when judged by eye or against a noise floor, without any change in LWP_95.
%
% Panel (a) varies r_e at fixed atmosphere; panel (b) varies the sub-cloud
% emissivity and cloud temperature at fixed r_e.  The markers show LWP_95.

clear; clc;
here = fileparts(mfilename('fullpath'));  if isempty(here), here = pwd; end
addpath(here);
fig_dir = fullfile(here, 'figures');  if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
P = seb_column_params();
sig = P.const.sigma_SB_W_m2_K4;
LWP = linspace(0, 80, 400);
T_sk = 248;  T_la = 251;  T_cl = 250;  T_ua = 238;  e_ua = P.atm.emis_ua;

% the cloud slab also carries its clear-air gases: eps_slab = 1 - (1-eps_gas)(1-eps_c)
e_gas = P.atm.emis_gas_cl;
dlr = @(eps_c, e_la, Tcl) arrayfun(@(e) getfield(lw_grey_slabs([T_la; Tcl; T_ua], [e_la; 1-(1-e_gas)*(1-e); e_ua], T_sk, P.ice.emis_sfc, sig), 'LWD_sfc'), eps_c); %#ok<GFLD>

hf = new_figure_clean('Position', [80 80 1200 480]);
tl = tiledlayout(hf, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) droplet size
ax = nexttile(tl); hold(ax, 'on');
re_list = [6 10 15 20];  cols = parula(numel(re_list)+1);
for i = 1:numel(re_list)
    [eps_c, k] = cloud_lw_emissivity(LWP, 'liquid', P, re_list(i));
    y = dlr(eps_c, P.atm.emis_la, T_cl);
    LWP95 = 3/k;
    plot(ax, LWP, y, '-', 'Color', cols(i,:), 'LineWidth', 2, 'DisplayName', sprintf('r_e = %g um, k_{abs} = %.3f m^2 g^{-1}, LWP_{95} = %.0f g m^{-2}', re_list(i), k, LWP95));
    plot(ax, LWP95, interp1(LWP, y, LWP95), 'o', 'Color', cols(i,:), 'MarkerFaceColor', cols(i,:), 'HandleVisibility', 'off');
end
xlabel(ax, 'liquid water path [g m^{-2}]'); ylabel(ax, 'surface DLR [W m^{-2}]'); grid(ax, 'on');
title(ax, '(a) droplet size moves the saturation LWP (3/k_{abs})'); legend(ax, 'Location', 'southeast', 'FontSize', 8);

% (b) sub-cloud emissivity and cloud temperature
ax = nexttile(tl); hold(ax, 'on');
e0 = P.atm.emis_la;                                  % the default sub-cloud emissivity (0.40), not a hard-coded 0.45
cases = {0.30, 250, 'dry sub-cloud air (\epsilon_{la} = 0.30), T_c = 250 K'; ...
         e0,   250, sprintf('default (\\epsilon_{la} = %.2f), T_c = 250 K', e0); ...
         0.65, 250, 'moist sub-cloud air (\epsilon_{la} = 0.65), T_c = 250 K'; ...
         e0,   258, 'default, warmer cloud T_c = 258 K'};
cols = lines(size(cases,1));
[eps_c, k] = cloud_lw_emissivity(LWP, 'liquid', P);
for i = 1:size(cases,1)
    y = dlr(eps_c, cases{i,1}, cases{i,2});
    plot(ax, LWP, y, '-', 'Color', cols(i,:), 'LineWidth', 2, 'DisplayName', sprintf('%s: clear %.0f -> black %.0f W m^{-2}', cases{i,3}, y(1), y(end)));
    plot(ax, 3/k, interp1(LWP, y, 3/k), 'o', 'Color', cols(i,:), 'MarkerFaceColor', cols(i,:), 'HandleVisibility', 'off');
end
xlabel(ax, 'liquid water path [g m^{-2}]'); ylabel(ax, 'surface DLR [W m^{-2}]'); grid(ax, 'on');
title(ax, '(b) atmosphere changes the amplitude, not LWP_{95} (markers)'); legend(ax, 'Location', 'southeast', 'FontSize', 8);
title(tl, 'Surface DLR vs cloud liquid water path in the grey-slab model (T_{skin} = 248 K)');
exportgraphics(hf, fullfile(fig_dir, 'dlr_vs_lwp.png'), 'Resolution', 150); close(hf);
fprintf('LWP_95 = 3/k_abs: r_e = 10 um -> %.0f g/m2; 15 um -> %.0f; 20 um -> %.0f; 6 um -> %.0f\n', ...
    3./cloud_lw_emissivity_k(P, [10 15 20 6]));
fprintf('figure written to %s\n', fullfile(fig_dir, 'dlr_vs_lwp.png'));

function k = cloud_lw_emissivity_k(P, re)
k = zeros(size(re));
for i = 1:numel(re), [~, k(i)] = cloud_lw_emissivity(1, 'liquid', P, re(i)); end
end
