% TEST_COLUMN_MODEL  Validation of the column model against things we know exactly.
%
%   >> test_column_model
%
% Five checks, each printing PASS/FAIL (or a comparison table):
%
%   T1  Longwave grey-slab solver conserves energy exactly.
%   T2  Ice growth reproduces the ANALYTIC Stefan problem [CJ59 ch. XI]:
%       with the skin held at T_sk and no ocean heat flux, the thickness obeys
%          quasi-steady:  h^2 = h0^2 + 2 k (T_f - T_sk) t / (rho L)      (small Stefan number)
%          exact Neumann: h = 2 lambda sqrt(kappa t),  lambda e^{lambda^2} erf(lambda) = St/sqrt(pi)
%       with St = c (T_f - T_sk)/L the Stefan number.  The quasi-steady law is the
%       classical "ice grows as sqrt(time)" result (Stefan 1891; Maykut 1986;
%       Leppäranta 1993) and is why the growth RATE falls as 1/h -- Matt Shupe's
%       "growth is fastest for the thinnest ice", OV meeting 22 Oct 2026.
%       Run with 3 layers (the default) and 12 to show the discretisation error.
%   T3  The coupled run conserves energy (residual relative to the flux scale).
%   T4  The two winter Arctic states: clear-sky vs cloudy equilibrium skin
%       temperature and LWD (Stramler et al. 2011; [S25] Fig. 1 clusters).
%   T5  The instantaneous response dT_sk/dLWD against the linearised four-term
%       skin balance of [LSM], with lambda_G the frequency-dependent ice admittance.
%
% Figures are written next to this file as PNGs (invisible figures, exported).

clear; clc;
here = fileparts(mfilename('fullpath'));  if isempty(here), here = pwd; end
addpath(here);
fig_dir = fullfile(here, 'figures');  if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
fprintf('\n==================== test_column_model ====================\n');

% ========================================================================
% T1  grey-slab longwave conservation
% ========================================================================
P = seb_column_params();
rng(42);                                              % seeded, reproducible
worst = 0;
for trial = 1:200
    T = 200 + 80*rand(3,1);  e = rand(3,1);  Ts = 200 + 80*rand;  es = 0.9 + 0.1*rand;
    LW = lw_grey_slabs(T, e, Ts, es, P.const.sigma_SB_W_m2_K4);
    worst = max(worst, abs(LW.residual));
end
fprintf('T1  LW conservation: max |Q_sfc + sum(Q_slab) + OLR| = %.2e W/m2  -> %s\n', ...
    worst, pass_str(worst < 1e-9));

% ========================================================================
% T2  analytic Stefan problem
% ========================================================================
dT_K   = 20;                               % T_f - T_sk
h0_m   = 0.05;                             % start thin so h0 matters little
days   = 30;
rows = {};
hfig = new_figure_clean('Position', [100 100 1100 420]);
tl = tiledlayout(hfig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax1 = nexttile(tl); hold(ax1, 'on'); ax2 = nexttile(tl); hold(ax2, 'on');
lay = {struct('N', 3, 'fr', [0.1 0.8 0.1], 'lab', '3 layers [0.1 0.8 0.1]'), ...
       struct('N', 3, 'fr', [],            'lab', '3 equal layers'), ...
       struct('N', 12, 'fr', [],           'lab', '12 equal layers')};
cols = lines(numel(lay));
for iN = 1:numel(lay)
    N = lay{iN}.N;
    P = seb_column_params('N_layers', N, 'layer_fractions', lay{iN}.fr, 'h_init_m', h0_m, ...
                          't_end_days', days, 'F_ocean_fixed_W_m2', 0, 'mode', 'prescribed', 'dt_s', 120);
    P.ice.T_sk_prescribed_K = P.ocean.T_f_K - dT_K;
    R = seb_column_run(P, 'verbose', false);
    t_s = R.t_days*86400;
    k = P.ice.k_W_m_K;  rho = P.ice.rho_kg_m3;  L = P.const.L_fus_J_kg;
    cI = P.ice.c_J_kg_K;  kappa = P.ice.kappa_m2_s;
    % quasi-steady Stefan law (linear profile, no heat capacity)
    h_qs = sqrt(h0_m^2 + 2*k*dT_K*t_s/(rho*L));
    % exact Neumann (one-phase, h0 = 0): solve lambda e^{lambda^2} erf(lambda) = St/sqrt(pi)
    St  = cI*dT_K/L;
    lam = fzero(@(x) x.*exp(x.^2).*erf(x) - St/sqrt(pi), [1e-6 3]);
    % offset the Neumann solution so it also starts at h0 (shift the time origin)
    t0  = (h0_m/(2*lam))^2/kappa;
    h_nm = 2*lam*sqrt(kappa*(t_s + t0));
    err_qs = 100*(R.h_m(end) - h_qs(end))/h_qs(end);
    err_nm = 100*(R.h_m(end) - h_nm(end))/h_nm(end);
    rows(end+1,:) = {lay{iN}.lab, R.h_m(end), h_qs(end), h_nm(end), err_qs, err_nm}; %#ok<SAGROW>
    plot(ax1, R.t_days, 100*R.h_m, '-', 'Color', cols(iN,:), 'LineWidth', 1.8, ...
        'DisplayName', ['model, ' lay{iN}.lab]);
    plot(ax2, R.t_days, 100*(R.h_m - h_nm), '-', 'Color', cols(iN,:), 'LineWidth', 1.8, ...
        'DisplayName', lay{iN}.lab);
end
plot(ax1, R.t_days, 100*h_qs, 'k--', 'LineWidth', 1.2, 'DisplayName', 'quasi-steady Stefan law');
plot(ax1, R.t_days, 100*h_nm, 'k:',  'LineWidth', 1.6, 'DisplayName', 'exact Neumann solution');
xlabel(ax1, 'time [days]'); ylabel(ax1, 'ice thickness [cm]'); grid(ax1, 'on');
title(ax1, sprintf('Stefan problem: T_{sk} = T_f - %g K, F_w = 0', dT_K));
legend(ax1, 'Location', 'southeast'); 
xlabel(ax2, 'time [days]'); ylabel(ax2, 'model - Neumann [cm]'); grid(ax2, 'on');
title(ax2, 'discretisation error'); legend(ax2, 'Location', 'northeast');
exportgraphics(hfig, fullfile(fig_dir, 'test_T2_stefan_analytic.png'), 'Resolution', 150);
close(hfig);
fprintf('T2  Stefan problem after %d days (Stefan number St = %.3f):\n', days, St);
fprintf('    %-26s %12s %14s %14s %10s %10s\n', 'layers', 'model h [m]', 'quasi-steady', 'Neumann', 'err_qs %', 'err_Nm %');
for i = 1:size(rows,1)
    fprintf('    %-26s %12.4f %14.4f %14.4f %10.2f %10.2f\n', rows{i,:});
end
fprintf('    -> %s (every layering within 1 %% of the exact Neumann solution)\n', pass_str(all(abs([rows{:,6}]) < 1)));
fprintf('    growth rate at the end: %.2f cm/day for h = %.2f m; the law says rate = k dT/(rho L h)\n', ...
    100*86400*R.dhdt_bot_m_s(end), R.h_m(end));
fprintf('    = %.2f cm/day\n', 100*86400*k*dT_K/(rho*L*R.h_m(end)));

% ========================================================================
% T3  energy conservation of the coupled run
% ========================================================================
P = seb_column_params('t_end_days', 10);
R = seb_column_run(P, 'verbose', false);
flux_scale = mean(abs(R.OLR_W_m2));
fprintf('T3  coupled energy residual: mean %.2e, max %.2e W/m2 (OLR ~ %.0f W/m2) -> %s\n', ...
    R.energy.residual_mean_W_m2, R.energy.residual_max_W_m2, flux_scale, ...
    pass_str(R.energy.residual_max_W_m2 < 1e-3*flux_scale));
fprintf('    max Newton iterations per step: %d\n', max(R.newton_iter));

% ========================================================================
% T4  the two winter states
% ========================================================================
P_cl = seb_column_params('t_end_days', 40, 'phase', 'none');
P_ov = seb_column_params('t_end_days', 40, 'phase', 'liquid', 'WP_g_m2', 30);
R_cl = seb_column_run(P_cl, 'verbose', false);
R_ov = seb_column_run(P_ov, 'verbose', false);
fprintf('T4  equilibrium winter states (day 40):\n');
fprintf('    %-14s %8s %8s %8s %8s %8s %8s %8s %9s\n', 'state', 'T_sk', 'T_la', 'T_cl', 'LWD', 'LWU', 'SH_dn', 'F_top', 'dh/dt cm/d');
for R0 = [R_cl R_ov]
    fprintf('    %-14s %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %9.2f\n', ...
        R0.P.cloud.phase, R0.T_sk_K(end), R0.T_la_K(end), R0.T_cl_K(end), R0.LWD_sfc_W_m2(end), ...
        R0.LWU_sfc_W_m2(end), R0.SH_down_W_m2(end), R0.F_cond_top_up_W_m2(end), 100*86400*R0.dhdt_bot_m_s(end));
end
fprintf('    cloudy - clear: dT_sk = %.1f K, dLWD = %.1f W/m2\n', ...
    R_ov.T_sk_K(end) - R_cl.T_sk_K(end), R_ov.LWD_sfc_W_m2(end) - R_cl.LWD_sfc_W_m2(end));
fprintf('    (MOSAiC winter: clear/cloudy skin ~ -30/-20 C, LWD ~ 160/220 W/m2 [S25 Fig. 1-2])\n');

% ========================================================================
% T5  step response vs the linearised skin balance of [LSM]
% ========================================================================
% Spin up the cloudy state with a PRESCRIBED atmosphere (alpha = 0), add +10 W/m2
% to the skin (F_extra) and compare dT_sk(t) with
%   dT_sk = dF / (lambda_LW + lambda_SH + lambda_LH + lambda_G(P)),
% lambda_G(P) = Re[k q coth(q h)], q = sqrt(i 2 pi rho c /(k P)): the in-phase
% admittance of the ice slab at forcing period P [LSM; CJ59 sec. 2.6].
% For a semi-infinite conductor the surface temperature after a STEP in flux is
% dT = (2 dF/k) sqrt(kappa t/pi) [CJ59 sec. 2.9], which equals dF / Re Y(P) at
% P = 4t exactly; so the model at time t should match the P = 4t column once
% the top layer is thin enough to resolve the penetration depth sqrt(kappa t/pi).
% Equal 3-layer splitting of 1.5 m cannot (top layer 50 cm, diffusion time
% 2.7 days); the default [0.1 0.8 0.1] split (15 cm, 6 h) does from ~6 h on;
% 30 equal layers (5 cm) from ~1 h on.  This is the limitation to remember
% when reading hourly skin temperatures from a 3-layer model.
lay5 = {struct('N', 3, 'fr', [],           'lab', '3 equal'), ...
        struct('N', 3, 'fr', [0.1 0.8 0.1], 'lab', '3 [0.1 0.8 0.1]'), ...
        struct('N', 12, 'fr', [],          'lab', '12 equal'), ...
        struct('N', 30, 'fr', [],          'lab', '30 equal')};
dF = 10;  th_list = [1 3 6 12 24 48 72 168];
tab = zeros(numel(th_list), numel(lay5));
for il = 1:numel(lay5)
    P0 = seb_column_params('t_end_days', 30, 'mode', 'prescribed', 'WP_g_m2', 30, ...
                           'N_layers', lay5{il}.N, 'layer_fractions', lay5{il}.fr);
    R0 = seb_column_run(P0, 'verbose', false);
    P1 = P0;  P1.num.t_end_days = 8;
    Rc = seb_column_run(P1, 'init_state', R0.final_state, 'verbose', false);
    P2 = P1;  P2.atm.F_extra_down_W_m2 = dF;
    Rp = seb_column_run(P2, 'init_state', R0.final_state, 'verbose', false);
    for it = 1:numel(th_list)
        [~, i] = min(abs(Rc.t_days*24 - th_list(it)));
        tab(it, il) = Rp.T_sk_K(i) - Rc.T_sk_K(i);
    end
    if il == 1
        lamLW = Rc.lambda_LW_W_m2_K(1); lamSH = Rc.lambda_SH_W_m2_K(1); lamLH = Rc.lambda_LH_W_m2_K(1);
        k = P0.ice.k_W_m_K; rc = P0.ice.rho_c_J_m3_K; h = Rc.h_m(1);
    end
end
lamG = @(Pd_s) real( k*sqrt(1i*2*pi./Pd_s*rc/k) .* coth(sqrt(1i*2*pi./Pd_s*rc/k)*h) );
fprintf('T5  step dF = +%g W/m2 at the skin, prescribed atmosphere; lambda_LW = %.2f, lambda_SH = %.2f, lambda_LH = %.2f W/m2/K; h = %.2f m\n', ...
    dF, lamLW, lamSH, lamLH, h);
labs = cellfun(@(x) x.lab, lay5, 'UniformOutput', false);
fprintf('    dT_sk [K]:   %-8s', 't [h]'); fprintf('%16s', labs{:}); fprintf('%18s %18s\n', 'linear, P = 4t', 'linear, P = 2t');
for it = 1:numel(th_list)
    th = th_list(it);
    fprintf('                 %-8d', th); fprintf('%16.3f', tab(it,:));
    fprintf('%18.3f %18.3f\n', dF/(lamLW+lamSH+lamLH+lamG(4*th*3600)), dF/(lamLW+lamSH+lamLH+lamG(2*th*3600)));
end
fprintf('    steady state (lambda_G = k/h = %.2f): dT = %.3f K\n', k/h, dF/(lamLW+lamSH+lamLH+k/h));
fprintf('    -> the 30-layer model should track the P = 4t column; 3 equal layers over-respond for the first ~2 days\n');

fprintf('\n==================== done ====================\n');

function s = pass_str(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
