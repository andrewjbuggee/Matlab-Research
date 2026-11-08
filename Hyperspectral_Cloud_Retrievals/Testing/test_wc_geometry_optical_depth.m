%% test_wc_geometry_optical_depth.m
% Why does the synthetic-profile tau_c label need Q_ext = 2.32?
% Two libRadtran runs per cloud — the CURRENT wc_file geometry and a FIXED one —
% with libRadtran's own optical depth read back from its verbose output.
%
% By Andrew John Buggee (drafted with Claude), 8 Nov 2026
%
% =========================================================================
% PHYSICS (every equation checked against Hansen & Travis 1974, Space Sci.
% Rev. 16, 527-610; equation numbers below are theirs)
% =========================================================================
%
% (1) Extinction coefficient of a polydispersion, H&T (2.48):
%         k_ext = int pi r^2 Q_ext(r) n(r) dr                       [m^-1]
%
% (2) Effective radius, H&T (2.53):
%         r_eff = int r pi r^2 n dr / int pi r^2 n dr = <r^3>/<r^2>
%     and the size-distribution-averaged efficiency, H&T (2.58) (written there
%     for scattering; extinction is identical):
%         Qbar_ext = k_ext / G,   G = int pi r^2 n dr  (geometric cross
%                                                       section per volume)
%
% (3) Liquid water content: LWC = rho_w (4/3) pi int r^3 n dr. Dividing (1) by
%     LWC and using (2):
%         k_ext / LWC = Qbar G / (rho_w (4/3) pi <r^3> N)
%                     = 3 Qbar / (4 rho_w r_eff)
%     so
%         beta(z) = (3/4) Qbar(r_eff, alpha, lambda) LWC(z) / (rho_w r_eff(z))
%     This is EXACT for any n(r). Qbar is the area-weighted efficiency of
%     H&T (2.58), which here is correctly paired with LWC/r_eff (= (4/3) rho_w
%     G). In mixed units (LWC g m^-3, r_e um, rho_w = 1e6 g m^-3):
%         beta [km^-1] = 750 * Qbar * LWC / r_e
%     The libRadtran Mie tables store ext = beta/LWC [km^-1 / (g m^-3)], so
%     Qbar = ext * r_e / 750.
%
% (4) libRadtran's gamma distribution n(r) ~ r^alpha exp(-(alpha+3) r / r_e)
%     is H&T's standard distribution (2.56) with a = r_eff and
%     b = v_eff = 1/(alpha+3), because (1 - 3b)/b = alpha (H&T 2.57). So each
%     custom Mie table (one alpha) is H&T's standard distribution.
%
% (5) Cloud optical depth is int beta dz over the cloud. What differs below
%     is ONLY which vertical quadrature of beta libRadtran actually receives.
%
% =========================================================================
% THE TWO GEOMETRIES (levels z_1 = z_top ... z_N = z_base, N = 7 here)
% =========================================================================
%
% libRadtran >= 1.4 reads wc_file lines as LAYER quantities: a line at
% altitude z_i applies to the layer from z_i UP TO the next line above
% (libRadtran 2.0.6 user guide, wc_file: "interpreted as layer quantities").
% This script checks that by reading the per-layer tau that uvspec prints.
%
% CURRENT  (write_wc_file_from_in_situ.m, used by
%           hysics_refl_from_synthetic_NN_inputs.m):
%     lines at z_top + dz (zeros), z_1 ... z_N (level values), 0 km (zeros).
%     -> N cloudy layers. The top-level values fill [z_top, z_top + dz], an
%        extra layer ABOVE the profile, so the cloud spans z_base ... z_top+dz:
%            tau_cur = beta_1 dz + sum_{i=2..N} beta_i (z_{i-1} - z_i)
%
% FIXED    (this script, write_wc_file_mixed_layers below):
%     lines at z_top (zeros), z_2 ... z_N (one value per interval), 0 km.
%     -> N-1 cloudy layers, each filling [z_{k+1}, z_k] exactly, so the cloud
%        spans z_base ... z_top, the same as the profile and the LWP label.
%     Each layer holds the well-mixed contents of its two bounding levels:
%            LWC_k  = (LWC_k + LWC_{k+1}) / 2                  (mass is additive)
%            r_e,k  = LWC_k / [ (LWC_k/r_e,k + LWC_{k+1}/r_e,k+1) / 2 ]
%     The r_e,k line follows from H&T (2.53): mixing two equal volumes adds
%     their <r^3> (proportional to LWC) and their <r^2> (proportional to
%     LWC/r_e). With Qbar treated as constant across a layer, this layer gives
%     beta = (beta_k + beta_{k+1})/2, the trapezoid rule. LWP = trapz(LWC) and
%     the Q = 2 tau = 1500 trapz(LWC/r_e) are reproduced exactly.
%     (Why not 'interpret_as_level wc'? libRadtran 2.0.6 refuses it for 1D
%     clouds: "species_properties wc requires optical properties defined per
%     layer for 1D clouds".)
%
% For a uniform grid the two differ by the trapezoid end-point correction
%     tau_cur - tau_trapz = dz (beta_top + beta_base) / 2
%     tau_cur / tau_trapz = 1 + (beta_top + beta_base) / (2 betabar (N-1))
% with betabar = tau_trapz / (z_top - z_base). A vertically uniform cloud gives
% N/(N-1) = 7/6 = 1.167. Real profiles have small beta at cloud base, so the
% typical factor is smaller (first run, 30 clouds: median 1.07, 5-95% 1.02-1.27).
% Q_ext = 2.32 = 2 x 1.16 is therefore roughly Mie (Qbar/2 ~ 1.055) times
% geometry (~1.07-1.10), and no constant can follow the per-cloud geometry term.
%
% OUTPUTS (OUT_DIR): wc files, INP/OUT/ERR files, results .mat, and
%   fig1_showcase_cloud.png    one cloud: profiles, per-layer extinction in
%                              each geometry, cumulative tau, tau_c by method
%   fig2_all_clouds.png        label vs libRadtran tau for every cloud, which
%                              formula reproduces which geometry, and the
%                              end-point identity explaining the 7/6
%   fig3_reflectance.png       how much the reflectance changes (550, 2130 nm)
%
% Runtime: about 4 uvspec calls per cloud (~1 s each).
%%

clear; close all;

%% ---------------------------- Configuration -----------------------------

% Synthetic profiles from 06_synthetic_profile_generator.py (a numpy .npz).
% Arrays are (N_clouds, 7) with levels ordered cloud top -> base.
NPZ_FILE = ['/Users/andrewbuggee/Documents/VS_CODE/Python-Research/lasp-CU-paper-3/', ...
    'synthetic_profiles/synthetic_profiles_jointMVN_N1000_L7_2026-11-08.npz'];

N_CLOUDS  = 30;            % how many synthetic clouds to run (2 geometries each)
RNG_SEED  = 20261108;      % seed for picking which clouds; printed and saved

TAU_WAVELENGTH_NM   = 550;          % wavelength at which tau is defined [nm]
REFL_WAVELENGTHS_NM = [550, 2130];  % wavelengths for the reflectance comparison [nm]
QEXT_LABEL          = 2.32;         % constant in 06_synthetic_profile_generator.compute_tau_c

% Viewing geometry and surface: generic, fixed for every cloud
SZA_DEG        = 30;       % solar zenith angle [deg]
VZA_DEG        = 0;        % nadir view [deg]
SURFACE_ALBEDO = 0.05;     % dark ocean-like Lambertian surface [-]
ATM_FILE       = 'afglms.dat';   % AFGL mid-latitude summer
MOL_ABS_PARAM  = 'reptran coarse';
NUM_STREAMS    = 16;

% Where everything is written (beside this script)
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
OUT_DIR = fullfile(script_dir, 'wc_geometry_test_output');

%% ------------------------- Paths for this computer ----------------------

addpath(fullfile(script_dir, '..', '..', 'startup'));   % whatComputer.m
which_computer = whatComputer();
switch which_computer
    case 'andrewbuggee'
        LRT_ROOT = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/';
    case 'anbu8374'
        LRT_ROOT = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/';
    case 'curc'
        LRT_ROOT = '/projects/anbu8374/software/libRadtran-2.0.5/';
    otherwise
        error('Unknown computer "%s": set LRT_ROOT by hand.', which_computer);
end
UVSPEC_BIN = fullfile(LRT_ROOT, 'bin');
MIE_DIR    = fullfile(LRT_ROOT, 'data', 'wc', 'custom_mieTables');

% write_wc_file_from_in_situ.m (the production writer for the CURRENT geometry)
addpath(fullfile(script_dir, '..', 'LibRadTran-wrapper'));

% Sub-folders for the RT files
WC_DIR  = fullfile(OUT_DIR, 'wc_files');
INP_DIR = fullfile(OUT_DIR, 'inp_out');
FIG_DIR = fullfile(OUT_DIR, 'figures');
for d = {OUT_DIR, WC_DIR, INP_DIR, FIG_DIR}
    if ~isfolder(d{1}), mkdir(d{1}); end
end

%% ----------------------- Load the synthetic clouds ----------------------
% An .npz is a zip of .npy files. Unzip to a temporary folder and read the
% arrays with read_npy (local function at the bottom).

tmp_npz = tempname;
unzip(NPZ_FILE, tmp_npz);
re_all_um    = double(read_npy(fullfile(tmp_npz, 're.npy')));      % (N, 7) um, top -> base
lwc_all_gm3  = double(read_npy(fullfile(tmp_npz, 'lwc.npy')));     % (N, 7) g m^-3
z_all_km     = double(read_npy(fullfile(tmp_npz, 'z.npy')));       % (N, 7) km
tau_npz      = double(read_npy(fullfile(tmp_npz, 'tau_c.npy')));   % (N,) label stored by 06
if isfile(fullfile(tmp_npz, 'mean_alpha.npy'))
    alpha_all = double(read_npy(fullfile(tmp_npz, 'mean_alpha.npy')));   % (N,) per-cloud alpha
else
    a = double(read_npy(fullfile(tmp_npz, 'alpha.npy')));               % (N, 7) broadcast scalar
    alpha_all = a(:, 1);
end
rmdir(tmp_npz, 's');

n_total = size(re_all_um, 1);
n_lev   = size(re_all_um, 2);
fprintf('Loaded %d synthetic clouds with %d levels from\n  %s\n', n_total, n_lev, NPZ_FILE);

% Pick N_CLOUDS at random with a recorded seed
rng(RNG_SEED, 'twister');
cloud_idx = sort(randperm(n_total, min(N_CLOUDS, n_total)));
n_c = numel(cloud_idx);
fprintf('Running %d clouds (RNG seed %d).\n', n_c, RNG_SEED);

%% ------------------ Mie tables: ext(r_e) at 550 nm per alpha ------------
% The custom tables exist for these alpha values (same list as
% find_custom_mieTable_closest_to_alpha_profile.m). Production picks the
% table closest to the cloud's mean alpha; so does this script.
MIE_ALPHAS = [1:40, 45, 50, 55, 60, 80, 100, 125, 150];
mie_name   = @(a) sprintf('wc_mieTable_gamma_rEff_1-50microns_gammaDist_alpha_%03d.cdf', a);

mie_cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
alpha_table = zeros(n_c, 1);
for c = 1:n_c
    [~, k] = min(abs(alpha_all(cloud_idx(c)) - MIE_ALPHAS));
    alpha_table(c) = MIE_ALPHAS(k);
    if ~isKey(mie_cache, alpha_table(c))
        mie_cache(alpha_table(c)) = load_mie_ext(fullfile(MIE_DIR, mie_name(alpha_table(c))), ...
            TAU_WAVELENGTH_NM);
    end
end

%% ------------------- Analytic optical depths per cloud ------------------
% Each from the UNROUNDED profile; see the header for the equations.
tau_HT        = zeros(n_c, 1);   % Q = 2 (geometric optics limit), trapezoid
tau_label     = zeros(n_c, 1);   % Q = 2.32, trapezoid  (= compute_tau_c)
tau_mie_trapz = zeros(n_c, 1);   % Mie Qbar, trapezoid over z_base..z_top
tau_mie_cur   = zeros(n_c, 1);   % Mie Qbar, CURRENT layer geometry
tau_mie_fix   = zeros(n_c, 1);   % Mie Qbar, FIXED mixed-layer geometry
endpoint_x    = zeros(n_c, 1);   % (beta_top + beta_base) / (2 betabar)
qbar_levels   = zeros(n_c, n_lev);
dz_spread     = zeros(n_c, 1);   % max/min level spacing, to check the grid is uniform

for c = 1:n_c
    i   = cloud_idx(c);
    re  = re_all_um(i, :)';   lwc = lwc_all_gm3(i, :)';   z = z_all_km(i, :)';
    mie = mie_cache(alpha_table(c));

    % beta at every level, eq. (3) via the table: beta = ext(r_e) * LWC   [km^-1]
    beta_km = interp1(mie.reff_um, mie.ext, re, 'linear') .* lwc;
    qbar_levels(c, :) = (interp1(mie.reff_um, mie.ext, re, 'linear') .* re / 750)';

    tau_HT(c)        = 750 * 2          * abs(trapz(z, lwc ./ re));
    tau_label(c)     = 750 * QEXT_LABEL * abs(trapz(z, lwc ./ re));
    tau_mie_trapz(c) = abs(trapz(z, beta_km));

    % CURRENT: line z_i fills [z_i, z_{i-1}], with z_0 = z_top + dz
    dz_km  = -mean(diff(z));                        % as in write_wc_file_from_in_situ.m
    thick  = [dz_km; -diff(z)];                     % (N,) layer thicknesses [km]
    tau_mie_cur(c) = sum(beta_km .* thick);

    % FIXED: mixed layers, eqs. in the header
    [lwc_mix, re_mix] = mix_adjacent_levels(lwc, re);
    beta_mix = interp1(mie.reff_um, mie.ext, re_mix, 'linear') .* lwc_mix;
    tau_mie_fix(c) = sum(beta_mix .* (-diff(z)));

    betabar_km     = tau_mie_trapz(c) / (z(1) - z(end));
    endpoint_x(c)  = (beta_km(1) + beta_km(end)) / (2 * betabar_km);
    dz_spread(c)   = max(-diff(z)) / min(-diff(z));
end

fprintf('\nCheck: recomputed label vs tau_c stored in the .npz: max |rel diff| = %.2e\n', ...
    max(abs(tau_label ./ tau_npz(cloud_idx) - 1)));
fprintf('Level spacing max/min over all clouds: %.4f (1 = uniform grid)\n', max(dz_spread));
fprintf('Mie Qbar at %d nm over all levels: median %.3f, range [%.3f, %.3f]\n', ...
    TAU_WAVELENGTH_NM, median(qbar_levels(:)), min(qbar_levels(:)), max(qbar_levels(:)));

%% --------------------------- Radiative transfer -------------------------
% For each cloud and geometry:
%   * 550 nm with 'verbose' -> libRadtran's own per-layer water-cloud tau
%   * 550 and 2130 nm reflectance (quiet)

geom_names = {'current', 'fixed'};
tau_rt     = zeros(n_c, 2);            % libRadtran total water-cloud tau [cur, fix]
tau_file   = zeros(n_c, 2);            % tau predicted from the wc file as written
refl       = zeros(n_c, 2, numel(REFL_WAVELENGTHS_NM));   % (cloud, geometry, wavelength)
layers     = cell(n_c, 2);             % per-layer verbose output for each run

t_start = tic;
for c = 1:n_c
    i   = cloud_idx(c);
    re  = re_all_um(i, :)';   lwc = lwc_all_gm3(i, :)';   z = z_all_km(i, :)';
    mie_file = fullfile(MIE_DIR, mie_name(alpha_table(c)));
    mie = mie_cache(alpha_table(c));

    % ---- write the two wc files ----
    % CURRENT: the production function itself, called as in
    % hysics_refl_from_synthetic_NN_inputs.m (campaign/date/time only form the name)
    wc_cur = write_wc_file_from_in_situ(re, lwc, z, 'geomTestCurrent', '08-Nov-2026', 0, ...
        false, which_computer, i, [WC_DIR, filesep]);
    wc_cur = fullfile(WC_DIR, wc_cur);
    % FIXED: same number format as production, so only the geometry differs
    wc_fix = fullfile(WC_DIR, sprintf('WC_geomTestFixed_nn%d.DAT', i));
    write_wc_file_mixed_layers(wc_fix, re, lwc, z);
    wc_files = {wc_cur, wc_fix};

    for g = 1:2
        % tau that libRadtran SHOULD compute from this file, using the layer
        % convention and the same Mie table (checks our reading of libRadtran)
        tau_file(c, g) = tau_from_wc_file(wc_files{g}, mie);

        % verbose run at the tau wavelength -> per-layer optical depth
        base = fullfile(INP_DIR, sprintf('cloud%04d_%s_%dnm_verbose', i, geom_names{g}, TAU_WAVELENGTH_NM));
        write_inp([base '.INP'], LRT_ROOT, wc_files{g}, mie_file, TAU_WAVELENGTH_NM, ...
            SZA_DEG, VZA_DEG, SURFACE_ALBEDO, ATM_FILE, MOL_ABS_PARAM, NUM_STREAMS, true);
        run_uvspec(UVSPEC_BIN, base);
        layers{c, g}  = parse_verbose_wc_layers([base '.ERR']);
        tau_rt(c, g)  = sum(layers{c, g}.dtau);

        % reflectance at each wavelength
        for w = 1:numel(REFL_WAVELENGTHS_NM)
            base = fullfile(INP_DIR, sprintf('cloud%04d_%s_%dnm', i, geom_names{g}, REFL_WAVELENGTHS_NM(w)));
            write_inp([base '.INP'], LRT_ROOT, wc_files{g}, mie_file, REFL_WAVELENGTHS_NM(w), ...
                SZA_DEG, VZA_DEG, SURFACE_ALBEDO, ATM_FILE, MOL_ABS_PARAM, NUM_STREAMS, false);
            run_uvspec(UVSPEC_BIN, base);
            refl(c, g, w) = read_reflectance([base '.OUT']);
        end
    end
    fprintf('  cloud %4d (%2d/%d): tau_label %6.2f | libRadtran current %6.2f, fixed %6.2f\n', ...
        i, c, n_c, tau_label(c), tau_rt(c, 1), tau_rt(c, 2));
end
fprintf('RT finished in %.0f s.\n', toc(t_start));

%% ------------------------------ Summary ---------------------------------
pct = @(x) prctile(x, [50 5 95]);
fprintf('\n=== Ratios (median [5th, 95th percentile]) over %d clouds ===\n', n_c);
fprintf('Our reading of libRadtran (wc file -> tau, should be 1):\n');
fprintf('  libRadtran / predicted, current : %.4f [%.4f, %.4f]\n', pct(tau_rt(:,1) ./ tau_file(:,1)));
fprintf('  libRadtran / predicted, fixed   : %.4f [%.4f, %.4f]\n', pct(tau_rt(:,2) ./ tau_file(:,2)));
fprintf('Analytic formulas vs what libRadtran simulated:\n');
fprintf('  label (Q=2.32, trapz) / RT current : %.4f [%.4f, %.4f]\n', pct(tau_label ./ tau_rt(:,1)));
fprintf('  label (Q=2.32, trapz) / RT fixed   : %.4f [%.4f, %.4f]\n', pct(tau_label ./ tau_rt(:,2)));
fprintf('  Mie Qbar, layer sum   / RT current : %.4f [%.4f, %.4f]\n', pct(tau_mie_cur ./ tau_rt(:,1)));
fprintf('  Mie Qbar, trapz       / RT fixed   : %.4f [%.4f, %.4f]\n', pct(tau_mie_trapz ./ tau_rt(:,2)));
fprintf('Geometry effect alone:\n');
fprintf('  RT current / RT fixed              : %.4f [%.4f, %.4f]\n', pct(tau_rt(:,1) ./ tau_rt(:,2)));
fprintf('  end-point identity 1 + x/(N-1)     : %.4f [%.4f, %.4f]\n', pct(1 + endpoint_x/(n_lev-1)));
for w = 1:numel(REFL_WAVELENGTHS_NM)
    fprintf('  reflectance %4d nm, current/fixed - 1 : %+.2f%% [%+.2f%%, %+.2f%%]\n', ...
        REFL_WAVELENGTHS_NM(w), 100*(pct(refl(:,1,w) ./ refl(:,2,w)) - 1));
end

results = struct('cloud_idx', cloud_idx, 'rng_seed', RNG_SEED, 'npz_file', NPZ_FILE, ...
    'alpha_table', alpha_table, 'tau_HT', tau_HT, 'tau_label', tau_label, ...
    'tau_mie_trapz', tau_mie_trapz, 'tau_mie_cur', tau_mie_cur, 'tau_mie_fix', tau_mie_fix, ...
    'tau_rt', tau_rt, 'tau_file', tau_file, 'refl', refl, ...
    'refl_wavelengths_nm', REFL_WAVELENGTHS_NM, 'tau_wavelength_nm', TAU_WAVELENGTH_NM, ...
    'endpoint_x', endpoint_x, 'qbar_levels', qbar_levels, 'geom_names', {geom_names}, ...
    'sza_deg', SZA_DEG, 'vza_deg', VZA_DEG, 'albedo', SURFACE_ALBEDO);
results.layers = layers;
save(fullfile(OUT_DIR, 'wc_geometry_test_results.mat'), 'results');

%% ------------------------------- Figures --------------------------------
COL_CUR = [0.835 0.369 0.000];   % orange  - current geometry
COL_FIX = [0.000 0.447 0.698];   % blue    - fixed geometry
COL_ANA = [0.35 0.35 0.35];      % grey    - analytic formulas

% Showcase cloud: the one whose geometry ratio is closest to the sample median
geom_ratio = tau_rt(:,1) ./ tau_rt(:,2);
[~, cs] = min(abs(geom_ratio - median(geom_ratio)));
i   = cloud_idx(cs);
re  = re_all_um(i, :)';   lwc = lwc_all_gm3(i, :)';   z = z_all_km(i, :)';
mie = mie_cache(alpha_table(cs));
beta_km = interp1(mie.reff_um, mie.ext, re, 'linear') .* lwc;
dz_km   = -mean(diff(z));
z_top = z(1);  z_base = z(end);

% ---------------- Figure 1: one cloud, step by step ----------------
fig = new_figure('Visible', 'off', 'Position', [50 50 1500 900]);
tl  = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Synthetic cloud %d: what libRadtran is given  ' ...
    '(\\alpha = %d, %d nm)'], i, alpha_table(cs), TAU_WAVELENGTH_NM), 'FontWeight', 'bold');

% (a) r_e and (b) LWC at the 7 levels
ax = nexttile(tl);
plot(ax, re, z, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1.2); hold(ax, 'on');
mark_cloud_bounds(ax, z_top, z_base, dz_km);
xlabel(ax, 'r_e (\mum)'); ylabel(ax, 'Altitude (km)'); title(ax, '(a) Effective radius');
grid(ax, 'on');
ax = nexttile(tl);
plot(ax, lwc, z, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1.2); hold(ax, 'on');
mark_cloud_bounds(ax, z_top, z_base, dz_km);
xlabel(ax, 'LWC (g m^{-3})'); title(ax, '(b) Liquid water content'); grid(ax, 'on');

% (c) cumulative tau from cloud top down
ax = nexttile(tl); hold(ax, 'on');
for g = 1:2
    L = layers{cs, g};
    [zt, order] = sort(L.z_top_km, 'descend');
    dt = L.dtau(order);  zb = L.z_bot_km(order);
    z_edges   = [zt(1); zb];                   % km, top of first layer then each bottom
    tau_edges = [0; cumsum(dt)];
    col = COL_CUR; if g == 2, col = COL_FIX; end
    plot(ax, tau_edges, z_edges, '-', 'Color', col, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('libRadtran, %s: \\tau_c = %.2f', geom_names{g}, tau_rt(cs, g)));
end
tau_lab_cum = 750 * QEXT_LABEL * cumtrapz(-z, lwc ./ re);   % from the top down
plot(ax, tau_lab_cum, z, '--', 'Color', COL_ANA, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('label, Q = %.2f trapz: \\tau_c = %.2f', QEXT_LABEL, tau_label(cs)));
mark_cloud_bounds(ax, z_top, z_base, dz_km);
xlabel(ax, '\tau accumulated from cloud top'); title(ax, '(c) Optical depth vs altitude');
legend(ax, 'Location', 'southwest', 'FontSize', 9); grid(ax, 'on');

% (d) CURRENT and (e) FIXED: per-layer extinction from libRadtran's verbose
% output (rectangles; area = layer tau) over the level values (dots) and the
% trapezoid interpolant (dashed; its area = Mie trapezoid tau).
for g = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    L = layers{cs, g};
    col = COL_CUR; if g == 2, col = COL_FIX; end
    for k = 1:numel(L.dtau)
        b = L.dtau(k) / (L.z_top_km(k) - L.z_bot_km(k));     % km^-1
        is_extra = L.z_bot_km(k) >= z_top - 1e-6;            % layer above the profile top
        fc = col; if is_extra, fc = [0.95 0.75 0.10]; end
        patch(ax, [0 b b 0], [L.z_bot_km(k) L.z_bot_km(k) L.z_top_km(k) L.z_top_km(k)], fc, ...
            'FaceAlpha', 0.55, 'EdgeColor', col, 'LineWidth', 1);
        if is_extra
            text(ax, b*1.02, (L.z_bot_km(k)+L.z_top_km(k))/2, ...
                sprintf(' extra layer above z_{top}\n \\Delta\\tau = %.2f (%.0f%% of \\tau_c)', ...
                L.dtau(k), 100*L.dtau(k)/tau_rt(cs, g)), 'FontSize', 9, 'VerticalAlignment', 'middle');
        end
    end
    plot(ax, beta_km, z, 'k--', 'LineWidth', 1.2);
    plot(ax, beta_km, z, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    mark_cloud_bounds(ax, z_top, z_base, dz_km);
    xlim(ax, [0, 1.6*max(beta_km)]);
    xlabel(ax, '\beta_{ext} (km^{-1})');
    if g == 1
        ylabel(ax, 'Altitude (km)');
        title(ax, sprintf('(d) CURRENT: %d layers, z_{base} to z_{top}+dz', numel(L.dtau)));
    else
        title(ax, sprintf('(e) FIXED: %d mixed layers, z_{base} to z_{top}', numel(L.dtau)));
    end
    grid(ax, 'on');
end

% (f) tau_c by method
ax = nexttile(tl);
vals = [tau_HT(cs), tau_label(cs), tau_mie_trapz(cs), tau_mie_cur(cs), tau_rt(cs,1), tau_rt(cs,2)];
names = {'H&T Q=2, trapz', sprintf('label Q=%.2f, trapz', QEXT_LABEL), 'Mie Q, trapz', ...
    'Mie Q, current layers', 'libRadtran current', 'libRadtran fixed'};
cols = [COL_ANA; COL_ANA; COL_FIX*0.6+0.4; COL_CUR*0.6+0.4; COL_CUR; COL_FIX];
bh = barh(ax, vals, 'FaceColor', 'flat'); bh.CData = cols;
set(ax, 'YTick', 1:numel(vals), 'YTickLabel', names, 'YDir', 'reverse', 'FontSize', 9);
for k = 1:numel(vals)
    text(ax, vals(k), k, sprintf(' %.2f', vals(k)), 'VerticalAlignment', 'middle', 'FontSize', 9);
end
xlim(ax, [0, 1.2*max(vals)]); xlabel(ax, '\tau_c');
title(ax, '(f) \tau_c by method'); grid(ax, 'on');

exportgraphics(fig, fullfile(FIG_DIR, 'fig1_showcase_cloud.png'), 'Resolution', 200);
close(fig);

% ---------------- Figure 2: every cloud ----------------
fig = new_figure('Visible', 'off', 'Position', [50 50 1600 520]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%d synthetic clouds, \\tau at %d nm', n_c, TAU_WAVELENGTH_NM), 'FontWeight', 'bold');

% (a) training label vs libRadtran
ax = nexttile(tl); hold(ax, 'on');
lim = [0, 1.05*max([tau_label; tau_rt(:)])];
plot(ax, lim, lim, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
scatter(ax, tau_rt(:,1), tau_label, 40, COL_CUR, 'filled', 'MarkerFaceAlpha', 0.8, ...
    'DisplayName', 'vs libRadtran, current');
scatter(ax, tau_rt(:,2), tau_label, 40, COL_FIX, 'filled', 'MarkerFaceAlpha', 0.8, ...
    'DisplayName', 'vs libRadtran, fixed');
scatter(ax, tau_rt(:,2), tau_mie_trapz, 40, COL_FIX, 'LineWidth', 1.2, ...
    'DisplayName', 'Mie Q trapz label vs fixed');
axis(ax, [lim lim]); axis(ax, 'square');
xlabel(ax, '\tau_c computed by libRadtran'); ylabel(ax, '\tau_c label');
title(ax, sprintf('(a) Label (Q = %.2f) vs libRadtran', QEXT_LABEL));
legend(ax, 'Location', 'northwest', 'FontSize', 9); grid(ax, 'on');

% (b) which formula reproduces which geometry
ax = nexttile(tl); hold(ax, 'on');
R = {tau_label ./ tau_rt(:,1), tau_label ./ tau_rt(:,2), tau_mie_cur ./ tau_rt(:,1), tau_mie_trapz ./ tau_rt(:,2)};
labs = {'label / current', 'label / fixed', 'Mie layers / current', 'Mie trapz / fixed'};
rc = [COL_CUR; COL_FIX; COL_CUR; COL_FIX];
yline(ax, 1, 'k-', 'LineWidth', 1);
for k = 1:4
    swarmchart(ax, k*ones(n_c,1), R{k}, 30, rc(k,:), 'filled', 'MarkerFaceAlpha', 0.7, 'XJitterWidth', 0.5);
    text(ax, k, max([R{:}],[],'all') + 0.02, sprintf('%.3f', median(R{k})), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end
set(ax, 'XTick', 1:4, 'XTickLabel', labs, 'FontSize', 9); xlim(ax, [0.4 4.6]);
ylabel(ax, 'analytic \tau_c / libRadtran \tau_c');
title(ax, '(b) Which formula matches which geometry'); grid(ax, 'on');

% (c) end-point identity: the size of the geometry effect is set by the
% extinction at the top and base relative to the cloud mean
ax = nexttile(tl); hold(ax, 'on');
xx = linspace(0, max(endpoint_x)*1.1, 50);
plot(ax, xx, 1 + xx/(n_lev-1), 'k-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('1 + x/(N-1),  N = %d', n_lev));
scatter(ax, endpoint_x, geom_ratio, 40, COL_CUR, 'filled', 'MarkerFaceAlpha', 0.8, ...
    'DisplayName', 'libRadtran current / fixed');
plot(ax, 1, 1 + 1/(n_lev-1), 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', ...
    'DisplayName', sprintf('uniform cloud: %d/%d = %.3f', n_lev, n_lev-1, n_lev/(n_lev-1)));
xlabel(ax, 'x = (\beta_{top} + \beta_{base}) / (2 \beta_{mean})');
ylabel(ax, '\tau_{current} / \tau_{fixed}');
title(ax, '(c) Geometry effect = trapezoid end-point term'); legend(ax, 'Location', 'northwest', 'FontSize', 9);
grid(ax, 'on');

exportgraphics(fig, fullfile(FIG_DIR, 'fig2_all_clouds.png'), 'Resolution', 200);
close(fig);

% ---------------- Figure 3: reflectance consequence ----------------
fig = new_figure('Visible', 'off', 'Position', [50 50 1100 450]);
tl  = tiledlayout(fig, 1, numel(REFL_WAVELENGTHS_NM), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Reflectance change from the extra layer (SZA %g°, nadir view)', SZA_DEG), ...
    'FontWeight', 'bold');
for w = 1:numel(REFL_WAVELENGTHS_NM)
    ax = nexttile(tl); hold(ax, 'on');
    d = 100 * (refl(:,1,w) ./ refl(:,2,w) - 1);
    yline(ax, 0, 'k-');
    scatter(ax, tau_rt(:,2), d, 40, COL_CUR, 'filled', 'MarkerFaceAlpha', 0.8);
    xlabel(ax, '\tau_c, libRadtran fixed geometry');
    if w == 1, ylabel(ax, '100 (R_{current} / R_{fixed} - 1)  (%)'); end
    title(ax, sprintf('%d nm', REFL_WAVELENGTHS_NM(w))); grid(ax, 'on');
end
exportgraphics(fig, fullfile(FIG_DIR, 'fig3_reflectance.png'), 'Resolution', 200);
close(fig);

fprintf('\nFigures and results written to\n  %s\n', OUT_DIR);


%% =========================================================================
%                              Local functions
% =========================================================================

function [lwc_mix, re_mix] = mix_adjacent_levels(lwc, re)
% Well-mixed contents of each interval [z_{k+1}, z_k] (see header).
%   LWC adds as mass, (4/3) pi rho_w <r^3>.
%   LWC/r_e adds as cross section, (4/3) rho_w G, by H&T (2.53).
% So r_e of the mixture = total LWC / total (LWC/r_e).
lwc_mix = (lwc(1:end-1) + lwc(2:end)) / 2;                         % g m^-3
area    = (lwc(1:end-1)./re(1:end-1) + lwc(2:end)./re(2:end)) / 2;  % g m^-3 um^-1
re_mix  = lwc_mix ./ area;                                          % um
end

function write_wc_file_mixed_layers(file_path, re, lwc, z)
% FIXED geometry: one layer per interval, so the cloud spans exactly
% z_base ... z_top. Inputs ordered top -> base. libRadtran reads each line as
% the layer from that altitude up to the next line, so the line at z_{k+1}
% carries interval k, and a zero line at z_top closes the cloud.
[lwc_mix, re_mix] = mix_adjacent_levels(lwc, re);
z_lines   = [z(1); z(2:end)];           % z_top (closing zero line), then interval bottoms
lwc_lines = [0; lwc_mix];
re_lines  = [0; re_mix];
if z(end) > 0                           % clear layer from the surface to cloud base
    z_lines = [z_lines; 0]; lwc_lines = [lwc_lines; 0]; re_lines = [re_lines; 0];
end
fid = fopen(file_path, 'w');
fprintf(fid, '%s %10s %7s %8s \n', '#', 'z', 'LWC', 'R_eff');
fprintf(fid, '%s %10s %7s %8s \n', '#', '(km)', '(g/m^3)', '(micron)');
fprintf(fid, '%12.3f %7.4f %8.3f \n', [z_lines'; lwc_lines'; re_lines']);   % same format as production
fclose(fid);
end

function tau = tau_from_wc_file(file_path, mie)
% The tau libRadtran should compute from a wc file as WRITTEN (rounded):
% each line fills [z_i, z_{i+1}] (next line up) with beta = ext(r_e) * LWC.
d = readmatrix(file_path, 'FileType', 'text', 'CommentStyle', '#');
d = sortrows(d, 1);                     % ascending altitude
tau = 0;
for k = 1:size(d, 1) - 1
    if d(k, 2) > 0
        tau = tau + interp1(mie.reff_um, mie.ext, d(k, 3), 'linear') * d(k, 2) * (d(k+1, 1) - d(k, 1));
    end
end
end

function mie = load_mie_ext(file_path, wl_nm)
% ext(r_e) [km^-1 / (g m^-3)] from a libRadtran Mie table, interpolated
% linearly in wavelength to wl_nm. ncdump order is ext(nlam, nreff);
% ncread returns it transposed as (nreff, nlam).
ext = double(ncread(file_path, 'ext'));
wl  = double(ncread(file_path, 'wavelen'));
if max(wl) < 100, wl = wl * 1e3; end    % table stores um -> nm
if size(ext, 1) ~= numel(double(ncread(file_path, 'reff'))), ext = ext'; end
mie.reff_um = double(ncread(file_path, 'reff'));
mie.ext     = interp1(wl, ext', wl_nm, 'linear')';      % (nreff, 1)
end

function write_inp(inp_path, lrt_root, wc_file, mie_file, wl_nm, sza, vza, albedo, ...
    atm_file, mol_abs, n_streams, verbose)
% Minimal uvspec input: one wavelength, DISORT, TOA nadir reflectance.
fid = fopen(inp_path, 'w');
fprintf(fid, 'data_files_path %s\n', fullfile(lrt_root, 'data', filesep));
fprintf(fid, 'atmosphere_file %s\n', fullfile(lrt_root, 'data', 'atmmod', atm_file));
fprintf(fid, 'source solar %s\n', fullfile(lrt_root, 'data', 'solar_flux', 'kurudz_1.0nm.dat'));
fprintf(fid, 'mol_abs_param %s\n', mol_abs);
fprintf(fid, 'wavelength %g %g\n', wl_nm, wl_nm);
fprintf(fid, 'rte_solver disort\n');
fprintf(fid, 'number_of_streams %d\n', n_streams);
fprintf(fid, 'wc_file 1D %s\n', wc_file);
fprintf(fid, 'wc_properties %s interpolate\n', mie_file);
fprintf(fid, 'sza %g\nphi0 0\numu %.4f\nphi 0\n', sza, cosd(vza));
fprintf(fid, 'zout toa\nalbedo %g\n', albedo);
fprintf(fid, 'output_quantity reflectivity\n');
if verbose, fprintf(fid, 'verbose\n'); else, fprintf(fid, 'quiet\n'); end
fclose(fid);
end

function run_uvspec(bin_dir, base)
% uvspec must run from its bin folder; stdout -> .OUT, stderr (verbose) -> .ERR
cmd = sprintf('cd "%s" && ./uvspec < "%s.INP" > "%s.OUT" 2> "%s.ERR"', bin_dir, base, base, base);
status = system(cmd);
if status ~= 0
    error('uvspec failed (%d) for %s.INP:\n%s', status, base, fileread([base '.ERR']));
end
end

function R = read_reflectance(out_path)
% With output_quantity reflectivity and one umu/phi, the third line is
% "umu  u0u  uu(phi)"; the last number is the reflectance in the view direction.
lines = splitlines(strtrim(fileread(out_path)));
v = sscanf(lines{3}, '%f');
R = v(end);
end

function L = parse_verbose_wc_layers(err_path)
% Read the water-cloud columns of uvspec's optical_properties() table.
% Each row is one model layer; its z column is the layer BOTTOM and the row
% above gives the layer top (checked against Rayleigh dtau, which scales with
% layer thickness). Water cloud dtau = scattering + absorption columns.
txt   = fileread(err_path);
k0    = strfind(txt, '*** optical_properties()');
if isempty(k0), error('No optical_properties table in %s', err_path); end
lines = splitlines(txt(k0(1):end));
z = []; dsca = []; dabs = [];
for k = 1:numel(lines)
    tok = regexp(lines{k}, '^\s*(\d+)\s*\|\s*([-\d.]+)\s*\|(.*)$', 'tokens', 'once');
    if isempty(tok)
        if ~isempty(z) && startsWith(strtrim(lines{k}), 'sum'), break; end
        continue
    end
    f  = strsplit(tok{3}, '|');           % {rayleigh, aerosol, water cloud, ice cloud, mol}
    wc = sscanf(f{3}, '%f');              % [scatter, absorption, asymmetry]
    z(end+1, 1)    = str2double(tok{2});  %#ok<AGROW>
    dsca(end+1, 1) = wc(1);               %#ok<AGROW>
    dabs(end+1, 1) = wc(2);               %#ok<AGROW>
end
z_top = [NaN; z(1:end-1)];                % layer top = previous row's z
keep  = (dsca + dabs) > 0;
L.z_bot_km = z(keep);
L.z_top_km = z_top(keep);
L.dtau     = dsca(keep) + dabs(keep);
end

function mark_cloud_bounds(ax, z_top, z_base, dz)
% Dashed: profile top and base. Dotted: z_top + dz, where the CURRENT
% geometry's cloud actually ends.
yline(ax, z_top,      'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax, z_base,     'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax, z_top + dz, 'k:',  'LineWidth', 1.2, 'HandleVisibility', 'off');
ylim(ax, [z_base - 0.6*dz, z_top + 1.6*dz]);
end

function A = read_npy(file_path)
% Minimal reader for numpy .npy (format v1/v2, little-endian, C or Fortran order).
fid = fopen(file_path, 'r', 'l');
cleaner = onCleanup(@() fclose(fid));
magic = fread(fid, 6, 'uint8=>char')';
if ~strcmp(magic(2:end), 'NUMPY'), error('%s is not a .npy file', file_path); end
ver = fread(fid, 2, 'uint8');
if ver(1) == 1, hlen = fread(fid, 1, 'uint16'); else, hlen = fread(fid, 1, 'uint32'); end
hdr = fread(fid, hlen, 'uint8=>char')';
descr = regexp(hdr, '''descr'':\s*''([^'']+)''', 'tokens', 'once');  descr = descr{1};
shp   = regexp(hdr, '''shape'':\s*\(([^\)]*)\)', 'tokens', 'once');  shp = shp{1};
shape = str2double(strsplit(strtrim(shp), ','));  shape = shape(~isnan(shape));
if isempty(shape), shape = 1; end
fortran = contains(hdr, '''fortran_order'': True');
types = struct('f4', 'single', 'f8', 'double', 'i4', 'int32', 'i8', 'int64');
if descr(1) ~= '<' && descr(1) ~= '|', error('Only little-endian .npy supported'); end
prec = types.(descr(2:end));
data = fread(fid, prod(shape), ['*' prec]);
if isscalar(shape)
    A = data(:);
elseif fortran
    A = reshape(data, shape);
else
    A = permute(reshape(data, fliplr(shape)), numel(shape):-1:1);
end
end

function fig = new_figure(varargin)
% A figure with MATLAB's FACTORY graphics defaults, regardless of the dark
% root defaults in startup.m (same recipe as EPCAPE/supersaturation/functions/new_figure.m).
fig = figure(varargin{:});
defs = fieldnames(get(groot, 'default'));
for k = 1:numel(defs)
    if startsWith(defs{k}, 'defaultFigure', 'IgnoreCase', true), continue; end
    set(fig, defs{k}, get(groot, regexprep(defs{k}, '^default', 'factory', 'ignorecase')));
end
fig.Color = 'w';
end
