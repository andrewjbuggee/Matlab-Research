%% Isolate the difference between the in-situ and the r_e/LWC optical thickness
%
% DESCRIPTION:
%   The optical thickness stored on the ORACLES ensemble profiles (computed in
%   find_verticalProfiles_ORACLES.m, line ~492) and the one recomputed from the
%   measured r_e and LWC profiles (06_synthetic_profile_generator.py,
%   compute_tau_c) disagree by ~24% in the median. This script computes the
%   optical thickness five ways for the same cloud and pulls that disagreement
%   apart into its two independent causes.
%
% THE REFERENCE CALCULATION
%   The most accurate optical thickness available from the in-situ data makes no
%   assumption about the shape of the droplet size distribution and no
%   substitution for its moments. It is Hansen & Travis (1974) Equation (2.48),
%   integrated over height:
%
%       tau = int_z [ int_r  pi * r^2 * Q_ext(r, lambda) * n(r, z)  dr ] dz
%
%   All equation numbers below refer to Hansen, J. E. and Travis, L. D. (1974),
%   'Light Scattering in Planetary Atmospheres', Space Sci. Rev. 16, 527-610.
%
%   The ORACLES probes report n(r, z) directly as Nd: droplet concentration per
%   bin (cm^-3), on 172 radius bins spanning 1.6 um to 5.25 mm. Because Nd is
%   already a per-bin concentration rather than a density dN/dr, the inner
%   integral is a plain sum over bins with no bin-width factor:
%
%       beta(z) = pi * sum_i [ Q_ext(r_i) * r_i^2 * Nd_i(z) ]
%
%   Q_ext(r) comes from one of two sources, set by QEXT_SOURCE:
%     'miev0'  - a single monodisperse MIEV0 run on a fine radius grid at the
%                wavelength of interest (a few seconds, serves every profile)
%     'lookup' - the saved table
%                Q_ext_1nm_sampling_monodispersed_radii_1_to_500_microns_reduced.txt
%                (340 radii x 2901 wavelengths, 1 um radius steps)
%
%   Either way Q is AVERAGED WITHIN EACH MEASURED BIN, weighted by r^2, rather
%   than point-sampled at the bin center. Monodisperse Q_ext oscillates with a
%   period of about 0.27 um in radius near r = 2 um, which is finer than the
%   narrowest bin, so a point sample there is arbitrary in phase: MIEV0 gives
%   2.19 at r = 6.25 um where the lookup table gives 2.05 at r = 6 um. Both are
%   correct; neither is representative on its own. Averaging over the bin
%   removes the ambiguity, and the script reports both sources so the residual
%   sensitivity is visible rather than assumed.
%
%   Bins above the Q grid take the geometric-optics limit Q = 2; they carry less
%   than 0.01% of the area integral in every profile checked.
%
% THE TWO APPROXIMATE FORMULAS
%   find_verticalProfiles_ORACLES.m:
%       tau = pi * integral( Q_ext_avg * r_e^2 * N_c  dz )
%
%   06_synthetic_profile_generator.py:
%       tau = (3/4) * integral( Q * LWC / (rho_w * r_e) dz )
%           = 750 * Q * integral( LWC[g/m^3] / r_e[um] dz[km] )
%
%   Substituting LWC = (4/3)*pi*rho_w*N_c*<r^3> and r_e = <r^3>/<r^2>, the
%   second reduces to pi * integral( Q * N_c * <r^2> dz ), so the two differ by
%
%       tau_insitu / tau_reLWC = (Q_ext_avg / Q) * (r_e^2 / <r^2>)
%                                \___________/   \_____________/
%                                 extinction      distribution
%                                 efficiency          width
%
%   The width factor is the one that matters. H&T (2.58) DEFINES the
%   distribution-averaged efficiency factor as the coefficient divided by the
%   geometric cross-sectional area per unit volume, G = int(pi r^2 n(r) dr):
%
%       Q_ext_avg = k_ext / G = int( pi r^2 Q_ext n dr ) / int( pi r^2 n dr )
%
%   and average_mie_over_size_distribution computes exactly that (line 297). So
%   k_ext = Q_ext_avg * G = pi * Q_ext_avg * N * <r^2>: the averaged efficiency
%   must be paired with the second moment, not with r_e^2.
%
%   That the substitution biases tau HIGH is guaranteed, not merely observed.
%   By H&T (2.53), r_e = <r^3>/<r^2>, so
%
%       <r^2> <= r_e^2   <=>   <r^2>^3 <= <r^3>^2
%
%   which is Lyapunov's moment inequality, true for every size distribution,
%   with equality only for a monodisperse one. For the standard distribution
%   (2.56)-(2.57) the factor takes the closed form
%
%       <r^2> / r_e^2 = (1 - v_eff) * (1 - 2*v_eff)
%
%   but that form is specific to the gamma family, and it is only as accurate as
%   the v_eff put into it: with v_eff measured from the distribution via (2.54)
%   it is good to ~2%, while with v_eff = 1/(alpha+3) from the fitted gamma it is
%   off by ~6%, because the ORACLES distributions carry a separate drizzle mode.
%   Nothing in this script relies on either closed form; the measured <r^2> is
%   used throughout, and the predictions are printed only for comparison.
%
% WHY THE r_e/LWC FORM IS EXACT, NOT AN APPROXIMATION
%   With LWC = (4/3)*pi*rho_w*N*<r^3> and r_e = <r^3>/<r^2>,
%
%       (3/4)*Q*LWC/(rho_w*r_e) = (3/4)*Q*(4/3)*pi*N*<r^3>*<r^2>/<r^3>
%                               = pi*Q*N*<r^2>  =  k_ext
%
%   identically. The Hansen & Travis 3/2 coefficient is this expression with
%   Q = 2; carrying Q explicitly, the coefficient is 3Q/4. The VERIFICATION
%   section below checks each of these statements numerically on the data.
%
% OUTPUT:
%   A per-profile breakdown printed to the command window, an ensemble summary,
%   and a two-panel figure.
%
% SEE ALSO:
%   find_verticalProfiles_ORACLES, average_mie_over_size_distribution
%
% By Andrew John Buggee

clear variables

%% --- Configuration -----------------------------------------------------------

profile_to_show   = 42;      % profile examined in detail
n_profiles_ens    = 200;      % profiles used for the ensemble summary
wavelength_nm     = 550;     % wavelength of the in-situ optical thickness calc
QEXT_PYTHON       = 2.32;    % fixed Q assumed by compute_tau_c in the Python code
significance_lvl  = 0.1;     % passed to find_bestFitDist_dropDist

% Where Q_ext(r) comes from: 'miev0' runs Mie here, 'lookup' reads the saved
% table. Both are computed when the table is available, so the two can be
% compared; QEXT_SOURCE only picks which one the reported numbers use.
QEXT_SOURCE       = 'miev0';

% Monodisperse Mie grid for the reference calculation. 200 um is well past the
% radius holding 99.99% of the area integral in every ORACLES profile checked.
% The step must be fine enough to average the interference ripple within a bin.
mie_r_min_um      = 1.5;
mie_r_max_um      = 200;
mie_r_step_um     = 0.05;
QEXT_LARGE_DROP   = 2.0;     % geometric-optics limit, used above mie_r_max_um

which_computer = whatComputer();

if strcmp(which_computer, 'anbu8374')==true
    folderpath_root = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/';
    mie_folder      = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
    qext_lut_file   = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/',...
        'Q_ext_1nm_sampling_monodispersed_radii_1_to_500_microns_reduced.txt'];
elseif strcmp(which_computer, 'andrewbuggee')==true
    folderpath_root = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/';
    mie_folder      = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';
    % NOTE: interp_mie_computed_tables.m still points at the libRadtran-2.0.4
    % path for this table, which no longer holds it, so it is read directly here.
    qext_lut_file   = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/HySICS/Testing_Mie/',...
        'Q_ext_1nm_sampling_monodispersed_radii_1_to_500_microns_reduced.txt'];
elseif strcmp(which_computer, 'curc')==true
    folderpath_root = '/projects/anbu8374/Matlab-Research/';
    mie_folder      = '/projects/anbu8374/Mie_Calculations/';
    qext_lut_file   = ['/projects/anbu8374/Mie_Calculations/',...
        'Q_ext_1nm_sampling_monodispersed_radii_1_to_500_microns_reduced.txt'];
end

folderpath_oracles = [folderpath_root, 'Hyperspectral_Cloud_Retrievals/ORACLES/ORACLES_data/'];

oracles_filename = ['ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_',...
    'Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'];

% Per-level gamma fits, cached by build_synthetic_profile_inputs_from_insitu_and_ERA5.m
alpha_cache = [folderpath_oracles, 'gammaFit_alpha_', erase(oracles_filename, '.mat'), '.mat'];

%% --- Load --------------------------------------------------------------------

ds = load([folderpath_oracles, oracles_filename], 'ensemble_profiles');
ensemble_profiles = ds.ensemble_profiles;

if isfile(alpha_cache)==true
    alpha_all = load(alpha_cache, 'alpha_all');
    alpha_all = alpha_all.alpha_all;
else
    alpha_all = [];
end

%% --- Mie extinction efficiency on a fine radius grid (computed ONCE) ---------
% Every profile shares the same radius bins, so a single monodisperse Mie run
% serves the whole ensemble. This is the only expensive step and it takes a few
% seconds, not minutes.

bin_edges_um = double(ensemble_profiles{1}.drop_radius_bin_edges(:));

fprintf('\nRunning MIEV0 at %d nm over r = %.2f to %.1f um (step %.2f um)...\n', ...
    wavelength_nm, mie_r_min_um, mie_r_max_um, mie_r_step_um)

mie_fine = compute_mie_Qext_grid(mie_r_min_um, mie_r_max_um, mie_r_step_um, ...
    wavelength_nm, mie_folder, which_computer);

Q_bin_miev0 = bin_average_Qext(mie_fine, bin_edges_um, QEXT_LARGE_DROP);

fprintf('  %d radii computed.  Q_ext(%.2f um) = %.4f,  Q_ext(%.0f um) = %.4f\n', ...
    numel(mie_fine.r_um), mie_fine.r_um(1), mie_fine.Qext(1), ...
    mie_fine.r_um(end), mie_fine.Qext(end))

% --- The saved lookup table, when it is available ---
Q_bin_lookup = [];

if isfile(qext_lut_file)==true

    fprintf('Reading the saved Q_ext table...\n')
    mie_lut = read_Qext_lookup_table(qext_lut_file, wavelength_nm);
    Q_bin_lookup = bin_average_Qext(mie_lut, bin_edges_um, QEXT_LARGE_DROP);

else
    fprintf('Saved Q_ext table not found at:\n  %s\nUsing MIEV0 only.\n', qext_lut_file)
end

if strcmp(QEXT_SOURCE, 'lookup')==true && ~isempty(Q_bin_lookup)
    Q_bin = Q_bin_lookup;
else
    Q_bin = Q_bin_miev0;
    QEXT_SOURCE = 'miev0';
end

fprintf('  Using Q_ext source: %s\n', QEXT_SOURCE)

%% --- Detailed breakdown for one profile --------------------------------------

p_show = ensemble_profiles{profile_to_show};

out_one = decompose_tau(p_show, ...
    get_alpha(alpha_all, ensemble_profiles, profile_to_show, significance_lvl), ...
    wavelength_nm, QEXT_PYTHON, which_computer, Q_bin);

%% --- Verification of the underlying identities --------------------------------
% Each check below is run on the measured data, not asserted.

alpha_show = get_alpha(alpha_all, ensemble_profiles, profile_to_show, significance_lvl);
keep_show  = double(p_show.re(:)) > 1 & double(p_show.lwc(:)) > 0 & sum(double(p_show.Nd), 1)' > 1;
if length(alpha_show) < length(keep_show)
    alpha_show(end+1:length(keep_show), 1) = median(alpha_show);
end
v = verify_identities(p_show, alpha_show(keep_show));

fprintf('\n=====================================================================\n')
fprintf(' VERIFICATION against Hansen & Travis (1974)\n')
fprintf('=====================================================================\n\n')

fprintf('  (a) Nd is a per-bin concentration, so the inner integral is a plain\n')
fprintf('      sum over bins with no bin-width factor. If that reading were\n')
fprintf('      wrong, neither of these would close:\n')
fprintf('        LWC  from (4/3)pi rho_w sum(r^3 Nd) vs provider lwc : %+.4f %%\n', v.pct_lwc)
fprintf('        r_e  from sum(r^3 Nd)/sum(r^2 Nd)   vs provider re  : %+.4f %%   [H&T 2.53]\n\n', v.pct_re)

fprintf('  (b) v_eff from H&T (2.54) implies <r^4>/<r^2> = r_e^2 (1 + v_eff)\n')
fprintf('      for ANY size distribution. Measured departure:\n')
fprintf('        max |rel. difference| over all levels               : %.3e\n\n', v.max_rel_m4)

fprintf('  (c) The r_e/LWC form is algebraically identical to H&T (2.48):\n')
fprintf('        (3/4) Q LWC/(rho_w r_e)  vs  pi Q N <r^2>\n')
fprintf('        max |rel. difference| over all levels               : %.3e\n\n', v.max_rel_identity)

fprintf('  (d) Lyapunov: <r^2> <= r_e^2 for every distribution, so the r_e^2\n')
fprintf('      substitution can only bias tau high.\n')
fprintf('        levels violating <r^2> <= r_e^2                     : %d of %d\n\n', ...
    v.n_violations, v.n_levels)

fprintf('  (e) The gamma-family closed form <r^2>/r_e^2 = (1-v)(1-2v) is only\n')
fprintf('      as good as the v_eff put into it:\n')
fprintf('        measured <r^2>/r_e^2 (reciprocal of the width factor\n')
fprintf('        quoted later; note E[1/X] ~= 1/E[X] over levels)    : %.4f\n', v.width_measured_inv)
fprintf('        (1-v)(1-2v), v_eff measured via H&T (2.54)          : %.4f  (%+.2f %%)\n', ...
    v.width_gamma_inv, 100*(v.width_gamma_inv - v.width_measured_inv)/v.width_measured_inv)
fprintf('        (1-v)(1-2v), v_eff = 1/(alpha+3) from the gamma FIT : %.4f  (%+.2f %%)\n', ...
    v.width_gammaFit_inv, 100*(v.width_gammaFit_inv - v.width_measured_inv)/v.width_measured_inv)
fprintf('        The measured second moment is used everywhere below, so\n')
fprintf('        nothing in the comparison depends on either closed form.\n\n')

fprintf('\n=====================================================================\n')
fprintf(' PROFILE %d   (%s, %d levels)\n', profile_to_show, ...
    char(p_show.dateOfFlight, 'yyyy-MM-dd'), out_one.n_levels)
fprintf('=====================================================================\n\n')

ref = out_one.tau_full_mie;
pct = @(x) 100*(x - ref)/ref;

fprintf('  CLOUD OPTICAL THICKNESS, FIVE WAYS\n\n')
fprintf('    %-52s %8s   %11s\n', 'method', 'tau_c', '% diff')
fprintf('    %s\n', repmat('-', 1, 75))
fprintf('    %-52s %8.3f   %11s\n', ...
    '(0) full Mie over measured n(r,z)   [REFERENCE]', ref, 'reference')
fprintf('    %-52s %8.3f   %+10.2f %%\n', ...
    '(1) moment form, gamma-fit Q:  pi Q N <r^2>', out_one.tau_moment, pct(out_one.tau_moment))
fprintf('    %-52s %8.3f   %+10.2f %%\n', ...
    '(2) in-situ formula:           pi Q N r_e^2', out_one.tau_insitu, pct(out_one.tau_insitu))
fprintf('    %-52s %8.3f   %+10.2f %%\n', ...
    sprintf('(3) r_e/LWC, Q = %.2f:          750 Q LWC / r_e', QEXT_PYTHON), ...
    out_one.tau_reLWC, pct(out_one.tau_reLWC))
fprintf('    %-52s %8.3f   %+10.2f %%\n', ...
    '(4) r_e/LWC, Q = measured area-weighted Q', out_one.tau_reLWC_Qmeas, pct(out_one.tau_reLWC_Qmeas))
fprintf('    %s\n', repmat('-', 1, 75))
fprintf('    %-52s %8.3f   %+10.2f %%\n\n', ...
    'stored on the profile (max of tau)', out_one.tau_stored, pct(out_one.tau_stored))

fprintf('  THE TWO FACTORS SEPARATING (2) FROM (3)\n')
fprintf('    Q_ext_avg over the fitted gamma                    : %6.3f\n', out_one.Q_fit)
fprintf('    Q_ext_avg over the MEASURED distribution           : %6.3f\n', out_one.Q_measured)
fprintf('    Q assumed by the Python code                       : %6.3f   (ratio %.3f)\n', ...
    QEXT_PYTHON, out_one.Q_fit / QEXT_PYTHON)
fprintf('    r_e^2 / <r^2>, from the measured Nd                : %6.3f\n', out_one.width_measured)
fprintf('    r_e^2 / <r^2>, predicted from the fitted alpha     : %6.3f   (Hansen & Travis)\n\n', ...
    out_one.width_gammaFit)

fprintf('    (2)/(3) measured                                   : %6.3f\n', ...
    out_one.tau_insitu / out_one.tau_reLWC)
fprintf('    (Q_fit/Q_python) * (r_e^2/<r^2>)                   : %6.3f  =  %.3f * %.3f\n\n', ...
    (out_one.Q_fit/QEXT_PYTHON)*out_one.width_measured, ...
    out_one.Q_fit/QEXT_PYTHON, out_one.width_measured)

% How much does the choice of Q_ext source change the reference?
if ~isempty(Q_bin_lookup)

    out_lut = decompose_tau(p_show, ...
        get_alpha(alpha_all, ensemble_profiles, profile_to_show, significance_lvl), ...
        wavelength_nm, QEXT_PYTHON, which_computer, Q_bin_lookup);

    fprintf('  SENSITIVITY TO THE Q_ext SOURCE\n')
    fprintf('    reference tau, MIEV0 run here          : %7.3f\n', out_one.tau_full_mie)
    fprintf('    reference tau, saved lookup table      : %7.3f   (%+.2f %%)\n\n', ...
        out_lut.tau_full_mie, pct(out_lut.tau_full_mie))

end

fprintf('    Line (4) landing on line (0) is the point: with the correct Q and\n')
fprintf('    the second moment in place of r_e^2, the r_e/LWC formula reproduces\n')
fprintf('    the full Mie integral. Only the r_e^2 substitution separates them.\n\n')

%% --- Ensemble summary --------------------------------------------------------

n_ens = min(n_profiles_ens, length(ensemble_profiles));

r_moment = nan(n_ens, 1);   % (1) / (0)
r_insitu = nan(n_ens, 1);   % (2) / (0)
r_reLWC  = nan(n_ens, 1);   % (3) / (0)
r_reQmea = nan(n_ens, 1);   % (4) / (0)
r_Q      = nan(n_ens, 1);   % Q_fit / Q_python
r_width  = nan(n_ens, 1);   % r_e^2 / <r^2>
Q_meas   = nan(n_ens, 1);   % area-weighted Q over the measured distribution

fprintf('Running the same decomposition over %d profiles...\n', n_ens)

for nn = 1:n_ens

    try
        o = decompose_tau(ensemble_profiles{nn}, ...
            get_alpha(alpha_all, ensemble_profiles, nn, significance_lvl), ...
            wavelength_nm, QEXT_PYTHON, which_computer, Q_bin);
    catch
        continue
    end

    r_moment(nn) = o.tau_moment      / o.tau_full_mie;
    r_insitu(nn) = o.tau_insitu      / o.tau_full_mie;
    r_reLWC(nn)  = o.tau_reLWC       / o.tau_full_mie;
    r_reQmea(nn) = o.tau_reLWC_Qmeas / o.tau_full_mie;
    r_Q(nn)      = o.Q_fit / QEXT_PYTHON;
    r_width(nn)  = o.width_measured;
    Q_meas(nn)   = o.Q_measured;

end

good = ~isnan(r_insitu);

fprintf('\n  ENSEMBLE (n = %d profiles)\n\n', sum(good))
fprintf('    percent difference from the full-Mie reference\n\n')
fprintf('    %-46s %10s   %20s\n', 'method', 'median', '5th - 95th pctile')
fprintf('    %s\n', repmat('-', 1, 80))
print_pct('    (1) moment form, gamma-fit Q                  ', r_moment(good))
print_pct('    (2) in-situ formula                           ', r_insitu(good))
print_pct('    (3) r_e/LWC, Q = 2.32                         ', r_reLWC(good))
print_pct('    (4) r_e/LWC, Q = measured                     ', r_reQmea(good))
fprintf('\n    the two factors, as ratios\n\n')
print_stat('      Q_fit / Q_python             ', r_Q(good))
print_stat('      r_e^2 / <r^2>                ', r_width(good))
fprintf('\n    area-weighted Q over the measured distribution at %d nm\n', wavelength_nm)
print_stat('      Q_measured                   ', Q_meas(good))
fprintf('\n')

%% --- Figure ------------------------------------------------------------------

figure('Position', [100 100 1050 430])

subplot(1,2,1)
vals = [out_one.tau_full_mie, out_one.tau_moment, out_one.tau_insitu, ...
    out_one.tau_reLWC, out_one.tau_reLWC_Qmeas];
b = bar(vals);
b.FaceColor = 'flat';
b.CData = [0.15 0.15 0.15; 0.45 0.45 0.45; 0.70 0.15 0.15; 0.15 0.35 0.70; 0.35 0.60 0.35];
set(gca, 'XTickLabel', {'(0) full Mie', '(1) moment', '(2) in-situ', ...
    '(3) r_e/LWC', '(4) r_e/LWC'})
ylabel('\tau_c')
title(sprintf('Profile %d: same cloud, five calculations', profile_to_show))
text(1:5, vals, compose('%.2f', vals(:)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
xlabel('(3) uses Q = 2.32, (4) uses the measured Q')
grid on
ylim([0, 1.15*max(vals)])

subplot(1,2,2)
hold on
histogram(100*(r_insitu(good)-1), 15, 'FaceColor', [0.70 0.15 0.15], 'FaceAlpha', 0.55)
histogram(100*(r_reLWC(good)-1),  15, 'FaceColor', [0.15 0.35 0.70], 'FaceAlpha', 0.55)
histogram(100*(r_reQmea(good)-1), 15, 'FaceColor', [0.35 0.60 0.35], 'FaceAlpha', 0.55)
xline(0, 'k--', 'LineWidth', 1.5)
legend({'(2) in-situ', '(3) r_e/LWC, Q=2.32', '(4) r_e/LWC, Q measured', ...
    'full Mie reference'}, 'Location', 'northeast')
xlabel('percent difference from \tau_{full Mie}  (%)')
ylabel('Number of profiles')
title('Bias against the full Mie integral')
grid on


%% =============================================================================
%  Local functions
%  =============================================================================

function print_pct(label, ratio)
% Print a ratio column as a percent difference from the reference.
d = 100*(ratio - 1);
fprintf('%s %+9.2f %%   %+8.2f to %+7.2f %%\n', label, median(d), prctile(d,5), prctile(d,95));
end



function print_stat(label, v)
fprintf('%s: %5.3f  [%5.3f - %5.3f]\n', label, median(v), prctile(v,5), prctile(v,95));
end



function mie = compute_mie_Qext_grid(r_min_um, r_max_um, r_step_um, wavelength_nm, ...
    mie_folder, which_computer)
% Monodisperse Mie extinction efficiency on an evenly spaced radius grid.
% One MIEV0 run serves every profile and every altitude.

r = r_min_um : r_step_um : r_max_um;        % microns

% MIEV0 takes [r_start, r_end, r_step] and rebuilds the vector itself; nudge the
% step until its vector has the same length as ours (see the same guard in
% average_mie_over_size_distribution).
radius_step = mean(diff(r));
while numel(r(1):radius_step:r(end)) ~= numel(r)
    radius_step = radius_step - 1e-15;
end

[input_filename, output_filename] = write_mie_file('MIEV0', 'water', ...
    [r(1), r(end), radius_step], wavelength_nm, {'mono', []}, 'verbose', ...
    991, randi([0, 1000], 1), mie_folder);

[~] = runMIE(mie_folder, input_filename, output_filename, which_computer);

[ds, ~, ~] = readMIE(mie_folder, output_filename);

mie.r_um = r(:);
mie.Qext = reshape(ds.Qext(1, :), [], 1);

end



function v = verify_identities(p, alpha_prof_for_v)
% Check, on the measured data, every algebraic step the comparison relies on.
% All equation numbers are Hansen & Travis (1974).

r_um = double(p.drop_radius_bin_center(:));       % um
Nd   = double(p.Nd);                              % cm^-3 per bin
re   = double(p.re(:));                           % um    (provider)
lwc  = double(p.lwc(:));                          % g/m^3 (provider)

N    = sum(Nd, 1)';                               % cm^-3
m2   = (sum((r_um.^2) .* Nd, 1)') ./ max(N, 1e-30);
m3   = (sum((r_um.^3) .* Nd, 1)') ./ max(N, 1e-30);
m4   = (sum((r_um.^4) .* Nd, 1)') ./ max(N, 1e-30);

keep = re > 1 & lwc > 0 & N > 1;

r_e_dist = m3(keep) ./ m2(keep);                  % H&T (2.53): r_e = <r^3>/<r^2>

% --- (a) Does the per-bin reading of Nd reproduce the provider's fields? ---
% LWC = (4/3) pi rho_w N <r^3>. r um -> cm (1e-4), Nd cm^-3, rho_w = 1 g/cm^3,
% then g/cm^3 -> g/m^3 (1e6).
lwc_dist = (4/3) * pi * sum(((r_um .* 1e-4).^3) .* Nd, 1)' .* 1e6;

v.pct_lwc = 100 * median((lwc_dist(keep) - lwc(keep)) ./ lwc(keep));
v.pct_re  = 100 * median((r_e_dist - re(keep)) ./ re(keep));

% --- (b) v_eff from H&T (2.54), and the identity it forces ---
% v_eff = int (r-r_e)^2 pi r^2 n dr / (r_e^2 int pi r^2 n dr)
% Expanding gives <r^4>/<r^2> = r_e^2 (1 + v_eff) for ANY distribution.
Nd_k  = Nd(:, keep);
num   = sum(((r_um - r_e_dist').^2) .* (r_um.^2) .* Nd_k, 1)';
den   = sum((r_um.^2) .* Nd_k, 1)';
v_eff = num ./ (den .* r_e_dist.^2);

lhs = m4(keep) ./ m2(keep);
rhs = (r_e_dist.^2) .* (1 + v_eff);
v.max_rel_m4 = max(abs(lhs - rhs) ./ rhs);

% --- (c) The r_e/LWC form against H&T (2.48), at fixed Q ---
% Any Q works: the identity is independent of it.
Q_test    = 2.1;
N_m3      = N(keep) .* 1e6;
m2_m2     = m2(keep) .* 1e-12;
m3_m3     = m3(keep) .* 1e-18;
r_e_m     = r_e_dist .* 1e-6;

k_ext     = pi * Q_test .* N_m3 .* m2_m2;                       % H&T (2.48)
lwc_kgm3  = (4/3) * pi * 1e3 .* N_m3 .* m3_m3;                  % kg/m^3
k_from_re = (3/4) * Q_test .* lwc_kgm3 ./ (1e3 .* r_e_m);       % 3Q/4 form

v.max_rel_identity = max(abs(k_from_re - k_ext) ./ k_ext);

% --- (d) Lyapunov: <r^2> <= r_e^2 at every level ---
v.n_violations = sum(m2(keep) > r_e_dist.^2 + eps(r_e_dist.^2));
v.n_levels     = sum(keep);

% --- (e) Measured width factor vs the gamma-family closed form ---
w = (r_um.^2)' * Nd_k;                        % area weight per level
w = w' ./ sum(w);

v.width_measured_inv = sum(w .* (m2(keep) ./ r_e_dist.^2));
v.width_gamma_inv    = sum(w .* ((1 - v_eff) .* (1 - 2.*v_eff)));

% Same closed form, but with v_eff from the fitted gamma shape parameter
% (v_eff = 1/(alpha+3)) rather than measured from the distribution.
v_fit = 1 ./ (alpha_prof_for_v(:) + 3);
v.width_gammaFit_inv = sum(w .* ((1 - v_fit) .* (1 - 2.*v_fit)));

end



function mie = read_Qext_lookup_table(lut_file, wavelength_nm)
% Read the saved monodisperse Q_ext table at one wavelength.
%
% The file holds Q_ext on a (340 radius x 2901 wavelength) grid: radii 1-300 um
% in 1 um steps then 305-500 um in 5 um steps, wavelengths 100-3000 nm in 1 nm
% steps, ordered radius-fastest within each wavelength.

file_id = fopen(lut_file);

if file_id == -1
    error([newline, 'Could not open the Q_ext lookup table: ', lut_file, newline])
end

data = textscan(file_id, '%f', 'CommentStyle', '#', 'Delimiter', {'\r', ' '}, ...
    'MultipleDelimsAsOne', 1);
fclose(file_id);

Q_all = reshape(data{1}, 340, []);            % (radius x wavelength)

wl_grid = 100:1:3000;                         % nm
r_grid  = [1:300, 305:5:500];                 % microns

[~, idx_wl] = min(abs(wl_grid - wavelength_nm));

mie.r_um = r_grid(:);
mie.Qext = Q_all(:, idx_wl);

end



function Q_bin = bin_average_Qext(mie, bin_edges_um, qext_large_drop)
% Average Q_ext over each measured size bin, weighted by r^2.
%
% Monodisperse Q_ext oscillates faster in r than the narrowest bins are wide, so
% point-sampling it at a bin center picks an arbitrary phase of the ripple. The
% physical quantity each bin contributes is int(Q(r) r^2 n(r) dr) over the bin,
% and n is only known as a bin total, so the representative efficiency is the
% r^2-weighted mean of Q across the bin.

n_bins = length(bin_edges_um) - 1;
Q_bin  = zeros(n_bins, 1);

r_fine = mie.r_um(:);
Q_fine = mie.Qext(:);
w_fine = r_fine.^2;

for ii = 1:n_bins

    in_bin = r_fine >= bin_edges_um(ii) & r_fine < bin_edges_um(ii+1);

    if any(in_bin)
        Q_bin(ii) = sum(Q_fine(in_bin) .* w_fine(in_bin)) / sum(w_fine(in_bin));

    elseif bin_edges_um(ii) > r_fine(end)
        % Above the grid: geometric optics. These bins hold < 0.01% of the area.
        Q_bin(ii) = qext_large_drop;

    else
        % Bin narrower than the grid spacing: fall back to interpolation
        r_center  = 0.5 * (bin_edges_um(ii) + bin_edges_um(ii+1));
        Q_bin(ii) = interp1(r_fine, Q_fine, r_center, 'linear', 'extrap');
    end

end

end



function alpha_prof = get_alpha(alpha_all, ensemble_profiles, nn, significance_lvl)
% Per-level gamma shape parameter: from the cache when available, else fit it.

if ~isempty(alpha_all) && nn <= length(alpha_all)
    alpha_prof = double(alpha_all{nn}(:));
else
    [~, ~, gammaFit] = find_bestFitDist_dropDist(ensemble_profiles{nn}.Nd, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, significance_lvl);
    alpha_prof = double(gammaFit.alpha(:));
end

alpha_prof(~isfinite(alpha_prof) | alpha_prof <= 0) = 7;

end



function out = decompose_tau(p, alpha_prof, wavelength_nm, qext_python, ...
    which_computer, Q_bin)
% Compute tau five ways for one profile and return the decomposition factors.

r_um  = double(p.drop_radius_bin_center(:));      % um, bin centers
Nd    = double(p.Nd);                             % cm^-3 per bin, (n_bins x n_lev)
re    = double(p.re(:));                          % um      (provider)
lwc   = double(p.lwc(:));                         % g/m^3   (provider)
Nc    = double(p.total_Nc(:));                    % cm^-3   (provider)
z_m   = double(p.altitude(:));                    % m

n_lev = length(z_m);

% Pad or trim the alpha profile to the profile length
if length(alpha_prof) < n_lev
    alpha_prof(end+1:n_lev, 1) = median(alpha_prof);
else
    alpha_prof = alpha_prof(1:n_lev);
end

% --- Moments straight from the measured size distribution ---
N_dist = sum(Nd, 1)';                                            % cm^-3
m2     = (sum((r_um.^2) .* Nd, 1)') ./ max(N_dist, 1e-30);       % <r^2>, um^2
m3     = (sum((r_um.^3) .* Nd, 1)') ./ max(N_dist, 1e-30);       % <r^3>, um^3

% --- Keep only levels where every quantity is usable ---
keep = re > 1 & lwc > 0 & N_dist > 1 & isfinite(m2) & isfinite(m3);
if sum(keep) < 5
    error('too few usable levels')
end

r_e     = re(keep);
lwc_k   = lwc(keep);
Nc_k    = Nc(keep);
m2_k    = m2(keep);
z_k     = z_m(keep);
Nd_k    = Nd(:, keep);
alpha_k = alpha_prof(keep);

% -------------------------------------------------------------------------
% (0) REFERENCE: full Mie integral over the measured size distribution
%
%     beta(z) = pi * sum_i Q_ext(r_i) * r_i^2 * Nd_i(z)
%
%     Nd is a per-bin concentration, so no bin-width factor enters. Q_bin is the
%     r^2-weighted average of Q_ext across each measured bin, computed once by
%     the caller.
% -------------------------------------------------------------------------
r_bin_m2 = (r_um .* 1e-6).^2;                     % um -> m, squared
Nd_k_m3  = Nd_k .* 1e6;                           % cm^-3 -> m^-3

beta_full_mie = pi .* sum((Q_bin .* r_bin_m2) .* Nd_k_m3, 1)';    % m^-1

out.tau_full_mie = abs(trapz(z_k, beta_full_mie));

% Area-weighted Q implied by the measured distribution, per level
area_per_level = sum(r_bin_m2 .* Nd_k_m3, 1)';
Q_measured_lev = sum((Q_bin .* r_bin_m2) .* Nd_k_m3, 1)' ./ max(area_per_level, 1e-30);

% -------------------------------------------------------------------------
% (1) and (2): the moment form and the in-situ formula, both using the
%     area-weighted Q over the FITTED gamma distribution, exactly as
%     find_verticalProfiles_ORACLES.m does.
% -------------------------------------------------------------------------
[~, Qe_fit, ~] = average_mie_over_size_distribution(r_e', alpha_k', wavelength_nm, ...
    'water', 'gamma', which_computer, 1);
Qe_fit = reshape(Qe_fit(1, :), [], 1);

Nc_m3 = Nc_k .* 1e6;              % cm^-3 -> m^-3
m2_m2 = m2_k .* 1e-12;            % um^2  -> m^2
re_m  = r_e  .* 1e-6;             % um    -> m

% Both lines use the provider's total_Nc, exactly as find_verticalProfiles_ORACLES
% does, so their RATIO isolates <r^2> vs r_e^2 cleanly. Note that total_Nc runs
% ~0.7% below sum(Nd), so line (1) sits ~0.7% above the reference for that reason
% rather than because of the gamma-fit Q, which is worth only about -0.2%.
beta_moment = pi .* Qe_fit .* Nc_m3 .* m2_m2;       % correct moment pairing
beta_insitu = pi .* Qe_fit .* Nc_m3 .* re_m.^2;     % find_verticalProfiles_ORACLES.m

out.tau_moment = abs(trapz(z_k, beta_moment));
out.tau_insitu = abs(trapz(z_k, beta_insitu));

% -------------------------------------------------------------------------
% (3) and (4): the r_e / LWC form. 750 * Q * LWC[g/m^3] / r_e[um] is tau per
%     km, so divide by 1e3 to get an extinction coefficient in m^-1.
% -------------------------------------------------------------------------
beta_reLWC       = 750 .* qext_python    .* (lwc_k ./ r_e) ./ 1e3;
beta_reLWC_Qmeas = 750 .* Q_measured_lev .* (lwc_k ./ r_e) ./ 1e3;

out.tau_reLWC       = abs(trapz(z_k, beta_reLWC));
out.tau_reLWC_Qmeas = abs(trapz(z_k, beta_reLWC_Qmeas));

out.tau_stored = max(double(p.tau));

% --- Decomposition factors, weighted by their contribution to the reference ---
w = beta_full_mie ./ sum(beta_full_mie);

out.Q_fit          = sum(w .* Qe_fit);
out.Q_measured     = sum(w .* Q_measured_lev);
out.width_measured = sum(w .* (re_m.^2 ./ m2_m2));

v_eff = 1 ./ (alpha_k + 3);
out.width_gammaFit = sum(w .* (1 ./ ((1 - v_eff) .* (1 - 2.*v_eff))));

out.n_levels = sum(keep);

end
