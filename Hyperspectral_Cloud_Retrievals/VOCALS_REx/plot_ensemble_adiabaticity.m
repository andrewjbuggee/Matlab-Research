%% Plot the Wood (2005) adiabaticity of an ensemble of in-situ vertical profiles
%
% DESCRIPTION:
%   One figure, three panels, built from the .adiabaticity structure that
%   compute_adiabaticity_Wood2005 attaches to every profile:
%
%   (a) The per-level adiabaticity f_ad(z) = LWC(z)/q_ad(z) of every profile
%       against normalized height z* = (z - z_cb)/h, with the ensemble
%       median and 10th-90th percentile band in z* bins. This is the
%       in-situ analogue of Wood (2005) Fig. 4c.
%
%   (b) Histograms of the two scalar adiabaticities: the path method
%       f_ad = LWP_obs / LWP_ad (Wood's Table 3 definition) and the profile
%       method (mean of f_ad(z) over z* >= zstar_min).
%
%   (c) The two scalars against each other with the 1:1 line. Points below
%       the line are clouds whose upper levels are more depleted than their
%       lower levels (the path method weights the top more heavily).
%
% INPUTS:
%   ensemble_profiles - the ensemble in any of the layouts used across this
%                       code base:
%                         (a) struct of cell arrays  ensemble_profiles.adiabaticity{nn}
%                         (b) cell array of structs  ensemble_profiles{nn}.adiabaticity
%                         (c) struct array           ensemble_profiles(nn).adiabaticity
%   campaign_label    - (optional) string for the title, e.g. 'VOCALS-REx'
%
% OUTPUT:
%   adiab_summary - structure of 1 x n_profiles vectors:
%       .f_ad_lwp              path-method adiabaticity [-]
%       .f_ad_profile          profile-method adiabaticity [-]
%       .Gamma_ad_g_m3_per_km  adiabatic LWC lapse rate used [g/m^3 per km]
%       .T_mean_C              mean in-cloud temperature used [deg C]
%       .h_m                   cloud thickness [m]
%   Useful for segmenting the ensemble by adiabaticity regime the same way
%   the scripts already segment by LWP and optical depth.
%
% SEE ALSO:
%   compute_adiabaticity_Wood2005, ensemble_vertical_statistics,
%   ensemble_vertical_statistics_ORACLES
%
% By Andrew John Buggee
%%

function adiab_summary = plot_ensemble_adiabaticity(ensemble_profiles, campaign_label)

if nargin < 2 || isempty(campaign_label)
    campaign_label = '';
end

%% Pull the adiabaticity structure out of whichever layout was handed in

if iscell(ensemble_profiles)
    % (b) cell array of per-profile structs
    profiles = ensemble_profiles;

elseif isstruct(ensemble_profiles) && isscalar(ensemble_profiles) && ...
        isfield(ensemble_profiles, 'adiabaticity') && iscell(ensemble_profiles.adiabaticity)
    % (a) struct of cell arrays -- rebuild a cell array of minimal structs
    n_prof = numel(ensemble_profiles.adiabaticity);
    profiles = cell(1, n_prof);
    for nn = 1:n_prof
        profiles{nn}.adiabaticity = ensemble_profiles.adiabaticity{nn};
        if isfield(ensemble_profiles, 'LWC_threshold')
            profiles{nn}.LWC_threshold = ensemble_profiles.LWC_threshold{nn};
        end
        if isfield(ensemble_profiles, 'Nc_threshold')
            profiles{nn}.Nc_threshold = ensemble_profiles.Nc_threshold{nn};
        end
    end

elseif isstruct(ensemble_profiles)
    % (c) struct array of profiles
    profiles = num2cell(ensemble_profiles);

else
    error([newline, 'Unrecognized ensemble layout.', newline])
end

n_prof = numel(profiles);

if n_prof == 0
    error([newline, 'No profiles to plot.', newline])
end

if ~isfield(profiles{1}, 'adiabaticity') || isempty(profiles{1}.adiabaticity)
    error([newline, 'Profiles carry no .adiabaticity field. Re-run the profile finder ', ...
        '(find_verticalProfiles_VOCALS_REx_ver2 / find_verticalProfiles_ORACLES) so that ', ...
        'compute_adiabaticity_Wood2005 is applied.', newline])
end


%% Collect the scalars and the per-level curves

f_ad_lwp     = nan(1, n_prof);
f_ad_profile = nan(1, n_prof);
Gamma_ad     = nan(1, n_prof);          % g/m^3 per km
T_mean_C     = nan(1, n_prof);
h_m          = nan(1, n_prof);

z_star_all = [];
f_ad_z_all = [];

for nn = 1:n_prof

    a = profiles{nn}.adiabaticity;

    f_ad_lwp(nn)     = a.f_ad_lwp;
    f_ad_profile(nn) = a.f_ad_profile;
    Gamma_ad(nn)     = a.Gamma_ad_g_m3_per_km;
    T_mean_C(nn)     = a.T_mean_C;
    h_m(nn)          = a.h_m;

    keep = isfinite(a.z_star) & isfinite(a.f_ad_z);
    z_star_all = [z_star_all, a.z_star(keep)];          %#ok<AGROW>
    f_ad_z_all = [f_ad_z_all, a.f_ad_z(keep)];          %#ok<AGROW>

end

zstar_min = profiles{1}.adiabaticity.zstar_min;

adiab_summary.f_ad_lwp             = f_ad_lwp;
adiab_summary.f_ad_profile         = f_ad_profile;
adiab_summary.Gamma_ad_g_m3_per_km = Gamma_ad;
adiab_summary.T_mean_C             = T_mean_C;
adiab_summary.h_m                  = h_m;

% thresholds for the title, if the profiles carry them
threshold_str = '';
if isfield(profiles{1}, 'LWC_threshold') && isfield(profiles{1}, 'Nc_threshold')
    threshold_str = ['  |  LWC $\geq$ ', num2str(profiles{1}.LWC_threshold), ' $g/m^3$  |  ', ...
        '$N_c \geq$ ', num2str(profiles{1}.Nc_threshold), ' $cm^{-3}$'];
end


%% Bin the per-level ratio in normalized height for the median curve

n_bins    = 20;
bin_edges = linspace(0, 1, n_bins + 1);
bin_center = 0.5 * (bin_edges(1:end-1) + bin_edges(2:end));

lowP  = 10;
highP = 90;

f_med  = nan(n_bins, 1);
f_low  = nan(n_bins, 1);
f_high = nan(n_bins, 1);

for bb = 1:n_bins

    if bb == 1
        in_bin = z_star_all >= bin_edges(bb) & z_star_all <= bin_edges(bb+1);
    else
        in_bin = z_star_all >  bin_edges(bb) & z_star_all <= bin_edges(bb+1);
    end

    if any(in_bin)
        f_med(bb)  = median(f_ad_z_all(in_bin));
        f_low(bb)  = prctile(f_ad_z_all(in_bin), lowP);
        f_high(bb) = prctile(f_ad_z_all(in_bin), highP);
    end

end


%% Figure

clr_profile = mySavedColors(2, 'fixed');      % sea foam green
clr_lwp     = mySavedColors(5, 'fixed');      % salmon red
clr_lines   = [0.75, 0.75, 0.75];

figure;

% ---- (a) f_ad(z) versus normalized height, Wood (2005) Fig. 4c ---------
subplot(1, 3, 1)
hold on

for nn = 1:n_prof
    a = profiles{nn}.adiabaticity;
    plot(a.f_ad_z, a.z_star, '-', 'Color', clr_lines, 'LineWidth', 0.5)
end

ok = isfinite(f_low) & isfinite(f_high);
x = [f_low(ok); flipud(f_high(ok))];
y = [bin_center(ok)'; flipud(bin_center(ok)')];
fill(x, y, clr_profile, 'EdgeAlpha', 0, 'FaceAlpha', 0.3)
plot(f_med, bin_center, '-', 'Color', clr_profile, 'LineWidth', 2)

xline(1, 'k--', 'LineWidth', 1)
yline(zstar_min, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1)

x_max = max(1.5, prctile(f_ad_z_all(z_star_all >= zstar_min), 97.5));
xlim([0, x_max])
ylim([0, 1])
grid on; grid minor
xlabel('$LWC(z) \, / \, q_{ad}(z)$', 'Interpreter', 'latex')
ylabel('$z_* = (z - z_{cb}) / h$', 'Interpreter', 'latex')
title(['Per level: median, ', num2str(lowP), '--', num2str(highP), 'th pct.'], ...
    'Interpreter', 'latex')


% ---- (b) histograms of the two scalar adiabaticities ------------------
subplot(1, 3, 2)
hold on

f_all = [f_ad_lwp(isfinite(f_ad_lwp)), f_ad_profile(isfinite(f_ad_profile))];
hist_edges = linspace(0, max(1.5, prctile(f_all, 99)), 26);

histogram(f_ad_lwp, hist_edges, 'FaceColor', clr_lwp, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
histogram(f_ad_profile, hist_edges, 'FaceColor', clr_profile, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
xline(1, 'k--', 'LineWidth', 1)

grid on; grid minor
xlabel('$f_{ad}$', 'Interpreter', 'latex')
ylabel('Counts')
legend({['$LWP_{obs}/LWP_{ad}$  (median ', num2str(median(f_ad_lwp, 'omitnan'), '%.2f'), ')'], ...
    ['$\langle LWC/q_{ad} \rangle_{z_* \geq ', num2str(zstar_min), '}$  (median ', ...
    num2str(median(f_ad_profile, 'omitnan'), '%.2f'), ')']}, ...
    'Interpreter', 'latex', 'Location', 'best')
title('Cloud-integrated $f_{ad}$', 'Interpreter', 'latex')


% ---- (c) path method versus profile method ---------------------------
subplot(1, 3, 3)
hold on

scatter(f_ad_profile, f_ad_lwp, 30, h_m, 'filled', 'MarkerFaceAlpha', 0.8)
cb = colorbar;
cb.Label.String = 'cloud thickness h (m)';

lim_max = max(1.5, prctile(f_all, 99));
plot([0, lim_max], [0, lim_max], 'k--', 'LineWidth', 1)
xlim([0, lim_max]); ylim([0, lim_max])
axis square
grid on; grid minor
xlabel('$\langle LWC / q_{ad} \rangle$  (profile method)', 'Interpreter', 'latex')
ylabel('$LWP_{obs} / LWP_{ad}$  (path method)', 'Interpreter', 'latex')
title('Path vs profile method', 'Interpreter', 'latex')


gamma_str = ['$\Gamma_{ad}$ = ', num2str(median(Gamma_ad, 'omitnan'), '%.2f'), ...
    ' [', num2str(min(Gamma_ad), '%.2f'), ', ', num2str(max(Gamma_ad), '%.2f'), ...
    '] $g\,m^{-3}\,km^{-1}$ at the mean in-cloud $(T, p)$'];

if ~isempty(campaign_label)
    line1 = ['Wood (2005) adiabaticity -- ', campaign_label, ', ', num2str(n_prof), ' profiles', threshold_str];
else
    line1 = ['Wood (2005) adiabaticity -- ', num2str(n_prof), ' profiles', threshold_str];
end

sgtitle({line1, gamma_str}, 'Interpreter', 'latex')

set(gcf, 'Position', [0 0 1500 500])

end
