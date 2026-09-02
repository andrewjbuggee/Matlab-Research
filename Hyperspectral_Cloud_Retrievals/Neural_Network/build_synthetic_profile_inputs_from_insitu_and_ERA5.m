%% Build a compact in-situ + ERA5 data set for the synthetic profile generator
%
% DESCRIPTION:
%   Reads the saved VOCALS-REx and ORACLES ensemble vertical profiles, extracts
%   only the variables needed by
%
%     Python-Research/lasp-CU-paper-3/06_synthetic_profile_generator.py
%
%   pairs every profile with the spatiotemporally closest ERA5 temperature /
%   water vapor profile, and writes the result to a small number of .mat files.
%
%   Unlike the per-profile files in
%   '/Volumes/My Passport/neural_network_training_data/saz0_allProfiles/',
%   these files contain no reflectances, and the VOCALS-REx profiles DO carry
%   their gamma-distribution shape parameter (alpha), so both campaigns can be
%   used to train the generative model.
%
% VARIABLES SAVED (per campaign, one row / cell entry per profile):
%   re                            - in-cloud effective radius profile [microns]
%   lwc                           - in-cloud liquid water content profile [g/m^3]
%   z                             - in-cloud altitude profile [km]
%   alpha                         - in-cloud gamma shape parameter profile [unitless]
%   tau                           - in-cloud cumulative optical depth profile [unitless]
%   z_top, z_base, thickness      - cloud top / base altitude and depth [km]
%   lwp_insitu, lwp_computed      - liquid water path [g/m^2]
%   tau_c                         - total cloud optical thickness, from the
%                                   stored in-situ tau profile [unitless]
%   tau_c_computed                - total cloud optical thickness recomputed
%                                   from r_e and LWC [unitless]
%   n_levels                      - number of in-cloud measurements per profile
%   era5_T_K                      - ERA5 temperature profile [K], surface -> TOA
%   era5_vapor_concentration_cm3  - ERA5 water vapor number density [molec/cm^3]
%   era5_q_kgkg                   - ERA5 specific humidity [kg/kg]
%   era5_vapor_massDensity_kgm3   - ERA5 water vapor mass density [kg/m^3]
%   era5_pressure_hPa             - ERA5 pressure levels [hPa]
%   era5_GP_height_m              - ERA5 geometric altitude [m]
%   ... plus per-profile geolocation / timing / ERA5 match metadata.
%
%   All in-cloud profile vectors are column vectors ordered cloud top (index 1)
%   to cloud base (index end). All ERA5 profiles are ordered surface -> TOA,
%   matching the era5.datProfiles convention of
%   write_ERA5_radiosonde_DAT_with_multiPixels.
%
% OUTPUT FILES (written to output_folder):
%   insitu_era5_synthetic_inputs_vocals-rex_<date>.mat
%   insitu_era5_synthetic_inputs_oracles_<date>.mat
%   insitu_era5_synthetic_inputs_combined_<date>.mat
%
%   Saved with '-v7' so scipy.io.loadmat can read them directly.
%
% NOTES:
%   - A handful of the stored ORACLES tau profiles are corrupt (tau ~ 1e-4 for
%     clouds with LWP > 50 g/m^2). The script prints which profiles these are
%     and saves tau_c_computed alongside tau_c so they can still be used.
%   - 22 of the 88 VOCALS-REx profiles have non-conforming 2DC data and so lack
%     the combined effective radius field '.re'. Those profiles fall back to the
%     CDP-only '.re_CDP', matching
%     hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver4.m. Which field was
%     used is recorded per profile in 're_source'.
%   - VOCALS-REx profiles store the per-altitude gamma fit at
%     ensemble_profiles{n}.gammaFit.alpha. ORACLES profiles do not, so alpha is
%     computed here with find_bestFitDist_dropDist (the same call used by
%     hysics_refl_from_oracles_and_era5_SZA_loopGeometry.m). That fit is slow,
%     so the result is cached to disk and reused on later runs.
%   - The ERA5 temperature / vapor profiles are extracted with
%     findClosestProfile_ERA5_VOCALS_REx.m and findClosestProfile_ERA5_ORACLES.m.
%     The specific-humidity -> vapor number density conversion replicates
%     write_ERA5_radiosonde_DAT_with_multiPixels.m without writing any
%     libRadtran radiosonde .dat files.
%   - Reading the VOCALS-REx ERA5 files is the slow part of this script
%     (each netCDF file is ~80 MB and is re-read for every profile).
%
% SEE ALSO:
%   findClosestProfile_ERA5_VOCALS_REx, findClosestProfile_ERA5_ORACLES,
%   write_ERA5_radiosonde_DAT_with_multiPixels, find_bestFitDist_dropDist
%
% By Andrew John Buggee

%% --- Setup -------------------------------------------------------------------

clear variables

which_computer = whatComputer();

% --- Which campaigns to process ---
process_vocals  = true;
process_oracles = true;

% --- Processing flags ---
% These mirror the cleaning steps applied in the HySICS reflectance scripts so
% the synthetic generator is trained on the same profiles that were used to
% simulate reflectances.
inputs.remove_levels_re_le_1um   = true;   % drop levels with r_e <= 1 micron
inputs.remove_levels_re_ge_50um  = true;   % drop levels with r_e >= 50 microns
inputs.min_num_levels            = 3;      % skip profiles thinner than this
inputs.alpha_default             = 7;      % used when the gamma fit fails outright
inputs.alpha_max                 = 150;    % libRadtran chokes above this
inputs.gammaFit_significance_lvl = 0.1;    % passed to find_bestFitDist_dropDist

% Bulk Mie extinction efficiency used for the recomputed optical thickness.
% 2 is the geometric optics limit (Hansen & Travis, 1974); 2.32 is the value
% 06_synthetic_profile_generator.py uses, calibrated against libRadtran.
inputs.qext_eff                  = 2.32;
inputs.recompute_oracles_alpha   = false;  % true forces a re-fit of the ORACLES alphas
inputs.print_era5_status         = false;

% Process only the first N profiles of each campaign. Set to [] to process them
% all. Useful for a quick smoke test, since reading the ~80 MB VOCALS-REx ERA5
% files and fitting the ORACLES droplet distributions are both slow.
inputs.max_profiles_per_campaign = [];

%% --- Define folder paths -----------------------------------------------------

if strcmp(which_computer, 'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folderpath_root = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/';

elseif strcmp(which_computer, 'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    folderpath_root = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/';

elseif strcmp(which_computer, 'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    folderpath_root = '/projects/anbu8374/Matlab-Research/';

end

folderpath_vocals = [folderpath_root, 'Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
    'vocals_rex_data/NCAR_C130/SPS_1/'];

folderpath_oracles = [folderpath_root, 'Hyperspectral_Cloud_Retrievals/ORACLES/',...
    'ORACLES_data/'];

output_folder = [folderpath_root, 'Hyperspectral_Cloud_Retrievals/Neural_Network/',...
    'Training_data_set/synthetic_profile_inputs/'];

if ~isfolder(output_folder)
    mkdir(output_folder)
end

% --- Saved ensemble profile files ---
vocals_filename = ['ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_',...
    'Nc-threshold_25_04-Dec-2025.mat'];

oracles_filename = ['ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_',...
    'Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'];

% --- Where the (slow) ORACLES gamma fits get cached ---
oracles_alpha_cache = [folderpath_oracles, 'gammaFit_alpha_', ...
    erase(oracles_filename, '.mat'), '.mat'];

% --- Date string appended to every output file ---
date_str = char(datetime('today', 'Format', 'dd-MMM-yyyy'));

%% --- VOCALS-REx --------------------------------------------------------------

if process_vocals==true

    fprintf('\n=============== VOCALS-REx ===============\n')

    ds_vocals = load([folderpath_vocals, vocals_filename], 'ensemble_profiles');
    ensemble_profiles = truncate_profiles(ds_vocals.ensemble_profiles, inputs);

    n_total = length(ensemble_profiles);
    fprintf('Loaded %d VOCALS-REx ensemble profiles\n', n_total)

    % --- Grab the raw in-cloud profiles ---
    raw = cell(n_total, 1);

    for nn = 1:n_total

        p = ensemble_profiles{nn};

        % Profiles whose 2DC data didn't conform lack the combined CDP + 2DC
        % effective radius. Fall back to the CDP-only effective radius, which is
        % what hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver4.m does.
        if isfield(p, 're')==true
            raw{nn}.re        = double(p.re(:));         % microns
            raw{nn}.re_source = 're';
        else
            raw{nn}.re        = double(p.re_CDP(:));     % microns
            raw{nn}.re_source = 're_CDP';
        end

        raw{nn}.lwc   = double(p.lwc(:));                % g/m^3
        raw{nn}.z     = double(p.altitude(:)) ./ 1e3;    % meters -> km
        raw{nn}.tau   = double(p.tau(:));                % unitless
        raw{nn}.alpha = double(p.gammaFit.alpha(:));     % unitless

        raw{nn}.lwp_insitu = double(p.lwp(1));           % g/m^2

        % mid-profile geolocation and time
        mid = round(length(p.latitude)/2);
        raw{nn}.latitude       = double(p.latitude(mid));
        raw{nn}.longitude      = double(p.longitude(mid));
        raw{nn}.time_utc_hours = double(p.time_utc(round(end/2)));   % decimal hours UTC
        raw{nn}.dateOfFlight   = p.dateOfFlight;

    end

    vocals_out = package_campaign(raw, ensemble_profiles, 'vocals-rex', ...
        [folderpath_vocals, vocals_filename], inputs, which_computer);

    save_campaign_file(vocals_out, output_folder, ...
        ['insitu_era5_synthetic_inputs_vocals-rex_', date_str, '.mat'])

end

%% --- ORACLES -----------------------------------------------------------------

if process_oracles==true

    fprintf('\n================ ORACLES =================\n')

    ds_oracles = load([folderpath_oracles, oracles_filename], 'ensemble_profiles');
    ensemble_profiles = truncate_profiles(ds_oracles.ensemble_profiles, inputs);

    n_total = length(ensemble_profiles);
    fprintf('Loaded %d ORACLES ensemble profiles\n', n_total)

    % ---------------------------------------------------------------------
    % --- Gamma shape parameter (alpha) -----------------------------------
    % ---------------------------------------------------------------------
    % ORACLES profiles don't store a gamma fit, so fit the measured droplet
    % size distribution at every in-cloud altitude. This takes a while, so the
    % result is cached and reused.

    alpha_all = load_or_compute_oracles_alpha(ensemble_profiles, ...
        oracles_alpha_cache, oracles_filename, inputs);

    % --- Grab the raw in-cloud profiles ---
    raw = cell(n_total, 1);

    for nn = 1:n_total

        p = ensemble_profiles{nn};

        raw{nn}.re        = double(p.re(:));             % microns
        raw{nn}.re_source = 're';

        raw{nn}.lwc   = double(p.lwc(:));                % g/m^3
        raw{nn}.z     = double(p.altitude(:)) ./ 1e3;    % meters -> km
        raw{nn}.tau   = double(p.tau(:));                % unitless
        raw{nn}.alpha = double(alpha_all{nn}(:));        % unitless

        raw{nn}.lwp_insitu = double(p.lwp(1));           % g/m^2

        % mid-profile geolocation and time
        mid = round(length(p.latitude)/2);
        raw{nn}.latitude       = double(p.latitude(mid));
        raw{nn}.longitude      = double(p.longitude(mid));
        raw{nn}.time_utc_hours = double(p.time(mid)) / 3600;   % seconds -> decimal hours UTC
        raw{nn}.dateOfFlight   = p.dateOfFlight;

    end

    oracles_out = package_campaign(raw, ensemble_profiles, 'oracles', ...
        [folderpath_oracles, oracles_filename], inputs, which_computer);

    save_campaign_file(oracles_out, output_folder, ...
        ['insitu_era5_synthetic_inputs_oracles_', date_str, '.mat'])

end

%% --- Combined file -----------------------------------------------------------
% One file holding both campaigns, with a 'campaign' label per profile. This is
% the file the synthetic profile generator should read.

if process_vocals==true && process_oracles==true

    combined = merge_campaigns(vocals_out, oracles_out);

    save_campaign_file(combined, output_folder, ...
        ['insitu_era5_synthetic_inputs_combined_', date_str, '.mat'])

end

fprintf('\nDone. Files written to:\n  %s\n', output_folder)

%% =============================================================================
%  Local functions
%  =============================================================================

function ensemble_profiles = truncate_profiles(ensemble_profiles, inputs)
% Keep only the first max_profiles_per_campaign profiles (used for smoke tests).

if ~isempty(inputs.max_profiles_per_campaign)

    n_keep = min(length(ensemble_profiles), inputs.max_profiles_per_campaign);
    ensemble_profiles = ensemble_profiles(1:n_keep);

    fprintf(['*** inputs.max_profiles_per_campaign is set: only the first %d ',...
        'profiles will be processed ***\n'], n_keep)

end

end



function out = package_campaign(raw, ensemble_profiles, campaign_name, ...
    source_file, inputs, which_computer)
% Clean every in-cloud profile, pair it with ERA5, and collect the results.

n_total = length(raw);

% --- Storage (cells first; trimmed to the profiles that survive) ---
keep_idx   = false(n_total, 1);
skip_idx   = [];
skip_why   = {};

re_c    = cell(n_total, 1);
lwc_c   = cell(n_total, 1);
z_c     = cell(n_total, 1);
alpha_c = cell(n_total, 1);
tau_c_p = cell(n_total, 1);

n_levels     = nan(n_total, 1);
z_top        = nan(n_total, 1);
z_base       = nan(n_total, 1);
thickness    = nan(n_total, 1);
tau_c          = nan(n_total, 1);
tau_c_computed = nan(n_total, 1);
lwp_insitu     = nan(n_total, 1);
lwp_computed   = nan(n_total, 1);

latitude       = nan(n_total, 1);
longitude      = nan(n_total, 1);
time_utc_hours = nan(n_total, 1);
day_of_year    = nan(n_total, 1);
date_of_flight = cell(n_total, 1);
re_source      = cell(n_total, 1);

era5_c   = cell(n_total, 1);
era5_lat = nan(n_total, 1);
era5_lon = nan(n_total, 1);
era5_filename    = cell(n_total, 1);
era5_utcTime     = cell(n_total, 1);
era5_timeDiff_min = nan(n_total, 1);
era5_dist_m       = nan(n_total, 1);

t_start = tic;

for nn = 1:n_total

    % -----------------------------------------------------------------
    % --- Clean and orient the in-cloud profile -----------------------
    % -----------------------------------------------------------------
    [prof, why] = clean_insitu_profile(raw{nn}, inputs);

    if isempty(prof)
        skip_idx(end+1, 1) = nn;                                     %#ok<AGROW>
        skip_why{end+1, 1} = why;                                    %#ok<AGROW>
        continue
    end

    % -----------------------------------------------------------------
    % --- Pair with the closest ERA5 profile in space and time --------
    % -----------------------------------------------------------------
    try

        if strcmp(campaign_name, 'vocals-rex')==true
            era5 = findClosestProfile_ERA5_VOCALS_REx(ensemble_profiles{nn}, ...
                inputs.print_era5_status, which_computer);
        else
            era5 = findClosestProfile_ERA5_ORACLES(ensemble_profiles{nn}, ...
                inputs.print_era5_status, which_computer);
        end

    catch ME
        skip_idx(end+1, 1) = nn;                                     %#ok<AGROW>
        skip_why{end+1, 1} = ['no ERA5 match: ', ME.message];        %#ok<AGROW>
        continue
    end

    era5_prof = era5_atmospheric_profile(era5);

    if any(isnan(era5_prof.T_K)) || any(isnan(era5_prof.vapor_concentration_cm3))
        skip_idx(end+1, 1) = nn;                                     %#ok<AGROW>
        skip_why{end+1, 1} = 'ERA5 profile contains NaNs';            %#ok<AGROW>
        continue
    end

    % -----------------------------------------------------------------
    % --- Store --------------------------------------------------------
    % -----------------------------------------------------------------
    keep_idx(nn) = true;

    re_c{nn}    = prof.re;
    lwc_c{nn}   = prof.lwc;
    z_c{nn}     = prof.z;
    alpha_c{nn} = prof.alpha;
    tau_c_p{nn} = prof.tau;

    n_levels(nn)     = numel(prof.z);
    z_top(nn)        = prof.z(1);
    z_base(nn)       = prof.z(end);
    thickness(nn)    = prof.z(1) - prof.z(end);
    tau_c(nn)        = max(prof.tau);
    lwp_insitu(nn)   = raw{nn}.lwp_insitu;
    lwp_computed(nn) = 1e3 * abs(trapz(prof.z, prof.lwc));    % g/m^3 * km -> g/m^2

    % Optical thickness recomputed from the measured r_e and LWC profiles:
    %   dtau = (3/2) * Q_ext * LWC[kg/m^3] * dz[m] / (rho_w[kg/m^3] * r_e[m])
    %        = 750 * Q_ext * LWC[g/m^3] * dz[km] / r_e[microns]
    % A few of the stored in-situ tau profiles are corrupt (tau ~ 1e-4 for a
    % cloud with LWP > 50 g/m^2), so this gives a self-consistent alternative.
    tau_c_computed(nn) = 750 * inputs.qext_eff * ...
        abs(trapz(prof.z, prof.lwc ./ prof.re));

    latitude(nn)       = raw{nn}.latitude;
    longitude(nn)      = raw{nn}.longitude;
    time_utc_hours(nn) = raw{nn}.time_utc_hours;
    date_of_flight{nn} = char(raw{nn}.dateOfFlight, 'yyyy-MM-dd');
    day_of_year(nn)    = day(raw{nn}.dateOfFlight, 'dayofyear');
    re_source{nn}      = raw{nn}.re_source;

    era5_c{nn}            = era5_prof;
    era5_lat(nn)          = double(era5.geo.Latitude);
    era5_lon(nn)          = double(era5.geo.Longitude);
    era5_filename{nn}     = era5.metadata.era5_filename;
    era5_utcTime{nn}      = char(era5.metadata.era5_utcTime, 'yyyy-MM-dd HH:mm:ss');
    era5_timeDiff_min(nn) = minutes(era5.metadata.min_timeDiff);

    if isfield(era5.metadata, 'dist_btwn_era5_and_vocalsRex_m')==true
        era5_dist_m(nn) = era5.metadata.dist_btwn_era5_and_vocalsRex_m;
    else
        era5_dist_m(nn) = era5.metadata.dist_btwn_era5_and_oracles_m;
    end

    if mod(nn, 10)==0 || nn==n_total
        fprintf('  processed %d / %d profiles (%.1f s elapsed)\n', nn, n_total, toc(t_start))
    end

end

% -----------------------------------------------------------------
% --- Trim to the profiles that survived --------------------------
% -----------------------------------------------------------------
kk = find(keep_idx);

if isempty(kk)
    error([newline, 'No ', campaign_name, ' profiles survived cleaning and ERA5 pairing.', newline])
end

out.campaign     = campaign_name;
out.source_file  = source_file;
out.date_created = char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss'));
out.inputs       = inputs;

out.profile_idx  = kk;                        % index into the original ensemble_profiles

out.re    = re_c(kk);
out.lwc   = lwc_c(kk);
out.z     = z_c(kk);
out.alpha = alpha_c(kk);
out.tau   = tau_c_p(kk);

out.n_levels     = n_levels(kk);
out.z_top        = z_top(kk);
out.z_base       = z_base(kk);
out.thickness    = thickness(kk);
out.tau_c          = tau_c(kk);          % max of the stored in-situ tau profile
out.tau_c_computed = tau_c_computed(kk); % recomputed from r_e and LWC
out.lwp_insitu     = lwp_insitu(kk);
out.lwp_computed   = lwp_computed(kk);

out.latitude       = latitude(kk);
out.longitude      = longitude(kk);
out.time_utc_hours = time_utc_hours(kk);
out.date_of_flight = date_of_flight(kk);
out.day_of_year    = day_of_year(kk);
out.re_source      = re_source(kk);      % 're' (CDP + 2DC / CAS + 2DS-HVPS) or 're_CDP'

% --- ERA5 arrays (n_profiles x n_pressure_levels) ---
n_keep = length(kk);
n_lev  = length(era5_c{kk(1)}.pressure_hPa);

out.era5_pressure_hPa            = nan(n_keep, n_lev);
out.era5_T_K                     = nan(n_keep, n_lev);
out.era5_q_kgkg                  = nan(n_keep, n_lev);
out.era5_vapor_concentration_cm3 = nan(n_keep, n_lev);
out.era5_vapor_massDensity_kgm3  = nan(n_keep, n_lev);
out.era5_GP_height_m             = nan(n_keep, n_lev);

for ii = 1:n_keep

    e = era5_c{kk(ii)};

    if length(e.pressure_hPa) ~= n_lev
        error([newline, 'ERA5 profile ', num2str(kk(ii)), ' has ', ...
            num2str(length(e.pressure_hPa)), ' levels, expected ', num2str(n_lev), newline])
    end

    out.era5_pressure_hPa(ii, :)            = e.pressure_hPa(:)';
    out.era5_T_K(ii, :)                     = e.T_K(:)';
    out.era5_q_kgkg(ii, :)                  = e.q_kgkg(:)';
    out.era5_vapor_concentration_cm3(ii, :) = e.vapor_concentration_cm3(:)';
    out.era5_vapor_massDensity_kgm3(ii, :)  = e.vapor_massDensity_kgm3(:)';
    out.era5_GP_height_m(ii, :)             = e.GP_height_m(:)';

end

out.era5_latitude          = era5_lat(kk);
out.era5_longitude         = era5_lon(kk);
out.era5_filename          = era5_filename(kk);
out.era5_utcTime           = era5_utcTime(kk);
out.era5_timeDiff_minutes  = era5_timeDiff_min(kk);
out.era5_distance_m        = era5_dist_m(kk);

% --- Record what was thrown out and why ---
out.skipped_profile_idx = skip_idx;
out.skipped_reason      = skip_why;

fprintf('\n  Kept %d of %d %s profiles (%d skipped)\n', ...
    n_keep, n_total, campaign_name, length(skip_idx))

if ~isempty(skip_idx)
    [uniq_why, ~, idx_why] = unique(skip_why);
    counts = accumarray(idx_why, 1);
    for ii = 1:length(uniq_why)
        fprintf('    %3d x  %s\n', counts(ii), uniq_why{ii})
    end
end

% --- Flag profiles whose stored tau disagrees with the recomputed tau ---
tau_ratio  = out.tau_c ./ out.tau_c_computed;
idx_suspect = find(tau_ratio < 0.5 | tau_ratio > 2 | ~isfinite(tau_ratio));

if ~isempty(idx_suspect)
    fprintf(['\n  Warning: %d profile(s) have a stored optical thickness that differs\n',...
        '  from the value recomputed from r_e and LWC by more than a factor of 2.\n',...
        '  Use tau_c_computed for these. Profile indices:\n'], length(idx_suspect))
    fprintf('    %s\n', num2str(out.profile_idx(idx_suspect)'))
end

end



function [prof, why] = clean_insitu_profile(raw, inputs)
% Orient a raw in-cloud profile cloud-top-first, remove redundant altitudes and
% out-of-range droplet sizes, and fill any gaps in the gamma shape parameter.
% Returns prof = [] if the profile can't be used.

prof = [];
why  = '';

re    = raw.re;
lwc   = raw.lwc;
z     = raw.z;
tau   = raw.tau;
alpha = raw.alpha;

% --- All five vectors must describe the same levels ---
n_use = min([length(re), length(lwc), length(z), length(tau), length(alpha)]);

if n_use < inputs.min_num_levels
    why = 'fewer levels than min_num_levels';
    return
end

re = re(1:n_use);   lwc = lwc(1:n_use);   z = z(1:n_use);
tau = tau(1:n_use); alpha = alpha(1:n_use);

% --- Order every vector so index 1 is cloud top ---
if (z(2) - z(1)) > 0
    % profile was sampled bottom-up
    z = flipud(z);   lwc = flipud(lwc);   re = flipud(re);
    tau = flipud(tau);   alpha = flipud(alpha);
end

% --- Drop repeated altitudes (they break the interpolation onto a common grid) ---
[~, idx_unique] = unique(z, 'stable');

if length(idx_unique) < length(z)
    idx_unique = sort(idx_unique);
    z = z(idx_unique);   lwc = lwc(idx_unique);   re = re(idx_unique);
    tau = tau(idx_unique);   alpha = alpha(idx_unique);
end

% --- Fill gaps in the gamma shape parameter ---
alpha = fill_alpha_gaps(alpha, z, inputs);

% --- Remove levels with droplets outside the usable size range ---
idx_remove = false(size(re));

if inputs.remove_levels_re_le_1um==true
    idx_remove = idx_remove | (re <= 1);
end

if inputs.remove_levels_re_ge_50um==true
    idx_remove = idx_remove | (re >= 50);
end

% also drop anything non-finite in the core profile variables
idx_remove = idx_remove | ~isfinite(re) | ~isfinite(lwc) | ~isfinite(z);

if any(idx_remove)
    z(idx_remove) = [];   lwc(idx_remove) = [];   re(idx_remove) = [];
    tau(idx_remove) = [];   alpha(idx_remove) = [];
end

if length(z) < inputs.min_num_levels
    why = 'too few levels after cleaning';
    return
end

if (z(1) - z(end)) <= 0
    why = 'zero or negative cloud thickness';
    return
end

if all(lwc <= 0)
    why = 'no positive liquid water content';
    return
end

prof.re    = re;
prof.lwc   = lwc;
prof.z     = z;
prof.tau   = tau;
prof.alpha = alpha;

end



function alpha = fill_alpha_gaps(alpha, z, inputs)
% Interpolate (and extrapolate) over NaN / non-positive gamma shape parameters,
% then clamp to the range libRadtran's Mie tables can handle.

bad = ~isfinite(alpha) | alpha <= 0;

if sum(~bad) == 0
    alpha = repmat(inputs.alpha_default, size(alpha));
    return

elseif sum(~bad) == 1
    % a single good value: interp1 needs at least two sample points
    alpha(bad) = alpha(~bad);
    return

end

if any(bad)
    alpha(bad) = interp1(z(~bad), alpha(~bad), z(bad), 'linear', 'extrap');
end

% extrapolation can still push alpha out of range
alpha(~isfinite(alpha) | alpha <= 0) = inputs.alpha_default;
alpha(alpha > inputs.alpha_max)      = inputs.alpha_max;

end



function e = era5_atmospheric_profile(era5)
% Convert the ERA5 profile returned by findClosestProfile_ERA5_* into the same
% quantities write_ERA5_radiosonde_DAT_with_multiPixels stores in
% era5.datProfiles, without writing a libRadtran radiosonde file.
%
% All output vectors are column vectors ordered surface -> TOA.

pressure_hPa = double(era5.pressure(:));      % hPa
T_K          = double(era5.temperature(:));   % K
q_kgkg       = double(era5.specificHumidity(:));   % kg/kg

% --- Geopotential [m^2/s^2] -> geopotential height [m] -> geometric height [m] ---
Z_geopotential_m = double(era5.geopotential(:)) ./ 9.80665;

con = physical_constants;
altitude_m = con.R_earth .* Z_geopotential_m ./ (con.R_earth - Z_geopotential_m);

% --- Specific humidity -> vapor mass density [kg/m^3] and number density [molec/cm^3] ---
use_virtual_temp = true;    % ~1% more accurate than the dry-air approximation

[rho_v_kgm3, vapor_concentration_cm3] = specificHumidity2waterVaporDensity(...
    q_kgkg, T_K, pressure_hPa .* 100, use_virtual_temp);

% --- Sort surface -> TOA (increasing altitude) ---
[altitude_m, idx_sort] = sort(altitude_m, 'ascend');

e.GP_height_m             = altitude_m;
e.T_K                     = T_K(idx_sort);
e.pressure_hPa            = pressure_hPa(idx_sort);
e.q_kgkg                  = q_kgkg(idx_sort);
e.vapor_concentration_cm3 = vapor_concentration_cm3(idx_sort);
e.vapor_massDensity_kgm3  = rho_v_kgm3(idx_sort);

% --- Total precipitable water [kg/m^2 = mm], for reference ---
e.total_precipitable_water_mm = trapz(e.GP_height_m, e.vapor_massDensity_kgm3);

end



function alpha_all = load_or_compute_oracles_alpha(ensemble_profiles, cache_file, ...
    source_filename, inputs)
% ORACLES profiles don't carry a stored gamma fit, so fit the measured droplet
% size distribution at each in-cloud altitude. Cached because it's slow.

n_total = length(ensemble_profiles);

% --- Reuse the cache when it matches this ensemble file ---
if inputs.recompute_oracles_alpha==false && isfile(cache_file)==true

    ds_cache = load(cache_file);

    if isfield(ds_cache, 'alpha_all')==true && length(ds_cache.alpha_all)==n_total

        fprintf('Using cached ORACLES gamma fits:\n  %s\n', cache_file)
        alpha_all = ds_cache.alpha_all;
        return

    else
        fprintf('Cached ORACLES gamma fits do not match this ensemble file. Refitting...\n')
    end

end

fprintf(['Fitting a gamma distribution to the droplet size distribution at every\n',...
    'in-cloud altitude of %d ORACLES profiles. This takes a while...\n'], n_total)

alpha_all = cell(n_total, 1);
t_start = tic;

for nn = 1:n_total

    p = ensemble_profiles{nn};

    try

        [~, ~, gammaFit] = find_bestFitDist_dropDist(p.Nd, ...
            p.drop_radius_bin_edges, p.drop_radius_bin_center, ...
            inputs.gammaFit_significance_lvl);

        alpha_all{nn} = double(gammaFit.alpha(:));

    catch ME

        warning(['Gamma fit failed for ORACLES profile ', num2str(nn), ': ', ...
            ME.message, '. Using alpha = ', num2str(inputs.alpha_default), '.'])

        alpha_all{nn} = repmat(inputs.alpha_default, length(p.altitude), 1);

    end

    if mod(nn, 10)==0 || nn==n_total
        fprintf('  fit %d / %d profiles (%.1f s elapsed)\n', nn, n_total, toc(t_start))
    end

end

% --- Cache for next time (only when the full ensemble was fit) ---
if isempty(inputs.max_profiles_per_campaign)
    save(cache_file, 'alpha_all', 'source_filename', '-v7')
    fprintf('Cached ORACLES gamma fits to:\n  %s\n', cache_file)
end

end



function combined = merge_campaigns(a, b)
% Stack two campaign structures into one, tagging each profile with its campaign.

combined.campaign     = 'vocals-rex + oracles';
combined.source_file  = {a.source_file; b.source_file};
combined.date_created = char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss'));
combined.inputs       = a.inputs;

n_a = length(a.profile_idx);
n_b = length(b.profile_idx);

combined.campaign_each = [repmat({a.campaign}, n_a, 1); repmat({b.campaign}, n_b, 1)];
combined.profile_idx   = [a.profile_idx; b.profile_idx];

% --- Cell array fields (one entry per profile) ---
cell_fields = {'re', 'lwc', 'z', 'alpha', 'tau', 'date_of_flight', 're_source', ...
    'era5_filename', 'era5_utcTime'};

for ff = 1:length(cell_fields)
    combined.(cell_fields{ff}) = [a.(cell_fields{ff}); b.(cell_fields{ff})];
end

% --- Column vector fields ---
vec_fields = {'n_levels', 'z_top', 'z_base', 'thickness', 'tau_c', 'tau_c_computed', ...
    'lwp_insitu', 'lwp_computed', 'latitude', 'longitude', 'time_utc_hours', ...
    'day_of_year', 'era5_latitude', 'era5_longitude', 'era5_timeDiff_minutes', ...
    'era5_distance_m'};

for ff = 1:length(vec_fields)
    combined.(vec_fields{ff}) = [a.(vec_fields{ff}); b.(vec_fields{ff})];
end

% --- ERA5 matrix fields (n_profiles x n_pressure_levels) ---
mat_fields = {'era5_pressure_hPa', 'era5_T_K', 'era5_q_kgkg', ...
    'era5_vapor_concentration_cm3', 'era5_vapor_massDensity_kgm3', 'era5_GP_height_m'};

for ff = 1:length(mat_fields)

    if size(a.(mat_fields{ff}), 2) ~= size(b.(mat_fields{ff}), 2)
        error([newline, 'The two campaigns have different numbers of ERA5 pressure ',...
            'levels. They cannot be stacked into one array.', newline])
    end

    combined.(mat_fields{ff}) = [a.(mat_fields{ff}); b.(mat_fields{ff})];

end

end



function save_campaign_file(out, output_folder, filename)
% Save every field of the campaign structure as a top-level variable, using the
% v7 MAT format so scipy.io.loadmat can read the file directly.

full_path = [output_folder, filename];

% '-struct' saves each field of 'out' as its own top-level variable
save(full_path, '-struct', 'out', '-v7')

d = dir(full_path);
fprintf('  saved %s (%.1f MB)\n', filename, d.bytes/1e6)

end
