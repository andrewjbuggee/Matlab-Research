function hysics_refl_smoothingTest_from_insitu(folder_paths, measurement_idx, sza, ...
    output_dir, campaign, smoothing_windows_m)
%% Test the retrievability of vertical droplet-size structure with HySICS
%
% For a single in-situ measured droplet profile, this computes the
% top-of-atmosphere HySICS reflectance spectrum (636 channels, 351-2297 nm)
% at a fixed nadir geometry for:
%
%   version 1  : the RAW measured effective-radius profile r_e(z)
%   version 2+ : vertically SMOOTHED versions of r_e(z), one per entry in
%                smoothing_windows_m (moving average, window width in meters)
%
% Everything else about the cloud is held identical across versions:
%   - the cloud optical depth tau_c is FORCED to the measured value for every
%     version (libRadtran 'wc_modify tau set'), so tau is held constant;
%   - the measured liquid-water-content SHAPE lwc(z) is kept (its magnitude is
%     rescaled internally by libRadtran to hit tau_c, so the per-level LWC is
%     free to follow the droplet size, as intended);
%   - the gamma-distribution width (Mie table), atmosphere (ERA5), aerosols,
%     surface, trace gases, and geometry are identical for all versions.
%
% The ONLY thing that changes between versions is the vertical shape of
% r_e(z). Comparing the resulting spectra against the HySICS measurement
% uncertainty (0.3%) shows whether vertical droplet-size fluctuations are
% detectable -- a sufficient (not necessary) condition for retrievability.
%
% INPUTS:
%   folder_paths         - struct from define_folderPaths_for_HySICS
%   measurement_idx      - integer index of the ensemble profile to process
%   sza                  - solar zenith angle [degrees]
%   output_dir           - full path (trailing slash) where the .mat is saved
%   campaign             - 'vocalsRex' or 'oracles'
%   smoothing_windows_m  - (optional) vector of moving-average window widths
%                          [meters]. Defaults to define_re_smoothing_windows().
%
% Note: the calling script is responsible for addLibRadTran_paths and
% define_folderPaths_for_HySICS (do not call them here).
%
% By Andrew John Buggee


%% Defaults / setup

if nargin < 6 || isempty(smoothing_windows_m)
    smoothing_windows_m = define_re_smoothing_windows();
end
smoothing_windows_m = smoothing_windows_m(:)';      % row vector

% number of profile versions = raw + one per smoothing window
num_versions = 1 + numel(smoothing_windows_m);

which_computer = folder_paths.which_computer;


%% Delete old temporary files

delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])
delete([folder_paths.atm_folder_path, '*-aboveCloud.DAT'])
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])


%% Start parallel pool

start_parallel_pool(which_computer)


%% ------------------------------------------------------------------------
%  Load the in-situ ensemble profile for the requested campaign
%  ------------------------------------------------------------------------

if strcmp(campaign, 'vocalsRex')

    campaign_name = 'vocalsRex';

    if strcmp(which_computer, 'anbu8374')
        folderpath_air = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
            'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];
    elseif strcmp(which_computer, 'andrewbuggee')
        folderpath_air = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
            'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];
    elseif strcmp(which_computer, 'curc')
        folderpath_air = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
            'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];
    end

    saved_profiles_filename = ['ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25', ...
        '_drizzleLWP-threshold_5_10-Nov-2025.mat'];

    ds = load([folderpath_air, saved_profiles_filename]);
    ensemble_profiles = ds.ensemble_profiles;

elseif strcmp(campaign, 'oracles')

    campaign_name = 'oracles';

    if strcmp(which_computer, 'anbu8374')
        folderpath_air = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
            'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];
    elseif strcmp(which_computer, 'andrewbuggee')
        folderpath_air = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
            'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];
    elseif strcmp(which_computer, 'curc')
        folderpath_air = ['/projects/anbu8374/Matlab-Research/', ...
            'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];
    end

    saved_profiles_filename = ['ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10', ...
        '_no_rEff_greaterThan50_microns_16-Mar-2026.mat'];

    ds = load([folderpath_air, saved_profiles_filename]);
    ensemble_profiles = ds.ensemble_profiles;

else
    error([newline, 'campaign must be ''vocalsRex'' or ''oracles''.', newline])
end

% validate the index
total_measurements = numel(ensemble_profiles);
if measurement_idx < 1 || measurement_idx > total_measurements
    error('measurement_idx must be between 1 and %d', total_measurements);
end
fprintf('\n Processing %s measurement %d of %d\n', campaign_name, measurement_idx, total_measurements);

profile = ensemble_profiles{measurement_idx};


%% ------------------------------------------------------------------------
%  Extract and clean the profile: re, lwc, z, tau, alpha (aligned by level)
%  ------------------------------------------------------------------------

% --- effective radius [microns] ---
if isfield(profile, 're')
    re_um = profile.re(:);
elseif isfield(profile, 're_CDP')
    re_um = profile.re_CDP(:);
else
    error([newline, 'No effective-radius field (re / re_CDP) on this profile.', newline])
end

% --- liquid water content [g/m^3], altitude [m], cumulative optical depth ---
lwc_gm3 = profile.lwc(:);
z_m     = profile.altitude(:);
tau_cum = profile.tau(:);

% --- gamma-distribution width (alpha) per level, campaign-specific source ---
if strcmp(campaign, 'vocalsRex')
    alpha_prof = profile.gammaFit.alpha(:);
else
    % ORACLES: fit the per-level droplet distribution to get alpha
    significance_lvl = 0.1;
    [~, ~, gammaFit] = find_bestFitDist_dropDist(profile.Nd, ...
        profile.drop_radius_bin_edges, profile.drop_radius_bin_center, significance_lvl);
    alpha_prof = gammaFit.alpha(:);
end

% make sure every per-level vector is the same length before we align them
n_lvl = numel(re_um);
if numel(lwc_gm3) ~= n_lvl || numel(z_m) ~= n_lvl || numel(tau_cum) ~= n_lvl
    error([newline, 're, lwc, z, and tau profiles are not the same length.', newline])
end
if numel(alpha_prof) ~= n_lvl
    % alpha occasionally differs in length; pad/trim conservatively
    alpha_prof = interp1((1:numel(alpha_prof))', alpha_prof, ...
        linspace(1, numel(alpha_prof), n_lvl)', 'linear', 'extrap');
end

% --- orient all vectors cloud-top-first (index 1 = cloud top) ---
if (z_m(2) - z_m(1)) > 0
    z_m        = flipud(z_m);
    re_um      = flipud(re_um);
    lwc_gm3    = flipud(lwc_gm3);
    tau_cum    = flipud(tau_cum);
    alpha_prof = flipud(alpha_prof);
end

% --- interpolate over any NaN alpha values, guard against non-positive ---
idx_nan = isnan(alpha_prof);
if any(~idx_nan)
    alpha_prof(idx_nan) = interp1(z_m(~idx_nan), alpha_prof(~idx_nan), ...
        z_m(idx_nan), 'linear', 'extrap');
end
alpha_prof(alpha_prof <= 0 | isnan(alpha_prof)) = 7;

% --- drop levels outside the custom Mie table / cloud-edge artifacts ---
%     custom Mie tables span 1-50 microns
drop = isnan(re_um) | re_um >= 50 | re_um < 1;
re_um(drop)      = [];
lwc_gm3(drop)    = [];
z_m(drop)        = [];
tau_cum(drop)    = [];
alpha_prof(drop) = [];

% --- remove duplicate altitude levels (libRadtran needs unique z) ---
[~, iu] = unique(z_m, 'stable');
if numel(iu) < numel(z_m)
    re_um      = re_um(iu);
    lwc_gm3    = lwc_gm3(iu);
    z_m        = z_m(iu);
    tau_cum    = tau_cum(iu);
    alpha_prof = alpha_prof(iu);
end

if numel(re_um) < 3
    warning([newline, 'Profile %d (%s) has fewer than 3 usable cloud levels; skipping.', newline], ...
        measurement_idx, campaign_name);
    return
end

% --- cloud optical depth held fixed for every version ---
tau_c = max(tau_cum);

% altitude in km for the wc-file writer
z_km = z_m ./ 1e3;

% flight metadata for file naming
date_of_flight = profile.dateOfFlight;
if strcmp(campaign, 'vocalsRex')
    t_vec = profile.time_utc(:);
    time_of_flight = t_vec(round(numel(t_vec)/2));            % UTC hours
else
    t_vec = profile.time(:);
    time_of_flight = t_vec(round(numel(t_vec)/2)) / 3600;     % s -> UTC hours
end


%% ------------------------------------------------------------------------
%  ERA5 atmosphere closest to this profile
%  ------------------------------------------------------------------------

if strcmp(campaign, 'vocalsRex')
    era5 = findClosestProfile_ERA5_VOCALS_REx(profile, true, which_computer);
else
    era5 = findClosestProfile_ERA5_ORACLES(profile, true, which_computer);
end

clear ds ensemble_profiles profile


%% ------------------------------------------------------------------------
%  Define the libRadtran / HySICS inputs (matches the in-situ pipelines)
%  ------------------------------------------------------------------------

inputs.which_computer = which_computer;
print_libRadtran_err  = false;

% vertical profile / calculation type
inputs.RT.vert_homogeneous_str    = 'vert-non-homogeneous';
inputs.RT.num_re_parameters       = 2;
inputs.calc_type                  = 'simulated_spectra';
inputs.compute_weighting_functions = false;

% RTE solver
inputs.RT.rte_solver              = 'disort';
inputs.RT.num_streams             = 16;
inputs.RT.use_nakajima_phaseCorrection = true;

% source file (0.1 nm sampling)
inputs.RT.source_file             = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
inputs.RT.source_file_resolution  = 0.1;     % nm

% all 636 HySICS channels (351 - 2297 nm)
inputs.bands2run = (1:1:636)';

inputs.RT.monochromatic_calc          = false;
inputs.RT.compute_reflectivity_uvSpec = false;

% spectral response functions
spec_response = create_HySICS_specResponse(inputs.bands2run, inputs.RT.source_file, ...
    inputs.which_computer);

inputs.RT.wavelengths2run = zeros(length(inputs.bands2run), 2);
for ww = 1:length(inputs.bands2run)
    inputs.RT.wavelengths2run(ww, :) = [spec_response.wavelength(ww, 1), ...
        spec_response.wavelength(ww, end)];
end

% band model, atmosphere
inputs.RT.band_parameterization = 'reptran coarse';
inputs.RT.atm_file              = 'afglus.dat';

% ERA5 radiosonde (T, P, RH)
inputs.RT.use_radiosonde_file = true;
inputs.RT.radiosonde_num_vars = 3;
[inputs.RT.radiosonde_file, era5] = write_ERA5_radiosonde_DAT_with_multiPixels(era5, ...
    folder_paths, 1, [], inputs.RT.radiosonde_num_vars, [], inputs.RT.atm_file, true);

% surface
inputs.RT.surface_albedo = 0.04;             % ocean

% day of year from the actual flight date (Earth-Sun distance; cancels in dR)
inputs.RT.day_of_year = day(date_of_flight, 'dayofyear');

% cloud switches
inputs.RT.yesCloud               = true;
inputs.RT.modify_wc_opticalDepth = true;     % *** forces tau = tau_c ***
inputs.RT.parameterization_str   = 'mie';
inputs.RT.lambda_forTau          = 500;      % nm

% atmospheric grid
inputs.RT.define_atm_grid = false;

% ---- geometry: SINGLE NADIR VIEW, fixed SZA ----
inputs.RT.sensor_altitude = 'toa';
vza  = 0;        % nadir viewing zenith angle [deg]
vaz  = 0;        % viewing azimuth [deg] (degenerate at nadir)
phi0 = 0;        % solar azimuth [deg]
% sza is supplied as an input argument

% Cox-Munk ocean surface
inputs.RT.use_coxMunk  = true;
inputs.RT.wind_speed   = 3;                  % m/s

% cross-section models
inputs.RT.specify_cross_section_model = true;
inputs.RT.crs_model_rayleigh          = 'Bodhaine29';

% aerosols
inputs.RT.yesAerosols         = true;
inputs.RT.aerosol_type         = 4;          % maritime
inputs.RT.aerosol_opticalDepth = 0.1;

% water vapor
inputs.RT.modify_total_columnWaterVapor      = false;
inputs.RT.modify_aboveCloud_columnWaterVapor = false;
inputs.RT.waterVapor_column                  = 20;   % mm

% trace gases
inputs.RT.modify_CO2       = true;
inputs.RT.CO2_mixing_ratio = 418;            % ppm
inputs.RT.modify_N2  = false;  inputs.RT.N2_mixing_ratio  = 0;
inputs.RT.modify_NO2 = false;  inputs.RT.NO2_mixing_ratio = 0;
inputs.RT.modify_O2  = false;  inputs.RT.O2_mixing_ratio  = 0;
inputs.RT.modify_O3  = false;  inputs.RT.O3_mixing_ratio  = 0;

% scattering switches
inputs.RT.no_molecular_abs  = false;
inputs.RT.no_scattering_mol = false;
inputs.RT.no_scattering_aer = false;

% verbosity
if print_libRadtran_err
    inputs.RT.errMsg = 'verbose';
else
    inputs.RT.errMsg = 'quiet';
end

% wc parameterization
inputs.RT.wc_parameterization = 'mie interpolate';


%% ------------------------------------------------------------------------
%  Custom Mie table from the (raw) alpha profile -- shared by all versions
%  ------------------------------------------------------------------------
% Smoothing r_e does NOT change the gamma-distribution width, so the same Mie
% table is correct for every version. The table is interpolated over r_e
% inside libRadtran.

use_35_or_50 = 50;
out = find_custom_mieTable_closest_to_alpha_profile(alpha_prof, use_35_or_50, which_computer);

inputs.RT.distribution_var                       = alpha_prof;
inputs.RT.mean_distribution_var                  = mean(alpha_prof);
inputs.RT.mean_distribution_var_closest_filename = out.mie_table_filename_closest_to_mean;
inputs.RT.use_custom_mie_calcs                   = true;


%% ------------------------------------------------------------------------
%  Build the raw + smoothed r_e profiles and write one wc file per version
%  ------------------------------------------------------------------------
% Only r_e changes between versions; lwc(z) and z are identical. libRadtran
% rescales the cloud so tau = tau_c for every version.

re_versions = cell(num_versions, 1);
wc_filename = cell(num_versions, 1);

for vv = 1:num_versions

    if vv == 1
        window_m = 0;                          % raw
    else
        window_m = smoothing_windows_m(vv-1);
    end

    re_v = smooth_re_profile(re_um, z_m, window_m);
    re_versions{vv} = re_v;

    % unique file index per version so the wc files don't overwrite
    wc_idx = 100 * measurement_idx + vv;

    wc_filename{vv} = write_wc_file_from_in_situ(re_v, lwc_gm3, z_km, campaign_name, ...
        date_of_flight, time_of_flight, inputs.compute_weighting_functions, ...
        which_computer, wc_idx, folder_paths.libRadtran_water_cloud_files);
end


%% ------------------------------------------------------------------------
%  Build the (wavelength x version) job table
%  ------------------------------------------------------------------------

num_wl        = length(inputs.bands2run);
num_INP_files = num_wl * num_versions;

[WW, VV] = ndgrid(1:num_wl, 1:num_versions);
wl_idx_col  = WW(:);          % which spectral channel
ver_idx_col = VV(:);          % which profile version

% columns: [wl_lower, wl_upper, spec_response_idx, version_idx]
changing_variables = [inputs.RT.wavelengths2run(wl_idx_col, :), wl_idx_col, ver_idx_col];
num_cols = size(changing_variables, 2);


%% ------------------------------------------------------------------------
%  Pre-slice the source flux and spectral response per iteration
%  ------------------------------------------------------------------------

wavelength_vec = [min(inputs.RT.wavelengths2run, [], 'all'), ...
    max(inputs.RT.wavelengths2run, [], 'all')];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);

wl_perturb = inputs.RT.source_file_resolution / 3;   % nm

spec_response_value = spec_response.value;

source_flux_perIter   = cell(num_INP_files, 1);
spec_response_perIter = cell(num_INP_files, 1);
for nn = 1:num_INP_files
    idx_wl = source_wavelength >= (changing_variables(nn, 1) - wl_perturb) & ...
             source_wavelength <= (changing_variables(nn, 2) + wl_perturb);
    source_flux_perIter{nn}   = source_flux(idx_wl);
    spec_response_perIter{nn} = spec_response_value(changing_variables(nn, 3), :);
end


%% ------------------------------------------------------------------------
%  Run radiative transfer (single parfor: write -> run -> read -> compute)
%  ------------------------------------------------------------------------

% minimal broadcast variables
libRadtran_inp       = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path       = folder_paths.libRadtran_water_cloud_files;

compute_reflectivity_uvSpec = inputs.RT.compute_reflectivity_uvSpec;
inputs_parfor.RT.rte_solver = inputs.RT.rte_solver;

date_str = date_of_flight;
time_str = time_of_flight;

inputFileName  = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);
Refl_model     = zeros(num_INP_files, 1);

tic

parfor nn = 1:num_INP_files

    wavelengths = changing_variables(nn, 1:2);
    ver_nn      = changing_variables(nn, num_cols);

    % unique INP / OUT filenames
    inputFileName{nn} = [num2str(mean(wavelengths)), 'nm_', campaign_name, '_', ...
        char(date_str), '_', num2str(time_str), '-UTC_meas', num2str(measurement_idx), ...
        '_ver', num2str(ver_nn), '_sza', num2str(sza), '.INP'];

    outputFileName{nn} = ['OUTPUT_', inputFileName{nn}(1:end-4)];

    % write the INP file (force tau = tau_c; fixed nadir geometry)
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, ...
        inputFileName{nn}, inputs, wavelengths, wc_filename{ver_nn}, [], tau_c, ...
        [], [], vza, sza, vaz, phi0);

    % run uvspec
    runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn}, which_computer);

    % read radiance (mW/nm/m^2/sr)
    [rad_calcs, ~, ~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, inputs_parfor, ...
        compute_reflectivity_uvSpec, vza, vaz);

    % compute reflectance
    [Refl_model(nn), ~] = reflectanceFunction_ver3(rad_calcs, ...
        source_flux_perIter{nn}, spec_response_perIter{nn}, sza, vza, vaz);

    % delete INP/OUT immediately to conserve scratch inodes
    inp_file_path = [libRadtran_inp, inputFileName{nn}];
    out_file_path = [libRadtran_inp, outputFileName{nn}, '.OUT'];
    if isfile(inp_file_path); delete(inp_file_path); end
    if isfile(out_file_path); delete(out_file_path); end

end

fprintf('\n%s profile %d: %d spectra computed in %.1f s\n', ...
    campaign_name, measurement_idx, num_versions, toc);


%% ------------------------------------------------------------------------
%  Reshape, compute differences, and the HySICS uncertainty envelope
%  ------------------------------------------------------------------------
% Refl(:, 1) = raw ; Refl(:, 1+w) = smoothed with smoothing_windows_m(w)

Refl = reshape(Refl_model, num_wl, num_versions);

% center wavelength of each HySICS channel [nm]
wl_nm = mean(spec_response.wavelength, 2);

% reflectance difference of each smoothed version relative to the raw profile
% dR(:, w) = R_raw - R_smoothed_w
dR = Refl(:, 1) - Refl(:, 2:end);

% HySICS measurement uncertainty (1-sigma) as an absolute reflectance envelope
inputs.measurement.uncert_hysics = 0.003;        % fraction of the reflectance
Refl_uncert_hysics = inputs.measurement.uncert_hysics .* Refl;   % per version

% EMIT for reference
inputs.measurement.uncert_emit = 0.04;
Refl_uncert_emit = inputs.measurement.uncert_emit .* Refl;


%% ------------------------------------------------------------------------
%  Save
%  ------------------------------------------------------------------------

inputs.folderpath_2save = output_dir;
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% record geometry and experiment metadata
geometry.sza  = sza;
geometry.vza  = vza;
geometry.vaz  = vaz;
geometry.phi0 = phi0;

filename = [output_dir, 'refl_smoothingTest_', campaign_name, '_meas', num2str(measurement_idx), ...
    '_sza', num2str(round(sza)), '_', char(date_of_flight), '_', num2str(time_of_flight), 'UTC', ...
    '_ran-', char(datetime('today')), '.mat'];

save(filename, ...
    'wl_nm', 'Refl', 'dR', 'Refl_uncert_hysics', 'Refl_uncert_emit', ...
    'smoothing_windows_m', 're_versions', 'z_m', 'lwc_gm3', 'tau_c', ...
    'alpha_prof', 'geometry', 'campaign_name', 'measurement_idx', ...
    'inputs', 'spec_response', 'era5');

fprintf('Saved %s profile %d -> %s\n', campaign_name, measurement_idx, filename);

end
