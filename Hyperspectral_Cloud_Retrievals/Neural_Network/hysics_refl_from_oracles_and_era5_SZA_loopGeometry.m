

function hysics_refl_from_oracles_and_era5_SZA_loopGeometry(folder_paths, measurement_idx, sza, output_dir)
%% Generate simulated HySICS measurements from ORACLES in-situ profiles
%
% Analogous to hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver2.m,
% adapted for the ORACLES P-3 dataset.
%
% The following ORACLES in-situ measured values are used to create simulated
% measurements:
%
%   (1) In-situ effective radius profile       - from CAS / 2DS / HVPS-3
%   (2) In-situ liquid water content profile   - combined probes
%   (3) In-situ cloud optical depth profile    - computed in find_verticalProfiles_ORACLES
%   (4) In-situ cloud top and base altitude
%   (5) ERA5 T, P, and RH - from the closest ERA5 grid point and time step
%
% Spectra are computed for a range of viewing zenith angles (vza), viewing
% azimuth angles (vaz), and solar azimuth angles (phi0) at the user-supplied
% solar zenith angle (sza).  One file is saved per SZA so that many SZA
% values can run in parallel on a supercomputer.
%
% INPUTS:
%   folder_paths    - structure from define_folderPaths_for_HySICS
%   measurement_idx - integer index specifying which ORACLES ensemble profile
%                     to process
%   sza             - solar zenith angle [degrees]
%   output_dir      - full path to the directory where output .mat files
%                     will be saved
%
% Note: Do not call addLibRadTran_paths or define_folderPaths_for_HySICS
% inside this function -- the calling script is responsible for those.

%%  Delete old temporary files

% INP and OUT files
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% water cloud files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% above-cloud water vapor profiles
delete([folder_paths.atm_folder_path, '*-aboveCloud.DAT'])

% Mie INP/OUT files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])


%% Start parallel pool

start_parallel_pool(folder_paths.which_computer)


%% Load the ORACLES ensemble profiles

which_computer = whatComputer();

if strcmp(which_computer, 'anbu8374')

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folderpath_oracles = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

elseif strcmp(which_computer, 'andrewbuggee')

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    folderpath_oracles = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

elseif strcmp(which_computer, 'curc')

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    folderpath_oracles = ['/projects/anbu8374/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

end


% Find the most recent saved ensemble profiles file (with precipitation, which
% includes all profiles).  If you want only non-precipitating profiles, change
% the search pattern accordingly.
ensemble_files = dir([folderpath_oracles, 'ensemble_profiles_with_precip_from_*.mat']);

if isempty(ensemble_files)
    % fall back: look for any ensemble profiles file
    ensemble_files = dir([folderpath_oracles, 'ensemble_profiles_*.mat']);
end

if isempty(ensemble_files)
    error([newline, 'No ORACLES ensemble profiles file found in: ', folderpath_oracles, newline])
end

% Use the most recently modified file
[~, idx_newest] = max([ensemble_files.datenum]);
saved_profiles_filename = ensemble_files(idx_newest).name;

ds_oracles = load([folderpath_oracles, saved_profiles_filename]);

% Field naming: ensemble_profiles is a cell array of profile structs
ensemble_profiles = ds_oracles.ensemble_profiles;

% Identify the campaign
campaign_name = 'oracles';


%% Validate measurement_idx

total_measurements = length(ensemble_profiles);

if measurement_idx < 1 || measurement_idx > total_measurements
    error('measurement_idx must be between 1 and %d', total_measurements);
end

fprintf('\n Processing ORACLES measurement %d of %d\n', measurement_idx, total_measurements);


%% Find the ERA5 vertical profile closest to this ORACLES profile

era5 = findClosestProfile_ERA5_ORACLES(ensemble_profiles{measurement_idx}, ...
    true, which_computer);


%% Unpack folder paths to avoid large broadcast variables in parfor

libRadtran_inp       = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path       = folder_paths.libRadtran_water_cloud_files;
which_computer       = folder_paths.which_computer;

inputs.which_computer = which_computer;


%% Define RT inputs

print_libRadtran_err = false;

% Vertical profile type
inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';
inputs.RT.num_re_parameters    = 2;

% Calculation type
inputs.calc_type               = 'simulated_spectra';
inputs.compute_weighting_functions = false;

% RTE solver
inputs.RT.rte_solver           = 'disort';
inputs.RT.num_streams           = 16;
inputs.RT.use_nakajima_phaseCorrection = true;


% ---- Source file ----
inputs.RT.source_file            = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
inputs.RT.source_file_resolution = 0.1;         % nm

% ---- Wavelength / band selection ----
% All 636 HySICS spectral channels (351 – 2297 nm)
inputs.bands2run = (1:1:636)';

inputs.RT.monochromatic_calc          = false;
inputs.RT.compute_reflectivity_uvSpec = false;


% ---- Spectral response functions ----
spec_response = create_HySICS_specResponse(inputs.bands2run, inputs.RT.source_file, ...
    inputs.which_computer);

inputs.RT.wavelengths2run = zeros(length(inputs.bands2run), 2);

for ww = 1:length(inputs.bands2run)
    inputs.RT.wavelengths2run(ww, :) = [spec_response.wavelength(ww, 1), ...
        spec_response.wavelength(ww, end)];
end

% ---- Band parameterization ----
inputs.RT.band_parameterization = 'reptran coarse';

% ---- Atmospheric data file ----
inputs.RT.atm_file = 'afglus.dat';

% ---- ERA5 radiosonde ----
inputs.RT.use_radiosonde_file  = true;
inputs.RT.radiosonde_num_vars  = 3;

[inputs.RT.radiosonde_file, era5] = write_ERA5_radiosonde_DAT_with_multiPixels(era5, ...
    folder_paths, 1, [], inputs.RT.radiosonde_num_vars, [], ...
    inputs.RT.atm_file, true);

% ---- Surface ----
inputs.RT.surface_albedo = 0.04;            % ocean water albedo

% ---- Day of year (computed from the actual ORACLES flight date) ----
inputs.RT.day_of_year = day(ensemble_profiles{measurement_idx}.dateOfFlight, 'dayofyear');

% ---- Cloud switch ----
inputs.RT.yesCloud             = true;
inputs.RT.modify_wc_opticalDepth = true;
inputs.RT.use_custom_mie_calcs   = false;       % will be overridden below after alpha fit
inputs.RT.parameterization_str   = 'mie';
inputs.RT.lambda_forTau          = 500;         % nm

% ---- Atmospheric grid ----
inputs.RT.define_atm_grid = false;

% ---- Sensor altitude ----
inputs.RT.sensor_altitude = 'toa';

% ---- Solar / viewing geometry ----
% Solar azimuth angles (0-360 deg clockwise from due south)
phi0 = [0, 45, 90, 180];       % degrees

% Viewing zenith angles sampled linearly in cos(vza)
mu_sample = linspace(cosd(0), cosd(65), 8);
vza = acosd(mu_sample);         % degrees

% Viewing azimuth angles (0-360 deg clockwise from due north)
vaz = [0, 45, 90, 180];         % degrees

% ---- Cox-Munk ocean surface ----
inputs.RT.use_coxMunk  = true;
inputs.RT.wind_speed   = 3;             % m/s

% ---- Cross-section models ----
inputs.RT.specify_cross_section_model = true;
inputs.RT.crs_model_rayleigh          = 'Bodhaine29';

% ---- Aerosols ----
inputs.RT.yesAerosols        = true;
inputs.RT.aerosol_type        = 4;       % maritime aerosols
inputs.RT.aerosol_opticalDepth = 0.1;

% ---- Water vapor ----
inputs.RT.modify_total_columnWaterVapor    = false;
inputs.RT.modify_aboveCloud_columnWaterVapor = false;
inputs.RT.waterVapor_column                = 20;   % mm

% ---- Trace gases ----
inputs.RT.modify_CO2         = true;
inputs.RT.CO2_mixing_ratio   = 418;     % ppm

inputs.RT.modify_N2          = false;
inputs.RT.N2_mixing_ratio    = 0;

inputs.RT.modify_NO2         = false;
inputs.RT.NO2_mixing_ratio   = 0;

inputs.RT.modify_O2          = false;
inputs.RT.O2_mixing_ratio    = 0;

inputs.RT.modify_O3          = false;
inputs.RT.O3_mixing_ratio    = 0;

% ---- Scattering switches ----
inputs.RT.no_molecular_abs    = false;
inputs.RT.no_scattering_mol   = false;
inputs.RT.no_scattering_aer   = false;

% ---- libRadtran error message verbosity ----
if print_libRadtran_err
    inputs.RT.errMsg = 'verbose';
else
    inputs.RT.errMsg = 'quiet';
end

% ---- Wc parameterization ----
inputs.RT.wc_parameterization = 'mie interpolate';


%% Count geometries and pre-allocate

num_meas = 1;
num_wl   = length(inputs.bands2run);
num_vza  = length(vza);
num_vaz  = length(vaz);
num_saz  = length(phi0);

num_INP_files = num_wl * num_vza * num_saz * num_vaz;

inputs.calc_type = 'simulated_spectra';

wc_filename    = cell(num_meas, 1);
tau_c          = zeros(num_meas, 1);
date_of_flight = cell(num_meas, 1);
time_of_flight = zeros(num_meas, 1);
re             = cell(num_meas, 1);
lwc            = cell(num_meas, 1);
z              = cell(num_meas, 1);
tau            = cell(num_meas, 1);
alpha_param    = cell(num_meas, 1);


%% Extract the cloud profile for this measurement

nn = 1;

% --- Core profile variables (column vectors) ---
lwc{nn} = ensemble_profiles{measurement_idx}.lwc';         % g/m^3
z{nn}   = ensemble_profiles{measurement_idx}.altitude' ./ 1e3;  % km
tau{nn} = ensemble_profiles{measurement_idx}.tau';

% --- Flight metadata ---
date_of_flight{nn} = ensemble_profiles{measurement_idx}.dateOfFlight;

% Time: seconds since midnight UTC → decimal hours UTC (for file naming)
time_of_flight(nn) = ensemble_profiles{measurement_idx}.time(round(end/2)) / 3600;

% --- Total cloud optical depth ---
tau_c(nn) = ensemble_profiles{measurement_idx}.tau(end);

% --- Effective radius profile ---
% ORACLES profiles store '.re' (combined CAS + 2DS + HVPS effective radius)
re{nn} = ensemble_profiles{measurement_idx}.re';    % microns

% --- Orient all profile vectors so index 1 = cloud top ---
if (z{nn}(2) - z{nn}(1)) > 0
    % Profile was sampled bottom-up; flip to cloud-top-first ordering
    z{nn}   = flipud(z{nn});
    lwc{nn} = flipud(lwc{nn});
    re{nn}  = flipud(re{nn});
    tau{nn} = flipud(tau{nn});
end


%% Fit the droplet size distribution to get the gamma shape parameter

% This is the ORACLES equivalent of the pre-stored gammaFit.alpha from
% VOCALS-REx profiles.  We compute it here because find_verticalProfiles_ORACLES
% does not store the fit result on the profile struct.

significance_lvl = 0.1;

[~, ~, gammaFit] = find_bestFitDist_dropDist( ...
    ensemble_profiles{measurement_idx}.Nd, ...
    ensemble_profiles{measurement_idx}.drop_radius_bin_edges, ...
    ensemble_profiles{measurement_idx}.drop_radius_bin_center, ...
    significance_lvl);

alpha_prof = gammaFit.alpha;            % per-altitude gamma shape parameter
z_prof     = ensemble_profiles{measurement_idx}.altitude;

% Interpolate / extrapolate over any NaN values in the alpha profile
idx_nan = isnan(alpha_prof);
alpha_prof_valid  = alpha_prof;  alpha_prof_valid(idx_nan) = [];
z_prof_valid      = z_prof;      z_prof_valid(idx_nan)     = [];

if isempty(alpha_prof_valid)
    % Edge case: all NaN -- use a sensible default
    new_alpha_prof = repmat(7, size(alpha_prof))';
    warning('All alpha values are NaN for ORACLES profile %d; using default alpha = 7.', measurement_idx);
else
    new_alpha_prof = interp1(z_prof_valid, alpha_prof_valid, z_prof, 'linear', 'extrap')';
end

% Guard against non-positive alpha values
new_alpha_prof(new_alpha_prof <= 0 | isnan(new_alpha_prof)) = 7;

% !! use custom mie tables that span 1-35 microns or 1-50 microns !!
use_35_or_50 = 50;

% Find the closest pre-computed Mie table
out = find_custom_mieTable_closest_to_alpha_profile(new_alpha_prof, use_35_or_50, which_computer);

inputs.RT.distribution_var                    = new_alpha_prof;
inputs.RT.mean_distribution_var               = mean(new_alpha_prof);
inputs.RT.mean_distribution_var_closest_filename = out.mie_table_filename_closest_to_mean;
inputs.RT.use_custom_mie_calcs                = true;

alpha_param{nn} = new_alpha_prof;


%% Handle re >= 50 µm (limit of the custom Mie look-up table)

if any(re{nn} >= use_35_or_50)

    idx_remove = re{nn} >= use_35_or_50;
    n_removed  = sum(idx_remove);

    re{nn}(idx_remove)  = [];
    lwc{nn}(idx_remove) = [];
    z{nn}(idx_remove)   = [];
    tau{nn}(idx_remove) = [];

    % Update tau_c to the maximum remaining optical depth
    tau_c(nn) = tau{nn}(end);

    fprintf('  Warning: removed %d level(s) with re >= %d um from ORACLES profile %d\n', ...
        n_removed, use_35_or_50, measurement_idx);

end


%% Remove levels with re <= 0.1 µm (thin gap / bad data)

idx_0 = re{nn} <= 0.1;

if any(idx_0)
    re{nn}(idx_0)  = [];
    lwc{nn}(idx_0) = [];
    z{nn}(idx_0)   = [];
    tau{nn}(idx_0) = [];
end


%% Write the water cloud file

wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name, ...
    date_of_flight{nn}, time_of_flight(nn), ...
    inputs.compute_weighting_functions, which_computer, ...
    measurement_idx, wc_folder_path);


%% Free the ensemble profiles from memory

clear ds_oracles ensemble_profiles


%% Build the state-vector index table (one row per INP file)

[WW, VZ, VA, SA] = ndgrid(1:num_wl, 1:num_vza, 1:num_vaz, 1:num_saz);

idx = [VZ(:), VA(:), SA(:), WW(:)];

changing_variables_allStateVectors = [ ...
    vza(idx(:,1))', ...
    vaz(idx(:,2))', ...
    phi0(idx(:,3))', ...
    inputs.RT.wavelengths2run(idx(:,4), :) ];

% Append a column with the spectral-response index (increases chronologically)
changing_variables_allStateVectors = [changing_variables_allStateVectors, ...
    repmat((1:num_wl)', num_vza * num_saz * num_vaz, 1)];

num_cols = size(changing_variables_allStateVectors, 2);

if num_INP_files ~= size(changing_variables_allStateVectors, 1)
    error([newline, 'Number of rows in changing_variables not equal to number of INP files', newline])
end


%% Pre-compute per-iteration source flux and spectral response slices

wavelength_vec = [min(inputs.RT.wavelengths2run, [], 'all'), ...
    max(inputs.RT.wavelengths2run, [], 'all')];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);

wl_perturb = inputs.RT.source_file_resolution / 3;   % nm

spec_response_value = spec_response.value;

source_flux_perIter = cell(num_INP_files, 1);
for nn = 1:num_INP_files
    idx_wl = source_wavelength >= (changing_variables_allStateVectors(nn, num_cols-2) - wl_perturb) & ...
             source_wavelength <= (changing_variables_allStateVectors(nn, num_cols-1) + wl_perturb);
    source_flux_perIter{nn} = source_flux(idx_wl);
end

spec_response_perIter = cell(num_INP_files, 1);
for nn = 1:num_INP_files
    spec_response_perIter{nn} = spec_response_value(changing_variables_allStateVectors(nn, num_cols), :);
end


%% Extract minimal inputs fields for the parfor broadcast variable

compute_reflectivity_uvSpec = inputs.RT.compute_reflectivity_uvSpec;
rte_solver_str              = inputs.RT.rte_solver;

inputs_parfor.RT.rte_solver = rte_solver_str;

% Pre-allocate output arrays
inputFileName  = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);
Refl_model_allStateVectors = zeros(num_INP_files, 1);

% Re-read the measurement index into local vars to avoid struct broadcast
nn_meas    = 1;
wc_file_1  = wc_filename{nn_meas};
tau_c_1    = tau_c(nn_meas);
date_str_1 = date_of_flight{nn_meas};
time_str_1 = time_of_flight(nn_meas);


%% Run radiative transfer (single parfor loop: write → run → read → compute)

tic

parfor nn = 1:num_INP_files
% for nn = 1:num_INP_files

    vza_nn  = changing_variables_allStateVectors(nn, 1);
    vaz_nn  = changing_variables_allStateVectors(nn, 2);
    phi0_nn = changing_variables_allStateVectors(nn, 3);
    wavelengths = changing_variables_allStateVectors(nn, num_cols-2:num_cols-1);

    % ---- Unique INP / OUT filenames ----
    inputFileName{nn} = [num2str(mean(wavelengths)), 'nm_', campaign_name, '_', ...
        char(date_str_1), '_', ...
        num2str(time_str_1), '-UTC_VR-meas_', ...
        num2str(measurement_idx), '_vza', num2str(vza_nn), ...
        '_sza', num2str(sza), '_vaz', num2str(vaz_nn), ...
        '_saz', num2str(phi0_nn), '.INP'];

    outputFileName{nn} = ['OUTPUT_', inputFileName{nn}(1:end-4)];

    % ---- Write INP file ----
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName{nn}, inputs, ...
        wavelengths, wc_file_1, [], tau_c_1, [], [], vza_nn, sza, vaz_nn, phi0_nn);

    % ---- Run uvspec ----
    runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn}, which_computer);

    % ---- Read output ----
    [rad_calcs, ~, ~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, inputs_parfor, ...
        compute_reflectivity_uvSpec, vza_nn, vaz_nn);

    % ---- Compute reflectance ----
    [Refl_model_allStateVectors(nn), ~] = reflectanceFunction_ver3(rad_calcs, ...
        source_flux_perIter{nn}, spec_response_perIter{nn}, ...
        sza, vza_nn, vaz_nn);

    % ---- Delete INP and OUT files immediately to conserve scratch inodes ----
    inp_file_path = [libRadtran_inp, inputFileName{nn}];
    out_file_path = [libRadtran_inp, outputFileName{nn}, '.OUT'];

    if isfile(inp_file_path)
        delete(inp_file_path)
    end
    if isfile(out_file_path)
        delete(out_file_path)
    end

end

disp([newline, 'Calculations took ', num2str(toc), ' seconds'])


%% Reshape reflectance array

% Rows span wavelengths; columns span unique viewing geometries
Refl_model_allStateVectors = reshape(Refl_model_allStateVectors, num_wl, []);


%% Add Gaussian noise

inputs.measurement.uncert       = 0.003;   % fraction of measurement
inputs.measurement.standard_dev = inputs.measurement.uncert / 3;

Refl_model_with_noise_allStateVectors = ...
    (inputs.measurement.standard_dev .* Refl_model_allStateVectors) ...
    .* randn(size(Refl_model_allStateVectors)) ...
    + Refl_model_allStateVectors;

Refl_model_uncert_allStateVectors = ...
    inputs.measurement.uncert .* Refl_model_with_noise_allStateVectors;


%% Save output

inputs.folderpath_2save = output_dir;

if ~exist(inputs.folderpath_2save, 'dir')
    mkdir(inputs.folderpath_2save)
end

filename = [inputs.folderpath_2save, 'simulated_spectra_HySICS_reflectance_', ...
    num2str(numel(inputs.bands2run)), 'bands_', ...
    num2str(100 * inputs.measurement.uncert), '%_uncert_', ...
    campaign_name, '_inSitu_re_lwc_tauC_z_', ...
    char(date_of_flight{1}), '_', num2str(time_of_flight), 'UTC_prof-nn_', ...
    num2str(measurement_idx), '_vzaRange_', num2str(round(vza(1))), ...
    '-', num2str(round(vza(end))), '_vazRange_', num2str(round(vaz(1))), ...
    '-', num2str(round(vaz(end))), '_sza_', num2str(round(sza)), ...
    '_sazRange_', num2str(round(phi0(1))), '-', num2str(round(phi0(end))), ...
    '_sim-ran-on-', char(datetime('today')), '.mat'];

save(filename, ...
    'Refl_model_allStateVectors', ...
    'Refl_model_with_noise_allStateVectors', ...
    'Refl_model_uncert_allStateVectors', ...
    'inputs', 'spec_response', ...
    'changing_variables_allStateVectors', ...
    're', 'tau', 'lwc', 'z', 'alpha_param', 'era5');

fprintf('Successfully completed and saved ORACLES measurement %d\n', measurement_idx);

end
