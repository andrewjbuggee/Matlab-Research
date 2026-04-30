function hysics_refl_from_synthetic_NN_inputs(input_file, idx, folder_paths, output_dir)
%% Generate HySICS measurements from a row of the synthetic-cloud training-input .nc
%
% Synthetic-cloud counterpart to hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver3.
% The VOCALS-REx ensemble loading, tau_c index adjustments, and find-closest-ERA5
% step are dropped — this function reads everything it needs straight from one
% row of the .nc file produced by 07_build_training_inputs.py.
%
% Each row is one synthetic cloud paired with one (SZA, SAZ, VZA, VAZ).
% No outer geometry sweep:  num_INP_files = num_wl  (just spectral channels).
%
% INPUT:
%   input_file    Full path to the training-input NetCDF.
%   idx           1-based row (cloud index) to process.
%   folder_paths  Same struct used by ver3 (.libRadtran_inp,
%                 .libRadtran_water_cloud_files, .atm_folder_path,
%                 .libRadtran_mie_folder, .libRadtran_data,
%                 .which_computer, etc.).
%   output_dir    Directory to write the per-cloud .mat output (with
%                 trailing '/').
%
% By Andrew John Buggee (synthetic-input variant, drafted with Claude)
%%

%% Tidy up old INP/OUT/wc/atmmod artefacts from a previous run

delete([folder_paths.libRadtran_inp,                '*.INP'])
delete([folder_paths.libRadtran_inp,                '*.OUT'])
delete([folder_paths.libRadtran_water_cloud_files,  '*.DAT'])
delete([folder_paths.atm_folder_path,               '*-aboveCloud.DAT'])
delete([folder_paths.libRadtran_mie_folder,         '*.INP'])
delete([folder_paths.libRadtran_mie_folder,         '*.OUT'])

%% Parallel pool

start_parallel_pool(folder_paths.which_computer)

which_computer = whatComputer();

campaign_name = 'syntheticORACLES';

%% Load the per-cloud row from the NetCDF input file

if ~isfile(input_file)
    error('Training-input file not found: %s', input_file);
end

n_clouds = ncreadatt(input_file, '/', 'title');                 %#ok<NASGU>
% number of clouds → length of the 'cloud' dimension
ncinfo_data = ncinfo(input_file);
dim_names   = {ncinfo_data.Dimensions.Name};
dim_lengths = [ncinfo_data.Dimensions.Length];
N_clouds    = dim_lengths(strcmp(dim_names, 'cloud'));
N_levels    = dim_lengths(strcmp(dim_names, 'level_cloud'));
N_atmos     = dim_lengths(strcmp(dim_names, 'level_atmos'));

if idx < 1 || idx > N_clouds
    error('idx must be in [1, %d]; got %d', N_clouds, idx);
end

fprintf('\nProcessing synthetic cloud %d of %d\n', idx, N_clouds);

% read this row only
re_um            = ncread(input_file, 're_um',            [1, idx], [N_levels, 1]);
lwc_g_per_m3     = ncread(input_file, 'lwc_g_per_m3',     [1, idx], [N_levels, 1]);
z_km             = ncread(input_file, 'z_km',             [1, idx], [N_levels, 1]);
alpha_scalar     = ncread(input_file, 'alpha',            idx,      1);
tau_c_scalar     = ncread(input_file, 'tau_c',            idx,      1);
LWP_g_per_m2     = ncread(input_file, 'LWP_g_per_m2',     idx,      1);             %#ok<NASGU>
sza              = ncread(input_file, 'sza_deg',          idx,      1);
saz              = ncread(input_file, 'saz_deg',          idx,      1);
vza              = ncread(input_file, 'vza_deg',          idx,      1);
vaz              = ncread(input_file, 'vaz_deg',          idx,      1);
T_K              = ncread(input_file, 'T_K',              [1, idx], [N_atmos, 1]);
vapor_cm3        = ncread(input_file, 'vapor_molec_per_cm3', [1, idx], [N_atmos, 1]);
pressure_hPa     = ncread(input_file, 'pressure_hPa',     1,        N_atmos);

% Profiles arrive top→base; column-vector form is what the rest of the chain expects.
re_um            = double(re_um(:));
lwc_g_per_m3     = double(lwc_g_per_m3(:));
z_km             = double(z_km(:));
T_K              = double(T_K(:));
vapor_cm3        = double(vapor_cm3(:));
pressure_hPa     = double(pressure_hPa(:));

%% Inputs needed to build INP files

num_meas = 1;
inputs.which_computer = which_computer;

libRadtran_inp       = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path       = folder_paths.libRadtran_water_cloud_files;

print_libRadtran_err = false;

% Vertical-profile + RT solver block ---------------------------------
inputs.RT.vert_homogeneous_str   = 'vert-non-homogeneous';
inputs.RT.num_re_parameters      = 2;
inputs.calc_type                 = 'simulated_spectra';
inputs.compute_weighting_functions = false;
inputs.RT.rte_solver             = 'disort';
inputs.RT.num_streams            = 16;
inputs.RT.use_nakajima_phaseCorrection = true;

% Source spectrum ----------------------------------------------------
inputs.RT.source_file            = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
inputs.RT.source_file_resolution = 0.1;     % nm

% HySICS spectral channels ------------------------------------------
inputs.bands2run = (1:1:636)';

inputs.RT.monochromatic_calc            = false;
inputs.RT.compute_reflectivity_uvSpec   = false;

% Spectral response functions ---------------------------------------
spec_response = create_HySICS_specResponse(inputs.bands2run, inputs.RT.source_file, ...
    inputs.which_computer);

inputs.RT.wavelengths2run = zeros(length(inputs.bands2run), 2);
for ww = 1:length(inputs.bands2run)
    inputs.RT.wavelengths2run(ww, :) = [spec_response.wavelength(ww, 1), ...
                                        spec_response.wavelength(ww, end)];
end

% Band model + base atmosphere file ---------------------------------
inputs.RT.band_parameterization = 'reptran coarse';
inputs.RT.atm_file              = 'afglus.dat';

% Surface, day-of-year (Earth-Sun distance only — sampled SZA already
% encodes the actual solar geometry) ---------------------------------
inputs.RT.surface_albedo  = 0.04;
inputs.RT.day_of_year     = 316;     % nominal; reflectance is insensitive at this precision

% Cox-Munk + Bodhaine cross sections --------------------------------
inputs.RT.use_coxMunk                  = true;
inputs.RT.wind_speed                   = 3;
inputs.RT.specify_cross_section_model  = true;
inputs.RT.crs_model_rayleigh           = 'Bodhaine29';

% Aerosols + gas concentration overrides ----------------------------
inputs.RT.yesAerosols          = true;
inputs.RT.aerosol_type         = 4;
inputs.RT.aerosol_opticalDepth = 0.1;

inputs.RT.modify_total_columnWaterVapor      = false;
inputs.RT.modify_aboveCloud_columnWaterVapor = false;
inputs.RT.modify_CO2          = true;     inputs.RT.CO2_mixing_ratio = 418;
inputs.RT.modify_N2           = false;    inputs.RT.N2_mixing_ratio  = 0;
inputs.RT.modify_NO2          = false;    inputs.RT.NO2_mixing_ratio = 0;
inputs.RT.modify_O2           = false;    inputs.RT.O2_mixing_ratio  = 0;
inputs.RT.modify_O3           = false;    inputs.RT.O3_mixing_ratio  = 0;

inputs.RT.no_molecular_abs   = false;
inputs.RT.no_scattering_mol  = false;
inputs.RT.no_scattering_aer  = false;

if print_libRadtran_err
    inputs.RT.errMsg = 'verbose';
else
    inputs.RT.errMsg = 'quiet';
end

% Cloud + sensor geometry ------------------------------------------
inputs.RT.yesCloud              = true;
inputs.RT.modify_wc_opticalDepth = true;
inputs.RT.sensor_altitude        = 'toa';

%% Write the synthetic ERA5-style radiosonde file

inputs.RT.use_radiosonde_file = true;
inputs.RT.radiosonde_num_vars = 3;

inputs.RT.radiosonde_file = write_synthetic_radiosonde_DAT(pressure_hPa, T_K, ...
    vapor_cm3, folder_paths, [], inputs.RT.radiosonde_num_vars, idx, false);

%% Mie table closest to this cloud's mean alpha

% the .nc stores a single mean alpha per cloud, so broadcast across levels
new_alpha_prof = repmat(alpha_scalar, length(re_um), 1);

use_35_or_50 = 50;        % custom Mie tables span 1–50 microns
out = find_custom_mieTable_closest_to_alpha_profile(new_alpha_prof, use_35_or_50, which_computer);

inputs.RT.distribution_var      = new_alpha_prof;
inputs.RT.mean_distribution_var = mean(new_alpha_prof);
inputs.RT.mean_distribution_var_closest_filename = out.mie_table_filename_closest_to_mean;
inputs.RT.use_custom_mie_calcs  = true;
inputs.RT.wc_parameterization   = 'mie interpolate';

%% Write the water-cloud .DAT file (per-cloud profile)

% write_wc_file_from_in_situ expects column vectors; profiles already top→base.
date_of_flight = {char(datetime('today'))};
time_of_flight = idx;        % numeric tag for filename uniqueness
tau_c          = tau_c_scalar;

% Drop near-zero r_e levels (libRadtran chokes on 0); should never trigger
% for synthetic data but kept for safety.
idx_zero = re_um <= 0.1;
if any(idx_zero)
    re_um(idx_zero)        = [];
    lwc_g_per_m3(idx_zero) = [];
    z_km(idx_zero)         = [];
end

if any(re_um >= 50)
    error('Synthetic cloud %d has r_e >= 50 µm; outside the Mie table range.', idx);
end

wc_filename = write_wc_file_from_in_situ(re_um, lwc_g_per_m3, z_km, campaign_name, ...
    date_of_flight{1}, time_of_flight, ...
    inputs.compute_weighting_functions, which_computer, ...
    idx, wc_folder_path);

%% Build per-INP geometry table — one row per spectral channel

num_wl        = length(inputs.bands2run);
num_INP_files = num_wl;     % no outer geometry sweep

% [sza, vza, vaz, saz, wl_lo, wl_hi, band_idx]
changing_variables_allStateVectors = [...
    repmat([sza, vza, vaz, saz], num_INP_files, 1), ...
    inputs.RT.wavelengths2run, ...
    (1:num_wl)' ];

num_cols = size(changing_variables_allStateVectors, 2);

%% Pre-compute per-iteration source flux + spectral response slices

wavelength_vec = [min(inputs.RT.wavelengths2run, [], 'all'), ...
                  max(inputs.RT.wavelengths2run, [], 'all')];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);
wl_perturb = inputs.RT.source_file_resolution / 3;

spec_response_value = spec_response.value;

source_flux_perIter   = cell(num_INP_files, 1);
spec_response_perIter = cell(num_INP_files, 1);
for nn = 1:num_INP_files
    idx_wl = source_wavelength >= (changing_variables_allStateVectors(nn, num_cols-2) - wl_perturb) & ...
             source_wavelength <= (changing_variables_allStateVectors(nn, num_cols-1) + wl_perturb);
    source_flux_perIter{nn}   = source_flux(idx_wl);
    spec_response_perIter{nn} = spec_response_value(changing_variables_allStateVectors(nn, num_cols), :);
end

% minimal fields needed inside parfor
compute_reflectivity_uvSpec = inputs.RT.compute_reflectivity_uvSpec;
inputs_parfor.RT.rte_solver = inputs.RT.rte_solver;

inputFileName  = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);
Refl_model_allStateVectors = zeros(num_INP_files, 1);

%% Pipelined parfor — write INP, run uvspec, read OUT, compute reflectance

tic
parfor nn = 1:num_INP_files

    sza_nn  = changing_variables_allStateVectors(nn, 1);
    vza_nn  = changing_variables_allStateVectors(nn, 2);
    vaz_nn  = changing_variables_allStateVectors(nn, 3);
    saz_nn  = changing_variables_allStateVectors(nn, 4);
    wavelengths = changing_variables_allStateVectors(nn, num_cols-2:num_cols-1);

    inputFileName{nn} = [num2str(mean(wavelengths)), 'nm_', campaign_name, '_cloud', ...
        num2str(idx), '_sza', num2str(sza_nn), '_saz', num2str(saz_nn), ...
        '_vza', num2str(vza_nn), '_vaz', num2str(vaz_nn), '.INP'];
    outputFileName{nn} = ['OUTPUT_', inputFileName{nn}(1:end-4)];

    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName{nn}, ...
        inputs, wavelengths, wc_filename, [], tau_c, [], [], ...
        vza_nn, sza_nn, vaz_nn, saz_nn);

    runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn}, which_computer);

    [rad_calcs, ~, ~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, ...
        inputs_parfor, compute_reflectivity_uvSpec, vza_nn, vaz_nn);

    [Refl_model_allStateVectors(nn), ~] = reflectanceFunction_ver3(rad_calcs, ...
        source_flux_perIter{nn}, spec_response_perIter{nn}, ...
        sza_nn, vza_nn, vaz_nn);

    % delete INP/OUT immediately to conserve scratch inodes
    inp_path = [libRadtran_inp, inputFileName{nn}];
    out_path = [libRadtran_inp, outputFileName{nn}, '.OUT'];
    if isfile(inp_path), delete(inp_path); end
    if isfile(out_path), delete(out_path); end
end
toc

%% Add Gaussian noise (HySICS + EMIT uncertainty)

inputs.measurement.uncert_hysics       = 0.003;
inputs.measurement.uncert_emit         = 0.04;
inputs.measurement.standard_dev_hysics = inputs.measurement.uncert_hysics;
inputs.measurement.standard_dev_emit   = inputs.measurement.uncert_emit;

Refl_model_with_noise_allStateVectors_hysics = ...
    (inputs.measurement.standard_dev_hysics .* Refl_model_allStateVectors) ...
    .* randn(size(Refl_model_allStateVectors)) + Refl_model_allStateVectors;
Refl_model_uncert_allStateVectors_hysics = ...
    inputs.measurement.uncert_hysics .* Refl_model_with_noise_allStateVectors_hysics;

Refl_model_with_noise_allStateVectors_emit = ...
    (inputs.measurement.standard_dev_emit .* Refl_model_allStateVectors) ...
    .* randn(size(Refl_model_allStateVectors)) + Refl_model_allStateVectors;
Refl_model_uncert_allStateVectors_emit = ...
    inputs.measurement.uncert_emit .* Refl_model_with_noise_allStateVectors_emit;

%% Save output

if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

% pack the per-cloud truth arrays in the same shape ver3 used: one cell per
% measurement — here just one — so downstream tooling that expects {nn} works.
re  = {re_um};
lwc = {lwc_g_per_m3};
z   = {z_km};
tau = {tau_c_scalar};

% pack the synthetic atmosphere into an "era5"-shaped struct so the saved
% file looks structurally similar to ver3 outputs.
era5.datProfiles.GP_height           = nan(N_atmos, 1);   % not needed downstream
era5.datProfiles.T                   = T_K;
era5.datProfiles.p                   = pressure_hPa;
era5.datProfiles.q                   = nan(N_atmos, 1);
era5.datProfiles.vapor_concentration = vapor_cm3;
era5.datProfiles.vapor_massDensity   = nan(N_atmos, 1);

inputs.folderpath_2save = output_dir;

filename = [inputs.folderpath_2save, ...
    'simulated_spectra_HySICS_reflectance_', num2str(numel(inputs.bands2run)), 'bands_', ...
    num2str(100 * inputs.measurement.uncert_hysics), '%_uncert_', campaign_name, ...
    '_cloud', num2str(idx), '_sza', num2str(round(sza)), ...
    '_saz', num2str(round(saz)), '_vza', num2str(round(vza)), ...
    '_vaz', num2str(round(vaz)), ...
    '_sim-ran-on-', char(datetime('today')), '.mat'];

save(filename, 'Refl_model_allStateVectors', ...
    'Refl_model_with_noise_allStateVectors_hysics', ...
    'Refl_model_uncert_allStateVectors_hysics', ...
    'Refl_model_with_noise_allStateVectors_emit', ...
    'Refl_model_uncert_allStateVectors_emit', ...
    'inputs', 'spec_response', 'changing_variables_allStateVectors', ...
    're', 'tau', 'lwc', 'z', 'era5');

fprintf('Successfully completed and saved synthetic cloud %d → %s\n', idx, filename);

end
