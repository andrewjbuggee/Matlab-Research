function out = hysics_refl_from_synthetic_NN_inputs_ver2(input_file, idx, folder_paths, spec_response, settings)
%% HySICS reflectance for one synthetic cloud: a few libRadtran runs, not 636
%
% ver2 of hysics_refl_from_synthetic_NN_inputs.m. The physics and the
% libRadtran inputs are identical; only the way the spectrum is computed and
% returned changes.
%
% What changed and why
% --------------------
% 1. ONE libRadtran run per spectral SEGMENT instead of one per HySICS band.
%    The original runs uvspec 636 times, each over its band's ±20 nm response
%    window. Band centers are only ~3.1 nm apart, so neighbouring windows
%    overlap and every wavelength is computed ~13 times. Here the 636 bands
%    are split into settings.n_segments contiguous groups; uvspec runs once
%    over the union of each group's windows, and each band's reflectance is
%    taken from its own window of that output with the SAME function as before
%    (reflectanceFunction_ver3). REPTRAN computes at the same representative
%    wavelengths whatever range is requested, so the result is identical to
%    the per-band runs (verified bit-for-bit on 17 bands of cloud 553,
%    13 Nov 2026, and to machine precision on all 636 bands of three clouds,
%    14 Nov 2026). Segments, rather than one run over 332-2317 nm, lower the
%    peak memory of each uvspec run: a single full-range run peaked at 5.0 GB
%    on the laptop, the 332-856 nm segment at 3.1 GB, against CURC's Alpine
%    cluster, which has 3.8 GB per core. Memory does not fall in proportion
%    to the range, so budget ~4 GB per MATLAB worker (uvspec + worker).
% 2. No parfor and no parallel pool inside. Parallelism is across clouds, in
%    run_synthetic_spectra_chunk.m.
% 3. Nothing is deleted except this cloud's own files. The original deleted
%    every *.INP / *.OUT / *.DAT in the shared folders, which would destroy
%    other clouds' files when clouds run concurrently in the same folders.
% 4. The spectral response functions are passed in (computed once per chunk with
%    create_HySICS_specResponse_ver2) rather than recomputed for every cloud.
% 5. Nothing is saved to disk. The function returns the clean reflectance and
%    this cloud's truth fields (that values used to define the cloud and that are
%    used to compare with the retrieved output); run_synthetic_spectra_chunk.m 
%    writes one file per chunk. No noise is added: the HDF converter regenerates
%    noise from the clean reflectance.
%
% INPUTS
%   input_file     full path to the training-input NetCDF (07_build_training_inputs.py)
%   idx            1 cloud index (row of the 'cloud' dimension)
%   folder_paths   struct from define_folderPaths_for_HySICS: .libRadtran_inp,
%                  .libRadtran_data, .libRadtran_water_cloud_files,
%                  .atm_folder_path
%   spec_response  struct from create_HySICS_specResponse_ver2 for bands 1:636
%   settings       (optional) struct:
%                    .n_segments       number of uvspec runs per cloud (default 4)
%                    .source_file      solar source file (default hybrid 1 nm file)
%                    .source_flux, .source_wavelength
%                                      (optional) precomputed solar flux on the
%                                      source grid covering all bands
%                    .keep_files       keep INP/OUT/wc/radiosonde (default false)
%                    .verbose          print per-segment timing (default true)
%
% OUTPUT
%   out  struct: .refl (636x1, clean reflectance, 1/sr), .idx, geometry
%        (.sza .saz .vza .vaz, deg), truth (.re_um .lwc_g_per_m3 .z_km
%        (7x1, top -> base), .tau_c, .LWP_g_per_m2, .alpha, .T_K and
%        .vapor_cm3 (37x1), .pressure_hPa), .mie_table, .segment_ranges_nm,
%        .t_segments_s, .t_total_s
%
% By Andrew John Buggee (synthetic-input variant)

%%

t_start = tic;

% ---------------- settings -------------------------------------------------
if nargin < 5 || isempty(settings), settings = struct(); end
if ~isfield(settings, 'n_segments'),  settings.n_segments  = 4;     end
if ~isfield(settings, 'source_file')
    settings.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
end
if ~isfield(settings, 'keep_files'),  settings.keep_files  = false; end
if ~isfield(settings, 'verbose'),     settings.verbose     = true;  end

which_computer = whatComputer();
campaign_name  = 'syntheticORACLES';

%% Load this cloud's row from the NetCDF input file (identical to ver1)

if ~isfile(input_file)
    error('Training-input file not found: %s', input_file);
end

ncinfo_data = ncinfo(input_file);
dim_names   = {ncinfo_data.Dimensions.Name};
dim_lengths = [ncinfo_data.Dimensions.Length];
N_clouds    = dim_lengths(strcmp(dim_names, 'cloud'));
N_levels    = dim_lengths(strcmp(dim_names, 'level_cloud'));
N_atmos     = dim_lengths(strcmp(dim_names, 'level_atmos'));

if idx < 1 || idx > N_clouds
    error('idx must be in [1, %d]; got %d', N_clouds, idx);
end

re_um            = ncread(input_file, 're_um',            [1, idx], [N_levels, 1]);
lwc_g_per_m3     = ncread(input_file, 'lwc_g_per_m3',     [1, idx], [N_levels, 1]);
z_km             = ncread(input_file, 'z_km',             [1, idx], [N_levels, 1]);
alpha_scalar     = ncread(input_file, 'alpha',            idx,      1);
tau_c_scalar     = ncread(input_file, 'tau_c',            idx,      1);
LWP_g_per_m2     = ncread(input_file, 'LWP_g_per_m2',     idx,      1);
sza              = ncread(input_file, 'sza_deg',          idx,      1);
saz              = ncread(input_file, 'saz_deg',          idx,      1);
vza              = ncread(input_file, 'vza_deg',          idx,      1);
vaz              = ncread(input_file, 'vaz_deg',          idx,      1);
T_K              = ncread(input_file, 'T_K',              [1, idx], [N_atmos, 1]);
vapor_cm3        = ncread(input_file, 'vapor_molec_per_cm3', [1, idx], [N_atmos, 1]);
pressure_hPa     = ncread(input_file, 'pressure_hPa',     1,        N_atmos);

% Profiles arrive top -> base; column vectors are what the rest of the chain expects.
re_um            = double(re_um(:));
lwc_g_per_m3     = double(lwc_g_per_m3(:));
z_km             = double(z_km(:));
T_K              = double(T_K(:));
vapor_cm3        = double(vapor_cm3(:));
pressure_hPa     = double(pressure_hPa(:));

% keep the truth as read, before any level is dropped below
out.idx = idx;
out.sza = double(sza);  out.saz = double(saz);  out.vza = double(vza);  out.vaz = double(vaz);
out.re_um = re_um;  out.lwc_g_per_m3 = lwc_g_per_m3;  out.z_km = z_km;
out.tau_c = double(tau_c_scalar);  out.LWP_g_per_m2 = double(LWP_g_per_m2);
out.alpha = double(alpha_scalar);
out.T_K = T_K;  out.vapor_cm3 = vapor_cm3;  out.pressure_hPa = pressure_hPa;

%% libRadtran inputs (identical to ver1)

inputs.which_computer = which_computer;

libRadtran_inp       = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path       = folder_paths.libRadtran_water_cloud_files;

% Vertical-profile + RT solver block
inputs.RT.vert_homogeneous_str   = 'vert-non-homogeneous';
inputs.RT.num_re_parameters      = 2;
inputs.calc_type                 = 'simulated_spectra';
inputs.compute_weighting_functions = false;
inputs.RT.rte_solver             = 'disort';
inputs.RT.num_streams            = 16;
inputs.RT.use_nakajima_phaseCorrection = true;

% Source spectrum
inputs.RT.source_file            = settings.source_file;
inputs.RT.source_file_resolution = 0.1;     % nm

% HySICS spectral channels and their response windows
inputs.bands2run = (1:1:636)';
if size(spec_response.value, 1) ~= numel(inputs.bands2run)
    error('spec_response must hold all %d bands (got %d)', numel(inputs.bands2run), ...
        size(spec_response.value, 1));
end
inputs.RT.monochromatic_calc            = false;
inputs.RT.compute_reflectivity_uvSpec   = false;
inputs.RT.wavelengths2run = [spec_response.wavelength(:, 1), spec_response.wavelength(:, end)];   % nm

% Band model + base atmosphere file
inputs.RT.band_parameterization = 'reptran coarse';
inputs.RT.atm_file              = 'afglus.dat';

% Surface, day-of-year (Earth-Sun distance only)
inputs.RT.surface_albedo  = 0.04;
inputs.RT.day_of_year     = 316;

% Cox-Munk + Bodhaine cross sections
inputs.RT.use_coxMunk                  = true;
inputs.RT.wind_speed                   = 3;
inputs.RT.specify_cross_section_model  = true;
inputs.RT.crs_model_rayleigh           = 'Bodhaine29';

% Aerosols + gas concentration overrides
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
inputs.RT.errMsg             = 'quiet';

% Cloud + sensor geometry
inputs.RT.yesCloud              = true;
inputs.RT.modify_wc_opticalDepth = false;
inputs.RT.sensor_altitude        = 'toa';

%% Synthetic ERA5-style radiosonde file (named by cloud index)

inputs.RT.use_radiosonde_file = true;
inputs.RT.radiosonde_num_vars = 3;
inputs.RT.radiosonde_file = write_synthetic_radiosonde_DAT(pressure_hPa, T_K, ...
    vapor_cm3, folder_paths, [], inputs.RT.radiosonde_num_vars, idx, false);

%% Mie table closest to this cloud's alpha

new_alpha_prof = repmat(alpha_scalar, length(re_um), 1);
use_35_or_50 = 50;        % custom Mie tables span 1-50 microns
mie = find_custom_mieTable_closest_to_alpha_profile(new_alpha_prof, use_35_or_50, which_computer);

inputs.RT.distribution_var      = new_alpha_prof;
inputs.RT.mean_distribution_var = mean(new_alpha_prof);
inputs.RT.mean_distribution_var_closest_filename = mie.mie_table_filename_closest_to_mean;
inputs.RT.use_custom_mie_calcs  = true;
inputs.RT.wc_parameterization   = 'mie interpolate';
out.mie_table = mie.mie_table_filename_closest_to_mean;

%% Water-cloud .DAT file (mixed-layer geometry, named by cloud index)

date_of_flight = char(datetime('today'));
time_of_flight = idx;        % numeric tag for filename uniqueness
tau_c          = tau_c_scalar;

idx_zero = re_um <= 0.1;     % libRadtran cannot take r_e = 0
if any(idx_zero)
    re_um(idx_zero)        = [];
    lwc_g_per_m3(idx_zero) = [];
    z_km(idx_zero)         = [];
end
if any(re_um >= 50)
    error('Synthetic cloud %d has r_e >= 50 um; outside the Mie table range.', idx);
end

wc_filename = write_wc_file_mixed_layers(re_um, lwc_g_per_m3, z_km, campaign_name, ...
    date_of_flight, time_of_flight, inputs.compute_weighting_functions, which_computer, ...
    idx, wc_folder_path);

%% Solar source flux on the source grid (precomputed per chunk if supplied)

wl_all = inputs.RT.wavelengths2run;                 % (636 x 2) nm, each band's window
if isfield(settings, 'source_flux') && isfield(settings, 'source_wavelength')
    source_flux       = settings.source_flux;
    source_wavelength = settings.source_wavelength;
else
    [source_flux, source_wavelength] = read_solar_flux_file( ...
        [min(wl_all, [], 'all'), max(wl_all, [], 'all')], inputs.RT.source_file);
end
% half a source-grid step of tolerance when selecting a window, as in ver1
wl_perturb = inputs.RT.source_file_resolution / 3;

%% One uvspec run per segment of bands

n_bands = numel(inputs.bands2run);
n_seg   = max(1, min(settings.n_segments, n_bands));
seg_edges = round(linspace(0, n_bands, n_seg + 1));   % band groups (contiguous)

inputs_read.RT.rte_solver = inputs.RT.rte_solver;    % fields readUVSPEC_ver2 needs
refl = nan(n_bands, 1);
out.segment_ranges_nm = zeros(n_seg, 2);
out.t_segments_s      = zeros(n_seg, 1);

for s = 1:n_seg

    bands_s = (seg_edges(s) + 1):seg_edges(s + 1);
    % the run must cover every band window in this group
    seg_range = [min(wl_all(bands_s, 1)), max(wl_all(bands_s, 2))];
    out.segment_ranges_nm(s, :) = seg_range;

    inp_name = sprintf('cloud%06d_seg%02d_%.1f-%.1fnm.INP', idx, s, seg_range(1), seg_range(2));
    out_name = ['OUTPUT_', inp_name(1:end-4)];

    t_seg = tic;
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inp_name, ...
        inputs, seg_range, wc_filename, [], tau_c, [], [], vza, sza, vaz, saz);
    runUVSPEC_ver2(libRadtran_inp, inp_name, out_name, which_computer);
    [rad, ~, ~] = readUVSPEC_ver2(libRadtran_inp, out_name, inputs_read, ...
        inputs.RT.compute_reflectivity_uvSpec, vza, vaz);

    % reflectance of each band from its own window of this run's output
    for nn = bands_s
        sel = rad.wavelength >= (wl_all(nn, 1) - wl_perturb) & ...
              rad.wavelength <= (wl_all(nn, 2) + wl_perturb);
        ds = rad;
        ds.wavelength     = rad.wavelength(sel);
        ds.radiance.value = rad.radiance.value(sel);
        idx_src = source_wavelength >= (wl_all(nn, 1) - wl_perturb) & ...
                  source_wavelength <= (wl_all(nn, 2) + wl_perturb);
        refl(nn) = reflectanceFunction_ver3(ds, source_flux(idx_src), ...
            spec_response.value(nn, :), sza, vza, vaz);
    end
    out.t_segments_s(s) = toc(t_seg);

    if ~settings.keep_files
        delete_if_present([libRadtran_inp, inp_name]);
        delete_if_present([libRadtran_inp, out_name, '.OUT']);
    end
    if settings.verbose
        fprintf('  cloud %d segment %d/%d (%.1f-%.1f nm, %d bands): %.1f s\n', idx, s, n_seg, ...
            seg_range(1), seg_range(2), numel(bands_s), out.t_segments_s(s));
    end
end

% this cloud's own input files
if ~settings.keep_files
    delete_if_present([wc_folder_path, wc_filename]);
    delete_if_present(inputs.RT.radiosonde_file);
end

out.refl = refl;
out.t_total_s = toc(t_start);

end


function delete_if_present(f)
% delete one file if it exists (never a wildcard)
if isfile(f)
    delete(f);
end
end
