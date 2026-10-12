%% Simulate HySICS top-of-atmosphere radiance and reflectance spectra for Peter
% Single-layer, vertically homogeneous liquid water clouds with one
% effective radius and several cloud optical thicknesses.
%
% Workflow:
%   (1) Define the scene settings (geometry, solar source, atmosphere,
%       surface) in inputs_measurement.
%   (2) Build the full libRadtran input structure with
%       create_HySICS_inputs_forPeter_8_Oct_2026.m. That function sets the
%       cloud (r_e, tau_c, cloud top/base), the RT solver, gases, aerosols,
%       and the HySICS bands to run (all 636 channels).
%   (3) Build the HySICS spectral response functions (SRFs) and define the
%       wavelength range of each uvspec run from them.
%   (4) Write one water cloud (wc) file per (r_e, tau_c) pair.
%   (5) Write one INP file per (r_e, tau_c, HySICS channel).
%   (6) Run uvspec, then convolve the monochromatic TOA radiance with each
%       channel's SRF to get the band radiance and band reflectance.
%   (7) Save and plot TOA radiance vs wavelength and TOA reflectance vs
%       wavelength, one line per optical thickness.
%
% Cost: number of uvspec runs = 636 * num_tauC * num_rEff. With the default
% inputs (r_e = 10 um, tau_c = 10:10:60) that is 3816 runs.
%
% By Andrew John Buggee


clear variables



%% Define folder paths and start a parallel pool

% folder_paths holds the libRadtran INP/OUT folder, data folder, water
% cloud folder and Mie folder for whichever computer this runs on.
% An input of 0 uses the base folder names (no numbered extension).
folder_paths = define_folderPaths_for_HySICS(0);

% Defining custom directory for Peter's calculations
folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    'For_Peter/HySICS_spectra/'];

% Define where figures should be saved
folder_paths.saved_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/For_Peter/saved_figures/';

which_computer = folder_paths.which_computer;

% Unpack folder paths into plain variables so the parfor loops below don't
% broadcast the whole folder_paths structure
libRadtran_inp = folder_paths.libRadtran_inp;                       % INP and OUT files are written here
libRadtran_data_path = folder_paths.libRadtran_data;                % libRadtran data directory (absolute path)
wc_folder_path = folder_paths.libRadtran_water_cloud_files;         % water cloud .DAT files are written here
mie_folder_path = folder_paths.libRadtran_mie_folder;               % write_wc_file runs Mie calcs here


% Make sure the INP/OUT and wc folders exist. If they don't, fopen returns
% -1 on the workers and write_INP_file fails with "Invalid fileID".
if ~exist(libRadtran_inp, 'dir')
    mkdir(libRadtran_inp)
end

if ~exist(wc_folder_path, 'dir')
    mkdir(wc_folder_path)
end


% Delete old INP and OUT files so stale outputs can't be read by mistake
delete([libRadtran_inp, '*.INP'])
delete([libRadtran_inp, '*.OUT'])


% Start the parallel pool (reuses an existing pool if one is running)
start_parallel_pool(which_computer)



%% Define the scene settings passed to the inputs function

% create_HySICS_inputs_forPeter_8_Oct_2026 copies these fields from
% inputs_measurement.RT, so they must all be defined here.

% ----- Band model -----
% REPTRAN coarse molecular absorption parameterization (Gasteiger et al.,
% 2014, JQSRT, doi:10.1016/j.jqsrt.2014.06.024)
inputs_measurement.RT.band_parameterization = 'reptran coarse';


% ----- Solar source -----
% TSIS-1 Hybrid Solar Reference Spectrum (Coddington et al., 2021, GRL,
% doi:10.1029/2020GL091709). Despite the file name, these data are sampled
% every 0.1 nm and are in mW/(m^2 nm). The libRadtran output radiance takes
% the units of the source file, so radiance is in mW/(m^2 nm sr).
inputs_measurement.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
inputs_measurement.RT.source_file_resolution = 0.1;         % nm


% ----- Atmosphere -----
% AFGL U.S. Standard Atmosphere (Anderson et al., 1986, AFGL-TR-86-0110)
inputs_measurement.RT.atm_file = 'afglus.dat';


% ----- Surface -----
% Ocean surface albedo. The inputs function also turns on the Cox-Munk
% ocean BRDF (brdf_cam u10 = 3 m/s), which sets the surface reflection.
inputs_measurement.RT.surface_albedo = 0.04;


% ----- Sensor altitude -----
% 'toa' tells uvspec to output radiance at the top of the atmosphere
inputs_measurement.RT.sensor_altitude = 'toa';


% ----- Sun and viewing geometry -----
% Same geometry as the EMIT pixel used in Paper 1 (Fig. 3a) and the
% Paper 2 HySICS simulations
inputs_measurement.RT.sza = 28;          % degrees - solar zenith angle
inputs_measurement.RT.phi0 = 0;         % degrees - solar azimuth angle (libRadtran convention)
inputs_measurement.RT.vza = 0;        % degrees - viewing zenith angle
inputs_measurement.RT.vaz = 0;      % degrees - viewing azimuth angle (libRadtran convention)


% ----- Cloud top and base -----
% Not defined here, so the inputs function uses its default of
% z_topBottom = [1.5, 1] km (cloud top, cloud base). Uncomment to override.
% inputs_measurement.RT.z_topBottom = [1.5, 1];     % km above surface



%% Build the full libRadtran input structure

% false = don't write the libRadtran error message file
print_libRadtran_err = true;

inputs = create_HySICS_inputs_forPeter_8_Oct_2026(inputs_measurement, print_libRadtran_err);

% store which computer this was run on
inputs.which_computer = which_computer;

% label the calculation type
inputs.calc_type = 'simulated_spectra';

% One droplet size parameter: a single r_e for the whole (homogeneous) cloud.
% write_wc_file needs this value.
inputs.RT.num_re_parameters = 1;

% Don't modify the total or above-cloud column water vapor. Use the
% water vapor profile in the atmosphere file as is.
inputs.RT.modify_total_columnWaterVapor = false;
inputs.RT.modify_aboveCloud_columnWaterVapor = false;



%% Create the HySICS spectral response functions

% spec_response.value(ww,:) is the SRF of channel ww and
% spec_response.wavelength(ww,:) is its wavelength grid (nm). Each grid
% spans +/- 20 nm around the channel center and sits on the solar source
% wavelength grid, because libRadtran computes radiance on the source grid.
% The SRFs are Gaussians whose FWHM is fit to the measured HySICS
% instrument line shapes (407-1813 nm) and extrapolated to the full range.
spec_response = create_HySICS_specResponse(inputs.bands2run, inputs.RT.source_file, ...
    which_computer);


% The wavelength range of each uvspec run is the lower and upper bound of
% that channel's SRF grid
num_wl = length(inputs.bands2run);
inputs.RT.wavelengths2run = zeros(num_wl, 2);     % nm - [lower, upper] bound of each channel

for ww = 1:num_wl

    inputs.RT.wavelengths2run(ww,:) = [spec_response.wavelength(ww, 1),...
        spec_response.wavelength(ww, end)];

end

% Use the SRF-weighted mean wavelength as the channel wavelength for
% plotting. This is the centroid of each Gaussian SRF.
wl_center_nm = sum(spec_response.value .* spec_response.wavelength, 2) ./ ...
    sum(spec_response.value, 2);     % nm



%% Define every combination of r_e, tau_c and wavelength

% length of each independent variable
num_rEff = length(inputs.RT.re);
num_tauC = length(inputs.RT.tau_c);

num_cloud_files = num_rEff*num_tauC;         % one wc file per (r_e, tau_c) pair
num_INP_files = num_cloud_files*num_wl;      % one INP file per (r_e, tau_c, channel)


% changing_variables steps through r_e, tau_c, and wavelength.
% In for-loop speak it would be:
% for rr = 1:num_rEff
%   for tt = 1:num_tauC
%       for ww = 1:num_wl
% Columns: [r_e (um), tau_c, wl_lower (nm), wl_upper (nm), channel index]
changing_variables = [reshape(repmat(inputs.RT.re, num_tauC*num_wl, 1), [], 1),...
    repmat(reshape(repmat(inputs.RT.tau_c, num_wl, 1), [], 1), num_rEff, 1),...
    repmat(inputs.RT.wavelengths2run, num_rEff*num_tauC, 1),...
    repmat((1:num_wl)', num_rEff*num_tauC, 1)];



%% Write one water cloud file for each (r_e, tau_c) pair

% Row of changing_variables where each unique (r_e, tau_c) pair starts
idx_unique_cloud = 1:num_wl:num_INP_files;

wc_filename_unique = cell(num_cloud_files, 1);

tic

parfor nn = 1:num_cloud_files

    % write_wc_file converts (r_e, tau_c) to liquid water content using Mie
    % calculations at lambda_forTau (500 nm), for a gamma droplet size
    % distribution. 'vert-homogeneous' gives a single-layer cloud between
    % z_topBottom(2) and z_topBottom(1). nn is the unique file index so
    % parallel workers don't write over each other's files.
    temp = write_wc_file(changing_variables(idx_unique_cloud(nn), 1),...
        changing_variables(idx_unique_cloud(nn), 2), inputs.RT.z_topBottom,...
        inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
        inputs.RT.distribution_var, inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
        inputs.RT.indVar, inputs.compute_weighting_functions, which_computer,...
        nn, inputs.RT.num_re_parameters, wc_folder_path, mie_folder_path);

    wc_filename_unique{nn} = temp{1};

end

% Every channel of a given (r_e, tau_c) pair uses the same wc file
wc_filename = reshape(repmat(wc_filename_unique', num_wl, 1), [], 1);

disp([newline, 'Finished writing water cloud files.', newline])
toc



%% Write all INP files

inputFileName = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);

tic

parfor nn = 1:num_INP_files

    % wavelength range of this uvspec run (nm)
    wavelengths = changing_variables(nn, 3:4);

    % ------------------------------------------------
    % ---- Define the input and output filenames! ----
    % ------------------------------------------------
    % each name contains the channel wavelength, r_e and tau_c so it is
    % unique and traceable
    inputFileName{nn} = ['forPeter_', num2str(mean(wavelengths)), '_nm_rEff_',...
        num2str(changing_variables(nn,1)), '_tauC_', num2str(changing_variables(nn,2)), '_',...
        inputs.RT.atm_file(1:end-4), '.INP'];

    outputFileName{nn} = ['OUTPUT_', inputFileName{nn}(1:end-4)];


    % ------------------ Write the INP File --------------------
    % The last argument is the optical depth passed to 'wc_modify tau set'.
    % write_INP_file only writes that line when custom Mie tables are used.
    % With the default libRadtran Mie table (use_custom_mie_calcs = false)
    % tau_c is set by the LWC written in the wc file.
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName{nn},...
        inputs, wavelengths, wc_filename{nn}, [], changing_variables(nn,2));

end

disp([newline, 'Finished writing INP files.', newline])
toc



%% Prepare the solar flux and SRFs for each channel

% Read the solar flux over the full HySICS range
wavelength_vec = [min(inputs.RT.wavelengths2run, [], "all"), max(inputs.RT.wavelengths2run, [], "all")];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);   % W/(m^2 nm)

% Add and subtract a small fraction of the source resolution so rounding
% errors don't drop the end points when selecting wavelengths
wl_perturb = inputs.RT.source_file_resolution/3;   % nm


% Store the source flux and SRF of each channel in cell arrays indexed by
% channel number (column 5 of changing_variables). Inside the parfor loop
% these are indexed by the channel number, not the loop variable, so they
% are broadcast to every worker. They are small (636 x ~401 values each,
% a few MB), so the broadcast cost is negligible compared to uvspec.
source_flux_per_channel = cell(num_wl, 1);      % W/(m^2 nm)
srf_per_channel = cell(num_wl, 1);              % SRF weights on the source grid

for ww = 1:num_wl

    idx_wl = source_wavelength >= (inputs.RT.wavelengths2run(ww,1) - wl_perturb) &...
        source_wavelength <= (inputs.RT.wavelengths2run(ww,2) + wl_perturb);

    source_flux_per_channel{ww} = source_flux(idx_wl);
    srf_per_channel{ww} = reshape(spec_response.value(ww,:), [], 1);

end



%% Run uvspec, then compute the band radiance and band reflectance

% One parfor iteration per INP file (one row of changing_variables). Every
% iteration is a single uvspec run, so all workers stay busy no matter how
% the runs split between channels and optical thicknesses.
Rad_TOA_vec = zeros(num_INP_files, 1);         % mW/(m^2 nm sr) - SRF-weighted TOA radiance
Refl_TOA_vec = zeros(num_INP_files, 1);        % unitless - TOA reflectance (pi*I / (mu0*F0))

compute_reflectivity_uvSpec = inputs.RT.compute_reflectivity_uvSpec;

tic

parfor nn = 1:num_INP_files

    % channel index of this INP file -> its SRF and solar flux
    channel_idx = changing_variables(nn, 5);
    srf = srf_per_channel{channel_idx};            % SRF weights (column vector)
    flux = source_flux_per_channel{channel_idx};   % W/(m^2 nm)


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------
    runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn}, which_computer);


    % read .OUT file
    % radiance is in units of mW/(m^2 nm sr)
    [ds, ~, ~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, inputs,...
        compute_reflectivity_uvSpec);


    % ----- Band radiance -----
    % SRF-weighted mean of the spectral radiance, with the SRF normalized to unit area:
    %   L_band = int( L(lambda) * SRF(lambda) dlambda ) / int( SRF(lambda) dlambda )
    % This exact form is written for MODIS thermal-band radiance (MODIS L1B ATBD,
    % Xiong et al., v4 2013, updated 2020, Eq. 3.5). For MODIS reflected solar bands,
    % L = rho*cos(theta)*E_sun / (pi*d^2) (ATBD Eq. 4.5), where E_sun is the
    % RSR-weighted, band-averaged solar irradiance. If E_sun = int(SRF*F)/int(SRF),
    % that L equals L_band above (derived here, not stated in the ATBD).
    L_lambda = reshape(ds.radiance.value, [], 1);         % mW/(m^2 nm sr)
    Rad_TOA_vec(nn) = trapz(ds.wavelength, L_lambda .* srf) ./ trapz(ds.wavelength, srf);


    % ----- Band reflectance -----
    % reflectanceFunction_ver2 computes R(lambda) = pi*L / (mu0*F0) and
    % averages it weighted by SRF * F0 (King et al., 1992, IEEE TGRS,
    % doi:10.1109/36.124212; Platnick et al., 2017, IEEE TGRS,
    % doi:10.1109/TGRS.2016.2610522). This is the same reflectance
    % used in the Paper 1 and Paper 2 HySICS simulations.
    [Refl_TOA_vec(nn), ~] = reflectanceFunction_ver2(inputs, ds, flux, srf);

end

disp([newline, 'Finished running uvspec.', newline])
toc


% Rearrange into [channel, tau_c, r_e]. The rows of changing_variables step
% through channel fastest, then tau_c, then r_e, so a column-major reshape
% puts each index in the right place.
Rad_TOA = reshape(Rad_TOA_vec, num_wl, num_tauC, num_rEff);       % mW/(m^2 nm sr)
Refl_TOA = reshape(Refl_TOA_vec, num_wl, num_tauC, num_rEff);     % unitless



%% Save the output

% Save next to this script in the For_Peter folder


if ~exist(folder_paths.HySICS_simulated_spectra, 'dir')
    mkdir(folder_paths.HySICS_simulated_spectra)
end

filename = [folder_paths.HySICS_simulated_spectra, 'simulated_HySICS_TOA_rad_refl_homogeneous_cloud_rEff_',...
    num2str(inputs.RT.re(1)), 'um_', num2str(num_tauC), '-tauC_sza_', num2str(round(inputs.RT.sza)),...
    '_vza_', num2str(round(inputs.RT.vza)), '_sim-ran-on-', char(datetime("today")), '.mat'];

% Add a revision number so earlier runs aren't overwritten
rev = 1;
filename_rev = [filename(1:end-4), '_rev', num2str(rev), '.mat'];
while isfile(filename_rev)
    rev = rev + 1;
    filename_rev = [filename(1:end-4), '_rev', num2str(rev), '.mat'];
end

save(filename_rev, "Rad_TOA", "Refl_TOA", "wl_center_nm", "inputs", "spec_response",...
    "changing_variables");

disp([newline, 'Saved output to: ', filename_rev, newline])



%% Plot TOA radiance and TOA reflectance vs wavelength

% Plot the first effective radius. One line per optical thickness.
rr_2plot = 1;

% Optical thickness increases monotonically, so use a sequential colormap.
% Darker lines are thinner clouds and lighter lines are thicker clouds.
% line_colors = parula(num_tauC + 1);
% line_colors = line_colors(1:num_tauC, :);       % drop the final yellow, which is hard to see on white
line_colors = mySavedColors(62:(62+length(inputs.RT.tau_c)), 'fixed');

% Font sizes
ttl_fnt = 25;
sub_ttl_fnt = 20;
ax_fnt = 25;
lgnd_fnt = 20;

legend_str = cell(length(inputs.RT.tau_c), 1);
for tt = 1:length(inputs.RT.tau_c)
    legend_str{tt} = ['$\tau_c = $', num2str(inputs.RT.tau_c(tt))];
end

% Text describing the scene, used as the subtitle of both figures
scene_str = ['$r_e = $', num2str(inputs.RT.re(rr_2plot)), ' $\mu m$, ',...
    'cloud top $= $', num2str(inputs.RT.z_topBottom(1)), ' km, ',...
    'cloud base $= $', num2str(inputs.RT.z_topBottom(2)), ' km, ',...
    '$sza = $', num2str(inputs.RT.sza), '$^{\circ}$, ',...
    '$vza = $', num2str(inputs.RT.vza), '$^{\circ}$'];




% ------------------------------------
% ----- Figure 1: TOA radiance -----
% ------------------------------------
fig1 = figure;

% Convert radiance to W/m^2/nm/sr
for tt = 1:length(inputs.RT.tau_c)

    plot(wl_center_nm, 1e-3 .* Rad_TOA(:, tt, rr_2plot), '-', 'LineWidth', 2, 'Color', line_colors(tt,:))
    hold on

end

grid on; grid minor
xlim([min(wl_center_nm), max(wl_center_nm)])
xlabel('Wavelength ($nm$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel('TOA Radiance ($W \, m^{-2} \, nm^{-1} \, sr^{-1}$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
title('TOA Radiance over HySICS spectral range - homogeneous liquid water cloud', 'Interpreter', 'latex',...
    'FontSize',ttl_fnt)
subtitle(scene_str, 'Interpreter', 'latex', 'FontSize', sub_ttl_fnt)
legend(legend_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', lgnd_fnt)
set(gcf, 'Position', [0 0 1200 700])

% -------------------------------------
% ---------- Save figure --------------
% save .fig file
saveas(fig1,[folder_paths.saved_figs,'TOA Radiance over HySICS spectral range - homogeneous liquid water cloud.fig']);


% save .png with 500 DPI resolution
exportgraphics(fig1,[folder_paths.saved_figs,'TOA Radiance over HySICS spectral range - homogeneous liquid water cloud.jpg'],...
    'Resolution', 500);
% -------------------------------------
% -------------------------------------








% ------------------------------------
% ---- Figure 2: TOA reflectance ----
% ------------------------------------
fig2 = figure;

for tt = 1:length(inputs.RT.tau_c)

    plot(wl_center_nm, Refl_TOA(:, tt, rr_2plot), '-', 'LineWidth', 2, 'Color', line_colors(tt,:))
    hold on

end

grid on; grid minor
xlim([min(wl_center_nm), max(wl_center_nm)])
xlabel('Wavelength ($nm$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
% R = pi*L / (mu0*F0) is unitless
ylabel('TOA Reflectance', 'Interpreter', 'latex', 'FontSize', ax_fnt)
title('TOA Reflectance over HySICS spectral range - homogeneous liquid water cloud', 'Interpreter', 'latex',...
    'FontSize', ttl_fnt)
subtitle(scene_str, 'Interpreter', 'latex', 'FontSize', sub_ttl_fnt)
legend(legend_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize',lgnd_fnt)
set(gcf, 'Position', [0 0 1200 700])

% -------------------------------------
% ---------- Save figure --------------
% save .fig file
saveas(fig2,[folder_paths.saved_figs,'TOA Reflectance over HySICS spectral range - homogeneous liquid water cloud.fig']);


% save .png with 500 DPI resolution
exportgraphics(fig2,[folder_paths.saved_figs,'TOA Reflectance over HySICS spectral range - homogeneous liquid water cloud.jpg'],...
    'Resolution', 500);
% -------------------------------------
% -------------------------------------
