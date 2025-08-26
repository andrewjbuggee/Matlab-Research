%% This script has been created to develop an above cloud column water vapor retrieval to 
% use as an a priori for my hyperspectral retrieval algorithm



% By Andrew John Buggee

%%

clear variables

% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


%%
% Would you like each function to print messages to the command terminal on
% the status of the computation?
print_status_updates = true;

% Would you like to print long libRadtran error messages?
print_libRadtran_err = true;


%% Define the HySICS folders for data and storage

[folder_paths, which_computer] = define_folderPaths_for_HySICS(0);



%% -- Start Parallel pool

start_parallel_pool(which_computer)


%%   Delete old files?

% First, delete files in the HySICS INP folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% delete old water vapor profiles
delete([folder_paths.atm_folder_path, '*-aboveCloud.DAT'])

% delete old MIE files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])


%% LOAD SIMULATED HYSICS DATA

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------




    %     filename = ['simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile',...
    %         '_66Bands_20mm-aboveCloud-WV_sim-ran-on-08-Jul-2025_rev1'];  % sza = 0, vza = 0

    % r_top = 10, r_bot = 5, tau_c = 11, tcwv = 14mm, 66 bands from first paper with 0.001%
    % uncertainty, viewing geometry based on EMIT measurement from 27
    % January, 2024
    filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11',...
        '_tcwv_14_vza_7_vaz_210_sza_10_saz_91_sim-ran-on-14-Aug-2025.mat'];



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev1.mat']); % sza = 10, vza = 0

    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-12-May-2025_rev1.mat']); % sza = 0, vza = 0

    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev2.mat']); % sza = 20, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev3.mat']); % sza = 30, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev4.mat']); % sza = 40, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev5.mat']); % sza = 50, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev6.mat']); % sza = 60, vza = 0

    % r_top = 9.5, r_bot = 4, tau_c = 6
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-17-Jun-2025_rev1.mat';

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 47 bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_47Bands_20mm-totalColumnWaterVapor_sim-ran-on-07-Jul-2025_rev1';

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 66
    % Bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_66Bands_20mm-totalColumnWaterVapor_sim-ran-on-08-Jul-2025_rev1';
  

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, ALL bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_allBands_20mm-totalColumnWaterVapor_sim-ran-on-08-Jul-2025_rev1';


    % r_top = 9.5, r_bot = 4, tau_c = 6, 66 bands from first paper with 1%
    % uncertainty
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_HySICS_reflectance_66bands_with_1%_uncertainty_sim-ran-on-12-Jul-2025_rev1.mat';


    % r_top = 10, r_bot = 3, tau_c = 29, tcwv = 15mm, 66 bands from first paper with 0.001%
    % uncertainty
    % filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_3_tauC_29_tcwv_15_vza_0',...
    %     '_vaz_210_sza_0_saz_0_sim-ran-on-12-Aug-2025.mat'];


    % r_top = 10, r_bot = 5, tau_c = 8, tcwv = 14mm, 66 bands from first paper with 0.001%
    % uncertainty - nadir viewing with overhead sun
    % filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_8_tcwv_14_vza_0',...
    %     '_vaz_210_sza_0_saz_0_sim-ran-on-12-Aug-2025.mat'];

    % r_top = 10, r_bot = 5, tau_c = 8, tcwv = 14mm, 66 bands from first paper with 0.001%
    % uncertainty, viewing geometry based on EMIT measurement from 27
    % January, 2024
    % filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_8_tcwv_14_vza_7',...
    %     '_vaz_210_sza_10_saz_91_sim-ran-on-14-Aug-2025.mat'];


    % r_top = 10, r_bot = 5, tau_c = 11, tcwv = 14mm, 66 bands from first paper with 0.001%
    % uncertainty, viewing geometry based on EMIT measurement from 27
    % January, 2024
    filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7',...
        '_vaz_210_sza_10_saz_91_sim-ran-on-14-Aug-2025.mat'];

    % test file with just 5 wavelengths
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_5wavelength_test_sim-ran-on-10-Jun-2025_rev1.mat';




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 47
    % bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    %filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-07-Jul-2025_rev1.mat';


    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 66
    % bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_66Bands_20mm-aboveCloud-WV_sim-ran-on-08-Jul-2025_rev1.mat';


    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, all
    % bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_allBands_20mm-totalColumnWaterVapor_sim-ran-on-08-Jul-2025_rev1.mat';


    % r_top = 10, r_bot = 5, tau_c = 8, total_column_waterVapor = 14, 66
    % bands
    % filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_8_tcwv_14_vza_0_vaz_210',...
    %     '_sza_0_saz_0_sim-ran-on-12-Aug-2025.mat'];

    % r_top = 10, r_bot = 5, tau_c = 10, total_column_waterVapor = 20, 66
    % bands
    filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_10_tauC_11_tcwv_11_vza_0',...
    '_vaz_210_sza_0_saz_0_sim-ran-on-12-Aug-2025.mat'];


    % r_top = 10, r_bot = 5, tau_c = 10, total_column_waterVapor = 20, 66
    % bands
    % filename = ['simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_10_tcwv_20_vza_4_vaz_257',...
    %     '_sza_31_saz_96_sim-ran-on-11-Aug-2025.mat'];


end


simulated_measurements = load([folder_paths.HySICS_simulated_spectra,filename]);



% *** Check to see if these measure have added uncertainty or not ***

if isfield(simulated_measurements, 'Refl_model_with_noise')==true

    disp([newline, 'Using measurements with added uncertianty...', newline])

    % Then we're using measurements with noise and we set this to be the
    % Reflectance measurements
    simulated_measurements.Refl_model = simulated_measurements.Refl_model_with_noise;

end


%% Create the name of the file to save all output to

rev = 1;

folder_paths.saveOutput_filename = [folder_paths.HySICS_retrievals,'acpw_HySICS_',...
    'bands_',num2str(100*simulated_measurements.inputs.measurement.uncert), '%_uncert',...
        '_rTop_', num2str(simulated_measurements.inputs.RT.r_top),...
        '_rBot_', num2str(simulated_measurements.inputs.RT.r_bot),...
        '_tauC_', num2str(simulated_measurements.inputs.RT.tau_c),...
        '_tcwv_', num2str(simulated_measurements.inputs.RT.waterVapor_column),...
        '_vza_', num2str(round(simulated_measurements.inputs.RT.vza)),...
        '_vaz_', num2str(round(simulated_measurements.inputs.RT.vaz)),...
        '_sza_', num2str(round(simulated_measurements.inputs.RT.sza)),...
        '_saz_', num2str(round(simulated_measurements.inputs.RT.phi0)),...
        '_sim-ran-on-',char(datetime("today")),'.mat'];



while isfile(filename)
    rev = rev+1;
    if rev<10
        filename = [filename(1:end-5), num2str(rev),'.mat'];
    elseif rev>10
        filename = [filename(1:end-6), num2str(rev),'.mat'];
    end
end


%% Grab the reflectance measurements at three spectral channels


% The window channel where water vapor absorption is negligible
R_window_872 = simulated_measurements.Refl_model(8);        % 1/sr - HySICS channel centered around 881 nm

% The weakly absorbing water vapor absorption channel
R_waterVap_900 = simulated_measurements.Refl_model(10);     % 1/sr - HySICS channel centered around 900 nm

% The strongly absorbing water vapor absorption channel
R_waterVap_1127 = simulated_measurements.Refl_model(21);     % 1/sr - HySICS channel centered around 1127 nm

R_measured = [R_window_872; R_waterVap_900; R_waterVap_1127];

%% Load simulated measurements

% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';

% We're modeling a homoegenous cloud layer using the TBLUT retrieval, so
% the number of free re parameters is 1

inputs.RT.num_re_parameters = 1;

% ---- First, let's simulate water clouds ----


% Define the parameters of the INP file

[inputs, spec_response] = create_uvSpec_DISORT_inputs_for_HySICS(inputs, true, simulated_measurements,...
                            'exact', print_libRadtran_err);


% only keep the three wavelenths used in this multispectral estimate
inputs.bands2run = [inputs.bands2run(8); inputs.bands2run(10); inputs.bands2run(21)];

inputs.RT.wavelengths2run = [inputs.RT.wavelengths2run(8, :); inputs.RT.wavelengths2run(10, :);...
                                inputs.RT.wavelengths2run(21, :)];

spec_response.value = [spec_response.value(8, :); spec_response.value(10, :); spec_response.value(21, :)];

spec_response.wavelength = [spec_response.wavelength(8, :); spec_response.wavelength(10, :); spec_response.wavelength(21, :)];

%% Set the total column water vapor?

inputs.RT.modify_total_columnWaterVapor = false;             % modify the full column

inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % don't modify the column above the cloud


%% Store the values used to generate the measurements

% Store the simulated state vector used to create the measurements
inputs.measurement.r_top = simulated_measurements.inputs.RT.r_top;      % microns
inputs.measurement.r_bot = simulated_measurements.inputs.RT.r_bot;      % microns
inputs.measurement.tau_c = simulated_measurements.inputs.RT.tau_c;      % optical depth
inputs.measurement.actpw = aboveCloud_CWV_simulated_hysics_spectra(simulated_measurements.inputs); % kg/m^2 (equivelant to mm)


%% Load the tblut retrieval and set the effective radius and optical depth accordingly

% Load tblut_retrieval
load(['dropletRetrieval_HySICS_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7',...
    '_vaz_210_sza_10_saz_91_sim-ran-on-23-Aug-2022.mat'], 'tblut_retrieval');

inputs.RT.re = tblut_retrieval.minRe;
inputs.RT.tau_c = tblut_retrieval.minTau;

%% Run simulated calculations using the same geometry and TBLUT estimates across many above cloud column water vapor amounts

acpw_sim = 3:0.25:30;    % mm
% num wavelengths
num_wl = length(inputs.bands2run);

% length of each independent variable
num_tcpw = length(acpw_sim);


num_INP_files = num_tcpw * num_wl;

inputFileName = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);


% changing variable steps through tcpw and wavelength
% in for loop speak, it would be:
% for pw = 1:num_tcpw
%   for ww = 1:num_wl
changing_variables_allStateVectors = [reshape(repmat(acpw_sim, num_wl,1), [],1),...
                                     repmat(inputs.RT.wavelengths2run, num_tcpw, 1)];


% Add a final column that includes the index for the spectral response
% function. These always increase chronologically
changing_variables_allStateVectors = [changing_variables_allStateVectors, repmat((1:num_wl)',  num_tcpw, 1)];

% Write the water cloud file


wc_filename = write_wc_file(tblut_retrieval.minRe, tblut_retrieval.minTau,...
    inputs.RT.z_topBottom,inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
    inputs.RT.distribution_var, inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
    inputs.RT.indVar, inputs.compute_weighting_functions, which_computer,...
    1, inputs.RT.num_re_parameters, folder_paths.libRadtran_water_cloud_files,...
    folder_paths.libRadtran_mie_folder);

wc_filename = wc_filename{1};






% Now write all the INP files
parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files


    % set the wavelengths for each file
    wavelengths = changing_variables_allStateVectors(nn, end-2:end-1);

    % create a custom water vapor profile
    custom_waterVapor_profile = alter_aboveCloud_columnWaterVapor_profile(inputs, changing_variables_allStateVectors(nn,1),...
        folder_paths.atm_folder_path);

    % ------------------------------------------------
    % ---- Define the input and output filenames! ----
    % ------------------------------------------------
    % input_names need a unique identifier. Let's give them the nn value so
    % they can be traced, and are writen over in memory


    inputFileName{nn} = [num2str(mean(wavelengths)), '_','nm_re_', num2str(tblut_retrieval.minRe),...
        '_tauC_', num2str(tblut_retrieval.minTau), '_acpw_',num2str(changing_variables_allStateVectors(nn,1)),'mm.INP'];



    outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


    % ------------------ Write the INP File --------------------
    write_INP_file(folder_paths.libRadtran_inp, folder_paths.libRadtran_data, folder_paths.libRadtran_water_cloud_files,...
        inputFileName{nn}, inputs, wavelengths, wc_filename, [], [], custom_waterVapor_profile, []);


end




%% Calculate Reflectance

% Read the solar flux file over the wavelength range specified
wavelength_vec = [min(inputs.RT.wavelengths2run,[],"all"), max(inputs.RT.wavelengths2run, [], "all")];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);   % W/nm/m^2

% we will add and subtract a small fraction of the source file resolution
% to ensure rounding errors don't cause an issue when selecting the
% wavelengths needed from the source file
wl_perturb = inputs.RT.source_file_resolution/3;   % nm



% define only the spec_response so the wavelengths are passed into the
% memory of the parallel for loop
spec_response_value = spec_response.value;


tic



% store the reflectances
Refl_model_allStateVectors = zeros(num_INP_files, 1);


parfor nn = 1:num_INP_files
% for nn = 1:num_INP_files


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------

    % compute INP file
    runUVSPEC_ver2(folder_paths.libRadtran_inp, inputFileName{nn}, outputFileName{nn},which_computer);


    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC_ver2(folder_paths.libRadtran_inp, outputFileName{nn}, inputs,...
        inputs.RT.compute_reflectivity_uvSpec);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    idx_wl = source_wavelength>=(changing_variables_allStateVectors(nn,end-2) - wl_perturb) &...
        source_wavelength<=(changing_variables_allStateVectors(nn, end-1) + wl_perturb);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    [Refl_model_allStateVectors(nn), ~] = reflectanceFunction_ver2(inputs, ds,...
        source_flux(idx_wl), spec_response_value(changing_variables_allStateVectors(nn,end),:));



end


toc

%% Which value is associated with the minimum RMS?

RMS = sqrt( mean( (repmat(R_measured, 1, num_tcpw) - reshape(Refl_model_allStateVectors, num_wl, [])).^2, 1) );

[~, idx_min] = min(RMS);

min_acpw = acpw_sim(idx_min)