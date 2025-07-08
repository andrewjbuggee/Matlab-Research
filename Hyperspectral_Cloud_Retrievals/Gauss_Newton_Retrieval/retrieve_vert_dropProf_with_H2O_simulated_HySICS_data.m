%% Retreive vertical profiles in a loop using simulated HySICS reflectance measurements
% The loop can be used to change pixels or to change retrieval settings

% This script uses my own TBLUT algorithm as the apriori values



% By Andrew John Buggee

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


%% Define the HySICS folders for the machine you're using



% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ***** Define the HySICS Folder with the simulated measurements *****
    folder_paths.HySICS_simulated_spectra = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];

    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Droplet_profile_retrievals/'];

    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/';


    % water cloud file location
    folder_paths.water_cloud_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';



elseif strcmp(which_computer,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % ***** Define the HySICS Folder with the simulated measurements *****
    folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];

    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/HySICS/'];

     % water cloud file location
    folder_paths.water_cloud_folder_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/'];




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


    % Define the HySICS simulated spectrum folder

    folder_paths.HySICS_simulated_spectra = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra';


    % water cloud file location
    folder_paths.water_cloud_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';

    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/scratch/alpine/anbu8374/HySICS/';

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    folder_paths.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';

end

% If the folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_inp, 'dir')

    mkdir(folder_paths.libRadtran_inp)

end





%%   Delete old files?
% First, delete files in the HySICS folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.water_cloud_folder_path, '*.DAT'])



%% LOAD SIMULATED HYSICS DATA

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-14-May-2025_rev1.mat';  % sza = 0, vza = 0
        


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

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 47
    % bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_47Bands_20mm-totalColumnWaterVapor_sim-ran-on-07-Jul-2025_rev1';
    
    % test file with just 5 wavelengths
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_5wavelength_test_sim-ran-on-10-Jun-2025_rev1.mat';

    


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 47
    % bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-07-Jul-2025_rev1.mat';



end


simulated_measurements = load([folder_paths.HySICS_simulated_spectra,filename]);


%% Compute the Two-Band Look-up Table retrieval of effective radius and optical depth

tic
%tblut_retrieval = TBLUT_for_HySICS(simulated_measurements, folder_paths);
tblut_retrieval = TBLUT_for_HySICS_ver2(simulated_measurements, folder_paths);
toc



%% CREATE GAUSS-NEWTON INPUTS

% We use the estimates calcualted by the TBLUT as our a priori
GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS_ver2(simulated_measurements);
disp('Dont forget to check the inputs and change if needed!!')


%% We're retrieving above cloud column water vapor. Make sure input settings are correct

GN_inputs.RT.modify_total_columnWaterVapor = false;             % don't modify the full column
GN_inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % modify the column above the cloud

%% CREATE MODEL PRIOR AND COVARIANCE MATRIX AND MEASUREMENT COVARIANCE

% I don't need anything but the covariance matrix and the expected values
%inputs = create_model_prior(inputs,data_inputs);

% -------------------------------------------------------
% do you want to use your estimates or the MODIS estimate?
% -------------------------------------------------------

use_TBLUT_estimates = true;

GN_inputs = create_model_prior_covariance_HySICS_ver2(GN_inputs, tblut_retrieval, use_TBLUT_estimates);

GN_inputs = create_HySICS_measurement_covariance(GN_inputs, simulated_measurements);


%% CALCULATE RETRIEVAL PARAMETERS

tic

% --------------------------------------------------------------
% ---------------- Retrieve Vertical Profile! ------------------
% --------------------------------------------------------------
[GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_HySICS(GN_inputs, simulated_measurements, folder_paths);
% --------------------------------------------------------------
% --------------------------------------------------------------

toc

%% Save the Outputs!
rev = 1;
if modisInputs.flags.useAdvection == true

    filename = [modisInputs.savedCalculations_folderName, 'GN_inputs_outputs_withAdvection_rt-cov_',num2str(r_top_apriori_percentage*100),...
        '_rb-cov_', num2str(r_bot_apriori_percentage_vector(rb)*100),'_tc-cov_', num2str(tau_c_apriori_percentage_vector(tc)*100),...
        '_',char(datetime("today")), '_rev', num2str(rev), '.mat'];

    % Check to see if this file name already exists
    while isfile(filename)==true

        rev = rev+1;

        filename = [modisInputs.savedCalculations_folderName, 'GN_inputs_outputs_withAdvection_rt-cov_',num2str(r_top_apriori_percentage*100),...
            '_rb-cov_', num2str(r_bot_apriori_percentage_vector(rb)*100),'_tc-cov_', num2str(tau_c_apriori_percentage_vector(tc)*100),...
            '_',char(datetime("today")), '_rev', num2str(rev), '.mat'];

    end

    save(filename,"GN_outputs","GN_inputs", "vocalsRex", "modisInputs",...
        "r_top_apriori_percentage", "r_bot_apriori_percentage", "tau_c_apriori_percentage");



else

    filename = [modisInputs.savedCalculations_folderName,'GN_inputs_outputs_withoutAdvection__rt-cov_',num2str(r_top_apriori_percentage*100),...
        '_rb-cov_', num2str(r_bot_apriori_percentage_vector(rb)*100),'_tc-cov_', num2str(tau_c_apriori_percentage_vector(tc)*100),...
        '_',char(datetime("today")), '_rev', num2str(rev),'.mat'];

    % Check to see if this file name already exists
    while isfile(filename)==true

        rev = rev+1;

        filename = [modisInputs.savedCalculations_folderName, 'GN_inputs_outputs_withAdvection_rt-cov_',num2str(r_top_apriori_percentage*100),...
            '_rb-cov_', num2str(r_bot_apriori_percentage_vector(rb)*100),'_tc-cov_', num2str(tau_c_apriori_percentage_vector(tc)*100),...
            '_',char(datetime("today")), '_rev', num2str(rev), '.mat'];

    end

    save(filename,"GN_outputs","GN_inputs", "vocalsRex", "modisInputs",...
        "r_top_apriori_percentage", "r_bot_apriori_percentage", "tau_c_apriori_percentage");



end



toc
%% PLOT RETRIEVED VERTICAL PROFILE WITH MODIS RETRIEVAL

modis_pixel_2_plot = 1;
plot_vocalsRex_with_MODIS_retrieved_re_and_vertProf_retrieval(vocalsRex, modis, modisInputs, GN_outputs, GN_inputs, modis_pixel_2_plot)


%% FIND ALL FILES WHERE R_TOP AND R_BOT COV VARY AND MAKE PLOTS

listing = dir([modisInputs.savedCalculations_folderName]);

% save all posterior covariance matrices
retreived_cov = [];

% which pixel would you like to plot?
% This is an index value
modis_pixel_2_plot = 2;

% compute the L2 norm value of the variance of each retrieved variable
L2_mag_total_var = nan(1, length(listing));

% create an empty array for the RMS residual
rms_residual = zeros(length(listing), length(vocalsRex.modisIndex_minDist));


% loop through and read covariance of retrieved variables
for nn = 1:length(listing)

    % We're looking for files with a unique covariance matrix. These files
    % have names longer than 60 characters
    if length(listing(nn).name)>=57


        % yes, it is a file that was run with a changing covariance
        % load the data
        d = load(listing(nn).name);


        % read the retrieval covaraince
        retreived_cov = cat(3, retreived_cov, d.GN_outputs.posterior_cov(:,:,modis_pixel_2_plot));


        % to determine which file had the lowest overall variance
        % between all of the retrieved variables, we need to compute
        % the L2 for each file. If no file, leave as zero.
        L2_mag_total_var(nn) = sqrt(retreived_cov(1,1,end).^2 + retreived_cov(2,2,end).^2 + retreived_cov(3,3,end).^2);


        % plot the retrieved profile
        plot_vocalsRex_with_MODIS_retrieved_re_and_vertProf_retrieval(vocalsRex, modis, modisInputs, d.GN_outputs, d.GN_inputs, modis_pixel_2_plot)


        % Plot the RMS residual versus the r_top and tau_c covariance
        %             r_top_apriori_percentage(nn) = d.r_top_apriori_percentage;
        %             r_bot_apriori_percentage(nn) = d.r_bot_apriori_percentage;
        %             tau_c_apriori_percentage(nn) = d.tau_c_apriori_percentage;

        % Store all values in an array for each file

        for pp = 1:length(d.GN_outputs.rms_residual)

            rms_residual(nn, pp) = d.GN_outputs.rms_residual{pp}(end);

        end



    else

        rms_residual(nn,:) = nan;

    end

end

% Find the minimum rms residual
[min_rms, min_rms_idx] = min(rms_residual, [], 'all');
% find the file and pixel associated with the minimum rms
[file_num, pixel_idx] = ind2sub(size(rms_residual), min_rms_idx);
disp([newline, 'File with smallest rms residual: ', listing(file_num).name,newline,...
    'Pixel with smallest rms residual: ', num2str(pixel_idx), newline,...
    'The minimum rms residual is: ', num2str(round(min_rms, 6))])


% find the collective minimum variance of the retrieved variables
% first set

[min_val, min_idx] = min(L2_mag_total_var);
