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


% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.


% ------------------------------------------------
% ------ Folders on the CU Super Computer --------
% ------------------------------------------------


% Define the HySICS folder name

folder_paths.HySICS_simulated_spectra = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
    'HySICS/Simulated_spectra/'];

% ---- Define where the retrievals will be stored ---
folder_paths.HySICS_retrievals = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
    'HySICS/Droplet_profile_retrievals/'];

% Define the folder path where all .INP files will be saved
folder_paths.libRadtran_inp = ['/scratch/alpine/anbu8374/HySICS/INP_OUT/'];


% water cloud file location
folder_paths.water_cloud_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';






%%   Delete old files?
% First, delete files in the HySICS folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.water_cloud_folder_path, '*.DAT'])



%% LOAD SIMULATED HYSICS DATA


% ------------------------------------------------
% ------ Folders on the CU Super Computer --------
% ------------------------------------------------

filename = 'simulated_measurement_HySICS_reflectance_35-nonWaterVaporBands_10mm_total-WV-content_sim-ran-on-11-Jul-2025.mat';



simulated_measurements = load([folder_paths.HySICS_simulated_spectra,filename]);


%%  *** Start parallel pool ***

% Is parpool running?
p = gcp('nocreate');
if isempty(p)==true

    % first read the local number of workers avilabile.
    p = parcluster('local');
    % start the cluster with the number of workers available
    if p.NumWorkers>64
        % Likely the amilan128c partition with 2.1 GB per core
        % Leave some cores for overhead
        parpool(p.NumWorkers - 8);

    elseif p.NumWorkers<=64 && p.NumWorkers>10

        parpool(p.NumWorkers);

    elseif p.NumWorkers<=10

        parpool(p.NumWorkers);

    end

end


%% Create the name of the file to save all output to

rev = 1;

folder_paths.saveOutput_filename = [folder_paths.HySICS_retrievals,'dropletRetrieval_HySICS_', num2str(numel(simulated_measurements.inputs.bands2run)),...
    'bands_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];



while isfile(filename)
    rev = rev+1;
    if rev<10
        filename = [filename(1:end-5), num2str(rev),'.mat'];
    elseif rev>10
        filename = [filename(1:end-6), num2str(rev),'.mat'];
    end
end



%% Compute the Two-Band Look-up Table retrieval of effective radius and optical depth

tic
tblut_retrieval = TBLUT_for_HySICS_ver2(simulated_measurements, folder_paths);
toc



%% CREATE GAUSS-NEWTON INPUTS

% We use the estimates calcualted by the TBLUT as our a priori
GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS(simulated_measurements);
disp('Dont forget to check the inputs and change if needed!!')



%% CREATE MODEL PRIOR AND COVARIANCE MATRIX AND MEASUREMENT COVARIANCE

% I don't need anything but the covariance matrix and the expected values
%inputs = create_model_prior(inputs,data_inputs);

% -------------------------------------------------------
% do you want to use your estimates or the MODIS estimate?
% -------------------------------------------------------

use_TBLUT_estimates = true;

GN_inputs = create_model_prior_covariance_HySICS(GN_inputs, tblut_retrieval, use_TBLUT_estimates);

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



%%
% ----------------------------------------------
% ------------ SAVE OUTPUT STRUCTURE -----------
% ----------------------------------------------

% Save the version without an measurement uncertainty. Then we can add
% uncertainty and save the new file


% If the folder path doesn't exit, create a new directory
if ~exist(folder_paths.HySICS_retrievals, 'dir')

    mkdir(folder_paths.HySICS_retrievals)

end

if exist(folder_paths.saveOutput_filename, 'file')==true
    % append
    save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "simulated_measurements", "folder_paths", '-append');

else
    save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "simulated_measurements", "folder_paths");

end





toc
