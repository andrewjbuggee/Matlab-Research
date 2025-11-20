%% This script is designed for running many retrievals simultaneously using SLURMs job array



% By Andrew John Buggee
%%

% clear variables
clear variables;



%% Define the folders for libRadtran calculation and define the computer in use

folder_paths = define_folderPaths_for_HySICS(2);
which_computer = folder_paths.which_computer;


%% Would you like to print status updates and/or the libRadtran error file?

print_status_updates = false;

print_libRadtran_err = false;


%% LOAD SIMULATED HYSICS DATA

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


    % % load all filenames in the folder defined above.
    % %     filenames = dir([folder_paths.HySICS_simulated_spectra, '*.mat']);
    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);




    % define the folder where the spectra are located
    folder_paths.HySICS_simulated_spectra = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-16-Sep-2025.mat']);



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the spectra are located
    % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/vza_7_subset_pt1percent/'];

    % define the folder where the spectra are located
    folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %      'simulated_HySICS_reflectance_35bands_with_1%_uncertainty_sim-ran-on-12-Jul-2025_rev1.mat']);

    % load all filenames in the folder defined above.
    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);
    % % filenames = dir([folder_paths.HySICS_simulated_spectra, '*.mat']);


    % Use folder on MATLAB Drive
    % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Simulated_spectra/',...
    %     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91/'];


    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_5_tcwv_14_vza_7*.mat']);

    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.1%_uncert_rTop_10_rBot_5_tauC_5_tcwv_14*.mat']);


    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);
    % filenames = dir([folder_paths.HySICS_simulated_spectra, '*.mat']);






    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14',...
        '_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-16-Sep-2025.mat']);



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


    % load measurements from curc storage
    % folder_paths.HySICS_simulated_spectra = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
    %     'HySICS/Simulated_spectra/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/'];

    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_8_vza_7*.mat']);



    folder_paths.HySICS_simulated_spectra = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/paper2_variableSweep/rTop_10/test_folder_mimic_localMachine_results/'];

    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_HySICS_reflectance_35bands*.mat']);





    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;

    % load all filenames in the folder defined above that start with a
    % specific block of text









end

%% Run the droplet profile retrieval

files{1} = filenames.name;


% --------------------------------------------
% --- Assume 10mm total column water vapor ---
% --------------------------------------------
% *** Retrieve r_top, r_bot, tau_c ***
[tblut_retrieval_10, GN_inputs_10, GN_outputs_10] = run_retrieval_dropletProfile_HySICS_ver4_logNewCov_noACPW_10(files,...
    folder_paths, print_status_updates, print_libRadtran_err);


% --------------------------------------------
% --- Assume 15mm total column water vapor ---
% --------------------------------------------
% *** Retrieve r_top, r_bot, tau_c ***
[tblut_retrieval_15, GN_inputs_15, GN_outputs_15] = run_retrieval_dropletProfile_HySICS_ver4_logNewCov_noACPW_15(files,...
    folder_paths, print_status_updates, print_libRadtran_err);



% --------------------------------------------
% --- Assume 20mm total column water vapor ---
% --------------------------------------------
% *** Retrieve r_top, r_bot, tau_c ***
[tblut_retrieval_20, GN_inputs_20, GN_outputs_20] = run_retrieval_dropletProfile_HySICS_ver4_logNewCov_noACPW_20(files,...
    folder_paths, print_status_updates, print_libRadtran_err);



% --------------------------------------------
% --- Assume 25mm total column water vapor ---
% --------------------------------------------
% *** Retrieve r_top, r_bot, tau_c ***
[tblut_retrieval_25, GN_inputs_25, GN_outputs_25] = run_retrieval_dropletProfile_HySICS_ver4_logNewCov_noACPW_25(files,...
    folder_paths, print_status_updates, print_libRadtran_err);


