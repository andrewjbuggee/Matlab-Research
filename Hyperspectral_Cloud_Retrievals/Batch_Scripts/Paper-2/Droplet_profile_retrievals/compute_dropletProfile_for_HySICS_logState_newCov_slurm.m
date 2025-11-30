%% This script is designed for running many retrievals simultaneously using SLURMs job array



% By Andrew John Buggee
%%

% clear variables
clear variables;



%% Define the folders for libRadtran calculation and define the computer in use

folder_paths = define_folderPaths_for_HySICS(3);
which_computer = folder_paths.which_computer;


%% Would you like to print status updates and/or the libRadtran error file?

print_status_updates = true;

print_libRadtran_err = true;


%% LOAD SIMULATED HYSICS DATA

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


    % load all filenames in the folder defined above.
    %     filenames = dir([folder_paths.HySICS_simulated_spectra, '*.mat']);
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

    % load all filenames in the folder defined above.
    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);
    % % filenames = dir([folder_paths.HySICS_simulated_spectra, '*.mat']);


    % Use folder on MATLAB Drive
    % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Simulated_spectra/',...
    %     'paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91/'];


    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_5_tcwv_20_vza_7*.mat']);




    % % define the folder where the spectra are located
    % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

    % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/vza_7_subset_1percent/'];

    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_HySICS_reflectance_66bands_with_1%_uncertainty_sim-ran-on-12-Jul-2025_rev1.mat']);


    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14',...
    %     '_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-16-Sep-2025.mat']);



    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/',...
        'log_newCov_subset_allBands_VR_inSitu_1/'];

    % all bands
    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_636bands_0.3%_uncert_vocalsRex_inSitu_re_lwc_tauC_z_',...
    %     '31-Oct-2008_9.0864UTC_prof-nn_24_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-28-Nov-2025.mat']);

    % all bands - testing failed file
    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_636bands_0.3%_uncert_vocalsRex_inSitu_re_lwc_tauC_z_',...
    %     '04-Nov-2008_6.09UTC_prof-nn_35_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-28-Nov-2025.mat']);


    % using only 66 bands 
    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.3%_uncert_vocalsRex_inSitu_re_lwc_tauC_z_',...
    %     '15-Oct-2008_18.4297UTC_prof-nn_1_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-29-Nov-2025.mat']);


    % 66 bands - ensemble profile 35
    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.3%_uncert_vocalsRex_inSitu_re_lwc_tauC_z_',...
        '04-Nov-2008_6.09UTC_prof-nn_35_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-30-Nov-2025.mat']); 



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


    % load measurements from curc storage
    folder_paths.HySICS_simulated_spectra = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;

    % load all filenames in the folder defined above that start with a
    % specific block of text


    % filenames = dir([folder_paths.HySICS_simulated_spectra,...
    %     'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_5_tcwv_20_vza_7*.mat']);



%     filenames = dir([folder_paths.HySICS_simulated_spectra,...
%         'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);


filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.3%_uncert_rTop_9.2516_rBot_5.3192_tauC_6.1312_tcwv_14',...
        '_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-16-Sep-2025.mat']);


end

%% Set the directory of the retrieved droplet profiles

if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Droplet_profile_retrievals/paper2_variableSweep/test_logSpace_newCov_with_VR_inSitu_meas/'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end 

%% Run the droplet profile retrieval

files{1} = filenames.name;

% % *** Retrieve r_top, r_bot, tau_c, and cwvs ***
% [tblut_retrieval, acpw_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropletProfile_HySICS_ver4_log_newCov(files,...
%     folder_paths, print_status_updates, print_libRadtran_err);


% *** Retrieve r_top, r_bot, tau_c, and cwvs ***
[tblut_retrieval, acpw_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropProf_HySICS_ver4_log_newCov_with_forMo_uncert(...
    files,folder_paths, print_status_updates, print_libRadtran_err);


