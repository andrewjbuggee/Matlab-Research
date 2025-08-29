%% This script is designed for running many retrievals simultaneously using SLURMs job array



% By Andrew John Buggee
%%

% clear variables
clear variables;

% add folders to the path
addpath(genpath('/projects/anbu8374/Matlab-Research'));
addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); 
addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); 
addLibRadTran_paths; 

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
  filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);
 


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % load all filenames in the folder defined above.
    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);
    % filenames = dir([folder_paths.HySICS_simulated_spectra, '*.mat']);



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % load all filenames in the folder defined above that start with a
    % specific block of text
    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);


end

%% Run the droplet profile retrieval


[tblut_retrieval, acpw_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropletProfile_HySICS_ver2(filenames.name, folder_paths, print_status_updates, print_libRadtran_err);


