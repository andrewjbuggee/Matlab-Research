%% SLURM wrapper for run_dropProf_EMIT with coincident EMIT and Terra data



print_status_updates = true;
print_libRadtran_err = true;
plot_figures = false; 
save_figures = false; 


folder_paths = define_EMIT_dataPath_and_saveFolders(3);



%% Define the folder of the coincident data set between EMIT and Terra

% ---------------------------------------------
% ---------- PICK COINCIDENT DATA SET  --------
% ---------------------------------------------

which_computer = folder_paths.which_computer;

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folder_paths.coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Terra_data/southEast_pacific/'];

    % folder_paths.coincident_dataFolder = '2024-09-12/';


    % 11 Pixels with H less than 1.6     ** Use this data set **
    folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Terra_data/southEast_pacific/'];



    % 10 Pixels with H less than 1.85      ** Use this data set **
    % But only 2 aren't masked out by EMIT cloud filter!
    folder_paths.coincident_dataFolder = '2024_5_17_T183930_1/';




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Terra_data/southEast_pacific/'];


    % 10 Pixels with H less than 1.85      ** Use this data set **
    % But only 2 aren't masked out by EMIT cloud filter!
    folder_paths.coincident_dataFolder = '2024_5_17_T183930_1/';

end

%%

[GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] =...
    run_dropProf_acpw_retrieval_EMIT_overlap_Terra_ver5(folder_paths, print_status_updates, print_libRadtran_err,...
    plot_figures, save_figures);