%% SLURM wrapper for run_dropProf_EMIT with coincident EMIT and Terra data



print_status_updates = true;
print_libRadtran_err = false;
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
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Terra_data/',...
        'southEast_pacific/likely_stratus/'];


    % 11 Pixels with H less than 1.6     ** Use this data set **
    folder_paths.coincident_dataFolder = '2023_3_6_T151922/';



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Terra_data/',...
        'southEast_pacific/likely_stratus/'];

    % 2 Pixels with H less than 0.85
    % folder_paths.coincident_dataFolder = '2023_3_6_T151922/';

    % 3 Pixels with H less than 1.1
    folder_paths.coincident_dataFolder = '2023_3_6_T152010/';




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end

%%

[GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] =...
    run_dropProf_acpw_retrieval_EMIT_overlap_Terra_ver5(folder_paths, print_status_updates, print_libRadtran_err,...
    plot_figures, save_figures);