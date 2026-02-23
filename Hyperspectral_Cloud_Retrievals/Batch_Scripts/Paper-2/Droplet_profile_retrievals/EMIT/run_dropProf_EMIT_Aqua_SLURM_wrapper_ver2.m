%% New SLURM Wrapper for EMIT-Aqua overlapping pixel retrieval with enhanced parallelization

clear variables


%%

print_status_updates = true;
print_libRadtran_err = true;
plot_figures = false; 
save_figures = false; 



%% Define the folder of the coincident data set between EMIT and Aqau

% ---------------------------------------------
% ---------- PICK COINCIDENT DATA SET  --------
% ---------------------------------------------

which_computer = whatComputer;

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    mat_file_path = [];



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    mat_file_path = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/',...
        'Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/overlap_data_per_pixel/',...
        'overlap_EMIT_pixel_005_2024_5_17_T183918_1.mat'];


    % *** Define output directory ***
    rev = 9;
    output_dir = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/',...
        'Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/Droplet_profile_retrievals/take_', num2str(rev)];

    if ~exist(output_dir, 'dir')
        mkdir(output_dir)

    else

        while exist(output_dir, 'dir')==7

            rev = rev + 1;
            output_dir = [output_dir(1:end-1), num2str(rev), '/'];

        end

        mkdir(output_dir)

    end


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;



end



%% Run retrieval!

% this should be the SLURM_ARRAY_TASK_ID number
folder_extension_number = 5;

% This is the final number in the directory name where the output will be
% stored
folder_rev_num = 9;


[GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = ...
    run_retrieval_singlePixel_EMIT_Aqua(mat_file_path, folder_extension_number, ...
    print_status_updates, print_libRadtran_err, folder_rev_num, output_dir);