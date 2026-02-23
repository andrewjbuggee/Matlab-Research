%% Find overlapping pixels between EMIT and MODIS, AIRS and AMSR-E on board Aqua

clear variables

print_status_updates = false;
print_libRadtran_err = false;
plot_figures = false;
save_figures = false;


folder_paths = define_EMIT_dataPath_and_saveFolders(5);


which_computer = folder_paths.which_computer;

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folder_paths.coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];






elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];


    %     % Grab filenames in drive
    %     sub_directories = dir(folder_paths.coincident_dataPath);
    %     idx_2delete = [];
    %     for nn = 1:length(sub_directories)
    %
    %         if strcmp(sub_directories(nn).name(1), '2')~=true
    %
    %             idx_2delete = [idx_2delete, nn];
    %
    %         end
    %
    %     end
    %
    % % delete rows that don't have retrieval filenames
    % sub_directories(idx_2delete) = [];

    % folder_paths.sub_directories = {'2023_9_16_T191118_1/', '2023_9_16_T191130_1/', '2023_9_16_T191142_1/',...
    %                                        '2024_1_13_T194658_1/', '2024_5_17_T183906_1/', '2024_5_17_T183918_1/',...
    %                                        '2024_5_17_T183930_1/'};

    sub_directories = {'2023_9_16_T191118_1/', '2023_9_16_T191130_1/', '2023_9_16_T191142_1/',...
        '2024_1_13_T194658_1/', '2024_5_17_T183906_1/', '2024_5_17_T183918_1/',...
        '2024_5_17_T183930_1/'};

    ver = 1;
    folder_paths.output_dir = [folder_paths.coincident_dataPath, 'overlap_data_per_pixel_ver', num2str(ver), '/'];


    if ~exist(folder_paths.output_dir, 'dir')
        mkdir(folder_paths.output_dir)

    else

        while exist(folder_paths.output_dir, 'dir')==7

            ver = ver + 1;
            folder_paths.output_dir = [folder_paths.coincident_dataPath, 'overlap_data_per_pixel_ver', num2str(ver), '/'];

        end

    end


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    % addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    % addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    sub_directories = {'2023_9_16_T191118_1/', '2023_9_16_T191130_1/', '2023_9_16_T191142_1/',...
        '2024_1_13_T194658_1/', '2024_5_17_T183906_1/', '2024_5_17_T183918_1/',...
        '2024_5_17_T183930_1/'};




    % IF running this on the CU supercomputer, store the data on the scratch
    % directory. Scratch has 10 TB of storage and is designed for high
    % performance reading/writing, but everything stored is scratch is deleted
    % each 90 days.
    ver = 1;

    % store on scratch
    folder_paths.output_dir = ['/scratch/alpine/anbu8374/EMIT_pix_overlap_with_Aqua_paper2_ver', num2str(ver), '/'];


    if ~exist(folder_paths.output_dir, 'dir')
        mkdir(folder_paths.output_dir)

    else

        while exist(folder_paths.output_dir, 'dir')==7

            ver = ver + 1;
            folder_paths.output_dir = ['/scratch/alpine/anbu8374/EMIT_pix_overlap_with_Aqua_paper2_ver', num2str(ver), '/'];

        end

    end


end








%%

% Define overlap criteria

criteria.cld_phase = 'water';
criteria.cld_cvr = 1;              % cloud fraction
criteria.cld_tau_min = 3;          % cloud optical depth minimum
criteria.cld_tau_max = 50;         % cloud optical depth maximum
criteria.H = 2;                    % horizontal inhomogeneity index
criteria.findN_smallest_H = true;
criteria.H_N_smallest = 1;         % keep the 5 MODIS pixels with smallest H

emit_pixels_per_modis = 2;        % select 30 EMIT pixels within each MODIS pixel


for nn = 1:length(sub_directories)


    n_saved = save_overlap_data_perPixel_EMIT_Aqua(criteria, folder_paths, true, emit_pixels_per_modis,...
        sub_directories{nn});

end


