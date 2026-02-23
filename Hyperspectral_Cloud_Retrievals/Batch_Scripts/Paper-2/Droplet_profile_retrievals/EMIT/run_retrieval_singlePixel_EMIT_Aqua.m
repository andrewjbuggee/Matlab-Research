%% Run EMIT droplet profile retrieval for a single pixel from a pre-saved .mat file
%
% This function loads a per-pixel .mat file (created by
% save_overlap_data_perPixel_EMIT_Aqua) and runs the droplet profile
% retrieval for that single EMIT pixel. It is designed to be called from a
% SLURM job array where each task processes one pixel independently.
%
% INPUTS:
%   mat_file_path           - Full path to the per-pixel .mat file
%   folder_extension_number - Number used by define_EMIT_dataPath_and_saveFolders
%                             to create unique INP/OUT/wc/mie directories
%                             (typically SLURM_ARRAY_TASK_ID)
%   print_status_updates   - (logical) Print progress to console
%   print_libRadtran_err   - (logical) Print libRadtran error messages
%         folder_rev_num   - Number representing the final character in the
%                            output directory

%
% OUTPUTS:
%   GN_inputs        - Gauss-Newton retrieval input structure
%   GN_outputs       - Gauss-Newton retrieval output structure
%   tblut_retrieval  - TBLUT retrieval results
%   acpw_retrieval   - ACPW retrieval results
%   folder_paths     - Updated folder paths structure
%
% By Andrew John Buggee
%%

function [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = ...
    run_retrieval_singlePixel_EMIT_Aqua(mat_file_path, folder_extension_number, ...
    print_status_updates, print_libRadtran_err, folder_rev_num)


%% Load the per-pixel .mat file

if print_status_updates == true
    disp([newline, 'Loading per-pixel data from: ', mat_file_path, newline])
end

loaded_data = load(mat_file_path);

% The per-pixel .mat files store single-pixel structures with the suffix
% '_pp' (created by save_overlap_data_perPixel_EMIT_Aqua using
% extract_pixel_from_struct). Rename them to the standard variable names
% expected by the retrieval functions.
emit = loaded_data.emit_pp;
modis = loaded_data.modis_pp;
airs = loaded_data.airs_pp;
amsr = loaded_data.amsr_pp;
overlap_pixels = loaded_data.overlap_pixels_pp;
pixel_num = loaded_data.pixel_num;        % always 1 for single-pixel files

% Use the folder_paths from define_EMIT_dataPath_and_saveFolders for
% libRadtran directories (unique per SLURM task), but keep the coincident
% data paths from the saved file
saved_folder_paths = loaded_data.folder_paths;

folder_paths = define_EMIT_dataPath_and_saveFolders(folder_extension_number);

% Restore the coincident data paths from the saved file
folder_paths.coincident_dataPath = saved_folder_paths.coincident_dataPath;
folder_paths.coincident_dataFolder = saved_folder_paths.coincident_dataFolder;
folder_paths.L1B_fileName_emit = saved_folder_paths.L1B_fileName_emit;
folder_paths.L1B_fileName_modis = saved_folder_paths.L1B_fileName_modis;


%% Add paths if running on the supercomputer

which_computer = folder_paths.which_computer;

if strcmp(which_computer, 'curc') == true

    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));

end


%% Check if the EMIT pixel is masked out

if all(isnan(emit.radiance.measurements(:, pixel_num)))

    warning([newline, 'EMIT pixel ', num2str(pixel_num), ' is masked out. Exiting.', newline])

    GN_inputs = [];
    GN_outputs = [];
    tblut_retrieval = [];
    acpw_retrieval = [];

    return

end

disp([newline, 'Retrieving Profile for pixel ', num2str(pixel_num), '...', newline])


%% Set libRadtran INP directory

folder_paths.libRadtran_inp = [folder_paths.libRadtran_inp, 'EMIT_', folder_paths.coincident_dataFolder(1:end-1), '_',...
    folder_paths.L1B_fileName_emit{1}(27:30), '/'];

% Create the directory if it doesn't exist
if ~exist(folder_paths.libRadtran_inp, 'dir')
    mkdir(folder_paths.libRadtran_inp)
end


%% Delete old files

% Delete files in the INP/OUT folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% Delete old wc files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% Delete old atm files
delete([folder_paths.atm_folder_path, '*.DAT'])

% Delete old mie files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])


%% Set output filename


folder_paths.saveOutput_directory = [folder_paths.coincident_dataPath,...
    'Droplet_profile_retrievals/take_', num2str(folder_rev_num), '/'];

% If the directory doesn't exist, create it
if ~exist(folder_paths.saveOutput_directory, 'dir') == true

    mkdir(folder_paths.saveOutput_directory)

% else
% 
%     while exist(folder_paths.saveOutput_directory, 'dir') == true
% 
%         folder_rev_num = folder_rev_num + 1;
% 
%         folder_paths.saveOutput_directory = [folder_paths.coincident_dataPath,...
%             folder_paths.coincident_dataFolder,...
%             'Droplet_profile_retrievals/take_', num2str(folder_rev_num), '/'];
% 
%     end
% 
%     mkdir(folder_paths.saveOutput_directory)

end

rev = 1;
folder_paths.saveOutput_filename = [folder_paths.saveOutput_directory,...
    'EMIT_dropRetrieval_', folder_paths.coincident_dataFolder(1:end-3),...
    '_EMIT-pixel-lin-idx_', num2str(overlap_pixels.emit.linear_idx),...
    '_ran-on-', char(datetime("today")), '_rev', num2str(rev), '.mat'];


while isfile(folder_paths.saveOutput_filename)
    rev = rev + 1;
    if rev < 10
        folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-5), num2str(rev), '.mat'];
    elseif rev > 10
        folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-6), num2str(rev), '.mat'];
    end
end


%% Run the single-pixel retrieval

[GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = retrieve_dropProf_acpw_EMIT_Aqua_singlePix_ver2(emit,...
    modis, airs, overlap_pixels, folder_paths, print_libRadtran_err, print_status_updates, pixel_num);


end
