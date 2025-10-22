%% Retreive vertical profiles just as the batch scripts are written for slurm

% This script retrieves 4 variables: r_top, r_bot, tau_c, and cwvs



% By Andrew John Buggee



%% Load paths

% addpath(genpath('/projects/anbu8374/Matlab-Research')); 
% addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); 
% addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); 

clear variables; 
addLibRadTran_paths; 

folder_paths = define_folderPaths_for_HySICS(1);

folder_paths.HySICS_simulated_spectra = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
    'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/testGRC_results/']; 

folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
    'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/testGRC_results/'];

print_status_updates = false; 
print_libRadtran_err = false; 

%% Determine the names of the files

% Grab filenames in drive
filenames = dir(folder_paths.HySICS_simulated_spectra);
idx_2delete = [];
for nn = 1:length(filenames)

    if strcmp(filenames(nn).name(1), 's')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames(idx_2delete) = [];

% now create a cell array
file_list = cell(length(filenames), 1);

for nn = 1:length(filenames)
    file_list{nn} = filenames(nn).name;
end

%%


[tblut_retrieval, acpw_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropletProfile_HySICS_ver3_lowUncertainty(file_list, folder_paths, print_status_updates, print_libRadtran_err);