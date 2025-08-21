%% Define all the folder paths necessary to read HySICS data, run libRadtran using HySICS data, and save output

% By Andrew John Buggee

%%

function [folder_paths, which_computer] = define_folderPaths_for_HySICS()



% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ***** Define the HySICS Folder with the simulated measurements *****
    folder_paths.HySICS_simulated_spectra = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/paper2_variableSweep/'];

    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Droplet_profile_retrievals/'];

    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/';


    % libRadtran data folder
    folder_paths.libRadtran_data = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';


    % libRadtran atmosphere folder location
    folder_paths.atm_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/';


    % water cloud file location
    folder_paths.libRadtran_water_cloud_files = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';

    % mie folder location
    folder_paths.libRadtran_mie_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';




elseif strcmp(which_computer,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % ***** Define the HySICS Folder with the simulated measurements *****
    % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

    folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/'];


    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];


    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/HySICS/INP_OUT_2/'];


    % libRadtran data folder
    folder_paths.libRadtran_data = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/'];


    % water cloud file location
    folder_paths.libRadtran_water_cloud_files = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/'];

    % libRadtran atmosphere folder location
    folder_paths.atm_folder_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
                       'LibRadTran/libRadtran-2.0.4/data/atmmod/'];


    % mie folder location
    folder_paths.libRadtran_mie_folder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/',...
        'libRadtran-2.0.4/Mie_Calculations_2/'];




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


    % Define the HySICS simulated spectrum folder

    folder_paths.HySICS_simulated_spectra = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
        'Simulated_spectra/paper2_variableSweep/'];


    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/';

    
    % libRadtran data folder
    folder_paths.libRadtran_data = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


    % water cloud file location
    folder_paths.libRadtran_water_cloud_files = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';


    % libRadtran atmosphere folder location
    folder_paths.atm_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/atmmod/';


    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/scratch/alpine/anbu8374/HySICS/INP_OUT_2/';


    % mie folder location
    folder_paths.libRadtran_mie_folder = '/scratch/alpine/anbu8374/Mie_Calculations_2/';


   


end




% If the libRadtran INP folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_inp, 'dir')

    mkdir(folder_paths.libRadtran_inp)

end


% If the mie INP folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_mie_folder, 'dir')

    mkdir(folder_paths.libRadtran_mie_folder)

end



end

