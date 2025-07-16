%%


% By Andrew John Buggee

%%

function [folder_paths] = define_EMIT_dataPath_and_saveFolders()

% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the EMIT data folder path

    folder_paths.emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/EMIT';


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    folder_paths.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';

    % water cloud file location
    folder_paths.water_cloud_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';


    % Define the folder path where the mat files of reflectances will be
    % saved
    folder_paths.reflectance_calcs = folder_paths.emitDataPath;


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the EMIT data folder path

    folder_paths.emitDataPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP and .OUT files will be saved
    folder_paths.libRadtran_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/'];


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    folder_paths.libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];

    % water cloud file location
    folder_paths.water_cloud_folder_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/'];


    % Define the folder path where the mat files of reflectances will be
    % saved
    folder_paths.reflectance_calcs = folder_paths.emitDataPath;


elseif strcmp(whatComputer,'curc')==true


    % ------ Folders on the CU Supercomputer /projects folder --------


    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/scratch/alpine/anbu8374/EMIT/INP_OUT/';


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    folder_paths.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


    % Define the EMIT data folder path
    folder_paths.emitDataPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

    % water cloud file location
    folder_paths.water_cloud_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';

    % Define the folder path where the mat files of reflectances will be
    % saved
    folder_paths.reflectance_calcs = folder_paths.emitDataPath;


end


% If the folder path doesn't exit, create a new directory
if ~exist(folder_paths.reflectance_calcs, 'dir')

    mkdir(folder_paths.reflectance_calcs);

end



% If the folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_inp, 'dir')

    mkdir(folder_paths.libRadtran_inp);

end

