%%


% By Andrew John Buggee

%%

function [emitDataPath, folder2save] = define_EMIT_dataPath_and_saveFolders()

% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the EMIT data folder path

    emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

    % Define the folder path where all .INP files will be saved
    folder2save.libRadTran_INP_OUT = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/';


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';
    

    % Define the folder path where the mat files of reflectances will be
    % saved
    folder2save.reflectance_calcs = emitDataPath;


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the EMIT data folder path

    emitDataPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP and .OUT files will be saved
    folder2save.libRadTran_INP_OUT = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/'];


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];


    % Define the folder path where the mat files of reflectances will be
    % saved
    folder2save.reflectance_calcs = emitDataPath;


elseif strcmp(whatComputer,'curc')==true


    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    folder2save.reflectance_calcs = '/scratch/alpine/anbu8374/hyperspectral_retrieval/';



    % Define the folder path where all .INP files will be saved
    folder2save.libRadTran_INP_OUT = '/scratch/alpine/anbu8374/hyperspectral_retrieval/';
    % If the folder path doesn't exit, create a new directory
    if ~exist(folder2save.libRadTran_INP_OUT, 'dir')

        mkdir(folder2save.libRadTran_INP_OUT)

    end

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


    % Define the EMIT data folder path
    emitDataPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


end

