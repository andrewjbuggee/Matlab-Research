% --- Add Folders to Path specific to this computer -----

computer_name = whatComputer;

if strcmp(computer_name,'anbu8374')==true
    % ----- LASP MAC ------
    % LibRadTran Data Folders
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/albedo/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/altitude/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/correlated_k/');

    % LibRadTran Bin folder
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/bin/');

elseif strcmp(computer_name,'andrewbuggee')==true
    % ----- MACBOOK ------


    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/solar_flux/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/albedo/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/altitude/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/atmmod/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/correlated_k/');

    % LibRadTran Bin folder
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/bin/');


elseif strcmp(computer_name,'curc')==true
    % ----- CU SUPERCOMPUTER ------

    % Add paths to all functions needed
    addpath('/projects/anbu8374/Hyperspectral_Cloud_Retrieval/')
    addpath('/projects/anbu8374/Hyperspectral_Cloud_Retrieval/Bayes_Inverse_Functions/')
    addpath('/projects/anbu8374/LibRadTran-wrapper/')
    addpath('/projects/anbu8374/MODIS-Cloud-Retrieval/')
    addpath('/projects/anbu8374/MODIS_data/')
    addpath('/projects/anbu8374/VOCALS_REx/')
    addpath('/projects/anbu8374/Radiative_Transfer_Physics/')

end

