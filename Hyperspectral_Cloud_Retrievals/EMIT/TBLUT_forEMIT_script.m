% Compute Cloud Optical Properties from EMIT data using the two-wavelength
% look-up-table method

% By Andrew John Buggee

clear variables

%% Load EMIT Data

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitFolder = '17_Jan_2024_coast/';

% -------------------------------------



% Load modis data and create input structure

% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the EMIT data folder path

    emitPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

    % Define the LibRadTran path
    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];


    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the EMIT data folder path

    emitPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

    % Define the LibRadTran path
    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4'];

    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];




end



[emit,L1B_fileName] = retrieveEMIT_data([emitPath, emitFolder]);


%% Define the pixels to use for the retrieval


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% row = 1112;
% col = 974;

% 17_Jan_2024_coast - large optical depth
pixels2use.row = [912, 913];
pixels2use.col = 929;

% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));



%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_emit_inputs_TBLUT(folder2save, L1B_fileName, emit);

% *** Check Inputs ***

%% Define the spectral response function of EMIT for the desired Bands

% create the spectral response functions
emit.spec_response = create_EMIT_specResponse(emit, inputs);

%% Define the solar source file name and read in the solar source data

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);

%% Convert radiance measurements to TOA reflectance for the desired pixels

emit = convert_EMIT_radiance_2_reflectance(emit, inputs);


%% Check to see if the pixel in question is looking at a cloud made of liquid water


inputs = check_EMIT_therodynamic_phase(emit, pixels2use);














%%