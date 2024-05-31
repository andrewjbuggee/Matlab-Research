% Compute Cloud Optical Properties from EMIT data using the two-wavelength
% look-up-table method

% By Andrew John Buggee

clear variables

%% Load EMIT Data

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitDataFolder = '17_Jan_2024_coast/';

% -------------------------------------


% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the EMIT data folder path

    emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP files will be saved
    folder2save.libRadTran_INP_OUT = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/'];

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

    % Define the folder path where the mat files of reflectances will be
    % saved
    folder2save.reflectance_calcs = emitDataPath;




end



[emit,L1B_fileName] = retrieveEMIT_data([emitDataPath, emitDataFolder]);


%% Define the pixels to use for the retrieval


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% row = 1112;
% col = 974;

% 17_Jan_2024_coast - optical depth of 4 and 4.5
% pixels2use.row = [912, 913];
% pixels2use.col = [929, 929];

% 17_Jan_2024_coast - optical depth of 6.6
% pixels2use.row = [932];
% pixels2use.col = [960];

% 17_Jan_2024_coast - optical depth of 3.2 and 3.8
pixels2use.row = [932, 932];
pixels2use.col = [970, 969];

% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));



%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_emit_inputs_TBLUT(emitDataFolder, folder2save, emit);

% *** Check Inputs ***

%% Define the spectral response function of EMIT for the desired Bands

% create the spectral response functions
emit.spec_response = create_EMIT_specResponse(emit, inputs);

%% Define the solar source file name and read in the solar source data

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);


%% Convert radiance measurements to TOA reflectance for the desired pixels

emit = convert_EMIT_radiance_2_reflectance(emit, inputs);


%% Check to see if the pixel in question is looking at a cloud made of liquid water
tic
inputs = check_EMIT_therodynamic_phase(emit, inputs);
toc

%% ----- Create .INP files for EMTI TBLUT -----


if inputs.flags.writeINPfiles == true

    [names.inp, inputs] = write_INP_file_4EMIT_homogenous(inputs, pixels2use, emit);

    % now lets write the output names

    names.out = writeOutputNames(names.inp);
else

    % if the files already exist, just grab the names!
    [names.inp, inputs] = getMODIS_INPnames_withClouds(emit.solar, inputs, pixels2use);
    names.out = writeOutputNames(names.inp);
end


%% ----- Run uvspec and calculate Reflectance Function using LibRadTran -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different values of effective radius and optical depth

if inputs.flags.runUVSPEC == true

    % 1st output - R is the reflectance integrated over a bandwidth
    % 2nd output - Rl is the reflectance at each spectral bin
    tic
    [R,~, inputs] = runReflectanceFunction_4EMIT(inputs, names, emit.spec_response.value);
    toc

    % Save the pixels2use structure
    save([inputs.folder2save.reflectance_calcs, inputs.reflectance_calculations_fileName],...
        "pixels2use", "-append"); % save inputSettings to the same folder as the input and output file

elseif inputs.flags.runUVSPEC == false

    load([inputs.savedCalculations_folderName,inputs.saveCalculations_fileName] ,'inputs','R');

end


%% ----- Compare Reflectance Function of MODIS with Theoretical Calculations (Grid Search) -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

minVals = leastSquaresGridSearch_EMIT(emit.reflectance, R, inputs);


%% ------------ Make plots! ------------

plot2ReflectanceFuncBands_EMIT(emit, R, inputs, pixels2use, 'king')

