%% Sensitivity to detecting a change in the radius at cloud bottom

% By Andrew John Buggee

%% Compute when the change in reflectivity is greater than the measurement uncertainty for different cloud optical depths

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;

%% Define the EMIT data file to use by defining the data folder

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitDataFolder = '17_Jan_2024_coast/';

% -------------------------------------


%% Load EMIT data and define folders 

[emitDataPath, folder2save] = define_EMIT_dataPath_and_saveFolders();

[emit,L1B_fileName] = retrieveEMIT_data([emitDataPath, emitDataFolder]);


%% Define the pixels to use for the retrieval


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% row = 1112;
% col = 974;

% 17_Jan_2024_coast - small optical depth
% pixels2use.row = [912, 913];
% pixels2use.col = [929, 929];

% 17_Jan_2024_coast - optical depth of 6.6
% pixels2use.row = [932];
% pixels2use.col = [960];

% 17_Jan_2024_coast - optical depth of 3.2 and 3.8
% pixels2use.row = [932, 932];
% pixels2use.col = [970, 969];

% 17_Jan_2024_coast - optical depth of 3.2
pixels2use.row = [932];
pixels2use.col = [970];


% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));



%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_emit_inputs_hyperspectral_top_bottom(emitDataFolder, folder2save, L1B_fileName, emit);
%inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

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


%% Compute the radiance measurement uncertainty 

emit.radiance.uncertainty = compute_EMIT_radiance_uncertainty(emit);


%% Compute the reflectance uncertainty

emit.reflectance.uncertainty = compute_EMIT_reflectance_uncertainty(emit, inputs);


%% Compute the change in the reflectance

% define the current state vector [r_top, r_bottom, tau_c]
state_vector = [12, 6, 3];

% Compute the reflectances for the above state vector
measurement_estimate = compute_forward_model_4EMIT_top_bottom(emit, state_vector, inputs, pixels2use, 1)';

measurement_change = compute_reflectanceChange_due_to_rBottom_change(emit, state_vector, measurement_estimate, inputs,...
    pixels2use, pp, false);

