%% Create figures for IRS 2024 Presentation

% By Andrew John Buggee

%% Reflectance spectra example over cloudy pixel

clear variables

% --- Load EMIT Data ---

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


% --- Define the pixels to use for the retrieval ---

% 17_Jan_2024_coast - optical depth of 6.6
pixels2use.row = 932;
pixels2use.col = 960;

% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));



% --- Remove data that is not needed ---

emit = remove_unwanted_emit_data(emit, pixels2use);

inputs = create_emit_inputs_TBLUT(emitDataFolder, folder2save, emit);

% --- Define the spectral response function of EMIT for the desired Bands
% ---

emit.spec_response = create_EMIT_specResponse(emit, inputs);

% --- Define the solar source file name and read in the solar source data
% ---

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);


% --- Convert radiance measurements to TOA reflectance for the desired
% pixels ---

emit = convert_EMIT_radiance_2_reflectance(emit, inputs);

% --- Create plot ---

figure;
plot(emit.radiance.wavelength, emit.reflectance.value, '.-', 'MarkerSize', 20,...
    'LineWidth',1, 'Color', mySavedColors(3, 'fixed'))
xlabel('Wavelength ($nm$)', Interpreter='latex', FontSize=30)
ylabel('Reflectance ($1/sr$)', Interpreter='latex', FontSize=30)
grid on; grid minor

% set figure size
set(gcf, 'Position', [0 0 1200 650])
