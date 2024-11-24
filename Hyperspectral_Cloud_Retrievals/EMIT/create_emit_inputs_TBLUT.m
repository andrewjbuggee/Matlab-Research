%% ----- CREATE INPUTS NEEDED TO COMPUTE TBLUT METHOD ON EMIT DATA -----


% INPUTS:
%   (1) folderName - 

%   (2) L1B_fileName - 


% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_emit_inputs_TBLUT(emitDataFolder, folder2save, emit)


%% Which computer are we using?



if strcmp(whatComputer,'anbu8374')

    EMIT_dataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

elseif strcmp(whatComputer,'andrewbuggee')

    EMIT_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'EMIT/EMIT_data/'];

else
    error('I dont recognize this computer user name')
end

%%

% --- SAVE THE EMIT FILE NAME ----
inputs.emitDataFolder = emitDataFolder;

% read the contents of the EMIT data folder
folder_contents = dir([EMIT_dataPath, emitDataFolder]);

% ----- Save the L1B file name -----
for nn = 1:length(folder_contents)

    if length(folder_contents(nn).name)>5 && strcmp(folder_contents(nn).name(1:12), 'EMIT_L1B_RAD')==true

        inputs.L1B_filename = folder_contents(nn).name;

    end

end 
 


% Define which EMIT bands to run
% band 38 has a center wavelength of 656 nm
% band 235 has a center wavelength of 2123 nm
inputs.bands2run = [38, 235]; % these are the bands that we will run uvspec with
inputs.bands2plot = [38, 235]; % these are the EMIT bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte

% if interpGridScaleFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 150; % scale factor the will be used to increase the grid size for interpolation.


% --------------------------------------------
% Create a new folder to save all calculations
% --------------------------------------------

% Define the folder that stores the inputs and calculated reflectanes
% using todays date
data_date = datetime([inputs.L1B_filename(18:21), '-', inputs.L1B_filename(22:23), '-', inputs.L1B_filename(24:25)],...
    'InputFormat','yyyy-MM-dd');

% Store the file name for the libRadTran INP and OUT files
inputs.folder2save.libRadTran_INP_OUT = [folder2save.libRadTran_INP_OUT, 'EMIT_',char(data_date),...
    '_time_', inputs.L1B_filename(27:30), '/'];


% This is the folder where the reflectance calculations will be stored
inputs.folder2save.reflectance_calcs = [folder2save.reflectance_calcs, emitDataFolder]; 

% This is the name of the .mat file with the reflectance calcs
inputs.reflectance_calculations_fileName = ['TBLUT_reflectance_calculations_', char(datetime("today")),'.mat'];




% ------------------
% ----- FLAGS! -----
% ------------------

% define flags that tell the codes to either run certain things, or don't
% run certain things

inputs.flags.findSuitablePixels = false; % if true, this will search the modis data set for pixels to use

% if true, the code will load an older set of pixels that has already been used before, and 
% likely has INP files. If false, it tells the code to find a new random subset of pixels
inputs.flags.loadPixelSet = true; 
inputs.flags.writeINPfiles = true; % if true, this will create inp files for each the length of vector pixel.row
inputs.flags.runUVSPEC = true; % if true, this will run all of the inp files create from the above flag through uvspec
inputs.flags.plotMLS_figures = false; % this will tell the leasSquaresGridSearch code to plot 





% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------


% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
inputs.RT.use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
inputs.RT.band_parameterization = 'reptran coarse';
%band_parameterization = 'reptran_channel modis_terra_b07';
% ------------------------------------------------------------------------



% ---------------------------------------------------------
% ------ Define the Solar Flux file and it's resolution ---
% ---------------------------------------------------------
% resolution should match the value listed in the file name

%inputs.RT.source_file = 'kurudz_0.1nm.dat';
%inputs.RT.source_file_resolution = 0.1;         % nm

%inputs.RT.source_file = 'kurudz_1.0nm.dat';
%inputs.RT.source_file_resolution = 1;         % nm

%inputs.RT.source_file = 'hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.dat';
%inputs.RT.source_file_resolution = 0.001;         % nm

%inputs.RT.source_file = 'hybrid_reference_spectrum_p025nm_resolution_c2022-11-30_with_unc.dat';
%inputs.RT.source_file_resolution = 0.005;         % nm

% inputs.RT.source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
% inputs.RT.source_file_resolution = 0.025;         % nm

inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
inputs.RT.source_file_resolution = 0.1;         % nm





% define the atmospheric data file
inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
inputs.RT.surface_albedo = 0.05;

% day of the year
inputs.RT.day_of_year = emit.day_of_year;




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
inputs.RT.linear_cloudFraction = false;
% if false, define the cloud cover percentage
inputs.RT.cloud_cover = 1;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
inputs.RT.use_MODIS_cloudTopHeight = false;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS above cloud water vapor? ---------
inputs.RT.use_MODIS_aboveCloudWaterVapor = false;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
inputs.RT.use_custom_mie_calcs = false;
% ------------------------------------------------------------------------
% This string is used to compute the LWC from optical depth and effective radius
% can be 'hu' or 'mie interpolate'
inputs.RT.wc_parameterization = 'mie interpolate';        % use the hu and stamnes parameterization for converting cloud properties to optical properties
% define the type of droplet distribution
inputs.RT.drop_distribution_str = 'gamma';
% define the distribution varaince
% 7 is the value libRadTran uses for liquid water clouds
inputs.RT.drop_distribution_var = 7;
% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';
% define how liquid water content will be computed
% can either be 'mie' or '2limit'
inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius


% --------------------------------------------------------------
% --------------------------------------------------------------



% ------------------------------------------------------------------------
% -------- Do you want to modify the column water vapor amount? ----------
inputs.RT.modify_waterVapor = false;

inputs.RT.waterVapor_column = 30;       % mm (kg/m^2) - of water condensed in a column
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 410;       % ppm - concentration of CO2
% ------------------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
inputs.RT.use_coxMunk = true;
inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
inputs.RT.yesAerosols = true;

inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------


% ----- Do you want a long error message? -----
% if so, set error message to 'verbose'. Otherwise, set error message to
% 'quiet'
inputs.RT.err_msg = 'verbose';









% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
