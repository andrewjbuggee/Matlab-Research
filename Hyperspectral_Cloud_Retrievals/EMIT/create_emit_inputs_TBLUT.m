%% ----- CREATE INPUTS NEEDED TO COMPUTE TBLUT METHOD ON EMIT DATA -----


% INPUTS:
%   (1) folderName - 

%   (2) L1B_fileName - 


% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_emit_inputs_TBLUT(folderName, L1B_fileName, emit)


% --- SAVE THE EMIT FILE NAME ----
inputs.emitDataFolder = folderName;



% ----- Save the L1B file name -----
inputs.L1B_filename = L1B_fileName{1};
 


% Define which EMIT bands to run
inputs.bands2run = [38, 235]; % these are the bands that we will run uvspec with
inputs.bands2search = [38, 235]; % these are the EMIT bands that are used in the retrieval problem
inputs.bands2plot = [38, 235]; % these are the EMIT bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte

% if interpGridScaleFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 150; % scale factor the will be used to increase the grid size for interpolation.

% --------------------------------------------
% Create a new folder to save all calculations
% --------------------------------------------
inputs.savedCalculations_folderName = [folderName, 'Retrieval_outputs_', char(datetime("today")),'/']; % this is the folder that all the saved calculations will go



inputs.saveCalculations_fileName = ['uvspec_calculations_', char(datetime("today")),'.mat'];



% Define the folder to save all the INP files in using the month, day and
% year
data_date = datetime([L1B_fileName{1}(18:21), '-', L1B_fileName{1}(22:23), '-', L1B_fileName{1}(24:25)],...
    'InputFormat','yyyy-MM-dd');

inputs.INP_folderName = ['EMIT_',char(data_date),'_time_',L1B_fileName{1}(27:30),'/']; % this is the folder name that the INP files will be written to 








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
inputs.RT.source_file_resolution = 0.1;           % nm


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
inputs.RT.drop_distribution_var = 10;
% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';
% define how liquid water content will be computed
% can either be 'mie' or '2limit'
inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius


% --------------------------------------------------------------
% --------------------------------------------------------------



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
inputs.RT.err_msg = 'quiet';









% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end