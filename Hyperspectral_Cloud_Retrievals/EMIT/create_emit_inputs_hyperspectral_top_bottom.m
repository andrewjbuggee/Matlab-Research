%% ----- CREATE INPUTS NEEDED TO COMPUTE DROPLET RADIUS PROFILE USING HYPERSPECTRAL MEASUREMENTS FROM EMIT -----


% INPUTS:
%   (1) emitDataFolder - 

%   (2) folder2save - 

%   (3) L1B_filename - 

%   (4) emit - 


% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_emit_inputs_hyperspectral_top_bottom(emitDataFolder, folderpaths, L1B_fileName, emit)


% --- SAVE THE EMIT FILE NAME ----
inputs.emitDataFolder = emitDataFolder;



% ----- Save the L1B file name -----
inputs.L1B_filename = L1B_fileName{1};
 


% -------------- Define which EMIT bands to run ------------------
% ****************************************************************
% *-*-*-*- Only keep wavelengths that avoid water vapor -*-*-*-*-*

% The following indexes are for wavelengths that avoid water vapor
% absopriton according to figure 5 from King and Vaughan, which shows the
% information content for r_top, r_bot, tau_c, and water vapor across
% wavelengths from 500 to 2500 nm

% inputs.bands2run = [17, 24, 31, 40, 52, 65, 86, 92, 93, 115, 117, 118, 119, 121,...
%     158, 159, 164, 165, 166, 167, 168, 174, 175, 221, 222, 226, 232, 234, 238,...
%     248, 252, 259]';              % these are the bands that we will run uvspec with

% --- New indexs - tried to improve avoidance of water vapor ---
% inputs.bands2run = [17, 24, 32, 40, 53, 67, 86, 89, 90, 117, 118, 119, 120, 121,...
% 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 227, 236,...
% 249, 250, 251, 252, 253, 254]';

% --- New New indexs - avoid water vapor and other absorbing gasses! ---
% inputs.bands2run = [17, 20, 24, 32, 38, 66, 67, 68, 86, 90, 92, 115, 116, 117,...
% 156, 157, 158, 175, 176, 230, 231, 233, 234, 235, 236, 237,...
% 249, 250, 251, 252, 253, 254]';

% --- New New New indexs - using HiTran - avoid water vapor and other absorbing gasses! ---
% inputs.bands2run = [8, 12, 17, 22, 25, 32, 39, 65, 66, 67, 68, 86, 87, 88, 89, 90,...
%     94, 114, 115, 116, 117, 123, 124, 154, 155, 156, 157, 158, 172, 175, 176, 230,...
%     231, 233, 234, 235, 236, 237, 249, 250, 251, 252, 253, 254]';

% % --- New New New New indexs - using HiTran - avoid water vapor and other absorbing gasses! With Pilewskie input ---
% inputs.bands2run = [12, 17, 22, 25, 32, 39, 65, 66, 67, 68, 86, 87, 88, 89, 90,...
%     94, 115, 116, 117, 156, 157, 158, 172, 175, 176,...
%     231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';


% --- New New New New New indexs - using HiTran - avoid water vapor and other absorbing gasses! With Pilewskie input ---
% libRadtran estimates of reflectance below 500 nm consistently
% overestimate the measured values from EMIT. Let's ignore wavelengths
% below 500
inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 86, 87, 88, 89, 90,...
    94, 115, 116, 117, 156, 157, 158, 172, 175, 176,...
    231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';


inputs.bands2plot = [38, 235];    % these are the EMIT bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte




% -----------------------------------------------------------
% --------------- Iterations and Convergence ----------------
% -----------------------------------------------------------

% define the number of iterations for the gauss-newton solver
inputs.GN_iterations = 5;


% define a percent threshold of the difference between successive
% iterations. If the percent difference is below the percent threshold,
% than the iterative process is stopped.
inputs.percent_change_limit = 0.03;

% Define the convergence limit. Convergence is defined using the residual,
% which is the difference between the true and estimated measurements.
% We take the RMS of the residual using all spectral channels. This is how
% we define the convergence limit. If the residual is the difference
% between the true measurement and the estimated measurement, and the true
% measurement has an uncertainty of 10%, then our estimate measurement
% should be within this uncertainty. Using MODIS, we can compute the
% RMS uncertainty vector and set this as the convergence limit.

inputs.convergence_limit = 0.005;  % generic convergence limit




% define the type of model prior pdf
inputs.model.prior = 'gaussian';




% define the number of model parameters to solve for
inputs.num_model_parameters = 3;


% -------------------------------------------
% --- Stuff for the Model Parameter Prior ---
% -------------------------------------------

% Using the King and Vaughn (2012) method, we retireve 3 or 4 parameters
%   (1) r_top = effective droplet size at the cloud top
%   (2) r_bottom = effective droplet size at the cloud bottom
%   (3) tau_c = cloud optical depth
%   (4) cwv = above-cloud column water vapor
% a good starting place is to assume the droplet size at cloud top and
% bottom are the same value

if inputs.num_model_parameters==3

    inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
        'Cloud Optical Depth'};

elseif inputs.num_model_parameters==4

    inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
        'Cloud Optical Depth', 'Above Cloud Columnn Water Vapor'};

end





% ---------------------------------------
% --- Stuff for the Measurement Prior ---
% ---------------------------------------


inputs.measurement.prior = 'gaussian';
% covaraince_type can be:
%   (1) 'independent - thus all off diagonal elements are 0
%   (2) 'computed' - uses measured data to compute covaraince
inputs.measurement.covariance_type = 'independent';




% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

% we have to assume a vertical profile of droplet size with cloud optical
% depth exsists. And we only retrieve the droplet size at the top and
% bottom. This is the method of King and Vaughn (2012)

% the options for the vertical droplet profile are:
%       (a) 'adiabatic' - this assumption forces the liquid water content to
%       be proportionl to z, the altitude.
%       (b) 'subadiabatic_aloft' - this assumption assumes there is
%       increasing entrainment and drying towards the cloud top.
%       (c) 'linear_with_z' - this constraint forces the effective droplet profile
%       to behave linearly with z (re(z)~z). Physically we are forcing subadiabtatic
%       behavior at mid-levels.
%       (d) 'linear_with_tau' - this constraint forces the effective
%       droplet radius to have linearly with optical depth (re(z)~tau).
%       Physically, this too forces subadiabatic behavior at mid-levels.
% x is determined by the choice of droplet profile within the function
% create_droplet_profile.m

inputs.model.profile.type = 'adiabatic';
inputs.model.profile.r_top = 10; % microns - value for our model
inputs.model.profile.r_bottom = 5; % microns - value for our model





% -----------------------------------------------------
% ---------------------- FLAGS! -----------------------
% -----------------------------------------------------

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

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% -------- Do you want to use the effective water vapor depth? -----------

% This comes from the retrieval of thermodynamic phase
inputs.RT.use_phaseRetrieval_columnWaterVapor = false;
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
%inputs.RT.drop_distribution_var = 26.75;
inputs.RT.drop_distribution_var = 7;
% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';
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


% ------------------------------------------------------------------------
% -------- Do you want to modify the column water vapor amount? ----------
inputs.RT.modify_waterVapor = true;

inputs.RT.waterVapor_column = 40;       % mm (kg/m^2) - of water condensed in a column
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 416;       % ppm - concentration of CO2
% ------------------------------------------------------------------------



%----------------------------------------------------------
% ------- Define the Cloud Top Height and Depth -----------
%----------------------------------------------------------


% Define a fixed cloud top height
inputs.RT.cloudTop_height = 3;           % km

% Define a custom cloud depth
inputs.RT.cloudDepth = 1;            % km


% Define number of layers to use in libRadTran when defining
% vertically inhomogenous clouds
inputs.RT.cloud_layers = 10;


% ----- Do you want a long error message? -----
% if so, set error message to 'verbose'. Otherwise, set error message to
% 'quiet'
inputs.RT.err_msg = 'verbose';




% --------------------------------------------
% Create a new folder to save all calculations
% --------------------------------------------

% Define the folder that stores the inputs and calculated reflectanes
% using todays date
data_date = datetime([L1B_fileName{1}(18:21), '-', L1B_fileName{1}(22:23), '-', L1B_fileName{1}(24:25)],...
    'InputFormat','yyyy-MM-dd');

% Store the file name for the libRadTran INP and OUT files
inputs.folder2save.libRadTran_INP_OUT = [folderpaths.libRadTran_INP_OUT, 'EMIT_',char(data_date),...
    '_time_', L1B_fileName{1}(27:30), '/'];


% This is the folder where the reflectance calculations will be stored
inputs.folder2save.reflectance_calcs = [folderpaths.reflectance_calcs, emitDataFolder]; 

% This is the name of the .mat file with the reflectance calcs
inputs.reflectance_calculations_fileName = ['hyperspectral_reflectance_calculations_', char(datetime("today")),'.mat'];








% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
