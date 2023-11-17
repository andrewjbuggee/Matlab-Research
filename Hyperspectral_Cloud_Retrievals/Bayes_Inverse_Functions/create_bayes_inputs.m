function bayes_inputs = create_bayes_inputs(modisInputs)

% ----- Define the number of pixels to estimate a profile for -----
%bayes_inputs.numPixels2Calculate = 4;
% For now, let's only compute this calculation for the MODIS pixel in
% question
% We will only use what is in the truth table!
bayes_inputs.numPixels2Calculate = modisInputs.pixels.num_2calculate;
% -------------------------------------------------------------------


% define the number of iterations for the gauss-newton solver
bayes_inputs.GN_iterations = 5;



% Define the convergence limit. Convergence is defined using the residual,
% which is the true measurement subtracted from the estiamted measurement.
% We take the RMS of the residual using all spectral channels. This is how
% we define the convergence limit. If the residual is the difference
% between the true measurement and the estimated measurement, and the true
% measurement has an uncertainty of 10%, then our estimate measurement
% should be within this uncertainty. For MODIS, the measuremen uncertainty
% for the reflectance is between 3 and 7%. So lets meet in the middle and
% say 5 %
bayes_inputs.convergence_limit = 0;

% define a percent threshold of the difference between successive
% iterations. If the percent difference is below the percent threshold,
% than the iterative process is stopped.
bayes_inputs.percent_change_limit = 0.03;

% define the type of model prior pdf
bayes_inputs.model.prior = 'gaussian';


% define the number of model parameters to solve for
bayes_inputs.num_model_parameters = 3;

% Define the spectral channels to use in the gauss-newton solver
% The data from 11-11-2008 at 18:50 measured erroneous values in the 1.6
% micron channel. If using this data, lets ignore this measurement

if strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1850/')==true
    
    bayes_inputs.bands2use = [1:5,7];
    
else 

    bayes_inputs.bands2use = 1:7;  % number of spectral bands to use

end



% -------------------------------------------
% --- Stuff for the Model Parameter Prior ---
% -------------------------------------------

% Using the King and Vaughn (2012) method, we retireve 3 parameters
%   (1) r_top = effective droplet size at the cloud top
%   (2) r_bottom = effective droplet size at the cloud bottom
%   (3) tau_c = cloud optical depth
% a good starting place is to assume the droplet size at cloud top and
% bottom are the same value



    

bayes_inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
    'Cloud Optical Depth'};


% ---------------------------------------
% --- Stuff for the Measurement Prior ---
% ---------------------------------------


bayes_inputs.measurement.prior = 'gaussian';
% covaraince_type can be:
%   (1) 'independent - thus all off diagonal elements are 0
%   (2) 'computed' - uses measured data to compute covaraince
bayes_inputs.measurement.covariance_type = 'independent';

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

bayes_inputs.model.profile.type = 'adiabatic';
bayes_inputs.model.profile.r_top = 10; % microns - value for our model
bayes_inputs.model.profile.r_bottom = 5; % microns - value for our model



% -----------------------------------------------
% ------------- Folder Locations  ---------------
% -----------------------------------------------
bayes_inputs.save_calcs_fileName = ['uvspec_CALCS_4Bayes_',date,'.mat'];




% -----------------------------------------------
% --------------- Define flags  -----------------
% -----------------------------------------------










% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------


% Define the number of streams to use in your radiative transfer model
bayes_inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
bayes_inputs.RT.use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
bayes_inputs.RT.band_parameterization = 'reptran coarse';
%band_parameterization = 'reptran_channel modis_terra_b07';
% ------------------------------------------------------------------------


% ---------------------------------------------------------
% ------ Define the Solar Flux file and it's resolution ---
% ---------------------------------------------------------
% resolution should match the value listed in the file name
bayes_inputs.RT.sourceFile_resolution = 1;                  % nm
% Define the source file
bayes_inputs.RT.source_file = '../data/solar_flux/kurudz_1.0nm.dat';

% define the atmospheric data file
bayes_inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
bayes_inputs.RT.surface_albedo = 0.05;

% day of the year
bayes_inputs.RT.day_of_year = str2double(modisInputs.L1B_filename(15:17));




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
bayes_inputs.RT.yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
bayes_inputs.RT.linear_cloudFraction = false;
% if false, define the cloud cover percentage
bayes_inputs.RT.cloud_cover = 1;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
bayes_inputs.RT.use_MODIS_cloudTopHeight = false;

% --- Do you want to use the VOCALS-REx cloud top height measurement? ----
bayes_inputs.RT.use_VOCALS_cloudTopHeight = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS above cloud water vapor? ---------
bayes_inputs.RT.use_MODIS_aboveCloudWaterVapor = false;
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% -------- Do you want to use the VOCALS measured cloud depth? -----------
bayes_inputs.RT.use_VOCALS_cloudDepth = true;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
bayes_inputs.RT.use_custom_mie_calcs = false;
% ------------------------------------------------------------------------
% This string is used to compute the LWC from optical depth and effective radius
% can be 'hu' or 'mie interpolate'
bayes_inputs.RT.wc_parameterization = 'mie interpolate';        % use the hu and stamnes parameterization for converting cloud properties to optical properties
% define the type of droplet distribution
bayes_inputs.RT.drop_distribution_str = 'gamma';
% define the distribution varaince
% 7 is the value libRadTran uses for liquid water clouds
bayes_inputs.RT.drop_distribution_var = 7;
% define whether this is a vertically homogenous cloud or not
bayes_inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';
% define how liquid water content will be computed
% can either be 'mie' or '2limit'
bayes_inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius


% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
bayes_inputs.RT.use_coxMunk = true;
bayes_inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
bayes_inputs.RT.yesAerosols = true;

bayes_inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
bayes_inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------


% ----- Do you want a long error message? -----
% if so, set error message to 'verbose'. Otherwise, set error message to
% 'quiet'
bayes_inputs.RT.err_msg = 'quiet';










end