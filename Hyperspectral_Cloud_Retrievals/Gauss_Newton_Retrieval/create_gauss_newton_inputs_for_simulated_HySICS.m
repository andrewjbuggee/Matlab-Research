function GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS()



% define the number of iterations for the gauss-newton solver
GN_inputs.GN_iterations = 5;


% define a percent threshold of the difference between successive
% iterations. If the percent difference is below the percent threshold,
% than the iterative process is stopped.
GN_inputs.percent_change_limit = 0.03;

% define the type of model prior pdf
GN_inputs.model.prior = 'gaussian';


% define the number of model parameters to solve for
GN_inputs.num_model_parameters = 3;




% -------------------------------------------
% --- Stuff for the Model Parameter Prior ---
% -------------------------------------------

% Using the King and Vaughn (2012) method, we retireve 3 parameters
%   (1) r_top = effective droplet size at the cloud top
%   (2) r_bottom = effective droplet size at the cloud bottom
%   (3) tau_c = cloud optical depth
% a good starting place is to assume the droplet size at cloud top and
% bottom are the same value



    

GN_inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
    'Cloud Optical Depth'};


% ---------------------------------------
% --- Stuff for the Measurement Prior ---
% ---------------------------------------


GN_inputs.measurement.prior = 'gaussian';
% covaraince_type can be:
%   (1) 'independent - thus all off diagonal elements are 0
%   (2) 'computed' - uses measured data to compute covaraince
GN_inputs.measurement.covariance_type = 'independent';

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

GN_inputs.model.profile.type = 'adiabatic';
GN_inputs.model.profile.r_top = 10; % microns - value for our model
GN_inputs.model.profile.r_bottom = 5; % microns - value for our model



% -----------------------------------------------
% ------------- Folder Locations  ---------------
% -----------------------------------------------
GN_inputs.save_calcs_fileName = ['uvspec_GaussNewton_calcs_',date,'.mat'];




% -----------------------------------------------
% --------------- Define flags  -----------------
% -----------------------------------------------










% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------

% Define the parameters of the INP file


% Determine which computer this is being run on
GN_inputs.which_computer = which_computer;

% Define the RTE Solver
GN_inputs.RT.rte_solver = 'disort';


% Define the number of streams to use in your radiative transfer model
GN_inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------- Define the source file and resolution -------------------

% source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
% source_file_resolution = 0.025;         % nm

% these data have 0.1nm sampling resolution, despite what the file name
% suggests
GN_inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
GN_inputs.RT.source_file_resolution = 0.1;         % nm

% these data have 1nm sampling resolution
% GN_inputs.RT.source_file = 'kurudz_1.0nm.dat';
% GN_inputs.RT.source_file_resolution = 1;         % nm

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------------------- DEFINE THE WAVELENGTHS! -------------------------
% ------------------------------------------------------------------------
% define the wavelength range. If monochromatic, enter the same number
% twice


% ----------------- Simulating HySICS spectral channels ------------------
% number of channels = 636 ranging from center wavelengths: [351, 2297]
% GN_inputs.bands2run = (1:1:636)';

% Paper 1 - Figures 7 and 8 - 35 spectral channels that avoid water vapor
% and other gaseous absorbers
GN_inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 217, 220,...
    222, 224, 227, 237, 288, 290, 293, 388, 390, 393,...
    426, 434, 436, 570, 574, 577, 579, 582, 613, 616,...
    618, 620, 623, 625]';


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Create the spectral response functions ------------------
% ------------------------------------------------------------------------

% ------------------------------------
% modeling the HySICS instrument...
% ------------------------------------
% Define the HySICS spectral response functions
spec_response = create_HySICS_specResponse(GN_inputs.bands2run, GN_inputs.RT.source_file, ...
    GN_inputs.which_computer);

% now define the wavelength range of each spectral channel
GN_inputs.RT.wavelengths2run = zeros(length(GN_inputs.bands2run), 2);

for ww = 1:length(GN_inputs.bands2run)
    % The wavelength vector for libRadTran is simply the lower and upper
    % bounds
    GN_inputs.RT.wavelengths2run(ww,:) = [spec_response.wavelength(ww, 1),...
        spec_response.wavelength(ww, end)];

end


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
GN_inputs.RT.use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
GN_inputs.RT.band_parameterization = 'reptran coarse';
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% define the atmospheric data file
GN_inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
GN_inputs.RT.albedo = 0.04;

% day of the year
GN_inputs.RT.day_of_year = 100;
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
GN_inputs.RT.yesCloud = true;

% define the cloud geometric depth
GN_inputs.RT.cloud_depth = 500;                % meters

% define the geometric location of the cloud top and cloud bottom
GN_inputs.RT.z_topBottom = [1.5, 1];          % km above surface


% Water Cloud depth
GN_inputs.RT.H = GN_inputs.RT.z_topBottom(1) - GN_inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --------------------- Various Cloud modeling GN_inputs --------------------
% ------------------------------------------------------------------------
% Do you want use your custom mie calculation file?
% If false, the default file is the precomputed mie table from libRadtran
GN_inputs.RT.use_custom_mie_calcs = false;


% define how liquid water content will be computed in the write_wc_file
% function.
GN_inputs.RT.parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
GN_inputs.RT.lambda_forTau = 500;            % nm
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------------- Cloud optical properties --------------------------
% ------------------------------------------------------------------------

% define whether this is a vertically homogenous cloud or not
GN_inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';


if strcmp(GN_inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % --------- HOMOGENOUS CLOUD -------------



    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    GN_inputs.RT.dist_var = 7;              % distribution variance

    % define the type of droplet distribution
    GN_inputs.RT.distribution_str = 'gamma';

    % GN_inputs.RT.re = 10.79;      % microns
    % GN_inputs.RT.tau_c = 7;

    GN_inputs.RT.re = 1:2:51;      % microns
    GN_inputs.RT.tau_c = [1:15, 20:5:100];


elseif strcmp(GN_inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------- NON-HOMOGENOUS CLOUD -------------

    % For a cloud with a nonconstant droplet profile
    % constraint - the physical constraint (string) - there are four
    %       different string options for a physical constraint:
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

    GN_inputs.RT.profile_type = 'adiabatic'; % type of water droplet profile

    GN_inputs.RT.n_layers = 10;                          % number of layers to model within cloud

    GN_inputs.RT.z = linspace(GN_inputs.RT.z_topBottom(1), GN_inputs.RT.z_topBottom(2), GN_inputs.RT.n_layers);        % km - altitude above ground vector

    GN_inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used

    GN_inputs.RT.dist_var = linspace(10,10, GN_inputs.RT.n_layers);              % distribution variance

    GN_inputs.RT.r_top = 12.565;     % microns
    GN_inputs.RT.r_bot = 4.135;        % microns
    GN_inputs.RT.tau_c = 6.424;

    % GN_inputs.RT.r_top = 9:10;     % microns
    % GN_inputs.RT.r_bot = 4:5;        % microns
    % GN_inputs.RT.tau_c = [5,10];

    %     GN_inputs.RT.r_top = 3:20;       % microns
    %     GN_inputs.RT.r_bot = 2:14;        % microns
    %     GN_inputs.RT.tau_c = [5.5, 6, 6.5, 7, 7.5];


    % define the type of droplet distribution
    GN_inputs.RT.distribution_str = 'gamma';

    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    GN_inputs.RT.dist_var = 7;              % distribution variance


end

% ------------------------------------------------------------------------



% Define the parameterization scheme used to comptue the optical quantities
% within the INP file, i.e. what function call that tells libRadtran how to
% compute scattering and optical quantities
if GN_inputs.RT.use_custom_mie_calcs==false

    GN_inputs.RT.wc_parameterization = 'mie interpolate';
    %GN_inputs.RT.wc_parameterization = 'hu';

else
    GN_inputs.RT.wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
end

% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% Define the altitude of the sensor
GN_inputs.RT.sensor_altitude = 'toa';          % top-of-atmosphere

% define the solar zenith angle
GN_inputs.RT.sza = 0;           % degree

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
GN_inputs.RT.phi0 = 0;         % degree

% define the viewing zenith angle
GN_inputs.RT.vza = 0; % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

GN_inputs.RT.vaz = 0;     % degree
% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
GN_inputs.RT.use_coxMunk = true;
GN_inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
GN_inputs.RT.yesAerosols = true;

GN_inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
GN_inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------




% --------------------------------------------------------
% --------- What is column water vapor amount? -----------

% Use a custom H2O profile
GN_inputs.RT.H2O_profile = 'afglus_H2O_none_inside_cloud.dat';


% Using measurements from the AMSR2 instrument, a passive microwave
% radiometer for 17 Jan 2024
GN_inputs.RT.modify_waterVapor = false;

GN_inputs.RT.waterVapor_column = 20;              % mm - milimeters of water condensed in a column
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
GN_inputs.RT.modify_CO2 = true;

GN_inputs.RT.CO2_mixing_ratio = 416;       % ppm
% ------------------------------------------------------------------------


% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
GN_inputs.RT.compute_reflectivity_uvSpec = false;
% --------------------------------------------------------------



% --------------------------------------------------------------
% Do you want to print an error message?
GN_inputs.RT.errMsg = 'quiet';
% --------------------------------------------------------------









end
