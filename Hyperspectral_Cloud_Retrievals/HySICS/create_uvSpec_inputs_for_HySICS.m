%% Create input structure for libRadtran INP files for the HySICS instrument

% By Andrew John Buggee

%%

function [inputs, spec_response] = create_uvSpec_inputs_for_HySICS(inputs)


% Define the parameters of the INP file

% *** compute weighting functions! ***
% This will create n wc_files where n is equal to the number of layers in
% the cloud. Starting with the entire cloud, each file will have one less
% layer
inputs.compute_weighting_functions = true;


% Are you simulating a measurement, or making forward model calculations
% for the retrieval?
% inputs.calc_type = 'simulated_spectra';
inputs.calc_type = 'monte_carlo';


% ----- Define the RTE Solver -----
% inputs.RT.rte_solver = 'disort';
inputs.RT.rte_solver = 'montecarlo';


% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;


% ---------------------------------------------------
% --------- Define Monte Carlo Parameters -----------
% ---------------------------------------------------
inputs.RT.mc.photons = 10000000;      % number of photons to use in the simulation
inputs.RT.mc.vroom = 'on';        % helps speed up calculations for particles with strong forward scattering
inputs.RT.mc.escape = 'on';       % calculates radiances via escape probabilities - speeds up computation

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------- Define the source file and resolution -------------------

% source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
% source_file_resolution = 0.025;         % nm

% these data have 0.1nm sampling resolution, despite what the file name
% suggests
inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
inputs.RT.source_file_resolution = 0.1;         % nm

% these data have 1nm sampling resolution
% inputs.RT.source_file = 'kurudz_1.0nm.dat';
% inputs.RT.source_file_resolution = 1;         % nm

% ------------------------------------------------------------------------


% ----- Define the day of the year to account for Earth-Sun distance -----
% inputs.RT.day_of_year = 239;



% ------------------------------------------------------------------------
% ---------------------- DEFINE THE WAVELENGTHS! -------------------------
% ------------------------------------------------------------------------
% define the wavelength range. If monochromatic, enter the same number
% twice


% ----------------- Simulating HySICS spectral channels ------------------
% number of channels = 636 ranging from center wavelengths: [351, 2297]
% inputs.bands2run = (1:1:636)';

% Paper 1 - Figures 7 and 8 - 35 spectral channels that avoid water vapor
% and other gaseous absorbers
% inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 217, 220,...
%     222, 224, 227, 237, 288, 290, 293, 388, 390, 393,...
%     426, 434, 436, 570, 574, 577, 579, 582, 613, 616,...
%     618, 620, 623, 625]';

% inputs.bands2run = [49, 426, 613]';

% test bands
% 500 nm 
% inputs.bands2run = 49;

% 1652 nm 
% inputs.bands2run = 426;

% 2122 nm 
% inputs.bands2run = 580;

% 2236 nm 
inputs.bands2run = 613;


% ------------------------------------------------------------------------
% Do you want to compute radiance/reflectance over a spectral region, or at
% a single wavelength?
% ------------------------------------------------------------------------
inputs.RT.monochromatic_calc = true;

% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
inputs.RT.compute_reflectivity_uvSpec = true;
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Create the spectral response functions ------------------
% ------------------------------------------------------------------------

% ------------------------------------
% modeling the HySICS instrument...
% ------------------------------------
% Define the HySICS spectral response functions
spec_response = create_HySICS_specResponse(inputs.bands2run, inputs.RT.source_file, ...
    inputs.which_computer);

% now define the wavelength range of each spectral channel
inputs.RT.wavelengths2run = zeros(length(inputs.bands2run), 2);


if inputs.RT.monochromatic_calc==true

    for ww = 1:length(inputs.bands2run)
        % The wavelength vector for libRadTran is simply the lower and upper
        % bounds
        inputs.RT.wavelengths2run(ww,:) = [round(mean(spec_response.wavelength(ww,:)), 1),...
            round(mean(spec_response.wavelength(ww, :)) ,1)];

    end


else

    for ww = 1:length(inputs.bands2run)
        % The wavelength vector for libRadTran is simply the lower and upper
        % bounds
        inputs.RT.wavelengths2run(ww,:) = [spec_response.wavelength(ww, 1),...
        spec_response.wavelength(ww, end)];

    end

end




% ------------------------------------------------------------------------
% ------------------------------------------------------------------------





% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
inputs.RT.band_parameterization = 'reptran coarse';
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% define the atmospheric data file
inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
% inputs.RT.surface_albedo = 0.04;        % Ocean water albedo
inputs.RT.surface_albedo = 0;             % Use a value of 0 when creating weighting functions

% day of the year
%inputs.RT.day_of_year = 17;
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% define the cloud geometric depth
inputs.RT.cloud_depth = 500;                % meters

% define the geometric location of the cloud top and cloud bottom
inputs.RT.z_topBottom = [1.5, 1];          % km above surface


% Water Cloud depth
inputs.RT.H = inputs.RT.z_topBottom(1) - inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --------------------- Various Cloud modeling inputs --------------------
% ------------------------------------------------------------------------
% Do you want use your custom mie calculation file?
% If false, the default file is the precomputed mie table from libRadtran
inputs.RT.use_custom_mie_calcs = false;


% define how liquid water content will be computed in the write_wc_file
% function.
inputs.RT.parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 500;            % nm
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------------- Cloud optical properties --------------------------
% ------------------------------------------------------------------------

% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';


if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % --------- HOMOGENOUS CLOUD -------------



    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    inputs.RT.distribution_var = 7;              % distribution variance

    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    inputs.RT.re = 5:2:9;      % microns
    inputs.RT.tau_c = 5:3:20;

    % inputs.RT.re = 1:2:51;      % microns
    % inputs.RT.tau_c = [1:15, 20:5:100];


elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

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

    % ** Values used in Platnick (2000) **
    inputs.RT.r_top = 12;     % microns
    inputs.RT.r_bot = 5;        % microns
    inputs.RT.tau_c = 8;

    % inputs.RT.r_top = 9.5167;     % microns
    % inputs.RT.r_bot = 4.0192;        % microns
    % inputs.RT.tau_c = 6.0312;

    % inputs.RT.r_top = 8.6;     % microns
    % inputs.RT.r_bot = 3.6;        % microns
    % inputs.RT.tau_c = 9.6;

    % inputs.RT.r_top = [9,10];     % microns
    % inputs.RT.r_bot = [4:6];        % microns
    % inputs.RT.tau_c = [10:5:25];

    % inputs.RT.r_top = 3:20;       % microns
    % inputs.RT.r_bot = 2:14;        % microns
    % inputs.RT.tau_c = [5.5, 6, 6.5, 7, 7.5];

    % inputs.RT.r_top = 3:20;       % microns
    % inputs.RT.r_bot = 2:14;        % microns
    % inputs.RT.tau_c = 7.5:0.5:15;


    inputs.RT.profile_type = 'adiabatic'; % type of water droplet profile

    % *** Use 250 if creating weighting functions using DISORT ***
    inputs.RT.n_layers = 250;                          % number of layers to model within cloud
    % inputs.RT.n_layers = 10;                          % number of layers to model within cloud

    % -------------------------------------------------------------------
    % define the independent variable used to define the effective radius
    % profile within the cloud
    % -------------------------------------------------------------------
    % if using altitude, z should be a vector starting at the cloud bottom
    % and increasing
    inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used
    inputs.RT.z_edges = linspace(inputs.RT.z_topBottom(2), inputs.RT.z_topBottom(1), inputs.RT.n_layers+1);   % km - the edges of each layer
    inputs.RT.z = linspace(inputs.RT.z_topBottom(2), inputs.RT.z_topBottom(1), inputs.RT.n_layers);        % km - altitude above ground vector

    % Try sampling random points between z top and bottom
    % inputs.RT.z_edges = sort(unique(abs(diff(inputs.RT.z_topBottom))*rand(1, inputs.RT.n_layers +1) + inputs.RT.z_topBottom(2)));   % km - the edges of each layer
    % inputs.RT.z = inputs.RT.z_edges(1:end-1) + diff(inputs.RT.z_edges);        % km - altitude above ground vector

    % If using optical depth, this vector should start with 0 (cloud top)
    % and end with the total cloud optical thickness.
    % inputs.RT.indVar = 'optical_depth';                    % string that tells the code which independent variable we used
    % inputs.RT.tau_edges = linspace(0, inputs.RT.tau_c, inputs.RT.n_layers+1)'; % define the boundaries of each tau layer
    % inputs.RT.zT_cloud_indVar = linspace(0, inputs.RT.tau_c, inputs.RT.n_layers)';        % optical depth



    inputs.RT.distribution_var = linspace(10,10, inputs.RT.n_layers);              % distribution variance



    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    inputs.RT.distribution_var = 7;              % distribution variance


end

% ------------------------------------------------------------------------



% Define the parameterization scheme used to comptue the optical quantities
% within the INP file, i.e. what function call that tells libRadtran how to
% compute scattering and optical quantities
if inputs.RT.use_custom_mie_calcs==false

    inputs.RT.wc_parameterization = 'mie interpolate';
    %inputs.RT.wc_parameterization = 'hu';

else
    inputs.RT.wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
end

% --------------------------------------------------------------
% --------------------------------------------------------------


% --------------------------------------------------------------
% ----------- Define the vertical atmospheric grid -----------
% --------------------------------------------------------------
inputs.RT.define_atm_grid=false;

% Set the vertical grid to include the cloud layers
inputs.RT.atm_z_grid = [0:0.5:inputs.RT.z_topBottom(2), inputs.RT.z_edges(2:end),...
    inputs.RT.z_topBottom(1)+1:1:20, 22:2:30, 35:5:50];   % km



% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% Define the altitude of the sensor
% How many layers to model in the cloud?
% if strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous')==true
%
%     inputs.RT.sensor_altitude = [0, sort(linspace(inputs.RT.z_topBottom(1), inputs.RT.z_topBottom(2), inputs.RT.n_layers+1))];          % top-of-atmosphere
%
% elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous')==true
%
%     inputs.RT.sensor_altitude = 'toa';
%
% end

% I think the sensor altitude, for now, is the cloud top
% inputs.RT.sensor_altitude = inputs.RT.z_topBottom(1);      % km - sensor altitude at cloud top
% inputs.RT.sensor_altitude = [0.1, 0.5, 0.9, inputs.RT.z_edges'];
% inputs.RT.sensor_altitude = [inputs.RT.z_edges'];
inputs.RT.sensor_altitude = 'toa';      % km - sensor altitude at cloud top


% define the solar zenith angle
% inputs.RT.sza = 31;           % degree
inputs.RT.sza = acosd(0.65);           % degree - for Platnick (2000)
% inputs.RT.sza = 0;           % degree


% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
%inputs.RT.phi0 = mod(293.8140 + 180, 360);
%inputs.RT.phi0 = -84 + 180;         % degree
inputs.RT.phi0 = 0;         % degree

% define the viewing zenith angle
% inputs.RT.vza = 4; % values are in degrees;                        % degree
inputs.RT.vza = acosd(0.85);                                         % degree - for Platnick (2000)
% inputs.RT.vza = 0; % values are in degrees;


% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

% inputs.RT.vaz = -103+360;     % degree
inputs.RT.vaz = 180;     % degree

% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
inputs.RT.use_coxMunk = true;
inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- specify various cross-section models for  -----------
inputs.RT.specify_cross_section_model = false;

inputs.RT.crs_model_rayleigh = 'Bodhaine29';               %  Rayleigh scattering cross section using Bodhaine et al. (1999) equation 29
% inputs.RT.crs_model_rayleigh = 'Bodhaine';                   %  Rayleigh scattering cross section using Bodhaine et al. (1999) equations 22-23

% ------------------------------------------------------------------------ 


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
inputs.RT.yesAerosols = true;

inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------




% --------------------------------------------------------
% --------- What is column water vapor amount? -----------

% Use a custom H2O profile
inputs.RT.H2O_profile = 'afglus_H2O_none_inside_cloud.dat';


% Using measurements from the AMSR2 instrument, a passive microwave
% radiometer for 17 Jan 2024
inputs.RT.modify_waterVapor = false;

inputs.RT.waterVapor_column = 0;              % mm - milimeters of water condensed in a column
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 416;       % ppm
% inputs.RT.CO2_mixing_ratio = 0;       % ppm
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to modify concentration of molecular nitrogen (N2)? ----

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_N2 = false;

inputs.RT.N2_mixing_ratio = 0;       % ppm
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to modify concentration of NO2? ----

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_NO2 = false;

inputs.RT.NO2_mixing_ratio = 0;       % ppm
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to modify concentration of molecular oxygen (O2)? ----

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_O2 = false;

inputs.RT.O2_mixing_ratio = 0;       % ppm
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want to modify concentration of Ozone (O3)? -----------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_O3 = false;

inputs.RT.O3_mixing_ratio = 0;       % ppm
% ------------------------------------------------------------------------



% --------------------------------------------------------------
% ------- Do you want to turn off molecular absorption? --------
% Note, that thermal emission of molecules is also switched off.
inputs.RT.no_molecular_abs = false;
% --------------------------------------------------------------


% --------------------------------------------------------------
% ------------ Do you want to turn off scattering? -------------

% Possible choises for the optional argument name are:
%   mol - Switch off molecular scattering.
%   aer - Switch off scattering by aerosols.
%   wc - Switch off scattering by water clouds.
%   ic - Switch off scattering by ice clouds.
%   profile - Switch off scattering by any profile defined in profile typename.
inputs.RT.no_scattering_mol = false;
inputs.RT.no_scattering_aer = false;
% --------------------------------------------------------------



% --------------------------------------------------------------
% Do you want to print an error message?
inputs.RT.errMsg = 'verbose';
% --------------------------------------------------------------





end
