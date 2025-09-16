%% ----- CREATE INPUTS NEEDED TO COMPUTE ACPW RETRIEVAL ON HySICS DATA -----


% INPUTS:
%   (1) folder_paths - folder path where the files should be written and
%   the caluclations stored

%   (2) inputs_measurement - we need to set up the inputs of the TBLUT
%   retrieval to be the same as those defined by the inputs of the
%   measurement. When we have true measurements, all we know is the
%   geometry of the sun and sensor. So this is the only information that
%   should be passed along.

% (3) print_libRadtran_err - true or false that tells the function to
% write and save the libRadtran error message file



% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_HySICS_inputs_ACPW(inputs_measurement, tblut_retrieval, print_libRadtran_err)






%%

% We're not computing weighting functions
inputs.compute_weighting_functions = false;



% Define which HySICS bands to run
% number of channels = 636 ranging from center wavelengths: [351, 2297]
% band 125 has a center wavelength of 731.75 nm
% band 154 has a center wavelength of 820.55 nm
% band 171 has a center wavelength of 872.65 nm
% band 180 has a center wavelength of 900.25 nm
% band 198 has a center wavelength of 955.35 nm
% band 254 has a center wavelength of 1127 nm
% inputs.bands2run = [171, 180, 254]; % these are the bands that we will run uvspec with
inputs.bands2run = [180, 198, 254]; % these are the bands that we will run uvspec with
inputs.bands2plot = inputs.bands2run;

% We're running calculations over spectral bands
inputs.RT.monochromatic_calc = false;

% if interpGridScaleFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 150; % scale factor the will be used to increase the grid size for interpolation.





% ------------------
% ----- FLAGS! -----
% ------------------

% define flags that tell the codes to either run certain things, or don't
% run certain things

inputs.flags.writeINPfiles = true; % if true, this will create inp files for each the length of vector pixel.row
inputs.flags.runUVSPEC = true; % if true, this will run all of the inp files create from the above flag through uvspec
inputs.flags.plotMLS_figures = false; % this will tell the leasSquaresGridSearch code to plot





% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------

% Define the RTE Solver
inputs.RT.rte_solver = 'disort';

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
inputs.RT.band_parameterization = inputs_measurement.RT.band_parameterization;
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

inputs.RT.source_file = inputs_measurement.RT.source_file;
inputs.RT.source_file_resolution = inputs_measurement.RT.source_file_resolution;         % nm





% define the atmospheric data file
inputs.RT.atm_file = inputs_measurement.RT.atm_file;

% define the surface albedo
inputs.RT.surface_albedo = inputs_measurement.RT.surface_albedo;

% day of the year
if isfield(inputs_measurement.RT, 'day_of_year')
    % day of the year
    inputs.RT.day_of_year = inputs_measurement.RT.day_of_year;       % value for pixel used in Figure 3.a from paper 1

end



% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% *** Use the TBLUT retreival estimates ***
inputs.RT.re = tblut_retrieval.minRe;     % microns
inputs.RT.tau_c = tblut_retrieval.minTau;



% define the cloud geometric depth
inputs.RT.cloud_depth = inputs_measurement.RT.cloud_depth;                % meters


% define the geometric location of the cloud top and cloud bottom
inputs.RT.z_topBottom = inputs_measurement.RT.z_topBottom;          % km above surface


% Water Cloud depth
inputs.RT.H = inputs.RT.z_topBottom(2) - inputs.RT.z_topBottom(1);                                % km - geometric thickness of cloud


% -------------------------------------------------------------------
% define the independent variable used to define the effective radius
% profile within the cloud
% -------------------------------------------------------------------
% if using altitude, z should be a vector starting at the cloud bottom
% and increasing
inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used


% ------------------------------------------------------------------------


% --------------------------------------------------------------
% ----------- Define the vertical atmospheric grid -----------
% --------------------------------------------------------------
inputs.RT.define_atm_grid=false;
% Set the vertical grid to include the cloud layers






% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
inputs.RT.use_custom_mie_calcs = false;
% ------------------------------------------------------------------------
% This string is used to compute the LWC from optical depth and effective radius
% can be 'hu' or 'mie interpolate'
inputs.RT.wc_parameterization = 'mie interpolate';        % use the hu and stamnes parameterization for converting cloud properties to optical properties
% define the type of droplet distribution
inputs.RT.distribution_str = 'gamma';
% define the distribution varaince
% 7 is the value libRadTran uses for liquid water clouds
inputs.RT.distribution_var = 7;
% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';

% We're modeling a homoegenous cloud layer using the TBLUT retrieval, so
% the number of free re parameters is 1
inputs.RT.num_re_parameters = 1;

% define how liquid water content will be computed
% can either be 'mie' or '2limit'
inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius

% define the wavelength used for the optical depth as the 650 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 500;            % nm


% Do you want to manually set the optical depth?
inputs.RT.modify_wc_opticalDepth = true;

% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% Define the altitude of the sensor
inputs.RT.sensor_altitude = inputs_measurement.RT.sensor_altitude;          % top-of-atmosphere

% define the solar zenith angle
inputs.RT.sza = inputs_measurement.RT.sza;           % degree

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
inputs.RT.phi0 = inputs_measurement.RT.phi0;         % degree

% define the viewing zenith angle
inputs.RT.vza = inputs_measurement.RT.vza; % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

inputs.RT.vaz = inputs_measurement.RT.vaz;     % degree
% --------------------------------------------------------------




% ------------------------------------------------------------------------
% -------- Do you want to modify the column water vapor amount? ----------
inputs.RT.modify_waterVapor = false;

% default value is 14.295 mm
inputs.RT.waterVapor_column = 40;       % mm (kg/m^2) - of water condensed in a column
% ------------------------------------------------------------------------



% -----------------------------------------------------------------------
% -------- Write a custom water vapor profile for above cloud -----------

% Alter the above cloud column water vapor amount
inputs.RT.modify_aboveCloud_columnWaterVapor = true;
% -----------------------------------------------------------------------



% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 416;       % ppm - concentration of CO2
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
% --------- specify various cross-section models for  -----------
inputs.RT.specify_cross_section_model = false;

inputs.RT.crs_model_rayleigh = 'Bodhaine29';               %  Rayleigh scattering cross section using Bodhaine et al. (1999) equation 29
% inputs.RT.crs_model_rayleigh = 'Bodhaine';                   %  Rayleigh scattering cross section using Bodhaine et al. (1999) equations 22-23

% ------------------------------------------------------------------------


% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
inputs.RT.compute_reflectivity_uvSpec = false;
% --------------------------------------------------------------


% --------------------------------------------------------------

% ----- Do you want a long error message? -----
% if so, set error message to 'verbose'. Otherwise, set error message to
% 'quiet'

if print_libRadtran_err==true

    inputs.RT.errMsg = 'verbose';

else

    inputs.RT.errMsg = 'quiet';

end

% --------------------------------------------------------------

% --------------------------------------------------------------













% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
