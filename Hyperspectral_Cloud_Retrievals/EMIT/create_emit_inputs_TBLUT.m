%% ----- CREATE INPUTS NEEDED TO COMPUTE TBLUT METHOD ON EMIT DATA -----


% INPUTS:
%   (1) folderName -

%   (2) folderpaths - folder path where the EMIT data is located

%   (3) emit -


% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_emit_inputs_TBLUT(emitDataFolder, folder_paths, emit, spec_response)


%%

% Determine which computer this is being run on
inputs.which_computer = whatComputer;

% --- SAVE THE EMIT FILE NAME ----
inputs.emitDataFolder = emitDataFolder;

% read the contents of the EMIT data folder
folder_contents = dir([folder_paths.emitDataPath, emitDataFolder]);

% ----- Save the L1B file name -----
for nn = 1:length(folder_contents)

    if length(folder_contents(nn).name)>5 && strcmp(folder_contents(nn).name(1:12), 'EMIT_L1B_RAD')==true

        inputs.L1B_filename = folder_contents(nn).name;

    end

end


% We're not computing weighting functions
inputs.compute_weighting_functions = false;

% Define which EMIT bands to run
% band 38 has a center wavelength of 656 nm
% band 235 has a center wavelength of 2123 nm
inputs.bands2run = [38, 235]; % these are the bands that we will run uvspec with
inputs.bands2plot = [38, 235]; % these are the EMIT bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte

 % ---- Define the wavelengths ----
inputs.RT.wavelengths2run = spec_response.wavelength([inputs.bands2run],[1,end]);


% if interpGridScaleFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 150; % scale factor the will be used to increase the grid size for interpolation.



% We're running calculations over spectral bands
inputs.RT.monochromatic_calc = false;


% --------------------------------------------
% Create a new folder to save all calculations
% --------------------------------------------

% Define the folder that stores the inputs and calculated reflectanes
% using todays date
data_date = datetime([inputs.L1B_filename(18:21), '-', inputs.L1B_filename(22:23), '-', inputs.L1B_filename(24:25)],...
    'InputFormat','yyyy-MM-dd');

% Store the file name for the libRadTran INP and OUT files
inputs.folder2save.libRadTran_INP_OUT = folder_paths.libRadtran_inp;


% This is the folder where the reflectance calculations will be stored
inputs.folder2save.reflectance_calcs = [folder_paths.reflectance_calcs, emitDataFolder];






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
inputs.RT.surface_albedo = 0.04;

% day of the year
inputs.RT.day_of_year = emit.day_of_year;




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;


% define the cloud geometric depth
inputs.RT.cloud_depth = 500;                % meters


% define the geometric location of the cloud top and cloud bottom
inputs.RT.z_topBottom = [1.5, 1];          % km above surface


% Water Cloud depth
inputs.RT.H = inputs.RT.z_topBottom(2) - inputs.RT.z_topBottom(1);                                % km - geometric thickness of cloud
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
inputs.RT.distribution_str = 'gamma';
% define the distribution varaince
% 7 is the value libRadTran uses for liquid water clouds
inputs.RT.distribution_var = 7;
% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';
% define how liquid water content will be computed
% can either be 'mie' or '2limit'
inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius

% define the wavelength used for the optical depth as the 650 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 500;            % nm


% -------------------------------------------------------------------
% define the independent variable used to define the effective radius
% profile within the cloud
% -------------------------------------------------------------------
% if using altitude, z should be a vector starting at the cloud bottom
% and increasing
inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used


% --------------------------------------------------------------
% ----------- Define the vertical atmospheric grid -----------
% --------------------------------------------------------------
inputs.RT.define_atm_grid=false;
% Set the vertical grid to include the cloud layers


% -----------------------------------------------------------------------
% -------------- Define the grid of r_e and tau_c values ----------------
% -----------------------------------------------------------------------
% MODIS only considers homogenous plane parallel clouds. Lets construct the
% re matrix needed to create homogenous water clouds using write_wc_file
inputs.RT.re = 3:2:24;      % microns
inputs.RT.tau_c = [1:10, 12.5, 15, 17.5,  20:5:50, 60];


% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% Define the altitude of the sensor
inputs.RT.sensor_altitude = 'toa';          % top-of-atmosphere

% -----------------------------
% define the solar zenith angle
% -----------------------------
inputs.RT.sza = emit.obs.solar.zenith;           % degree



% -------------------------------
% define the solar azimuith angle
% -------------------------------
% this is how we map EMIT-defined solar azimuth to the LibRadTran
% definition.
% LibRadTran defines the solar azimuth clockwise from South.
% So at 0deg the Sun is due south, 90deg the Sun is due West,
% and so on. EMIT defines the solar azimuth clockwise from due North.
% So we need to add 180deg to the EMIT values, but modulo 360, since it
% needs to wrap around.
inputs.RT.phi0 = mod(emit.obs.solar.azimuth + 180, 360);         % degree


% --------------------------------
% define the viewing azimuth angle
% --------------------------------
% we need the cosine of the zenith viewing angle
% positive values solve for upwelling radiance, where the sensor is
% defined to be looking down towrads the Earth's surface. negative
% values solve for downwelling radiance, where the sensor is looking
% upwards towards the sky

% Define the cosine of the zenith viewing angle
% ------------------------------------------------
% define the viewing zenith angle
inputs.RT.vza = double(emit.obs.sensor.zenith); % values are in degrees;                        % degree


% --------------------------------
% define the viewing azimuth angle
% --------------------------------
% LibRadTran defines the viewing azimuth clockwise from North.
% So at 0deg the sensor is in the north looking south, at 90deg
% the sensor is in the east looking West, and so on.
% EMIT defines the sensor azimuth clockwise from due North.
% So we don't need to change the emit values

inputs.RT.vaz = emit.obs.sensor.azimuth;
% --------------------------------------------------------------






% ------------------------------------------------------------------------
% -------- Do you want to modify the column water vapor amount? ----------
inputs.RT.modify_total_columnWaterVapor = false;

% default value is 14.295 mm
inputs.RT.waterVapor_column = 40;       % mm (kg/m^2) - of water condensed in a column
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

% Using values measured by OCO-3 on 27 Jan 2024 off the coast of Chile
inputs.RT.CO2_mixing_ratio = 418;       % ppm - concentration of CO2
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
inputs.RT.errMsg = 'verbose';
% --------------------------------------------------------------

% --------------------------------------------------------------





% --------------------------------------------------------------
% --- Create a file name for the droplet profile retrieval -----
% --------------------------------------------------------------

rev = 1;



inputs.saveOutput_fileName = folder_paths.saveOutput_filename;




while isfile(inputs.saveOutput_fileName)
    rev = rev+1;
    inputs.saveOutput_fileName = [inputs.saveOutput_fileName(1:end-5), num2str(rev),'.mat'];
end








% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
