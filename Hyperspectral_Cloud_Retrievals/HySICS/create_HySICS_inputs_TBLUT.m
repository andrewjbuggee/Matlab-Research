%% ----- CREATE INPUTS NEEDED TO COMPUTE TBLUT METHOD ON EMIT DATA -----


% INPUTS:
%   (1) folder_paths - folder path where the files should be written and
%   the caluclations stored

%   (2) inputs_measurement - we need to set up the inputs of the TBLUT
%   retrieval to be the same as those defined by the inputs of the
%   measurement. When we have true measurements, all we know is the
%   geometry of the sun and sensor. So this is the only information that
%   should be passed along. 



% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_HySICS_inputs_TBLUT(folder_paths, inputs_measurement)


%% Find computer and folders

% Determine which computer this is being run on
inputs.which_computer = whatComputer;


% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(inputs.which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';




elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];






elseif strcmp(inputs.which_computer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


end





%%





% Define which HySICS bands to run
% number of channels = 636 ranging from center wavelengths: [351, 2297]
% band 98 has a center wavelength of 649 nm
% band 582 has a center wavelength of 2131 nm
inputs.bands2run = [98, 582]; % these are the bands that we will run uvspec with
inputs.bands2plot = inputs.bands2run;

% We're running calculations over spectral bands
inputs.RT.monochromatic_calc = false;

% if interpGridScaleFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 150; % scale factor the will be used to increase the grid size for interpolation.


% --------------------------------------------
% Create a new folder to save all calculations
% --------------------------------------------


% Store the file name for the libRadTran INP and OUT files
inputs.save_inp_files = [folder_paths.libRadtran_inp, 'TBLUT_retrieval_',char(datetime("today")),'/'];




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
%inputs.RT.day_of_year = simulated_measurements.day_of_year;




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

inputs.RT.re = 3:2:24;      % microns
inputs.RT.tau_c = [1:10, 12.5, 15, 17.5,  20:5:50, 60];

% inputs.RT.re = 3:2:11;      % microns
% inputs.RT.tau_c = [1:10];

% define the cloud geometric depth
inputs.RT.cloud_depth = 500;                % meters

% define the geometric location of the cloud top and cloud bottom
inputs.RT.z_topBottom = [1.5, 1];          % km above surface


% Water Cloud depth
inputs.RT.H = inputs.RT.z_topBottom(1) - inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud
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

% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% Define the altitude of the sensor
inputs.RT.sensor_altitude = 'toa';          % top-of-atmosphere

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



% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 416;       % ppm - concentration of CO2
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



    inputs.save_mat_filename = [folder_paths.HySICS_retrievals,'droplet_profile_retrieval_',...
        'sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];




while isfile(inputs.save_mat_filename)
    rev = rev+1;
    inputs.save_mat_filename = [folder_paths.HySICS_retrievals,'droplet_profile_retrieval_',...
        'sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];
end

% --------------------------------------------------------------
% --------------------------------------------------------------







% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
