function hysics_refl_pt3_percent_in_situ_prof_and_tau_func_array(folder_paths, measurement_idx)
%% Generate measurements with in-situ measured droplet profiles and optical depth
% Everything else is defined by the user
% 
% INPUT:
%   measurement_idx - integer index (1-73) specifying which measurement to process

% Note: Don't clear variables or call addLibRadTran_paths / define_folderPaths_for_HySICS
% These should already be set by the calling script

%%  Delete old files?

% First, delete files in the HySICS INP folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% delete old water vapor profiles
delete([folder_paths.atm_folder_path, '*-aboveCloud.DAT'])

% delete old MIE files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])

%% Start parallel pool

start_parallel_pool(folder_paths.which_computer)

%% Load the ensemble profiles from VOCALS-REx




% load the save in-situ profiles
% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % Location of ensemble data
    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];



    ds = load([folderpath,...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_02-Nov-2025.mat']);


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % Location of ensemble data
    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
        'vocals_rex_data/NCAR_C130/SPS_1/'];

    saved_profiles_filename = ['ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25',...
        '_drizzleLWP-threshold_5_10-Nov-2025.mat'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    ds = load([folderpath, saved_profiles_filename]);



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

     % Location of ensemble data
    folderpath = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
        'vocals_rex_data/NCAR_C130/SPS_1/'];

    saved_profiles_filename = ['ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25',...
        '_drizzleLWP-threshold_5_10-Nov-2025.mat'];


    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    ds = load([folderpath, saved_profiles_filename]);

end


% define the campaign name from which the data came
campaign_name = 'vocalsRex';

%% Validate measurement_idx input

% Check that measurement_idx is valid
total_measurements = length(ds.ensemble_profiles);

if measurement_idx < 1 || measurement_idx > total_measurements
    error('measurement_idx must be between 1 and %d', total_measurements);
end

fprintf('Processing measurement %d of %d\n', measurement_idx, total_measurements);



%% Inputs needed to create the measurement input structure

% define the number of measurements being created - now just processing one
num_meas = 1;




% ----- unpack parallel for loop variables ------
% We want to avoid large broadcast variables!
libRadtran_inp = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path = folder_paths.libRadtran_water_cloud_files;
which_computer = folder_paths.which_computer;
% store which_computer in inputs structure
inputs.which_computer = which_computer;






% ---- First, let's simulate water clouds ----


% Define the parameters of the INP file
print_libRadtran_err = false;

% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';

% we model two free parameters, r_top and r_bot
inputs.RT.num_re_parameters = 2;



%% Define the inputs uvspec DISORT inputs for HySICS




% Define the parameters of the INP file

% *** compute weighting functions! ***
% This will create n wc_files where n is equal to the number of layers in
% the cloud. Starting with the entire cloud, each file will have one less
% layer
inputs.compute_weighting_functions = false;


% Are you simulating a measurement, or making forward model calculations
% for the retrieval?
inputs.calc_type = 'simulated_spectra';
% inputs.calc_type = 'forward_model_calcs_forRetrieval';
% inputs.calc_type = 'weighting_functions';



% ----- Define the RTE Solver -----
inputs.RT.rte_solver = 'disort';



% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;

% --- Do you want to use the Nakajima and Tanka radiance correction? -----


inputs.RT.use_nakajima_phaseCorrection = true;




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





% ------------------------------------------------------------------------
% ---------------------- DEFINE THE WAVELENGTHS! -------------------------
% ------------------------------------------------------------------------
% define the wavelength range. If monochromatic, enter the same number
% twice


% ----------------- Simulating HySICS spectral channels ------------------
% number of channels = 636 ranging from center wavelengths: [351, 2297]
inputs.bands2run = (1:1:636)';

% Paper 1 - Figures 7 and 8 - 35 spectral channels that avoid water vapor
% and other gaseous absorbers
% inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 217, 220,...
%     222, 224, 227, 237, 288, 290, 293, 388, 390, 393,...
%     426, 434, 436, 570, 574, 577, 579, 582, 613, 616,...
%     618, 620, 623, 625]';


% Using all 35 spectral channels above that avoid water vapor and other
% gaseous absorbers, AND 12 bands in the wings of water vapor absorption
% features for a total of 47 bands
% inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 217, 220,...
%     222, 224, 227, 237, 245, 249, 254, 264, 288, 290, 293,...
%     346, 351, 354, 360, 365, 367, 372, 379, 388, 390, 393, 426, 434, 436,...
%     570, 574, 577, 579, 582, 613, 616, 618, 620, 623, 625]';


% Using almost all 35 spectral channels above that avoid water vapor and other
% gaseous absorbers, AND 31 bands in the wings of water vapor absorption
% features for a total of 66 bands
% inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 180, 188,...
%     198, 217, 220, 222, 224, 227, 237, 245, 249, 254, 264, 288, 290, 293,...
%     346, 351, 354, 360, 365, 367, 372, 379, 388, 390, 393, 426, 434, 436,...
%     462, 468, 469, 520, 524, 525, 526, 527, 530, 531, 533, 535, 537, 539,...
%     543, 547, 570, 574, 577, 579, 582, 613, 616, 618, 620, 623, 625]';



% inputs.bands2run = [49, 426, 613]';
% inputs.bands2run = [49, 57, 288, 426, 613]';

% test bands
% 500 nm
% inputs.bands2run = 49;

% 1598 nm
% inputs.bands2run = 408;

% 1652 nm
% inputs.bands2run = 426;

% 2122 nm
% inputs.bands2run = 580;

% 2236 nm
% inputs.bands2run = 613;



% ------------------------------------------------------------------------
% Do you want to compute radiance/reflectance over a spectral region, or at
% a single wavelength?
% ------------------------------------------------------------------------
inputs.RT.monochromatic_calc = false;


% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---

inputs.RT.compute_reflectivity_uvSpec = false;
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

inputs.RT.surface_albedo = 0.04;            % Ocean water albedo
% inputs.RT.surface_albedo = 0;             % Use a value of 0 when creating weighting functions





% ----- Define the day of the year to account for Earth-Sun distance -----

% day of the year
inputs.RT.day_of_year = 316;       % value for pixel used in Figure 3.a from paper 1

% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------


inputs.RT.yesCloud = true;



% Do you want to manually set the optical depth?
inputs.RT.modify_wc_opticalDepth = true;

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

% define the wavelength used for the optical depth as the 500 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 500;            % nm
% ------------------------------------------------------------------------






% Define the parameterization scheme used to comptue the optical quantities
% within the INP file, i.e. what function call that tells libRadtran how to
% compute scattering and optical quantities

% inputs.RT.wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
inputs.RT.wc_parameterization = 'mie interpolate';


% --------------------------------------------------------------
% --------------------------------------------------------------


% --------------------------------------------------------------
% ----------- Define the vertical atmospheric grid -----------
% --------------------------------------------------------------

inputs.RT.define_atm_grid=false;






% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% Define the altitude of the sensor

% I think the sensor altitude, for now, is the cloud top
% inputs.RT.sensor_altitude = inputs.RT.z_topBottom(1);      % km - sensor altitude at cloud top
% inputs.RT.sensor_altitude = [0.1, 0.5, 0.9, inputs.RT.z_edges'];
% inputs.RT.sensor_altitude = [inputs.RT.z_edges'];
inputs.RT.sensor_altitude = 'toa';      % km - sensor altitude at cloud top







% -----------------------------
% define the solar zenith angle
% -----------------------------
inputs.RT.sza = 31;               % degree - value for pixel used in Figure 3.a from paper 1
% inputs.RT.sza = acosd(0.65);           % degree - for Platnick (2000)
% inputs.RT.sza = 10;           % degree




% -------------------------------
% define the solar azimuith angle
% -------------------------------

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
%inputs.RT.phi0 = mod(293.8140 + 180, 360);
inputs.RT.phi0 = -84 + 180;         % degree - value for pixel used in Figure 3.a from paper 1
% inputs.RT.phi0 = 0;         % degree





% --------------------------------
% define the viewing azimuth angle
% --------------------------------


inputs.RT.vza = 4.29;                                   % degree - value for pixel used in Figure 3.a from paper 1
% inputs.RT.vza = acosd(0.75);                              % degree - for Platnick (2000)
% inputs.RT.vza = 20;                             % values are in degrees;




% --------------------------------
% define the viewing azimuth angle
% --------------------------------


% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

% inputs.RT.vaz = -103+360;     % degree
% inputs.RT.vaz = 210;            % degree

% to properly map the MODIS azimuth angle onto the reference plane used by
% libRadTran...
inputs.RT.vaz = 360 + -102.79;     % degree - value for pixel used in Figure 3.a from paper 1



% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
inputs.RT.use_coxMunk = true;
inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- specify various cross-section models for  -----------
inputs.RT.specify_cross_section_model = true;

inputs.RT.crs_model_rayleigh = 'Bodhaine29';               %  Rayleigh scattering cross section using Bodhaine et al. (1999) equation 29
% inputs.RT.crs_model_rayleigh = 'Bodhaine';                   %  Rayleigh scattering cross section using Bodhaine et al. (1999) equations 22-23

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
inputs.RT.yesAerosols = true;

inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------




% --------------------------------------------------------------------
% --------- What is total column water vapor amount? -----------------

% Using measurements from the AMSR2 instrument, a passive microwave
% radiometer for 17 Jan 2024
inputs.RT.modify_total_columnWaterVapor = false;

inputs.RT.waterVapor_column = 20;   % mm - milimeters of water condensed in a column
% ------------------------------------------------------------------------



% -----------------------------------------------------------------------
% -------- Write a custom water vapor profile for above cloud -----------

% Alter the above cloud column water vapor amount
inputs.RT.modify_aboveCloud_columnWaterVapor = true;

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

if print_libRadtran_err==true

    inputs.RT.errMsg = 'verbose';

else

    inputs.RT.errMsg = 'quiet';

end
% --------------------------------------------------------------
% --------------------------------------------------------------
% --------------------------------------------------------------
% --------------------------------------------------------------











%% Define the geometry

% Numbers based off NOAA solar calculator (27 Jan 2024 off coast of Chile)
% inputs.RT.sza = 10;         % degrees - Solar Zenith Angle
% inputs.RT.phi0 = 91.45;        % degrees - Solar Azimuth Angle
% inputs.RT.vza = 7;        % degrees - Viewing Zenith Angle
% inputs.RT.vaz = 210;       % degrees - Viewing Azimuth Angle

% Replicating EMIT values from paper 1
inputs.RT.sza = 31;         % degrees - Solar Zenith Angle
inputs.RT.phi0 = 96;        % degrees - Solar Azimuth Angle
inputs.RT.vza = 4.29;        % degrees - Viewing Zenith Angle
inputs.RT.vaz = 257.21;       % degrees - Viewing Azimuth Angle




%% Set the total column water vapor?

inputs.RT.modify_total_columnWaterVapor = false;             % modify the full column

inputs.RT.modify_aboveCloud_columnWaterVapor = false;         % don't modify the column above the cloud



%%

% num wavelengths
num_wl = length(inputs.bands2run);




num_INP_files = num_meas * num_wl;


inputs.calc_type = 'simulated_spectra';


% First, write all the wc files - now just for the single measurement
wc_filename = cell(num_meas, 1);

% for storage
tau_c = zeros(num_meas, 1);
date_of_flight = cell(num_meas, 1);
time_of_flight = zeros(num_meas, 1);
re = cell(num_meas, 1);
lwc = cell(num_meas, 1);
z = cell(num_meas, 1);


% --------------------------------------------------------
% ---- Write Water Cloud file for single measurement! ----
% --------------------------------------------------------

% Use nn=1 since we're only processing one measurement
nn = 1;

% define lwc, re, and z as column vectors
% grab the LWC vector
lwc{nn} = ds.ensemble_profiles{measurement_idx}.lwc';     % g/m^3
% grab the altitude vector
z{nn} = (ds.ensemble_profiles{measurement_idx}.altitude') ./ 1e3;   % kilometers
% grab the date of flight
date_of_flight{nn} = ds.ensemble_profiles{measurement_idx}.dateOfFlight;
% grab the time of flight
time_of_flight(nn) = ds.ensemble_profiles{measurement_idx}.time_utc(round(length(ds.ensemble_profiles{measurement_idx}.time_utc)/2));  % UTC time

% store the optical depth of each profile
tau_c(nn) = max(ds.ensemble_profiles{measurement_idx}.tau);


if isfield(ds.ensemble_profiles{measurement_idx}, 're') == true

    % define the effective radius profile
    re{nn} = ds.ensemble_profiles{measurement_idx}.re';        % microns

    % rearrange z, lwc, and re so that all profiles are stored from cloud
    % top to bottom. Therefore, the first value in the column vector is
    % cloud top
    if (z{nn}(2) - z{nn}(1))>0
        % The the profile starts from cloud bottom and needs to be flipped
        z{nn} = flipud(z{nn});
        lwc{nn} = flipud(lwc{nn});
        re{nn} = flipud(re{nn});
    end

    % Sometimes droplets will be larger than 25 microns. If there are droplets
    % larger than 25 microns, skip this in-situ measurement for now. libRadtran
    % cannot process droplets larger than 25 microns, which is the limit of the
    % pre-computed mie table

    if any(re{nn}>=25)==true

        if strcmp(saved_profiles_filename, ['ensemble_profiles_without_precip_from_14_files_LWC',...
                '-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat'])==true &&...
                measurement_idx==11

            % only the last two measurements exceed 25 microns. Remove
            % them and use this profile
            idx_remove = re{nn}>=25;
            re{nn}(idx_remove) = [];
            lwc{nn}(idx_remove) = [];
            z{nn}(idx_remove) = [];

            tau_c(nn) = ds.ensemble_profiles{measurement_idx}.tau(end-2);

            wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
                date_of_flight{nn}, time_of_flight(nn),...
                inputs.compute_weighting_functions, which_computer,...
                measurement_idx, wc_folder_path);

        elseif strcmp(saved_profiles_filename, ['ensemble_profiles_without_precip_from_14_files_LWC',...
                '-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat'])==true &&...
                measurement_idx==15


            % only the last measurement exceeds 25 microns. Remove
            % them and use this profile
            idx_remove = re{nn} >= 25;
            re{nn}(idx_remove) = [];
            lwc{nn}(idx_remove) = [];
            z{nn}(idx_remove) = [];

            tau_c(nn) = ds.ensemble_profiles{measurement_idx}.tau(end-1);

            wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
                date_of_flight{nn}, time_of_flight(nn),...
                inputs.compute_weighting_functions, which_computer,...
                measurement_idx, wc_folder_path);


        elseif strcmp(saved_profiles_filename, ['ensemble_profiles_without_precip_from_14_files_LWC',...
                '-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat'])==true &&...
                measurement_idx==22


            % only the last measurement exceeds 25 microns. Remove
            % them and use this profile
            idx_remove = re{nn} >= 25;
            re{nn}(idx_remove) = [];
            lwc{nn}(idx_remove) = [];
            z{nn}(idx_remove) = [];

            tau_c(nn) = ds.ensemble_profiles{measurement_idx}.tau(end-1);

            wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
                date_of_flight{nn}, time_of_flight(nn),...
                inputs.compute_weighting_functions, which_computer,...
                measurement_idx, wc_folder_path);


        elseif strcmp(saved_profiles_filename, ['ensemble_profiles_without_precip_from_14_files_LWC',...
                '-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat'])==true &&...
                measurement_idx==36


            % only the last 12 measurements exceed 25 microns. Remove
            % them and use this profile
            idx_remove = re{nn} >= 25;
            re{nn}(idx_remove) = [];
            lwc{nn}(idx_remove) = [];
            z{nn}(idx_remove) = [];

            tau_c(nn) = ds.ensemble_profiles{measurement_idx}.tau(end-12);

            wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
                date_of_flight{nn}, time_of_flight(nn),...
                inputs.compute_weighting_functions, which_computer,...
                measurement_idx, wc_folder_path);


        elseif strcmp(saved_profiles_filename, ['ensemble_profiles_without_precip_from_14_files_LWC',...
                '-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat'])==true &&...
                measurement_idx==50


            % only the exclude the 60th measurement, which jumps to 30 microns. Remove
            % and use this profile
            idx_remove = re{nn} >= 25;
            re{nn}(idx_remove) = [];
            lwc{nn}(idx_remove) = [];
            z{nn}(idx_remove) = [];

            wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
                date_of_flight{nn}, time_of_flight(nn),...
                inputs.compute_weighting_functions, which_computer,...
                measurement_idx, wc_folder_path);


        elseif strcmp(saved_profiles_filename, ['ensemble_profiles_without_precip_from_14_files_LWC',...
                '-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat'])==true &&...
                measurement_idx==54


            % only the last measurement exceeds 25 microns. Remove
            % them and use this profile
            idx_remove = re{nn} >= 25;
            re{nn}(idx_remove) = [];
            lwc{nn}(idx_remove) = [];
            z{nn}(idx_remove) = [];

            tau_c(nn) = ds.ensemble_profiles{measurement_idx}.tau(end-1);

            wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
                date_of_flight{nn}, time_of_flight(nn),...
                inputs.compute_weighting_functions, which_computer,...
                measurement_idx, wc_folder_path);


        else

            wc_filename{nn} = NaN;

        end

    else

        % At times there are very thin gaps in the clouds
        % as well. In these cases the droplet sizes go to 0. Get rid of
        % these levels

        idx_0 = re{nn}<=0.1;
        sum(idx_0)
        if sum(idx_0)>0
            re{nn}(idx_0) = [];
            lwc{nn}(idx_0) = [];
            z{nn}(idx_0) = [];
        end

        wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
            date_of_flight{nn}, time_of_flight(nn),...
            inputs.compute_weighting_functions, which_computer,...
            measurement_idx, wc_folder_path);


    end



elseif isfield(ds.ensemble_profiles{measurement_idx}, 're_CDP') == true


    % define the effective radius profile
    re{nn} = ds.ensemble_profiles{measurement_idx}.re_CDP';        % microns

    % rearrange z, lwc, and re so that all profiles are stored from cloud
    % top to bottom. Therefore, the first value in the column vector is
    % cloud top
    if (z{nn}(2) - z{nn}(1))>0
        % The the profile starts from cloud bottom and needs to be flipped
        z{nn} = flipud(z{nn});
        lwc{nn} = flipud(lwc{nn});
        re{nn} = flipud(re{nn});
    end

    % Sometimes droplets will be larger than 25 microns. If there are droplets
    % larger than 25 microns, skip this in-situ measurement for now. libRadtran
    % cannot process droplets larger than 25 microns, which is the limit of the
    % pre-computed mie table

    if any(re{nn}>=25)==true

        wc_filename{nn} = NaN;

    else

        % At times there are very thin gaps in the clouds
        % as well. In these cases the droplet sizes go to 0. Get rid of
        % these levels

        idx_0 = re{nn}<=0.1;
        sum(idx_0)
        if sum(idx_0)>0

            re{nn}(idx_0) = [];
            lwc{nn}(idx_0) = [];
            z{nn}(idx_0) = [];
        end


        wc_filename{nn} = write_wc_file_from_in_situ(re{nn}, lwc{nn}, z{nn}, campaign_name,...
            date_of_flight{nn}, time_of_flight(nn),...
            inputs.compute_weighting_functions, which_computer,...
            measurement_idx, wc_folder_path);

    end



end


% *** we dont' need the ensemble profiles anymore ***
clear ds

% changing variable steps through each wavelength for this single measurement
% in for loop speak, it would be:
%   for ww = 1:num_wl
changing_variables_allStateVectors = {reshape(repmat( (1:num_meas), num_wl, 1), [],1),...
    reshape(repmat(tau_c', num_wl, 1), [], 1),...
    reshape(repmat(date_of_flight', num_wl, 1), [], 1),...
    reshape(repmat(time_of_flight', num_wl, 1), [], 1),...
    reshape(repmat(wc_filename', num_wl, 1), [], 1),...
    repmat(inputs.RT.wavelengths2run, num_meas, 1)};


% Add a final column that includes the index for the spectral response
% function. These always increase chronologically
changing_variables_allStateVectors = [changing_variables_allStateVectors,...
    repmat((1:num_wl)', num_meas, 1)];


% Double check that the number of rows for 'changing_variables' is equal
% to the number of INP files
if num_INP_files ~= size(changing_variables_allStateVectors{1}, 1)

    error([newline, 'Number of rows in changing_variables not equal to number of INP files', newline])

end


% --------------------------------------
% ------- Write all INP/OUT Files ------
% --------------------------------------
inputFileName = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);

parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files


    % set the wavelengths for each file
    wavelengths = changing_variables_allStateVectors{end-1}(nn,:);

    % ------------------------------------------------
    % ---- Define the input and output filenames! ----
    % ------------------------------------------------
    % input_names need a unique identifier. Let's give them the measurement_idx value and nn value
    % so they can be traced, and are written over in memory


    inputFileName{nn} = [num2str(mean(wavelengths)), 'nm_', campaign_name,'_',...
        char(changing_variables_allStateVectors{3}{nn}),'_',...
        num2str(changing_variables_allStateVectors{4}(nn)),'-UTC_meas_',...
        num2str(measurement_idx), '_nn_', num2str(nn), '.INP'];



    outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


    % ------------------ Write the INP File --------------------
    % force modify tau
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName{nn}, inputs,...
        wavelengths, changing_variables_allStateVectors{5}{nn}, [], changing_variables_allStateVectors{2}(nn));


end










%% Calculate Reflectance

% Read the solar flux file over the wavelength range specified
wavelength_vec = [min(inputs.RT.wavelengths2run,[],"all"), max(inputs.RT.wavelengths2run, [], "all")];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);   % W/nm/m^2

% we will add and subtract a small fraction of the source file resolution
% to ensure rounding errors don't cause an issue when selecting the
% wavelengths needed from the source file
wl_perturb = inputs.RT.source_file_resolution/3;   % nm



% define only the spec_response so the wavelengths are passed into the
% memory of the parallel for loop
spec_response_value = spec_response.value;


tic



% store the reflectances
Refl_model_allStateVectors = zeros(num_INP_files, 1);


parfor nn = 1:num_INP_files


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------

    % compute INP file
    runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn},which_computer);


    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, inputs,...
        inputs.RT.compute_reflectivity_uvSpec);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    idx_wl = source_wavelength>=(changing_variables_allStateVectors{end-1}(nn,1) - wl_perturb) &...
        source_wavelength<=(changing_variables_allStateVectors{end-1}(nn,2) + wl_perturb);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    [Refl_model_allStateVectors(nn), ~] = reflectanceFunction_ver2(inputs, ds,...
        source_flux(idx_wl), spec_response_value(changing_variables_allStateVectors{end}(nn),:));



end


toc


%% Rearrange the reflectances

Refl_model_allStateVectors = reshape(Refl_model_allStateVectors, num_wl, []);



%% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
% inputs.measurement.uncert = [0.003, 0.01:0.01:0.1];
inputs.measurement.uncert = 0.003;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve

% Compute the new synethtic measurement with gaussian noise
% *** Gaussian noise can be either positive or negative. Meaning, an
% uncertainty of 5% implies the true value can lie anywhere between
% +/- 5% of the measured value

% To sample a normal distribution with mean mu, and standard deviation s,
% we compute the following: y = s * randn() + mu

% We define the standard deviation as the measurement uncertainty divided
% by three. Therefore, after sample a large number of times, 99% of
% measurements will be within +/- measurement uncertainy of the mean

if any(inputs.measurement.uncert > 0)

    inputs.measurement.standard_dev = inputs.measurement.uncert/3;       % this is still just a fraction

    for uu = 1:length(inputs.measurement.uncert)

        clear Refl_model_with_noise_allStateVectors Refl_model_uncert_allStateVectors


        Refl_model_with_noise_allStateVectors = (inputs.measurement.standard_dev(uu) .* Refl_model_allStateVectors) .* randn(size(Refl_model_allStateVectors))...
            + Refl_model_allStateVectors;


        % define the synthetic relfectance uncertainty
        Refl_model_uncert_allStateVectors = inputs.measurement.uncert(uu) .* Refl_model_with_noise_allStateVectors;    % 1/sr



    end

end



%%
% ----------------------------------------------
% ---------- SAVE REFLECTANCE OUTPUT! ----------
% ----------------------------------------------

% Save the version without an measurement uncertainty. Then we can add
% uncertainty and save the new file

if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    inputs.folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/',...
        'log_newCov_subset_allBands_VR_inSitu_1/'];



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    inputs.folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/paper2_variableSweep/log_newCov_subset_allBands_VR_inSitu_1/'];


elseif strcmp(which_computer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % inputs.folderpath_2save = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
    %     'Simulated_spectra/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/'];

    % inputs.folderpath_2save = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
    %     'Simulated_spectra/paper2_variableSweep/rTop_10/vza_4_vaz_257_sza_31_saz_96_subset_newRetrieval3/'];


    inputs.folderpath_2save = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
        'Simulated_spectra/paper2_variableSweep/log_newCov_subset_allBands_VR_inSitu_1/'];



end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_2save, 'dir')

    mkdir(inputs.folderpath_2save)

end


% *** save all spectral calculations for this single in-situ measurement as one ***

for nn = 1:(num_INP_files/num_wl)

    % Grab the reflectance group for the given state vector
    Refl_model = Refl_model_allStateVectors(:,nn);
    Refl_model_with_noise = Refl_model_with_noise_allStateVectors(:,nn);
    Refl_model_uncert = Refl_model_uncert_allStateVectors(:,nn);



    inputs.RT.re = re{nn};
    inputs.RT.tau_c = tau_c(nn);
    inputs.RT.z = z{nn};
    inputs.RT.lwc = lwc{nn};
    inputs.RT.z_topBottom = [z{nn}(1), z{nn}(end)]; % km

    filename = [inputs.folderpath_2save,'simulated_spectra_HySICS_reflectance_',...
        num2str(numel(inputs.bands2run)), 'bands_',num2str(100*inputs.measurement.uncert), '%_uncert_',...
        campaign_name, '_inSitu_re_lwc_tauC_z_', char(date_of_flight{nn}),'_', num2str(time_of_flight(nn)),...
        'UTC_prof-nn_', num2str(measurement_idx), '_vza_', num2str(round(inputs.RT.vza)),...
        '_vaz_', num2str(round(inputs.RT.vaz)), '_sza_', num2str(round(inputs.RT.sza)),...
        '_saz_', num2str(round(inputs.RT.phi0)),...
        '_sim-ran-on-',char(datetime("today")),'.mat'];


    save(filename, "Refl_model", "Refl_model_with_noise", "Refl_model_uncert",...
        "inputs", "spec_response");

end

fprintf('Successfully completed measurement %d\n', measurement_idx);

end
