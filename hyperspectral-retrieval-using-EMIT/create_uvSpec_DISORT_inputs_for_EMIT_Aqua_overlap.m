%% Create input structure for libRadtran INP files for real EMIT data

% (1) inputs - structure that contains the inputs to run radiative transfer
% files. Likely this only includes a few things like which_computer etc.
% This very function provides many more inputs.


% (3) emit - this is the structure containing information on the EMIT
% measurements, such as wavelengths, solar and viewing geometry, and
% spectral response information

% ** Don't include pixel level attributes in this function **



% By Andrew John Buggee

%%

function inputs = create_uvSpec_DISORT_inputs_for_EMIT_Aqua_overlap(inputs, emit, print_libRadtran_err)


% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 8 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin>=2

else

    error([newline,'Should be at least 2 inputs: an input strucutre with at least the computer name ,',...
        [' and a flag telling the code to use parameters from a simualted measurement. If this flag ,',...
        ' is true, then a third input is required that includes the input structure for the ',...
        'simulated measurements.'], newline])
end






% Define the parameters of the INP file





% *** compute weighting functions! ***
% This will create n wc_files where n is equal to the number of layers in
% the cloud. Starting with the entire cloud, each file will have one less
% layer
inputs.compute_weighting_functions = false;


% Are you simulating a measurement, or making forward model calculations
% for the retrieval?
% inputs.calc_type = 'simulated_spectra';
inputs.calc_type = 'forward_model_calcs_forRetrieval';
% inputs.calc_type = 'weighting_functions';


% ----- Define the RTE Solver -----
inputs.RT.rte_solver = 'disort';



% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;





% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----

inputs.RT.use_nakajima_phaseCorrection = true;



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





% ------------------------------------------------------------------------
% ---------------------- DEFINE THE WAVELENGTHS! -------------------------
% ------------------------------------------------------------------------
% define the wavelength range. If monochromatic, enter the same number
% twice



% -------------- Define which EMIT bands to run ------------------
% ****************************************************************
% *-*-*-*- Only keep wavelengths that avoid water vapor -*-*-*-*-*

% --- Use all 285 spectral channels -
% inputs.bands2run = (1:285)';



% --- First 7 MODIS spectral channels -
% inputs.bands2run = [17, 24, 36, 65, 116, 170, 235]';

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
% inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 86, 87, 88, 89, 90,...
%     94, 115, 116, 117, 156, 157, 158, 172, 175, 176,...
%     231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';



% --- Indexes using same 35 as above, in addition to 29 water vapor bands ---
% This set has a total of 64 bands. They are not exactly the same set as
% the 66 HySICS bands used to retrieve column water vapor because the
% HySICS channels are more narrow.
inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 71, 74, 78, 86, 87, 88, 89, 90,...
    94, 97, 99, 101, 105, 115, 116, 117, 139, 141, 142, 145, 147, 148, 149, 151, 156,...
    157, 158, 172, 175, 176, 187, 189, 190, 210, 212, 213, 214, 215, 216, 217, 218, 219,...
    220, 222, 231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------



inputs.bands2plot = [38, 235];    % these are the EMIT bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte






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
inputs.RT.day_of_year = emit.day_of_year;       % value for pixel used in Figure 3.a from paper 1


% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------

inputs.RT.yesCloud = true;





% define the geometric location of the cloud top and cloud bottom
% inputs.RT.z_topBottom = [1.5, 1];          % km above surface
inputs.RT.z_topBottom = [1.25, 0.75];          % km above surface  - value for pixel used in Figure 3.a from paper 1



% Water Cloud depth
inputs.RT.H = inputs.RT.z_topBottom(1) - inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud

% Do you want to manually set the optical depth?
inputs.RT.modify_wc_opticalDepth = true;

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --------------------- Various Cloud modeling inputs --------------------
% ------------------------------------------------------------------------
% Do you want use your custom mie tables to conver optical props?
% If false, the default file is the precomputed mie table from libRadtran
% If true, it will use the 'custom_mieTables' folder to convert cloud
% optical properties to mie properties in order to compute the radiative
% transfer
inputs.RT.use_custom_mie_calcs = true;


% define how liquid water content will be computed in the write_wc_file
% function.
inputs.RT.parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 500 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 500;            % nm
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------------- Cloud optical properties --------------------------
% ------------------------------------------------------------------------


if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % --------- HOMOGENOUS CLOUD -------------



    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    inputs.RT.distribution_var = 10;              % distribution variance

    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    inputs.RT.re = 5:2:9;      % microns
    inputs.RT.tau_c = 5:3:20;

    % inputs.RT.re = 1:2:51;      % microns
    % inputs.RT.tau_c = [1:15, 20:5:100];


elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------- NON-HOMOGENOUS CLOUD -------------

    if inputs.RT.num_re_parameters==2

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

        inputs.RT.profile_type = 'adiabatic'; % type of water droplet profile


        % ** Values used in Platnick (2000) **
        % inputs.RT.r_top = 12;     % microns
        % inputs.RT.r_bot = 5;        % microns
        % inputs.RT.tau_c = 8;

        % simulating measurement from value for pixel used in Figure 3.a from paper 1
        inputs.RT.r_top = 9.2516;     % microns
        inputs.RT.r_bot = 5.3192;        % microns
        inputs.RT.tau_c = 6.1312;

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




    elseif inputs.RT.num_re_parameters==3

        % We have three variables, r_top, r_middle, and r_bot

        inputs.RT.profile_type = 'linear_with_tau'; % type of water droplet profile

        % ** Values used in Platnick (2000) **
        % inputs.RT.r_top = 12;     % microns
        % inputs.RT.r_bot = 5;        % microns
        % inputs.RT.tau_c = 8;

        % simulating measurement from value for pixel used in Figure 3.a from paper 1
        inputs.RT.r_top = 9.2516;     % microns
        inputs.RT.r_middle = 7;       % microns
        inputs.RT.r_bot = 5.3192;     % microns
        inputs.RT.tau_c = 6.1312;




    else

        error([newline, 'I need to know the number of free parameters defining the droplet profile', newline])


    end

    % *** Use 250 if creating weighting functions using DISORT ***
    % inputs.RT.n_layers = 250;                          % number of layers to model within cloud
    inputs.RT.n_layers = 20;                          % number of layers to model within cloud

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





    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***

    inputs.RT.distribution_var = linspace(10,10, inputs.RT.n_layers);              % distribution variance



end

% ------------------------------------------------------------------------



% Define the parameterization scheme used to comptue the optical quantities
% within the INP file, i.e. what function call that tells libRadtran how to
% compute scattering and optical quantities
if inputs.RT.use_custom_mie_calcs==false

    inputs.RT.wc_parameterization = 'mie interpolate';
    %inputs.RT.wc_parameterization = 'hu';

else
    % inputs.RT.wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
    warning([newline, 'Set the wc_parameterization with the directory and file for the custom mie table desired.', newline])
    
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









% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----

inputs.RT.use_coxMunk = true;
inputs.RT.wind_speed = 3;             % m/s


% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- specify various cross-section models for  ------------------

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
inputs.RT.modify_aboveCloud_columnWaterVapor = false;

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% -------- Do you want to use the effective water vapor depth? -----------

% This comes from the retrieval of thermodynamic phase
inputs.RT.use_phaseRetrieval_columnWaterVapor = false;
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

% Using values measured by OCO-3 on 27 Jan 2024 off the coast of Chile
inputs.RT.CO2_mixing_ratio = 418;       % ppm
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



















end
