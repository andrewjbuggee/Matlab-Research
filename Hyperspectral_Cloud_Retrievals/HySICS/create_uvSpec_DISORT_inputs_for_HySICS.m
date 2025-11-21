%% Create input structure for libRadtran INP files for the HySICS instrument

% (1) inputs - structure that contains the inputs to run radiative transfer
% files. Likely this only includes a few things like which_computer etc.
% This very function provides many more inputs.

% (2) load_parameters_from_measurement - this is a flag that, if true, will
% use the same input as the defined measurement for the GN inputs. But we
% will only use inputs that we could know from a real measurement, such as
% the geometry

% (3) sim_meas - this is the input structure used to create
% the simulated measurements. This includes the input structure and the
% spectral response function

% (4) sim_meas_likeness - string that tells the function how similar you
% want the forward model to be to the simulated measurements. Options are:
%   (a) 'exact' - the settings are exactly the same
%   (b) 'subset'

% (5) print_libRadtran_err - true or false that tells the function to
% write and save the libRadtran error message file

% By Andrew John Buggee

%%

function [inputs, spec_response] = create_uvSpec_DISORT_inputs_for_HySICS(inputs, load_parameters_from_measurement, sim_meas,...
    sim_meas_likeness, print_libRadtran_err)


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


if strcmp(sim_meas_likeness, 'exact')==true



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
    if load_parameters_from_measurement==true

        % Load the RTE solver
        inputs.RT.rte_solver = sim_meas.inputs.RT.rte_solver;

    else

        inputs.RT.rte_solver = 'disort';

    end

    % Define the number of streams to use in your radiative transfer model
    if load_parameters_from_measurement==true

        % Load the number of streams
        inputs.RT.num_streams = sim_meas.inputs.RT.num_streams;

    else

        inputs.RT.num_streams = 16;

    end



    % ------------------------------------------------------------------------


    % ------------------------------------------------------------------------
    % --- Do you want to use the Nakajima and Tanka radiance correction? -----
    if load_parameters_from_measurement==true

        % Load nakajima phase correction
        if isfield(sim_meas.inputs.RT, 'use_nakajima_phaseCorrection')==true

            inputs.RT.use_nakajima_phaseCorrection = sim_meas.inputs.RT.use_nakajima_phaseCorrection;
        else
            inputs.RT.use_nakajima_phaseCorrection = false;
        end

    else

        inputs.RT.use_nakajima_phaseCorrection = true;

    end

    % ------------------------------------------------------------------------






    % ------------------------------------------------------------------------
    % -------------- Define the source file and resolution -------------------
    if load_parameters_from_measurement==true

        % Load source file and resolution used
        inputs.RT.source_file = sim_meas.inputs.RT.source_file;
        inputs.RT.source_file_resolution = sim_meas.inputs.RT.source_file_resolution;

    else

        % source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
        % source_file_resolution = 0.025;         % nm

        % these data have 0.1nm sampling resolution, despite what the file name
        % suggests
        inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
        inputs.RT.source_file_resolution = 0.1;         % nm

        % these data have 1nm sampling resolution
        % inputs.RT.source_file = 'kurudz_1.0nm.dat';
        % inputs.RT.source_file_resolution = 1;         % nm

    end


    % ------------------------------------------------------------------------





    % ------------------------------------------------------------------------
    % ---------------------- DEFINE THE WAVELENGTHS! -------------------------
    % ------------------------------------------------------------------------
    % define the wavelength range. If monochromatic, enter the same number
    % twice

    if load_parameters_from_measurement==true

        % Load the simulated bands
        inputs.bands2run = sim_meas.inputs.bands2run;

    else

        % ----------------- Simulating HySICS spectral channels ------------------
        % number of channels = 636 ranging from center wavelengths: [351, 2297]
        % *** Run All Bands ***
        % inputs.bands2run = (1:1:636)';

        % First 7 MODIS spectral channels
        % inputs.bands2run = [39, 67, 98, 169, 291, 421, 581];

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
        inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 180, 188,...
            198, 217, 220, 222, 224, 227, 237, 245, 249, 254, 264, 288, 290, 293,...
            346, 351, 354, 360, 365, 367, 372, 379, 388, 390, 393, 426, 434, 436,...
            462, 468, 469, 520, 524, 525, 526, 527, 530, 531, 533, 535, 537, 539,...
            543, 547, 570, 574, 577, 579, 582, 613, 616, 618, 620, 623, 625]';



        % inputs.bands2run = [49, 426, 613]';
        %         inputs.bands2run = [49, 57, 288, 408, 613]';

        % test bands
        % 500 nm
        % inputs.bands2run = 49;

        % 1598 nm
        %         inputs.bands2run = 408;

        % 1652 nm
        % inputs.bands2run = 426;

        % 2122 nm
        %         inputs.bands2run = 580;

        % 2236 nm
        % inputs.bands2run = 613;

    end

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

    if load_parameters_from_measurement==true

        % store the spectral response
        spec_response = sim_meas.spec_response;

        % store the wavelengths to run
        inputs.RT.wavelengths2run = sim_meas.inputs.RT.wavelengths2run;


    else

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


    end


    % ------------------------------------------------------------------------
    % ------------------------------------------------------------------------





    % ------------------------------------------------------------------------
    % ----------------- What band model do you want to use? ------------------

    % reptran coarse is the default
    % if using reptran, provide one of the following: coarse (default), medium
    % or fine
    if load_parameters_from_measurement==true

        % Load the band parameterization used
        inputs.RT.band_parameterization = sim_meas.inputs.RT.band_parameterization;

    else

        inputs.RT.band_parameterization = 'reptran coarse';

    end
    % ------------------------------------------------------------------------



    % ------------------------------------------------------------------------
    % define the atmospheric data file
    if load_parameters_from_measurement==true

        % Load the atm file used to generate the measurements
        inputs.RT.atm_file = sim_meas.inputs.RT.atm_file;

    else

        inputs.RT.atm_file = 'afglus.dat';

    end


    % define the surface albedo
    if load_parameters_from_measurement==true

        % Load the surface albedo used for the simulated measurements
        inputs.RT.surface_albedo = sim_meas.inputs.RT.surface_albedo;

    else

        inputs.RT.surface_albedo = 0.04;            % Ocean water albedo
        % inputs.RT.surface_albedo = 0;             % Use a value of 0 when creating weighting functions

    end



    % ----- Define the day of the year to account for Earth-Sun distance -----
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'day_of_year')
            % day of the year
            inputs.RT.day_of_year = sim_meas.inputs.RT.day_of_year;       % value for pixel used in Figure 3.a from paper 1

        end

    else

        % day of the year
        inputs.RT.day_of_year = 316;       % value for pixel used in Figure 3.a from paper 1

    end
    % ------------------------------------------------------------------------




    % ------------------------------------------------------------------------
    % -------------- Do you want a cloud in your model? ----------------------
    if load_parameters_from_measurement==true

        % Load
        inputs.RT.yesCloud = sim_meas.inputs.RT.yesCloud;

    else

        inputs.RT.yesCloud = true;

    end



    % define the cloud geometric depth
    if load_parameters_from_measurement==true

        inputs.RT.cloud_depth = sim_meas.inputs.RT.cloud_depth;                % meters

    else

        inputs.RT.cloud_depth = 500;                % meters

    end

    % define the geometric location of the cloud top and cloud bottom
    if load_parameters_from_measurement==true

        inputs.RT.z_topBottom = sim_meas.inputs.RT.z_topBottom;          % km above surface  - value for pixel used in Figure 3.a from paper 1


    else

        % inputs.RT.z_topBottom = [1.5, 1];          % km above surface
        inputs.RT.z_topBottom = [1.25, 0.75];          % km above surface  - value for pixel used in Figure 3.a from paper 1

    end

    % Water Cloud depth
    inputs.RT.H = inputs.RT.z_topBottom(1) - inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud

    % Do you want to manually set the optical depth?
    inputs.RT.modify_wc_opticalDepth = false;

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

        inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used


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
        if load_parameters_from_measurement==true

            inputs.RT.distribution_var = sim_meas.inputs.RT.distribution_var;              % distribution variance

        else

            inputs.RT.distribution_var = linspace(10,10, inputs.RT.n_layers);              % distribution variance

        end

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
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'define_atm_grid')
            % Load the measure atm grid inputs
            inputs.RT.define_atm_grid = sim_meas.inputs.RT.define_atm_grid;

        end

        if isfield(sim_meas.inputs.RT, 'define_atm_grid')
            % Set the vertical grid to include the cloud layers
            inputs.RT.atm_z_grid = sim_meas.inputs.RT.define_atm_grid;  % km

        end

    else

        inputs.RT.define_atm_grid=false;

        % Set the vertical grid to include the cloud layers
        inputs.RT.atm_z_grid = [0:0.5:inputs.RT.z_topBottom(2), inputs.RT.z_edges(2:end),...
            inputs.RT.z_topBottom(1)+1:1:20, 22:2:30, 35:5:50];   % km

    end




    % --------------------------------------------------------------
    % ----------- Define the Solar and Viewing Gemometry -----------
    % --------------------------------------------------------------
    if load_parameters_from_measurement==true

        % Load the defined sensor altitude for the simualte measurements
        inputs.RT.sensor_altitude = sim_meas.inputs.RT.sensor_altitude;

    else

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


    end




    % -----------------------------
    % define the solar zenith angle
    % -----------------------------
    if load_parameters_from_measurement==true

        inputs.RT.sza = sim_meas.inputs.RT.sza;      % degrees

    else

        inputs.RT.sza = 31;               % degree - value for pixel used in Figure 3.a from paper 1
        % inputs.RT.sza = acosd(0.65);           % degree - for Platnick (2000)
        % inputs.RT.sza = 0;         % degree

    end


    % -------------------------------
    % define the solar azimuith angle
    % -------------------------------
    if load_parameters_from_measurement==true

        inputs.RT.phi0 = sim_meas.inputs.RT.phi0;      % degrees

    else

        % Define the solar azimuth measurement between values 0 and 360
        % The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
        % due north. The libRadTran solar azimuth is defined as 0-360 degrees
        % clockwise from due south. So they are separated by 180 degrees. To map
        % the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
        % 360
        %inputs.RT.phi0 = mod(293.8140 + 180, 360);
        inputs.RT.phi0 = -84 + 180;         % degree - value for pixel used in Figure 3.a from paper 1
        % inputs.RT.phi0 = 0;         % degree

    end



    % --------------------------------
    % define the viewing azimuth angle
    % --------------------------------
    if load_parameters_from_measurement==true

        inputs.RT.vza = sim_meas.inputs.RT.vza;         % degrees

    else

        inputs.RT.vza = 4.29;                                   % degree - value for pixel used in Figure 3.a from paper 1
        % inputs.RT.vza = acosd(0.75);                              % degree - for Platnick (2000)
        % inputs.RT.vza = 20;                             % values are in degrees;

    end


    % --------------------------------
    % define the viewing azimuth angle
    % --------------------------------
    if load_parameters_from_measurement==true

        inputs.RT.vaz = sim_meas.inputs.RT.vaz;      % degrees

    else

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

    end

    % --------------------------------------------------------------



    % --------------------------------------------------------------
    % --- Do you want to use the Cox-Munk Ocean Surface Model? -----

    if load_parameters_from_measurement==true

        inputs.RT.use_coxMunk = sim_meas.inputs.RT.use_coxMunk;
        inputs.RT.wind_speed = sim_meas.inputs.RT.wind_speed;             % m/s

    else

        inputs.RT.use_coxMunk = true;
        inputs.RT.wind_speed = 3;             % m/s

    end
    % --------------------------------------------------------------


    % ------------------------------------------------------------------------
    % --------- specify various cross-section models for  ------------------
    if load_parameters_from_measurement==true

        inputs.RT.specify_cross_section_model = sim_meas.inputs.RT.specify_cross_section_model;

        inputs.RT.crs_model_rayleigh = sim_meas.inputs.RT.crs_model_rayleigh;             %  Rayleigh scattering cross section using Bodhaine et al. (1999) equation 29
        % inputs.RT.crs_model_rayleigh = 'Bodhaine';                   %  Rayleigh scattering cross section using Bodhaine et al. (1999) equations 22-23


    else

        inputs.RT.specify_cross_section_model = true;

        inputs.RT.crs_model_rayleigh = 'Bodhaine29';               %  Rayleigh scattering cross section using Bodhaine et al. (1999) equation 29
        % inputs.RT.crs_model_rayleigh = 'Bodhaine';                   %  Rayleigh scattering cross section using Bodhaine et al. (1999) equations 22-23

    end
    % ------------------------------------------------------------------------


    % ------------------------------------------------------------------------
    % --------- Do you want boundary layer aerosols in your model? -----------
    if load_parameters_from_measurement==true

        inputs.RT.yesAerosols = sim_meas.inputs.RT.yesAerosols;

        inputs.RT.aerosol_type = sim_meas.inputs.RT.aerosol_type;               % 4 = maritime aerosols
        inputs.RT.aerosol_opticalDepth = sim_meas.inputs.RT.aerosol_opticalDepth;     % MODIS algorithm always set to 0.1

    else


        inputs.RT.yesAerosols = true;

        inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
        inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1

    end
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

    if load_parameters_from_measurement==true

        inputs.RT.modify_O3 = sim_meas.inputs.RT.modify_O3;

        inputs.RT.O3_mixing_ratio = sim_meas.inputs.RT.O3_mixing_ratio;       % ppm

    else

        % 400 ppm = 1.0019 * 10^23 molecules/cm^2
        inputs.RT.modify_O3 = false;

        inputs.RT.O3_mixing_ratio = 0;       % ppm

    end
    % ------------------------------------------------------------------------



    % --------------------------------------------------------------
    % ------- Do you want to turn off molecular absorption? --------
    % Note, that thermal emission of molecules is also switched off.
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'no_molecular_abs')
            % Load the setting for molecular absorption
            inputs.RT.no_molecular_abs = sim_meas.inputs.RT.no_molecular_abs;

        else

            inputs.RT.no_molecular_abs = false;

        end

    else

        inputs.RT.no_molecular_abs = false;

    end
    % --------------------------------------------------------------


    % --------------------------------------------------------------
    % ------------ Do you want to turn off scattering? -------------
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'no_scattering_mol')
            % Load the settings for molecular scattering and aerosol scattering
            inputs.RT.no_scattering_mol = sim_meas.inputs.RT.no_scattering_mol;

        end

        if isfield(sim_meas.inputs.RT, 'no_scattering_aer')
            inputs.RT.no_scattering_aer = sim_meas.inputs.RT.no_scattering_aer;
        end

    else

        % Possible choises for the optional argument name are:
        %   mol - Switch off molecular scattering.
        %   aer - Switch off scattering by aerosols.
        %   wc - Switch off scattering by water clouds.
        %   ic - Switch off scattering by ice clouds.
        %   profile - Switch off scattering by any profile defined in profile typename.
        inputs.RT.no_scattering_mol = false;
        inputs.RT.no_scattering_aer = false;

    end
    % --------------------------------------------------------------



    % --------------------------------------------------------------
    % Do you want to print an error message?

    if print_libRadtran_err==true

        inputs.RT.errMsg = 'verbose';

    else

        inputs.RT.errMsg = 'quiet';

    end
    % --------------------------------------------------------------




    % **************************************************************
    % **************************************************************
    % **************************************************************























elseif strcmp(sim_meas_likeness, 'subset')==true


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
    if load_parameters_from_measurement==true

        % Load the RTE solver
        inputs.RT.rte_solver = sim_meas.inputs.RT.rte_solver;

    else

        inputs.RT.rte_solver = 'disort';

    end

    % Define the number of streams to use in your radiative transfer model
    if load_parameters_from_measurement==true

        % Load the number of streams
        inputs.RT.num_streams = sim_meas.inputs.RT.num_streams;

    else

        inputs.RT.num_streams = 16;

    end



    % ------------------------------------------------------------------------


    % ------------------------------------------------------------------------
    % --- Do you want to use the Nakajima and Tanka radiance correction? -----
    if load_parameters_from_measurement==true

        % Load nakajima phase correction
        if isfield(sim_meas.inputs.RT, 'use_nakajima_phaseCorrection')==true

            inputs.RT.use_nakajima_phaseCorrection = sim_meas.inputs.RT.use_nakajima_phaseCorrection;
        else
            inputs.RT.use_nakajima_phaseCorrection = false;
        end

    else

        inputs.RT.use_nakajima_phaseCorrection = true;

    end

    % ------------------------------------------------------------------------






    % ------------------------------------------------------------------------
    % -------------- Define the source file and resolution -------------------
    if load_parameters_from_measurement==true

        % Load source file and resolution used
        inputs.RT.source_file = sim_meas.inputs.RT.source_file;
        inputs.RT.source_file_resolution = sim_meas.inputs.RT.source_file_resolution;

    else

        % source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
        % source_file_resolution = 0.025;         % nm

        % these data have 0.1nm sampling resolution, despite what the file name
        % suggests
        inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
        inputs.RT.source_file_resolution = 0.1;         % nm

        % these data have 1nm sampling resolution
        % inputs.RT.source_file = 'kurudz_1.0nm.dat';
        % inputs.RT.source_file_resolution = 1;         % nm

    end


    % ------------------------------------------------------------------------





    % ------------------------------------------------------------------------
    % ---------------------- DEFINE THE WAVELENGTHS! -------------------------
    % ------------------------------------------------------------------------
    % define the wavelength range. If monochromatic, enter the same number
    % twice

    if load_parameters_from_measurement==true

        % Load the simulated bands
        inputs.bands2run = sim_meas.inputs.bands2run;

    else

        % ----------------- Simulating HySICS spectral channels ------------------
        % number of channels = 636 ranging from center wavelengths: [351, 2297]
        %     inputs.bands2run = (1:1:636)';

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
        inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 180, 188,...
            198, 217, 220, 222, 224, 227, 237, 245, 249, 254, 264, 288, 290, 293,...
            346, 351, 354, 360, 365, 367, 372, 379, 388, 390, 393, 426, 434, 436,...
            462, 468, 469, 520, 524, 525, 526, 527, 530, 531, 533, 535, 537, 539,...
            543, 547, 570, 574, 577, 579, 582, 613, 616, 618, 620, 623, 625]';



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

    end

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

    if load_parameters_from_measurement==true

        % store the spectral response
        spec_response = sim_meas.spec_response;

        % store the wavelengths to run
        inputs.RT.wavelengths2run = sim_meas.inputs.RT.wavelengths2run;


    else

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


    end


    % ------------------------------------------------------------------------
    % ------------------------------------------------------------------------





    % ------------------------------------------------------------------------
    % ----------------- What band model do you want to use? ------------------

    % reptran coarse is the default
    % if using reptran, provide one of the following: coarse (default), medium
    % or fine
    if load_parameters_from_measurement==true

        % Load the band parameterization used
        inputs.RT.band_parameterization = sim_meas.inputs.RT.band_parameterization;

    else

        inputs.RT.band_parameterization = 'reptran coarse';

    end
    % ------------------------------------------------------------------------



    % ------------------------------------------------------------------------
    % define the atmospheric data file
    if load_parameters_from_measurement==true

        % Load the atm file used to generate the measurements
        inputs.RT.atm_file = sim_meas.inputs.RT.atm_file;

    else

        inputs.RT.atm_file = 'afglus.dat';

    end


    % define the surface albedo
    if load_parameters_from_measurement==true

        % Load the surface albedo used for the simulated measurements
        inputs.RT.surface_albedo = sim_meas.inputs.RT.surface_albedo;

    else

        inputs.RT.surface_albedo = 0.04;            % Ocean water albedo
        % inputs.RT.surface_albedo = 0;             % Use a value of 0 when creating weighting functions

    end



    % ----- Define the day of the year to account for Earth-Sun distance -----
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'day_of_year')
            % day of the year
            inputs.RT.day_of_year = sim_meas.inputs.RT.day_of_year;       % value for pixel used in Figure 3.a from paper 1

        end

    else

        % day of the year
        inputs.RT.day_of_year = 316;       % value for pixel used in Figure 3.a from paper 1

    end
    % ------------------------------------------------------------------------




    % ------------------------------------------------------------------------
    % -------------- Do you want a cloud in your model? ----------------------
    if load_parameters_from_measurement==true

        % Load
        inputs.RT.yesCloud = sim_meas.inputs.RT.yesCloud;

    else

        inputs.RT.yesCloud = true;

    end


    % define the cloud geometric depth
    inputs.RT.cloud_depth = 500;                % meters

    % define the geometric location of the cloud top and cloud bottom
    % inputs.RT.z_topBottom = [1.5, 1];          % km above surface
    inputs.RT.z_topBottom = [1.25, 0.75];          % km above surface  - value for pixel used in Figure 3.a from paper 1


    % Water Cloud depth
    inputs.RT.H = inputs.RT.z_topBottom(1) - inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud

    % Do you want to manually set the optical depth?
    inputs.RT.modify_wc_opticalDepth = false;

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


        inputs.RT.profile_type = 'adiabatic'; % type of water droplet profile

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
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'define_atm_grid')
            % Load the measure atm grid inputs
            inputs.RT.define_atm_grid = sim_meas.inputs.RT.define_atm_grid;

        end

        if isfield(sim_meas.inputs.RT, 'define_atm_grid')
            % Set the vertical grid to include the cloud layers
            inputs.RT.atm_z_grid = sim_meas.inputs.RT.define_atm_grid;  % km

        end

    else

        inputs.RT.define_atm_grid=false;

        % Set the vertical grid to include the cloud layers
        inputs.RT.atm_z_grid = [0:0.5:inputs.RT.z_topBottom(2), inputs.RT.z_edges(2:end),...
            inputs.RT.z_topBottom(1)+1:1:20, 22:2:30, 35:5:50];   % km

    end




    % --------------------------------------------------------------
    % ----------- Define the Solar and Viewing Gemometry -----------
    % --------------------------------------------------------------
    if load_parameters_from_measurement==true

        % Load the defined sensor altitude for the simualte measurements
        inputs.RT.sensor_altitude = sim_meas.inputs.RT.sensor_altitude;

    else

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


    end




    % -----------------------------
    % define the solar zenith angle
    % -----------------------------
    if load_parameters_from_measurement==true

        inputs.RT.sza = sim_meas.inputs.RT.sza;      % degrees

    else

        inputs.RT.sza = 31;               % degree - value for pixel used in Figure 3.a from paper 1
        % inputs.RT.sza = acosd(0.65);           % degree - for Platnick (2000)
        % inputs.RT.sza = 10;           % degree

    end


    % -------------------------------
    % define the solar azimuith angle
    % -------------------------------
    if load_parameters_from_measurement==true

        inputs.RT.phi0 = sim_meas.inputs.RT.phi0;      % degrees

    else

        % Define the solar azimuth measurement between values 0 and 360
        % The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
        % due north. The libRadTran solar azimuth is defined as 0-360 degrees
        % clockwise from due south. So they are separated by 180 degrees. To map
        % the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
        % 360
        %inputs.RT.phi0 = mod(293.8140 + 180, 360);
        inputs.RT.phi0 = -84 + 180;         % degree - value for pixel used in Figure 3.a from paper 1
        % inputs.RT.phi0 = 0;         % degree

    end



    % --------------------------------
    % define the viewing azimuth angle
    % --------------------------------
    if load_parameters_from_measurement==true

        inputs.RT.vza = sim_meas.inputs.RT.vza;         % degrees

    else

        inputs.RT.vza = 4.29;                                   % degree - value for pixel used in Figure 3.a from paper 1
        % inputs.RT.vza = acosd(0.75);                              % degree - for Platnick (2000)
        % inputs.RT.vza = 20;                             % values are in degrees;

    end


    % --------------------------------
    % define the viewing azimuth angle
    % --------------------------------
    if load_parameters_from_measurement==true

        inputs.RT.vaz = sim_meas.inputs.RT.vaz;      % degrees

    else

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

    end

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
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'no_molecular_abs')
            % Load the setting for molecular absorption
            inputs.RT.no_molecular_abs = sim_meas.inputs.RT.no_molecular_abs;

        else

            inputs.RT.no_molecular_abs = false;

        end

    else

        inputs.RT.no_molecular_abs = false;

    end
    % --------------------------------------------------------------


    % --------------------------------------------------------------
    % ------------ Do you want to turn off scattering? -------------
    if load_parameters_from_measurement==true

        if isfield(sim_meas.inputs.RT, 'no_scattering_mol')
            % Load the settings for molecular scattering and aerosol scattering
            inputs.RT.no_scattering_mol = sim_meas.inputs.RT.no_scattering_mol;

        end

        if isfield(sim_meas.inputs.RT, 'no_scattering_aer')
            inputs.RT.no_scattering_aer = sim_meas.inputs.RT.no_scattering_aer;
        end

    else

        % Possible choises for the optional argument name are:
        %   mol - Switch off molecular scattering.
        %   aer - Switch off scattering by aerosols.
        %   wc - Switch off scattering by water clouds.
        %   ic - Switch off scattering by ice clouds.
        %   profile - Switch off scattering by any profile defined in profile typename.
        inputs.RT.no_scattering_mol = false;
        inputs.RT.no_scattering_aer = false;

    end
    % --------------------------------------------------------------



    % --------------------------------------------------------------
    % Do you want to print an error message?

    if print_libRadtran_err==true

        inputs.RT.errMsg = 'verbose';

    else

        inputs.RT.errMsg = 'quiet';

    end
    % --------------------------------------------------------------



else


    error([newline, 'The input for the similiarity to the simulated measurements is not recognized', newline])


end






end
