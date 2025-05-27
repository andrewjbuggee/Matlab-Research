%% Compute weighting functions using libRadtran's MYSTIC Monte Carlo Model

% By Andrew John Buggee

clear variables


%% Which computer are you using?


clear inputs

% Determine which computer this is being run on
which_computer = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Monte_Carlo/'];


    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/testing_MYSTIC/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';




elseif strcmp(which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Monte_Carlo/'];


    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/testing_MYSTIC/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];






elseif strcmp(which_computer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = '/scratch/alpine/anbu8374/Thermodynamic_phase/';



    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = '/scratch/alpine/anbu8374/Thermodynamic_phase/';

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_inp, 'dir')

    mkdir(inputs.folderpath_inp)

end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_reflectance, 'dir')

    mkdir(inputs.folderpath_reflectance)

end



%% set up the inputs to create an INP file for MYSTIC!


% Define the parameters of the INP file

% *** compute weighting functions! ***
inputs.compute_weighting_functions = true;


% Are you simulating a measurement, or making forward model calculations
% for the retrieval?
inputs.calc_type = 'monte_carlo';

% Determine which computer this is being run on
inputs.which_computer = which_computer;

% Define the RTE Solver
inputs.RT.rte_solver = 'montecarlo';

% ---------------------------------------------------
% --------- Define Monte Carlo Parameters -----------
% ---------------------------------------------------
inputs.RT.mc.photons = 10000;      % number of photons to use in the simulation
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

% test bands
inputs.bands2run = [49]';


% ------------------------------------------------------------------------
% Do you want to compute radiance/reflectance over a spectral region, or at
% a single wavelength?
% ------------------------------------------------------------------------
inputs.RT.monochromatic_calc = true;


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

for ww = 1:length(inputs.bands2run)
    % The wavelength vector for libRadTran is simply the lower and upper
    % bounds
    inputs.RT.wavelengths2run(ww,:) = [round(mean(spec_response.wavelength(ww,:)), 1),...
        round(mean(spec_response.wavelength(ww, :)) ,1)];

end


% num wavelengths
num_wl = size(inputs.RT.wavelengths2run,1);


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
inputs.RT.surface_albedo = 0;

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

    inputs.RT.r_top = 12;     % microns
    inputs.RT.r_bot = 5;        % microns
    inputs.RT.tau_c = 8;

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


    inputs.RT.n_layers = 10;                          % number of layers to model within cloud

    % -------------------------------------------------------------------
    % define the independent variable used to define the effective radius
    % profile within the cloud
    % -------------------------------------------------------------------
    % if using altitude, z should be a vector starting at the cloud bottom
    % and increasing
    inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used
    inputs.RT.z_edges = linspace(inputs.RT.z_topBottom(2), inputs.RT.z_topBottom(1), inputs.RT.n_layers+1)';   % km - the edges of each layer
    inputs.RT.zT_cloud_indVar = linspace(inputs.RT.z_topBottom(2), inputs.RT.z_topBottom(1), inputs.RT.n_layers)';        % km - altitude above ground vector

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
inputs.RT.sensor_altitude = inputs.RT.z_topBottom(1);      % km - sensor altitude at cloud top


% define the solar zenith angle
% inputs.RT.sza = 31;           % degree
inputs.RT.sza = acosd(0.65);           % degree


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
inputs.RT.vza = acosd(0.85); % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

% inputs.RT.vaz = -103+360;     % degree
inputs.RT.vaz = 0;     % degree

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




% --------------------------------------------------------
% --------- What is column water vapor amount? -----------

% Use a custom H2O profile
inputs.RT.H2O_profile = 'afglus_H2O_none_inside_cloud.dat';


% Using measurements from the AMSR2 instrument, a passive microwave
% radiometer for 17 Jan 2024
inputs.RT.modify_waterVapor = false;

inputs.RT.waterVapor_column = 20;              % mm - milimeters of water condensed in a column
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 416;       % ppm
% ------------------------------------------------------------------------


% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
inputs.RT.compute_reflectivity_uvSpec = true;
% --------------------------------------------------------------



% --------------------------------------------------------------
% Do you want to print an error message?
inputs.RT.errMsg = 'verbose';
% --------------------------------------------------------------




%% Write the INP files

% update the basename
inputs.RT.mc.basename = [inputs.folderpath_inp, 'mc_',num2str(inputs.RT.wavelengths2run(1)),...
    'nm'];


idx = 0;


tic

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true



    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------

    % length of each independent variable
    num_rEff = length(inputs.RT.re);
    num_tauC = length(inputs.RT.tau_c);
    num_wl = length(inputs.bands2run);

    num_INP_files = num_rEff*num_tauC*num_wl;

    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);


    % changing variable steps through reff, tauC, and wavelength
    % in for loop speak, it would be:
    % for rr = 1:num_rEff
    %   for tt = 1:num_tauC
    %       for ww = 1:num_wl
    changing_variables = [reshape(repmat(inputs.RT.re, num_tauC*num_wl,1), [],1),...
        repmat(reshape(repmat(inputs.RT.tau_c, num_wl,1), [],1), num_rEff, 1),...
        repmat(inputs.RT.wavelengths2run, num_rEff*num_tauC, 1)];

    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, repmat((1:num_wl)', num_rEff*num_tauC, 1)];


    % First, write all the wc files
    temp_names = cell(num_rEff*num_tauC, 1);
    wc_filename = cell(num_INP_files, 1);

    % only jump on indexes where there is a unique r and tau pair
    idx_unique_indVars = 1:num_wl:num_INP_files;

    parfor nn = 1:length(idx_unique_indVars)

        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------

        temp = write_wc_file(changing_variables(idx_unique_indVars(nn), 1),...
            changing_variables(idx_unique_indVars(nn), 2),inputs.RT.z_topBottom,...
            inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
            inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
            inputs.which_computer, idx_unique_indVars(nn));

        temp_names{nn} = temp{1};

    end

    % the wc_filenames should be the same for different wavelengths
    for ww = 0:num_wl-1
        wc_filename(idx_unique_indVars+ww) = temp_names;
    end





    % Now write all the INP files
    parfor nn = 1:num_INP_files


        % set the wavelengths for each file
        wavelengths = changing_variables(nn, end-2:end-1);

        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory


        inputFileName{nn} = [num2str(mean(wavelengths)), '_',...
            'nm_rEff_', num2str(changing_variables(nn,1)), '_tauC_', num2str(changing_variables(nn,2)), '_',...
            inputs.RT.atm_file(1:end-4),'.INP'];



        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        % ------------------ Write the INP File --------------------
        write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
            wavelengths, wc_filename{nn});


    end







elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------------------------------------------
    % --------- NON-HOMOGENOUS CLOUD -------------
    % --------------------------------------------



    % length of each independent variable
    % length of each independent variable
    % num_rTop = length(inputs.RT.r_top);
    % num_rBot = length(inputs.RT.r_bot);
    % num_tauC = length(inputs.RT.tau_c);
    num_wl = length(inputs.bands2run);
    num_tau_layers = inputs.RT.n_layers;

    % num_INP_files = num_rTop*num_rBot*num_tauC*num_wl;
    num_INP_files = num_wl*num_tau_layers;


    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);


    % changing variable steps through reff, tauC, and wavelength
    % in for loop speak, it would be:
    % for rt = 1:num_rTop
    %   for rb = 1:num_rBot
    %       for tt = 1:num_tauC
    %           for ww = 1:num_wl
    % changing_variables = [reshape(repmat(inputs.RT.r_top, num_rBot*num_tauC*num_wl,1), [],1),...
    %     repmat(reshape(repmat(inputs.RT.r_bot, num_tauC*num_wl,1), [],1), num_rTop, 1),...
    %     repmat(reshape(repmat(inputs.RT.tau_c, num_wl,1), [],1), num_rBot*num_rTop, 1),...
    %     repmat(inputs.RT.wavelengths2run, num_rTop*num_rBot*num_tauC, 1)];

    changing_variables = [repmat(inputs.RT.wavelengths2run(1), num_tau_layers,1),repmat(inputs.RT.zT_cloud_indVar, num_wl,1)];


    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, repmat((1:num_wl)', num_tau_layers, 1)];

    % First, write all the wc files
    temp_names = cell(num_tau_layers, 1);
    wc_filename = cell(num_INP_files, 1);

    % need an re and tau_c file for each layer, because we slowly build the
    % cloud but adding layer after layer
    idx_unique_indVars = 1:num_tau_layers;


    % create a droplet profile
    re = create_droplet_profile2([inputs.RT.r_top, inputs.RT.r_bot], inputs.RT.zT_cloud_indVar,...
        inputs.RT.indVar, inputs.RT.profile_type);     % microns - effective radius vector

    % -----------------------------------
    % ---- Write a Water Cloud file! ----
    % -----------------------------------
    % re must be defined from cloud bottom to cloud top
    % z_topBottom must be defined as the cloud top height first, then
    % cloud bottom
    [temp, lwc, ext_bulk_coeff_per_LWC] = write_wc_file(re, inputs.RT.tau_c,...
        inputs.RT.z_topBottom, inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
        inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
        inputs.RT.indVar, inputs.compute_weighting_functions, inputs.which_computer, 1);

    temp_names = temp{1};

    % Compute the optical depth of each layer
    inputs.RT.tau_layers = (lwc.*ext_bulk_coeff_per_LWC.*(inputs.RT.z_edges(2) - inputs.RT.z_edges(1)))';  % the optical depth of each layer, starting from cloud top

    %parfor nn = 1:length(idx_unique_indVars)
    for nn = 1:length(idx_unique_indVars)




        temp = write_wc_file(re(end-(nn-1):end), sum(inputs.RT.tau_layers(end-(nn-1):end)),...
            flipud(inputs.RT.z_edges([end-nn,end]))', inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
            inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
            inputs.RT.indVar, inputs.which_computer, idx_unique_indVars(nn));

        temp_names{nn} = temp{1};

    end

    % the wc_filenames should be the same for different wavelengths
    for ww = 0:num_wl-1
        wc_filename(idx_unique_indVars+ww) = temp_names;
    end


    % Now write all the INP files
    % parfor nn = 1:num_INP_files
    for nn = 1:num_INP_files

        % set the wavelengths for each file
        wavelengths = changing_variables(nn);

        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory


        inputFileName{nn} = ['monteCarlo_',num2str(mean(wavelengths)), '_','nm_rTop_', num2str(inputs.RT.r_top),...
            '_rBot_', num2str(inputs.RT.r_bot),'_tauC_', num2str(changing_variables(nn,2)), '_',...
            inputs.RT.atm_file(1:end-4),'.INP'];



        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        % ------------------ Write the INP File --------------------
        write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
            wavelengths, wc_filename{nn});


    end



end

toc



%% Run the MYSTIC program

% define only the spec_response so the wavelengths are passed into the
% memory of the parallel for loop
spec_response_value = spec_response.value;


tic



% store the reflectances
Refl_model = zeros(num_INP_files, 1);


parfor nn = 1:num_INP_files
    % for ww = 1:size(inputs.RT.wavelengths2run, 1)


    disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------


    % compute INP file
    [inputSettings] = runUVSPEC(inputs.folderpath_inp, inputFileName{nn}, outputFileName{nn},...
        inputs.which_computer);

    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC(inputs.folderpath_inp, outputFileName{nn},inputSettings(2,:),...
        inputs.RT.compute_reflectivity_uvSpec);

    % Store the reflectance
    Refl_model(nn, :) = ds.radiance.value;       % radiance is in units of mW/nm/m^2/sr

    % % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    % [Refl_model(nn), ~] = reflectanceFunction(inputSettings(2,:), ds,...
    %     spec_response_value(changing_variables(nn,end),:));



end






toc



%% Compute the weighting functions using Platnick (2000)
% Equation 4 defines the weighting function as the normalized derivative of
% reflectance with respect to optical depth.


%% test code

z = [];
tau_vector = [];
R = [];

wm = diff(R)./diff(tau_vector) ./ R(1);