%% Testing the retrieval of column water vapor from SWIR hyperspectral measurements
% See how radiance estimates change when the water vapor profile changes


% By Andrew John Buggee


clear variables


%% Define the path location where INP files will be stored, and where Reflectances will be stored

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/SWIR_water Vapor_retrieval/column_water_vapor_retrieval/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/column_water_vapor_retrieval/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';




elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'SWIR_water Vapor_retrieval/column_water_vapor_retrieval/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/column_water_vapor_retrieval/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];






elseif strcmp(whatComputer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = '/scratch/alpine/anbu8374/Thermodynamic_phase/';



    % Define the folder path where all .INP files will be saved
    folderpath_inp = '/scratch/alpine/anbu8374/Thermodynamic_phase/';

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


end


% If the folder path doesn't exit, create a new directory
if ~exist(folderpath_inp, 'dir')

    mkdir(folderpath_inp)

end


% If the folder path doesn't exit, create a new directory
if ~exist(folderpath_reflectance, 'dir')

    mkdir(folderpath_reflectance)

end



%% ---- First, let's simulate water clouds ----


% Define the parameters of the INP file
clear inputs


% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------- Define the source file and resolution -------------------

% source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
% source_file_resolution = 0.025;         % nm

% these data have 0.1nm sampling resolution, despite what the file name
% suggests
% inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
% inputs.RT.source_file_resolution = 0.1;         % nm

% these data have 0.1nm sampling resolution, despite what the file name
% suggests
inputs.RT.source_file = 'kurudz_1.0nm.dat';
inputs.RT.source_file_resolution = 1;         % nm

% ------------------------------------------------------------------------


% ------- Define the instrument being modeled ---------
inputs.RT.instrument = 'HySICS';
% inputs.RT.instrument = 'EMIT';
% inputs.RT.instrument = 'arbitrary';


% ------------------------------------------------------------------------
% ---------------------- DEFINE THE WAVELENGTHS! -------------------------
% ------------------------------------------------------------------------
% define the wavelength range. If monochromatic, enter the same number
% twice

if strcmp(inputs.RT.instrument, 'arbitrary')==true

    % --- Compute reflectance from from 350 to 2500 nanometers ---
    inputs.RT.wavelength = 350:2500;


elseif strcmp(inputs.RT.instrument, 'EMIT')==true
    % ----------------- Simulating EMIT spectral channels ------------------

    % define the wavelength channels that cover the range between 1615 and 1730
    % microns
    % inputs.bands2run = find(emit.radiance.wavelength>=1615 & emit.radiance.wavelength<=1730)';


    % --- New New New New New indexs - using HiTran - avoid water vapor and other absorbing gasses! With Pilewskie input ---
    % libRadtran estimates of reflectance below 500 nm consistently
    % overestimate the measured values from EMIT. Let's ignore wavelengths
    % below 500
    inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 86, 87, 88, 89, 90,...
        94, 115, 116, 117, 156, 157, 158, 172, 175, 176,...
        231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';


elseif strcmp(inputs.RT.instrument, 'HySICS')==true
    % ----------------- Simulating HySICS spectral channels ------------------
    % number of channels = 460
    inputs.bands2run = (1:1:460)';

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Create the spectral response functions ------------------
% ------------------------------------------------------------------------

if strcmp(inputs.RT.instrument, 'arbitrary')==true
    % ---------------------------------------
    % If you're not modeling an instrument...
    % ---------------------------------------
    spec_response = ones(1, length(inputs.RT.wavelength));

elseif strcmp(inputs.RT.instrument, 'EMIT')==true
    % ------------------------------------
    % If modeling the EMIT instrument...
    % ------------------------------------
    % Define the EMIT spectral response functions
    spec_response = create_EMIT_specResponse(emit, inputs);
    % keep only the response functions for the wavelengths we care about
    spec_response.value = spec_response.value(inputs.bands2run, :);
    spec_response.wavelength = spec_response.wavelength(inputs.bands2run, :);

    % now define the wavelength range of each spectral channel
    inputs.RT.wavelength = zeros(length(inputs.bands2run), 2);

    for ww = 1:length(inputs.bands2run)

        % The wavelength vector for libRadTran is simply the lower and upper
        % bounds
        inputs.RT.wavelength(ww,:) = [spec_response.wavelength(ww, 1),...
            spec_response.wavelength(ww, end)];

    end

elseif strcmp(inputs.RT.instrument, 'HySICS')==true

    % ------------------------------------
    % if modeling the HySICS instrument...
    % ------------------------------------
    % Define the HySICS spectral response functions
    spec_response = create_HySICS_specResponse(inputs.bands2run, inputs.RT.source_file_resolution);

    % now define the wavelength range of each spectral channel
    inputs.RT.wavelength = zeros(length(inputs.bands2run), 2);

    for ww = 1:length(inputs.bands2run)
        % The wavelength vector for libRadTran is simply the lower and upper
        % bounds
        inputs.RT.wavelength(ww,:) = [spec_response{ww}(1, 1),...
            spec_response{ww}(end, 1)];

    end

end


% ------------------------------------------------------------------------
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
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% define the atmospheric data file
inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
inputs.RT.albedo = 0.04;

% day of the year
inputs.RT.day_of_year = 100;
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
    inputs.RT.dist_var = 7;              % distribution variance

    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    inputs.RT.re = 10.79;      % microns
    inputs.RT.tau_c = 7;


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

    inputs.RT.profile_type = 'adiabatic'; % type of water droplet profile

    inputs.RT.n_layers = 10;                          % number of layers to model within cloud

    inputs.RT.z = linspace(inputs.RT.z_topBottom(1), inputs.RT.z_topBottom(2), inputs.RT.n_layers);        % km - altitude above ground vector

    inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used

    inputs.RT.dist_var = linspace(10,10, inputs.RT.n_layers);              % distribution variance

    inputs.RT.r_top = 12.565;     % microns
    inputs.RT.r_bot = 4.135;        % microns
    inputs.RT.tau_c = 6.424;

    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    inputs.RT.dist_var = 7;              % distribution variance


end

% ------------------------------------------------------------------------



% Define the parameterization scheme used to comptue the optical quantities
% within the INP file, i.e. what function call that tells libRadtran how to
% compute scattering and optical quantities
if inputs.RT.use_custom_mie_calcs==false

    inputs.RT.wc_parameterization = 'mie interpolate';
    %inputs.RT.wc_parameterization = 'hu';

else
    %wc_parameterization = '../data/wc/mie/wc.mie_test.cdf interpolate';
    inputs.RT.wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
end

% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% define the solar zenith angle
inputs.RT.sza = 19.5688;           % degree

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
inputs.RT.phi0 = 113.8140;         % degree

% define the viewing zenith angle
inputs.RT.vza = 8.3134; % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

inputs.RT.vaz = 70.0849;     % degree
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

inputs.RT.CO2_mixing_ratio = 410;       % ppm
% ------------------------------------------------------------------------


% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
inputs.RT.compute_reflectivity_uvSpec = false;
% --------------------------------------------------------------



% --------------------------------------------------------------
% Do you want to print an error message?
inputs.RT.errMsg = 'verbose';
% --------------------------------------------------------------



%% Write each INP file and Calculate Reflectance



% num wavelengths
if size(inputs.RT.wavelength,1)>1 && size(inputs.RT.wavelength,2)>1

    wl_flag_range = true;
    num_wl = size(inputs.RT.wavelength,1);

else

    wl_flag_range = false;
    num_wl = length(inputs.RT.wavelength);
end


idx = 0;


tic


if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------


    % create a legend string
    lgnd_str = cell(1, length(inputs.RT.waterVapor_column));

    % store the reflectances
    Refl_model = zeros(num_wl, length(inputs.RT.re), length(inputs.RT.tau_c));


    for rr = 1:length(inputs.RT.re)



        for tc = 1:length(inputs.RT.tau_c)


            idx = idx + 1;
            % -----------------------------------
            % ---- Write a Water Cloud file! ----
            % -----------------------------------

            % ------------------------------------------------------
            % --------------------VERY IMPORTANT ------------------
            % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
            % ------------------------------------------------------
            wc_filename = write_wc_file(inputs.RT.re(rr), inputs.RT.tau_c(tc), inputs.RT.z_topBottom,...
                inputs.RT.lambda_forTau, inputs.RT.distribution_str, inputs.RT.dist_var,...
                inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str, idx);

            wc_filename = wc_filename{1};


            disp(['Iteration: [re, tc] = [', num2str(rr), '/', num2str(length(inputs.RT.re)),', ',...
                num2str(tc), '/', num2str(length(inputs.RT.tau_c)), ']', newline])


            parfor ww = 1:size(inputs.RT.wavelength, 1)
            % for ww = 1:size(inputs.RT.wavelength, 1)


                disp(['ww = ', num2str(ww),'/',num2str(num_wl), newline])


                % ------------------------------------------------
                % ---- Define the input and output filenames! ----
                % ------------------------------------------------
                % input_names need a unique identifier. Let's give them the nn value so
                % they can be traced, and are writen over in memory

                inputName = [num2str(inputs.RT.wavelength(ww, 1)), '-', num2str(inputs.RT.wavelength(ww, 2)), 'nm_',...
                    inputs.RT.atm_file(1:end-4),'.INP'];



                outputName = ['OUTPUT_',inputName(1:end-4)];



                % ----------------- ******************** ---------------------
                % ------------------ Write the INP File --------------------
                % ----------------- ******************** ---------------------

                % Open the old file for writing
                fileID = fopen([folderpath_inp,inputName], 'w');

                % Define which RTE solver to use
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n';
                fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');


                % Define the number of streams to keep track of when solving the equation
                % of radiative transfer
                % ------------------------------------------------
                formatSpec = '%s %u %5s %s \n\n';
                fprintf(fileID, formatSpec,'number_of_streams', inputs.RT.num_streams,' ', '# Number of streams');


                % Use phase function correction?
                % ------------------------------------------------
                if inputs.RT.use_nakajima_phaseCorrection==true
                    % define the pahse correction to be true
                    % ------------------------------------------------
                    formatSpec = '%s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'disort_intcor phase', ' ', '# Apply the Nakajima and Tanka radiance correction');
                end

                % Define location of the data files
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'data_files_path', libRadtran_data_path, ' ', '# Location of libRadtran data files');


                % Define the band model to use
                % of radiative transfer
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'mol_abs_param', inputs.RT.band_parameterization,' ', '# Band parameterization');


                % Define the location and filename of the atmopsheric profile to use
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'atmosphere_file ', [libRadtran_data_path, 'atmmod/', inputs.RT.atm_file],...
                    ' ', '# Location of atmospheric profile');

                % Define the location and filename of the extraterrestrial solar source
                % ---------------------------------------------------------------------
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'source solar', inputs.RT.source_file, ' ', '# Bounds between 250 and 10000 nm');



                % Define the surface albedo
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'albedo', inputs.RT.albedo, ' ', '# Surface albedo of the ocean');


                % Define the Water Cloud properties, if you want a cloud in your model
                % --------------------------------------------------------------------
                if inputs.RT.yesCloud==true

                    % Define the water cloud file
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n';
                    fprintf(fileID, formatSpec,'wc_file 1D', [libRadtran_data_path,'wc/',wc_filename], ' ', '# Location of water cloud file');
                    %fprintf(fileID, formatSpec,'wc_file 1D', [libRadtran_data_path,'wc/', wc_filename{rr,tc}{1}], ' ', '# Location of water cloud file');



                    % Define the technique or parameterization used to convert liquid cloud
                    % properties of r_eff and LWC to optical depth
                    % ----------------------------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'wc_properties', inputs.RT.wc_parameterization, ' ', '# optical properties parameterization technique');

                end


                if wl_flag_range==true

                    % Define the wavelengths for which the equation of radiative transfer will
                    % be solve
                    % -------------------------------------------------------------------------
                    formatSpec = '%s %f %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'wavelength', inputs.RT.wavelength(ww, 1), inputs.RT.wavelength(ww, 2), ' ', '# Wavelength range');

                else

                    % Define the wavelengths for which the equation of radiative transfer will
                    % be solve
                    % -------------------------------------------------------------------------
                    formatSpec = '%s %f %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'wavelength', inputs.RT.wavelength(ww), inputs.RT.wavelength(ww), ' ', '# Wavelength range');

                end

                % Spline interpolate the calculated spectrum between wavelengths lambda_0
                % and lambda_1 in steps of lambda_step, in nm.
                % -------------------------------------------------------------------------
                % formatSpec = '%s %f %f %5s %s \n\n';
                % fprintf(fileID, formatSpec,'spline', inputs.RT.wavelength(ww, 1), inputs.RT.wavelength(ww, 2), ' ', '# Wavelength range');




                if inputs.RT.use_coxMunk==true

                    % Define the wind speed for the Cox-Munk ocean surface bi-directional reflectance model
                    % be solve
                    % -------------------------------------------------------------------------
                    formatSpec = '%s %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'brdf_cam u10', inputs.RT.wind_speed, ' ', '# (m/s) Ocean Surface wind speed');

                end



                % Define the column water vapor amount
                % --------------------------------------------------------------------
                if inputs.RT.modify_waterVapor==true

                    % If true, modify the amount of column water vapor
                    % --------------------------------------------------------------
                    formatSpec = '%s %f %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'mol_modify H2O ', inputs.RT.waterVapor_column, ' MM', ' ', '# Column water vapor amount');


                end


                % Define the concentration of carbon dioxide
                % --------------------------------------------------------------------
                if inputs.RT.modify_CO2==true

                    % If true, modify the mixing ratio of carbon dioxide
                    % --------------------------------------------------------------
                    formatSpec = '%s %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'mixing_ratio CO2 ', inputs.RT.CO2_mixing_ratio, ' ', '# ppm of CO2');


                end



                % Define the Aerosol Layer properties, if you want a cloud in your model
                % --------------------------------------------------------------------
                if inputs.RT.yesAerosols==true

                    % Turn on default aersol layer, which occupies lower 2km of model
                    % --------------------------------------------------------------
                    formatSpec = '%s %5s %s \n';
                    fprintf(fileID, formatSpec,'aerosol_default', ' ', '# turn on Shettle (1989) boundary layer aerosols');


                    % Specify the Aerosl type
                    % 1=rural aersols,  4=maritime aersols,  5=Urban aerosols,
                    % 6=Tropospheric aerosols
                    % ------------------------------------------------
                    formatSpec = '%s %u %5s %s \n';
                    fprintf(fileID, formatSpec,'aerosol_haze', inputs.RT.aerosol_type, ' ', '# Aerosol type');


                    % Define aerosol layer optical depth
                    % ----------------------------------------------------------------------
                    formatSpec = '%s %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'aerosol_modify tau set', inputs.RT.aerosol_opticalDepth, ' ',...
                        '# Optical Depth of aerosol layer');

                end




                % Define the sensor altitude
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n';
                fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');

                % Define the solar zenith angle
                % ------------------------------------------------
                formatSpec = '%s %f %5s %s \n';
                fprintf(fileID, formatSpec,'sza', inputs.RT.sza, ' ', '# Solar zenith angle');

                % Define the solar azimuth angle
                % -------------------------------------------------------
                formatSpec = '%s %f %5s %s \n';
                fprintf(fileID, formatSpec,'phi0', inputs.RT.phi0, ' ', '# Solar azimuth angle');

                % Define the cosine of the zenith viewing angle
                % ------------------------------------------------
                formatSpec = '%s %f %5s %s \n';
                fprintf(fileID, formatSpec,'umu', round(cosd(inputs.RT.vza),4), ' ', '# Cosine of the zenith viewing angle');

                % Define the azimuth viewing angle
                % ------------------------------------------------
                formatSpec = '%s %f %5s %s \n\n';
                fprintf(fileID, formatSpec,'phi', inputs.RT.vaz, ' ', '# Azimuthal viewing angle');



                if inputs.RT.compute_reflectivity_uvSpec==true
                    % Set the output quantity to be reflectivity
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'output_quantity', 'reflectivity', ' ', '# Output is reflectance');
                end


                %     % Set the outputs
                %     % ------------------------------------------------
                %     formatSpec = '%s %s %5s %s \n\n';
                %     fprintf(fileID, formatSpec,'output_user', 'lambda edir edn eup uavgdir uavgdn uavgup uu', ' ', '# Output quantities');





                % Set the error message to quiet of verbose
                % ------------------------------------------------
                formatSpec = '%s';
                fprintf(fileID, formatSpec, inputs.RT.errMsg);


                % Close the file!
                fclose(fileID);
                % ----------------------------------------------------
                % ----------------------------------------------------




                % ----------------------------------------------------
                % --------------- RUN RADIATIVE TRANSFER -------------
                % ----------------------------------------------------


                % compute INP file
                [inputSettings] = runUVSPEC(folderpath_inp,inputName, outputName);

                % read .OUT file
                % radiance is in units of mW/nm/m^2/sr
                [ds,~,~] = readUVSPEC(folderpath_inp,outputName,inputSettings(2,:),...
                    inputs.RT.compute_reflectivity_uvSpec);

                % Store the Radiance
                %            Rad_model(rr, tc, ww, :) = ds.radiance.value;       % radiance is in units of mW/nm/m^2/sr

                % compute the reflectance
                % Refl_model(ww, rr, tc) = reflectanceFunction_4EMIT(inputSettings(2,:), ds,...
                %     spec_response.value(ww, :)');

                [Refl_model(ww, rr, tc), ~] = reflectanceFunction(inputSettings(2,:), ds, spec_response{ww}(:,2));








            end



        end

    end





elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------------------------------------------
    % --------- NON-HOMOGENOUS CLOUD -------------
    % --------------------------------------------

    % create a legend string
    lgnd_str = cell(1, length(inputs.RT.waterVapor_column));

    % store the reflectances
    Refl_model = zeros(num_wl, length(inputs.RT.r_top), length(inputs.RT.r_bot),...
        length(inputs.RT.tau_c));


    for rt = 1:length(inputs.RT.r_top)


        for rb = 1:length(inputs.RT.r_bot)


            for tc = 1:length(inputs.RT.tau_c)


                disp(['Iteration: [rt, rb, tc] = [', [num2str(rt),', ', num2str(rb), ', ', num2str(tc)], ']...', newline])


                idx = idx + 1;
                % -----------------------------------
                % ---- Write a Water Cloud file! ----
                % -----------------------------------

                % ------------------------------------------------------
                % --------------------VERY IMPORTANT ------------------
                % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
                % ------------------------------------------------------

                % -----------------------------------
                % ---- Write a Water Cloud file! ----
                % -----------------------------------
                % most uncertainties for the modis optical retrieval are between 2
                % and 10 percent. So lets round off all re values to the 1000th decimal
                % place

                re = create_droplet_profile2([inputs.RT.r_top(rt), inputs.RT.r_bot(rb)],...
                    inputs.RT.z, inputs.RT.indVar, inputs.RT.profile_type);     % microns - effective radius vector


                % ------------------------------------------------------
                % --------------------VERY IMPORTANT ------------------
                % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
                % ------------------------------------------------------
                wc_filename = write_wc_file(re, inputs.RT.tau_c(tc), inputs.RT.z_topBottom,...
                    inputs.RT.lambda_forTau, inputs.RT.distribution_str, inputs.RT.dist_var,...
                    inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str, rt*rb);

                wc_filename = wc_filename{1};



                parfor ww = 1:length(inputs.RT.wavelength)
                % for ww = 1:length(inputs.RT.wavelength)



                    disp(['ww = ', num2str(ww),'/',num2str(num_wl), newline])



                    % ------------------------------------------------
                    % ---- Define the input and output filenames! ----
                    % ------------------------------------------------
                    % input_names need a unique identifier. Let's give them the nn value so
                    % they can be traced, and are writen over in memory

                    inputName = [num2str(inputs.RT.wavelength(ww)), 'nm_',...
                        inputs.RT.atm_file(1:end-4),'.INP'];



                    outputName = ['OUTPUT_',inputName(1:end-4)];



                    % ----------------- ******************** ---------------------
                    % ------------------ Write the INP File --------------------
                    % ----------------- ******************** ---------------------

                    % Open the old file for writing
                    fileID = fopen([folderpath_inp,inputName], 'w');

                    % Define which RTE solver to use
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n';
                    fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');


                    % Define the number of streams to keep track of when solving the equation
                    % of radiative transfer
                    % ------------------------------------------------
                    formatSpec = '%s %u %5s %s \n\n';
                    fprintf(fileID, formatSpec,'number_of_streams', inputs.RT.num_streams,' ', '# Number of streams');


                    % Use phase function correction?
                    % ------------------------------------------------
                    if inputs.RT.use_nakajima_phaseCorrection==true
                        % define the pahse correction to be true
                        % ------------------------------------------------
                        formatSpec = '%s %5s %s \n\n';
                        fprintf(fileID, formatSpec,'disort_intcor phase', ' ', '# Apply the Nakajima and Tanka radiance correction');
                    end

                    % Define location of the data files
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'data_files_path', libRadtran_data_path, ' ', '# Location of libRadtran data files');


                    % Define the band model to use
                    % of radiative transfer
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'mol_abs_param', inputs.RT.band_parameterization,' ', '# Band parameterization');


                    % Define the location and filename of the atmopsheric profile to use
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'atmosphere_file ', [libRadtran_data_path, 'atmmod/', inputs.RT.atm_file],...
                        ' ', '# Location of atmospheric profile');

                    % Define the location and filename of the extraterrestrial solar source
                    % ---------------------------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'source solar', inputs.RT.source_file, ' ', '# Bounds between 250 and 10000 nm');



                    % Define the surface albedo
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'albedo', inputs.RT.albedo, ' ', '# Surface albedo of the ocean');


                    % Define the Water Cloud properties, if you want a cloud in your model
                    % --------------------------------------------------------------------
                    if inputs.RT.yesCloud==true

                        % Define the water cloud file
                        % ------------------------------------------------
                        formatSpec = '%s %s %5s %s \n';
                        fprintf(fileID, formatSpec,'wc_file 1D', [libRadtran_data_path,'wc/',wc_filename], ' ', '# Location of water cloud file');
                        %fprintf(fileID, formatSpec,'wc_file 1D', [libRadtran_data_path,'wc/', wc_filename{rr,tc}{1}], ' ', '# Location of water cloud file');



                        % Define the technique or parameterization used to convert liquid cloud
                        % properties of r_eff and LWC to optical depth
                        % ----------------------------------------------------------------------
                        formatSpec = '%s %s %5s %s \n\n';
                        fprintf(fileID, formatSpec,'wc_properties', inputs.RT.wc_parameterization, ' ', '# optical properties parameterization technique');

                    end


                    if wl_flag_range==true

                        % Define the wavelengths for which the equation of radiative transfer will
                        % be solve
                        % -------------------------------------------------------------------------
                        formatSpec = '%s %f %f %5s %s \n\n';
                        fprintf(fileID, formatSpec,'wavelength', inputs.RT.wavelength(ww, 1), inputs.RT.wavelength(ww, 2), ' ', '# Wavelength range');

                    else

                        % Define the wavelengths for which the equation of radiative transfer will
                        % be solve
                        % -------------------------------------------------------------------------
                        formatSpec = '%s %f %f %5s %s \n\n';
                        fprintf(fileID, formatSpec,'wavelength', inputs.RT.wavelength(ww), inputs.RT.wavelength(ww), ' ', '# Wavelength range');

                    end

                    % Spline interpolate the calculated spectrum between wavelengths lambda_0
                    % and lambda_1 in steps of lambda_step, in nm.
                    % -------------------------------------------------------------------------
                    % formatSpec = '%s %f %f %5s %s \n\n';
                    % fprintf(fileID, formatSpec,'spline', inputs.RT.wavelength(ww, 1), inputs.RT.wavelength(ww, 2), ' ', '# Wavelength range');






                    if inputs.RT.use_coxMunk==true

                        % Define the wind speed for the Cox-Munk ocean surface bi-directional reflectance model
                        % be solve
                        % -------------------------------------------------------------------------
                        formatSpec = '%s %f %5s %s \n\n';
                        fprintf(fileID, formatSpec,'brdf_cam u10', inputs.RT.wind_speed, ' ', '# (m/s) Ocean Surface wind speed');

                    end



                    % Define the column water vapor amount
                    % --------------------------------------------------------------------
                    if inputs.RT.modify_waterVapor==true

                        % If true, modify the amount of column water vapor
                        % --------------------------------------------------------------
                        formatSpec = '%s %f %s %5s %s \n\n';
                        fprintf(fileID, formatSpec,'mol_modify H2O ', inputs.RT.waterVapor_column, ' MM', ' ', '# Column water vapor amount');


                    end


                    % Define the concentration of carbon dioxide
                    % --------------------------------------------------------------------
                    if inputs.RT.modify_CO2==true

                        % If true, modify the mixing ratio of carbon dioxide
                        % --------------------------------------------------------------
                        formatSpec = '%s %f %5s %s \n\n';
                        fprintf(fileID, formatSpec,'mixing_ratio CO2 ', inputs.RT.CO2_mixing_ratio, ' ', '# ppm of CO2');


                    end



                    % Define the Aerosol Layer properties, if you want a cloud in your model
                    % --------------------------------------------------------------------
                    if inputs.RT.yesAerosols==true

                        % Turn on default aersol layer, which occupies lower 2km of model
                        % --------------------------------------------------------------
                        formatSpec = '%s %5s %s \n';
                        fprintf(fileID, formatSpec,'aerosol_default', ' ', '# turn on Shettle (1989) boundary layer aerosols');


                        % Specify the Aerosl type
                        % 1=rural aersols,  4=maritime aersols,  5=Urban aerosols,
                        % 6=Tropospheric aerosols
                        % ------------------------------------------------
                        formatSpec = '%s %u %5s %s \n';
                        fprintf(fileID, formatSpec,'aerosol_haze', inputs.RT.aerosol_type, ' ', '# Aerosol type');


                        % Define aerosol layer optical depth
                        % ----------------------------------------------------------------------
                        formatSpec = '%s %f %5s %s \n\n';
                        fprintf(fileID, formatSpec,'aerosol_modify tau set', inputs.RT.aerosol_opticalDepth, ' ',...
                            '# Optical Depth of aerosol layer');

                    end




                    % Define the sensor altitude
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n';
                    fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');

                    % Define the solar zenith angle
                    % ------------------------------------------------
                    formatSpec = '%s %f %5s %s \n';
                    fprintf(fileID, formatSpec,'sza', inputs.RT.sza, ' ', '# Solar zenith angle');

                    % Define the solar azimuth angle
                    % -------------------------------------------------------
                    formatSpec = '%s %f %5s %s \n';
                    fprintf(fileID, formatSpec,'phi0', inputs.RT.phi0, ' ', '# Solar azimuth angle');

                    % Define the cosine of the zenith viewing angle
                    % ------------------------------------------------
                    formatSpec = '%s %f %5s %s \n';
                    fprintf(fileID, formatSpec,'umu', round(cosd(inputs.RT.vza),4), ' ', '# Cosine of the zenith viewing angle');

                    % Define the azimuth viewing angle
                    % ------------------------------------------------
                    formatSpec = '%s %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'phi', inputs.RT.vaz, ' ', '# Azimuthal viewing angle');



                    if inputs.RT.compute_reflectivity_uvSpec==true
                        % Set the output quantity to be reflectivity
                        % ------------------------------------------------
                        formatSpec = '%s %s %5s %s \n\n';
                        fprintf(fileID, formatSpec,'output_quantity', 'reflectivity', ' ', '# Output is reflectance');
                    end


                    %     % Set the outputs
                    %     % ------------------------------------------------
                    %     formatSpec = '%s %s %5s %s \n\n';
                    %     fprintf(fileID, formatSpec,'output_user', 'lambda edir edn eup uavgdir uavgdn uavgup uu', ' ', '# Output quantities');





                    % Set the error message to quiet of verbose
                    % ------------------------------------------------
                    formatSpec = '%s';
                    fprintf(fileID, formatSpec, inputs.RT.errMsg);


                    % Close the file!
                    fclose(fileID);
                    % ----------------------------------------------------
                    % ----------------------------------------------------




                    % ----------------------------------------------------
                    % --------------- RUN RADIATIVE TRANSFER -------------
                    % ----------------------------------------------------


                    % compute INP file
                    [inputSettings] = runUVSPEC(folderpath_inp,inputName, outputName);

                    % read .OUT file
                    % radiance is in units of mW/nm/m^2/sr
                    [ds,~,~] = readUVSPEC(folderpath_inp,outputName,inputSettings(2,:),...
                        inputs.RT.compute_reflectivity_uvSpec);

                    % Store the Radiance
                    %            Rad_model(rr, tc, ww, :) = ds.radiance.value;       % radiance is in units of mW/nm/m^2/sr



                    if strcmp(inputs.RT.instrument, 'HySICS')

                        % compute the reflectance
                        [Refl_model(ww, rt, rb, tc), ~] = reflectanceFunction(inputSettings(2,:), ds, spec_response{ww}(:,2));

                    elseif strcmp(inputs.RT.instrument, 'EMIT')

                        % compute the reflectance
                        Refl_model(ww, rt, rb, tc) = reflectanceFunction_4EMIT(inputSettings(2,:), ds,...
                            spec_response.value(ww, :)');

                    elseif strcmp(inputs.RT.instrument, 'arbitrary')


                    end








                end



            end

        end

    end


end



toc


%% Plot the results!

figure;
if size(inputs.RT.wavelength,1)>1 && size(inputs.RT.wavelength,2)>1

    plot(mean(inputs.RT.wavelength,2), Refl_model, '-', 'linewidth', 5, 'markersize', 27, 'Color', mySavedColors(1, 'fixed'),...
        'linewidth', 3)

else

    plot(inputs.RT.wavelength, Refl_model, '-', 'linewidth', 5, 'markersize', 27, 'Color', mySavedColors(1, 'fixed'),...
        'linewidth', 3)

end

hold on

grid on; grid minor
xlabel('Wavelength (nm)','Interpreter', 'latex')
ylabel('Reflectance (1/sr)','Interpreter', 'latex')
set(gcf, 'Position', [0 0 1000 1000])
legend('Simulated Reflectance',  'Interpreter', 'latex', 'Fontsize', 30', 'location', 'best')
% title(['Simulated Reflectance - liquid water cloud - $r_e = $', num2str(re_2plot), ' $\mu m$, $\tau_c = $',...
%     num2str(tau_2plot)], 'Interpreter', 'latex')





