%% Simulate EMIT Measurements with Gaussian Noise


% By Andrew John Buggee


clear variables

%% Define the path location where INP files will be stored, and where Reflectances will be stored

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/EMIT/Thermodynamic_phase/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Thermodynamic_phase/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';


    % Define the EMIT data folder path
    emitPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/Thermodynamic_phase/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/Thermodynamic_phase/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];


    % Define the EMIT data folder path
    emitPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

elseif strcmp(whatComputer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = '/scratch/alpine/anbu8374/Thermodynamic_phase/';



    % Define the folder path where all .INP files will be saved
    folderpath_inp = '/scratch/alpine/anbu8374/Thermodynamic_phase/';
    % If the folder path doesn't exit, create a new directory
    if ~exist(folderpath_inp, 'dir')

        mkdir(folderpath_inp)

    end

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


    % Define the EMIT data folder path
    emitPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';







end



%% Simulate EMIT Measurements using the EMIT wavelength channels

% Pick any emit data set and pixel. We don't care about the data

% -------------------------------------
% ------- PICK EMIT DATA SET  ---------
% -------------------------------------

emitFolder = '17_Jan_2024_coast/';
%emitFolder = '17_Jan_2024_ocean/';


[emit,L1B_fileName] = retrieveEMIT_data([emitPath, emitFolder]);

% 17_Jan_2024_coast - my TBLUT algorithm found an optical depth of 6.57 and
% an effective radius of 10.79
% Good spectral agreement with:
%   r_e = 10.79
%   tau_c = 7
%   CWV = 40mm
%   CO_2 = 416ppm
pixels2use.row = 932;
pixels2use.col = 960;

% 17_Jan_2024_coast - TBLUT algorithm found optical depths of 6.79 9.22, 11.6, 14.53 19.8
% pixels2use.row = [932, 969, 969, 969, 969];
% pixels2use.col = [960, 989, 984, 980, 957];


% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));

% Remove excess data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% ---- First, let's simulate water clouds ----


% Define the parameters of the INP file
clear inputs

% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% ------- Define the source file and resolution ------

% source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
% source_file_resolution = 0.025;         % nm

% these data have 0.1nm sampling resolution, despite what the file name
% suggests
inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
inputs.RT.source_file_resolution = 0.1;         % nm

% ------------------------------------------------------------------------




% define the wavelength range. If monochromatic, enter the same number
% twice
% ------------------------------------------------------------------------

% --- Use wavelengths from 1.4 to 1.8 microns ---
% The wavelength vector for libRadTran is simply the lower and upper
% bounds
%inputs.RT.wavelength = [1400, 1800];  % nm

% ---------------------------- WAVELENGTHS! -------------------------------
% define the wavelength channels that cover the range between 1615 and 1730
% microns
%inputs.bands2run = find(emit.radiance.wavelength>=1615 & emit.radiance.wavelength<=1730)';


% define the wavelength channels that cover the range between 1615 and 1730
% nm and 2140 to 2260 nm
% inputs.bands2run = [find(emit.radiance.wavelength>=1615 & emit.radiance.wavelength<=1730)',...
%     find(emit.radiance.wavelength>=2140 & emit.radiance.wavelength<=2260)'];


% define the wavelength channels that cover the range between 1600 and 1750
% nm and 2100 to 2300 nm
% inputs.bands2run = [find(emit.radiance.wavelength>=1600 & emit.radiance.wavelength<=1750)',...
%     find(emit.radiance.wavelength>=2100 & emit.radiance.wavelength<=2300)'];


% define the wavelength channels that cover the range between 1000 and 1100, 1600 and 1750
% nm and 2100 to 2300 nm
% inputs.bands2run = [find(emit.radiance.wavelength>=1000 & emit.radiance.wavelength<=1100)',...
%     find(emit.radiance.wavelength>=1600 & emit.radiance.wavelength<=1750)',...
%     find(emit.radiance.wavelength>=2100 & emit.radiance.wavelength<=2300)'];

% Only compute reflectance at 500 nm
%[~, inputs.bands2run] = min(abs(emit.radiance.wavelength - 500));

% Compute all wavelengths above 1000 nm
%inputs.bands2run = find(emit.radiance.wavelength>=1000)';

% Compute all wavelengths below 650 nm
%inputs.bands2run = find(emit.radiance.wavelength<=650)';

% plot all EMIT wavelengths
% inputs.bands2run = find(emit.radiance.wavelength>=300 & emit.radiance.wavelength<=2600)';


% --- New New New New New indexs - using HiTran - avoid water vapor and other absorbing gasses! With Pilewskie input ---
% libRadtran estimates of reflectance below 500 nm consistently
% overestimate the measured values from EMIT. Let's ignore wavelengths
% below 500
inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 86, 87, 88, 89, 90,...
    94, 115, 116, 117, 156, 157, 158, 172, 175, 176,...
    231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';
% ------------------------------------------------------------------------

% Define the EMIT spectral response functions
spec_response = create_EMIT_specResponse(emit, inputs);
% keep only the response functions for the wavelengths we care about
spec_response_2run.value = spec_response.value(inputs.bands2run, :);
spec_response_2run.wavelength = spec_response.wavelength(inputs.bands2run, :);

% now define the wavelength range of each spectral channel
inputs.RT.wavelength = zeros(length(inputs.bands2run), 2);

for ww = 1:length(inputs.bands2run)

    % The wavelength vector for libRadTran is simply the lower and upper
    % bounds
    inputs.RT.wavelength(ww,:) = [spec_response_2run.wavelength(ww, 1),...
        spec_response_2run.wavelength(ww, end)];

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% ----**** Using custom spectral response and wavelength sampling ****----

% inputs.RT.wavelength_center = 350:(emit.radiance.wavelength(2) - emit.radiance.wavelength(1))/3:...
%                               emit.radiance.wavelength(end);     % nm

% inputs.RT.wavelength_center = emit.radiance.wavelength;     % nm
%
% inputs.RT.fwhm = linspace(emit.radiance.fwhm(1), emit.radiance.fwhm(1), length(inputs.RT.wavelength_center));
%
% spec_response_2run = create_gaussian_specResponse(inputs.RT.wavelength_center, inputs.RT.fwhm, inputs);
%
%
% % now define the wavelength range of each spectral channel
% inputs.RT.wavelength = zeros(length(inputs.RT.wavelength_center), 2);
%
% for ww = 1:length(inputs.RT.wavelength_center)
%
%     % The wavelength vector for libRadTran is simply the lower and upper
%     % bounds
%     inputs.RT.wavelength(ww,:) = [spec_response_2run.wavelength(ww, 1),...
%         spec_response_2run.wavelength(ww, end)];
%
% end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ----**** Computing radiance at 0.1 nm resolution ****----
% inputs.RT.wavelength_center = 350:inputs.RT.source_file_resolution:2500;     % nm
%
% % now define each monochromatic calculation as a seperate file
% inputs.RT.wavelength = zeros(length(inputs.RT.wavelength_center), 2);
%
% for ww = 1:length(inputs.RT.wavelength_center)
%
%     % The wavelength vector for libRadTran is simply the lower and upper
%     % bounds
%     inputs.RT.wavelength(ww,:) = [inputs.RT.wavelength_center(ww), inputs.RT.wavelength_center(ww)];
%
% end
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




% define the atmospheric data file
inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
inputs.RT.albedo = 0.05;

% day of the year
inputs.RT.day_of_year = str2double(emitFolder(1:2));

% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
%inputs.RT.linear_cloudFraction = false;
% if false, define the cloud cover percentage
%inputs.RT.percent_cloud_cover = 1;

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

% define the type of droplet distribution
inputs.RT.distribution_str = 'gamma';


% define how liquid water content will be computed in the write_wc_file
% function.
inputs.RT.parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 500;            % nm


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


% define the solar zenith angle
inputs.RT.sza = double(emit.obs.solar.zenith);           % degree

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
inputs.RT.phi0 = mod(double(emit.obs.solar.azimuth + 180), 360);         % degree

% define the viewing zenith angle
inputs.RT.vza = double(emit.obs.sensor.zenith); % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

inputs.RT.vaz = emit.obs.sensor.azimuth;     % degree


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

% default value is 14.295 mm
inputs.RT.waterVapor_column = 40;              % mm - milimeters of water condensed in a column
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 416;       % ppm
% ------------------------------------------------------------------------


% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
inputs.RT.compute_reflectivity_uvSpec = false;
% --------------------------------------------------------------


% --------------------------------------------------------------
% Do you want to print an error message?
inputs.RT.errMsg = 'quiet';
% --------------------------------------------------------------


% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);


% Convert radiance measurements to TOA reflectance for the desired pixels

emit = convert_EMIT_radiance_2_reflectance(emit, inputs);





%% Write each INP file and Calculate Reflectance




idx = 0;




tic


if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % --------- HOMOGENOUS CLOUD -------------
    % store the reflectances
    Refl_model = zeros(length(inputs.RT.re), length(inputs.RT.tau_c), size(inputs.RT.wavelength, 1));

    for rr = 1:length(inputs.RT.re)

        %lgnd_str{rr} = ['$r_e = $', num2str(inputs.RT.re(rr)), ' $\mu m$'];




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


            parfor ww = 1:size(inputs.RT.wavelength, 1)
                % for ww = 1:size(inputs.RT.wavelength, 1)


                disp(['Iteration: [re, tc] = [', num2str(rr), '/', num2str(length(inputs.RT.re)),', ',...
                    num2str(tc), '/', num2str(length(inputs.RT.tau_c)), ', ',...
                    num2str(ww), '/', num2str(size(inputs.RT.wavelength, 1)),']...', newline])






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



                % Define the wavelengths for which the equation of radiative transfer will
                % be solve
                % -------------------------------------------------------------------------
                formatSpec = '%s %f %f %5s %s \n\n';
                fprintf(fileID, formatSpec,'wavelength', inputs.RT.wavelength(ww, 1), inputs.RT.wavelength(ww, 2), ' ', '# Wavelength range');




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
                Refl_model(rr, tc, ww) = reflectanceFunction_4EMIT(inputSettings(2,:), ds,...
                    spec_response_2run.value(ww, :)');








            end



        end

    end





elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------- NON-HOMOGENOUS CLOUD -------------
    % store the Radiances
    Rad_model = zeros(length(inputs.RT.r_top), length(inputs.RT.r_bot), length(inputs.RT.tau_c), size(inputs.RT.wavelength, 1));


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
                    inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str, ww);

                wc_filename = wc_filename{1};



                parfor ww = 1:size(inputs.RT.wavelength, 1)


                    disp(['Iteration: ww = [',num2str(ww), '/', num2str(size(inputs.RT.wavelength, 1)),']...', newline])






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



                    % Define the wavelengths for which the equation of radiative transfer will
                    % be solve
                    % -------------------------------------------------------------------------
                    formatSpec = '%s %f %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'wavelength', inputs.RT.wavelength(ww, 1), inputs.RT.wavelength(ww, 2), ' ', '# Wavelength range');




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
                    Refl_model(rt, rb, tc, ww) = reflectanceFunction_4EMIT(inputSettings(2,:), ds,...
                        spec_response_2run.value(ww, :)');








                end



            end

        end

    end


end



toc



%% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
uncert = 0;

% Define the mean magnitude of the noise
mean_magnitude = emit.reflectance.value(inputs.bands2run) .* uncert;

% define the standard deviation
std = mean_magnitude./3;

% compute gaussian noise for a custom mean and standard deviation
noise = mean_magnitude + std.*randn(length(inputs.bands2run), 1);

Refl_model_with_noise = reshape(Refl_model, [], 1) + noise;



%%
% ----------------------------------------------
% ---------- SAVE REFLECTANCE OUTPUT! ----------
% ----------------------------------------------

rev = 1;




filename = [folderpath_reflectance,'reflectance_calcs_EMIT_water_cloud_sim-ran-on-',char(datetime("today")),...
    '_rev', num2str(rev),'.mat'];

while isfile(filename)
    rev = rev+1;
    filename = [folderpath_reflectance,'reflectance_calcs_EMIT_water_cloud_sim-ran-on-',char(datetime("today")),...
        '_rev', num2str(rev),'.mat'];
end


save(filename,"inputs", "Refl_model", "Refl_model_with_noise");





%%  Plot Reflectance Spectrum


figure;

plot(emit.radiance.wavelength, emit.reflectance.value, '-', 'linewidth', 5, 'markersize', 27, 'Color', mySavedColors(1, 'fixed'),...
    'linewidth', 3)

hold on



if isfield(inputs.RT, 're')

    re_2plot = inputs.RT.re(1); % microns
    tau_2plot = inputs.RT.tau_c(1);


    plot(emit.radiance.wavelength(inputs.bands2run), reshape(Refl_model(inputs.RT.re==re_2plot,...
        inputs.RT.tau_c==tau_2plot, :), [], 1),...
    '.-', 'linewidth', 2, 'markersize', 27, 'Color', mySavedColors(2, 'fixed'))

    title(['Simulated EMIT Reflectance - liquid water cloud - $r_e = $', num2str(re_2plot), ' $\mu m$, $\tau_c = $',...
    num2str(tau_2plot)], 'Interpreter', 'latex')

else

    rTop_2plot = inputs.RT.r_top(1); % microns
    rBot_2plot = inputs.RT.r_bot(1); % microns
    tau_2plot = inputs.RT.tau_c(1);

    plot(emit.radiance.wavelength(inputs.bands2run), reshape(Refl_model(inputs.RT.r_top==rTop_2plot,...
        inputs.RT.r_bot==rBot_2plot, inputs.RT.tau_c==tau_2plot, :), [], 1),...
    '.-', 'linewidth', 2, 'markersize', 27, 'Color', mySavedColors(2, 'fixed'))

    title(['Simulated EMIT Reflectance - liquid water cloud - $r_{top} = $', num2str(rTop_2plot),...
        ' $\mu m$, $r_{bot} = $', num2str(rBot_2plot),' $\mu m$, $\tau_c = $',...
    num2str(tau_2plot)], 'Interpreter', 'latex')

end







hold on




grid on; grid minor
xlabel('Wavelength (nm)','Interpreter', 'latex')
ylabel('Reflectance (1/sr)','Interpreter', 'latex')
set(gcf, 'Position', [0 0 1000 1000])
legend('EMIT Measurements', 'Calculated',  'Interpreter', 'latex', 'Fontsize', 30', 'location', 'best')


