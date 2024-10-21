%% Thermodynamic phase discrimination from EMIT simulated measurements


% By Andrew John Buggee

%% Let's simulate EMIT measurements by varrying the following:

% 1) cloud phase
% 2) cloud particle size
% 3) cloud optical depth
% 4) solar geometry
% 5) viewing geometry


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


    % Define the EMIT data folder path
    emitPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/EMIT/Thermodynamic_phase/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];

    % Define the EMIT data folder path
    emitPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

elseif strcmp(whatComputer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/EMIT/thermodynamic_phase/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = '/scratch/alpine/anbu8374/EMIT_reflectance_uniqueness/';


    % Define the EMIT data folder path
    emitPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

    if ~exist(folderpath_inp, 'dir')

        mkdir(folderpath_inp)
    end





end



%% Simulate EMIT Measurements using the EMIT wavelength channels

% Pick any emit data set and pixel. We don't care about the data

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitFolder = '17_Jan_2024_coast/';
%emitFolder = '17_Jan_2024_ocean/';


[emit,L1B_fileName] = retrieveEMIT_data([emitPath, emitFolder]);

% 17_Jan_2024_coast - my TBLUT algorithm found an optical depth of 6.6
pixels2use.row = 932;
pixels2use.col = 960;


% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));

% Remove excess data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% ---- First, let's simulate water clouds ----

% Define the parameters of the INP file


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

% define the wavelength channels that cover the range between 1550 and 1750
% microns. Leave room for the spectral channel width (~30nm)
inputs.RT.wavelength = [1550 - 15, 1750 + 15];


% create the spectral response functions
inputs.bands2run = find(emit.radiance.wavelength>=1550 & emit.radiance.wavelength<=1750)';

spec_response = create_EMIT_specResponse(emit, inputs);
% keep only the response functions for the wavelengths we care about
spec_response_2run.value = spec_response.value(inputs.bands2run, :);
spec_response_2run.wavelength = spec_response.wavelength(inputs.bands2run, :);

% % compute indices for each spectral channel
% for ww = 1:length(inputs.bands2run)
%
% end
%
% % clear the full spec_response
% clear spec_response




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



% ------------------------------------------------------------------------
inputs.RT.cloud_depth = 1000;                % meters

% define the geometric location of the cloud top and cloud bottom
inputs.RT.z_topBottom = [3, 2];          % km above surface



% Water Cloud depth
inputs.RT.H = inputs.RT.z_topBottom(1) - inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
inputs.RT.linear_cloudFraction = false;
% if false, define the cloud cover percentage
inputs.RT.percent_cloud_cover = 1;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
inputs.RT.use_custom_mie_calcs = false;
% ------------------------------------------------------------------------

% define the type of droplet distribution
inputs.RT.distribution_str = 'gamma';

% define the spread of the droplet distribution
inputs.RT.dist_var = 10;              % distribution variance

% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';

% define how liquid water content will be computed
inputs.RT.parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 1400;            % nm


% ------------------------------------------------------------------------
% -------------------- Cloud optical properties --------------------------
% ------------------------------------------------------------------------

inputs.RT.re = 10;      % microns
inputs.RT.tau_c = 10;
% ------------------------------------------------------------------------



% Define the parameterization scheme used to comptue the optical quantities
if inputs.RT.use_custom_mie_calcs==false
    inputs.RT.wc_parameterization = 'mie interpolate';
else
    %wc_parameterization = '../data/wc/mie/wc.mie_test.cdf interpolate';
    inputs.RT.wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
end

% --------------------------------------------------------------
% --------------------------------------------------------------


% define the solar zenith angle
inputs.RT.sza = 15;           % degree

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
inputs.RT.phi0 = 0;         % degree

% define the viewing zenith angle
inputs.RT.vza = 0; % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

inputs.RT.vaz = 0;     % degree


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

% Using measurements from the AMSR2 instrument, a passive microwave
% radiometer for 17 Jan 2024
inputs.RT.modify_waterVapor = false;

inputs.RT.waterVapor_column = 30;              % mm - milimeters of water condensed in a column
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


%% Write each INP file and Calculate Reflectance


Refl_model = zeros(length(inputs.RT.re), length(inputs.RT.tau_c), size(inputs.RT.wavelength, 1));

idx = 0;


tic
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
            num2str(tc), '/', num2str(length(inputs.RT.tau_c)),']...', newline])






        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory

        inputName = [num2str(inputs.RT.wavelength(1)), '-', num2str(inputs.RT.wavelength(2)), 'nm_',...
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


        % Define the band model to use
        % of radiative transfer
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'mol_abs_param', inputs.RT.band_parameterization,' ', '# Band model');


        % Define the location and filename of the atmopsheric profile to use
        % ------------------------------------------------
        formatSpec = '%s %5s %s \n';
        fprintf(fileID, formatSpec,['atmosphere_file ','../data/atmmod/', inputs.RT.atm_file],...
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
            fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename], ' ', '# Location of water cloud file');

            % Define the percentage of horizontal cloud cover
            % This is a number between 0 and 1
            % ------------------------------------------------
            formatSpec = '%s %f %5s %s \n';
            fprintf(fileID, formatSpec,'cloudcover wc', inputs.RT.percent_cloud_cover, ' ', '# Cloud cover percentage');


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
        fprintf(fileID, formatSpec,'wavelength', inputs.RT.wavelength(1), inputs.RT.wavelength(2), ' ', '# Wavelength range');




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
        fprintf(fileID, formatSpec,'verbose');


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

        % define the inputSettings used to run in the reflectance
        % calculation
        inputSettings_2run = inputSettings(2,:);
        
        % grab the source function
        source = inputSettings_2run{7};


        parfor ww = 1:size(inputs.bands2run, 1)

            ds2run = struct();
            % create a data structure for each spectral channel
            ds2run.wavelength = ds.wavelength(ds.wavelength>=spec_response_2run.wavelength(ww,1) & ...
                ds.wavelength<=spec_response_2run.wavelength(ww,end));

            ds2run.radiance.value = ds.radiance.value(ds.wavelength>=spec_response_2run.wavelength(ww,1) & ...
                ds.wavelength<=spec_response_2run.wavelength(ww,end));
            

            % change the source function so that it matches the spectral
            % channel
            inputSettings_2run{7} = source(1:100, :);



            % compute the reflectance
            Refl_model(rr, tc, ww) = reflectanceFunction_4EMIT(inputSettings(2,:), ds2run,...
                spec_response_2run.value(ww, :)');




        end


    end

end





% ----------------------------------------------
% ---------- SAVE REFLECTANCE OUTPUT! ----------
% ----------------------------------------------

rev = 1;




filename = [folderpath_reflectance,'reflectance_calcs_EMIT_sim-ran-on-',char(datetime("today")),...
    '_rev', num2str(rev),'.mat'];

while isfile(filename)
    rev = rev+1;
    filename = [folderpath_reflectance,'reflectance_calcs_EMIT_sim-ran-on-',char(datetime("today")),...
        '_rev', num2str(rev),'.mat'];
end

save(filename,"inputs", "Refl_model");

toc

