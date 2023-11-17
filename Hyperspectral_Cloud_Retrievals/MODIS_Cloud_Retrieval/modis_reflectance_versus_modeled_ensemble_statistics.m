%% Comparing MODIS reflectances with Modeled reflectances fro the first 7 channels of MODIS, and for 3 different days (3 different scenes)

% By Andrew John Buggee

%% Define the MODIS scenes to loop through for this analysis

clear variables

% Defien MODIS data folder

modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/';

% define the days to test in this analysis
data2test = {'2008_10_18/', '2008_10_25/', '2008_11_02/', '2008_11_09/', '2008_11_13/'};


% Define the folder path where all .INP files will be saved
folder2save = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/comparing_modis_libRadTran/';

%% Define the constraints


% define the number of pixels to use from each data set
inputs.pixels.n_pixels = 100;

% --- Define the droplet radius constraint ---
inputs.pixels.re_min_threshold = 3;               % microns
inputs.pixels.re_max_threshold = 25;                % microns

% ----- Define the cloud optical thickness limits -----
inputs.pixels.tau_min_threshold = 3;
inputs.pixels.tau_max_threshold = 100;

% ----- Define the uncertainty limits -----
inputs.pixels.retrieval_uncertainty_re = 10;               % percentage
inputs.pixels.retrieval_uncertainty_tau = 10;              % percentage
inputs.pixels.reflectance_uncertainty = 5;                 % percentage


% ----- Define the MODIS bands to evaluate in this analysis -----
inputs.bands2run = 1:7; % these are the bands that we will run uvspec with



%% Define the radiative transfer inputs

% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------


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
inputs.RT.sourceFile_resolution = 1;                  % nm
% Define the source file
inputs.RT.source_file = '../data/solar_flux/kurudz_1.0nm.dat';

% define the atmospheric data file
inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
inputs.RT.surface_albedo = 0.05;




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
inputs.RT.linear_cloudFraction = false;
% if false, define the cloud cover percentage
inputs.RT.cloud_cover = 1;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
inputs.RT.use_MODIS_cloudTopHeight = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS above cloud water vapor? ---------
inputs.RT.use_MODIS_aboveCloudWaterVapor = true;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
inputs.RT.use_custom_mie_calcs = false;
% ------------------------------------------------------------------------
% This string is used to compute the LWC from optical depth and effective radius
% can be 'hu' or 'mie interpolate'
inputs.RT.wc_parameterization = 'mie interpolate';        % use the hu and stamnes parameterization for converting cloud properties to optical properties
% define the type of droplet distribution
inputs.RT.drop_distribution_str = 'gamma';
% define the distribution varaince
% 7 is the value libRadTran uses for liquid water clouds
inputs.RT.drop_distribution_var = 7;
% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';
% define how liquid water content will be computed
% can either be 'mie' or '2limit'
inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius

% define the wavelength used for the optical depth as the 650 nm
inputs.RT.lambda_for_tau = modisBands(1);
inputs.RT.lambda_for_tau = inputs.RT.lambda_for_tau(1);            % nm


% --------------------------------------------------------------
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


% ----- Do you want a long error message? -----
% if so, set error message to 'verbose'. Otherwise, set error message to
% 'quiet'
inputs.RT.err_msg = 'quiet';




%% Step through each data set, grab 100 random pixels that match certain requirements

% preallocate a few cells and arrays
wc_filename = cell(inputs.pixels.n_pixels, length(data2test));
inputName = cell(inputs.pixels.n_pixels, length(data2test), length(inputs.bands2run));
outputName = cell(inputs.pixels.n_pixels, length(data2test), length(inputs.bands2run));
pixels2use = cell(1, length(data2test));





for dd = 1:length(data2test)

    % print the data set currently being evaluated
    disp([newline, 'Analyzing data set ', num2str(dd), ' ...', newline])



    % Load MODIS data set
    [modis, inputs.L1B_filename] = retrieveMODIS_data([modisFolder, data2test{dd}]);




    % -----------------------------------------------------------------------
    % --------------- Define the spectral response function -----------------
    % -----------------------------------------------------------------------

    % using the input 'bands2run', this defines the spectral bands that will be
    % written into INP files.

    % check to see if the MODIS instrument is aboard Terra or Aqua
    if strcmp(inputs.L1B_filename{1}(1:3), 'MOD')==true

        % Then read in the spectral response functions for the terra instrument
        inputs.spec_response = modis_terra_specResponse_func(inputs.bands2run, inputs.RT.sourceFile_resolution);

    elseif strcmp(inputs.L1B_filename{1}(1:3), 'MYD')==true

        % Then read in the spectral response functions for the Aqua instrument
        inputs.spec_response = modis_aqua_specResponse_func(inputs.bands2run, inputs.RT.sourceFile_resolution);

    end



    % ------------------------------------------------------------------------


    % day of the year
    inputs.RT.day_of_year = str2double(inputs.L1B_filename{1}(15:17));


    % ----------------------------------
    % ----- Find Suitable Pixels -------
    % ----------------------------------

    % check to see if a suitable pixels file already exists
    inputs.pixels_file_flag = isfile([modisFolder,data2test{dd},'suitablePixels.mat']);

    if inputs.pixels_file_flag == false

        % if no suitable pixels file exists, make one
        disp([newline, 'Searching for suitable pixels ...', newline])

        pixels = findSuitablePixel(modis,inputs);


        save([modisFolder,data2test{dd},'suitablePixels.mat'],'pixels','inputs');

        % save these pixels in a .mat file, along with the inputs


    else

        % load a previously computed suitbale pixels file
        load([modisFolder,data2test{dd},'suitablePixels.mat'],'pixels');

    end




    % sample from the indexes above
    pixels2use{dd}.idx = randsample(pixels.res1km.index, inputs.pixels.n_pixels, false);

    % determine the rows and columns
    [pixels2use{dd}.row, pixels2use{dd}.col] = ind2sub(size(modis.cloud.effRadius17), pixels2use{dd}.idx);

    % clear the pixels strucutre, we don't need it anymore
    clear('pixels');



    % ----- Step through each pixel -----
    for nn = 1:inputs.pixels.n_pixels

        % --- Store the MODIS reflectances -----
        data2compare.modis_refl{dd}(nn,:) = reshape(modis.EV1km.reflectance(pixels2use{dd}.row(nn), pixels2use{dd}.col(nn), :),[], length(inputs.bands2run));
        data2compare.modis_refl_uncert{dd}(nn,:) = 0.01 * reshape(modis.EV1km.reflectanceUncert(pixels2use{dd}.row(nn), pixels2use{dd}.col(nn), :),[], length(inputs.bands2run));


        % most uncertainties for the modis optical retrieval are between 2
        % and 10 percent. So lets round off all re and tau values to the 1000th decimal
        % place


        % ---- Store the effective radius ----
        data2compare.modis_re{dd}(nn) = round(modis.cloud.effRadius17(pixels2use{dd}.idx(nn)), 3);

        % ---- Store the effective radius uncertainty ----
        data2compare.modis_re_uncert{dd}(nn) = modis.cloud.effRad_uncert_17(pixels2use{dd}.idx(nn));

        % ---- Store the optical thickness ----
        data2compare.modis_opt_thickness{dd}(nn) = round(modis.cloud.optThickness17(pixels2use{dd}.idx(nn)), 3);

        % ---- Store the optical thickness uncertainty ----
        data2compare.modis_opt_thickness_uncert{dd}(nn) = modis.cloud.optThickness_uncert_17(pixels2use{dd}.idx(nn));

        % ---- Store the cloud top height ----
        data2compare.modis_cloud_top_height{dd}(nn) = modis.cloud.topHeight(pixels2use{dd}.idx(nn));            % m agl

        % ---- Store the column water vapor ----
        data2compare.modis_column_water_vapor{dd}(nn) = modis.cloud.aboveWaterVaporCol(pixels2use{dd}.idx(nn)) * 10;   % mm


        % ----- Store the solar Zenith angle -----
        data2compare.modis_sza{dd}(nn) = modis.solar.zenith(pixels2use{dd}.idx(nn));           % degree

        % ----- Store the solar azimuth angle angle -----
        % this is how we map MODIS azimuth of the sun to the LibRadTran measurement
        data2compare.modis_phi0{dd}(nn) = modis.solar.azimuth(pixels2use{dd}.idx(nn)) + 180;         % degree

        % ----- Store the viewing azimuth angle angle -----
        % define the viewing azimuth angle
        % to properly map the azimuth angle onto the reference plane used by
        % libRadTran, we need an if statement
        if modis.sensor.azimuth(pixels2use{dd}.idx(nn))<0
            data2compare.modis_phi{dd}(nn) = 360 + modis.sensor.azimuth(pixels2use{dd}.idx(nn));
        else
            data2compare.modis_phi{dd}(nn) = modis.sensor.azimuth(pixels2use{dd}.idx(nn));
        end


        % ----- Store the viewing zenith angle -----
        % define the viewing zenith angle
        data2compare.modis_vza{dd}(nn) = double(modis.sensor.zenith(pixels2use{dd}.idx(nn))); % values are in degrees;                        % degree






        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------


        % define the geometric location of the cloud top and cloud bottom
        if inputs.RT.use_MODIS_cloudTopHeight==false

            z_topBottom = [9,8];          % km above surface

        else

            % if the cloud top height is below 1 km, make the lower altitude 0
            % ----- FIX THE 2% SYSTEMATIC BIAS OF MY REFLECTANCE ----
            % My reflectances, when using MODIs retrieved cloud top height,
            % effective radius and clou doptical depth, are consistantly 2%
            % less than the MODIS measurements.

            if data2compare.modis_cloud_top_height{dd}(nn)>1000

                %                z_topBottom = [data2compare.modis_cloud_top_height{dd}(nn), data2compare.modis_cloud_top_height{dd}(nn) - 1000]./1000;      % km above surface

                % ---- TESTING CLOUD TOP HEIGHT UNCERTAINTY ------
                % arbitrarily add or subtract 1 km to each MODIS retrieval to simulate an
                % uncertainty of 1 km
                rand_sign = 0;
                while rand_sign==0
                    rand_sign = randi([-1, 1]);
                end

                zTop = data2compare.modis_cloud_top_height{dd}(nn)*(1 + rand_sign*0.5);         % add or subtract 1 km to the MODIS retrieval
                zBot = zTop - 1000;
                if zBot<0
                    zBot = 0;
                end

                z_topBottom = [zTop, zBot]./1000;      % km above surface
                % -------------------------------------------------------

            elseif data2compare.modis_cloud_top_height{dd}(nn)<=1000

                %                 z_topBottom = [data2compare.modis_cloud_top_height{dd}(nn), 0]./1000;      % km above surface

                % ---- TESTING CLOUD TOP HEIGHT UNCERTAINTY ------
                % arbitrarily add or subtract 1 km to each MODIS retrieval to simulate an
                % uncertainty of 1 km
                rand_sign = 0;
                while rand_sign==0
                    rand_sign = randi([-1, 1]);
                end

                zTop = data2compare.modis_cloud_top_height{dd}(nn)*(1 + rand_sign*0.5);         % add or subtract 1 km to the MODIS retrieval
                zBot = zTop - 1000;
                if zBot<0
                    zBot = 0;
                end

                z_topBottom = [zTop, zBot]./1000;      % km above surface
                % -------------------------------------------------------

            end

        end

        % ------------------------------------------------------
        % --------------------VERY IMPORTANT ------------------
        % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
        % ------------------------------------------------------
        wc_filename{nn,dd} = write_wc_file(data2compare.modis_re{dd}(nn),data2compare.modis_opt_thickness{dd}(nn),...
            z_topBottom, inputs.RT.lambda_for_tau, inputs.RT.drop_distribution_str, inputs.RT.drop_distribution_var, ...
            inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str, dd*nn);

        wc_filename{nn,dd} = wc_filename{nn,dd}{1};

        % ------------------------------------------------------------------------
        % ------------------------------------------------------------------------



        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory
        if inputs.RT.yesCloud==true

            % write an input file for all wavelenghts desire
            for ww = 1:length(inputs.bands2run)

                inputName{nn,dd,ww} = [num2str(floor((inputs.spec_response{ww}(end,1)-inputs.spec_response{ww}(1,1))/2 + inputs.spec_response{ww}(1,1))),...
                    'nm_withCloudLayer_',num2str(dd*nn),'nn_row',num2str(pixels2use{dd}.row(nn)), '_col', num2str(pixels2use{dd}.col(nn)), '_',...
                    inputs.RT.atm_file(1:end-4),'.INP'];

                outputName{nn, dd, ww} = ['OUTPUT_',inputName{nn,dd,ww}(1:end-4)];






                % ----------------- ******************** ---------------------
                % ------------------ Write the INP File --------------------
                % ----------------- ******************** ---------------------

                % Open the old file for writing
                fileID = fopen([folder2save,inputName{nn, dd, ww}], 'w');

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
                fprintf(fileID, formatSpec,['atmosphere_file ','../data/atmmod/', inputs.RT.atm_file],' ', '# Location of atmospheric profile');

                % Define the location and filename of the extraterrestrial solar source
                % ---------------------------------------------------------------------
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'source solar', inputs.RT.source_file, ' ', '# Bounds between 250 and 10000 nm');


                % Define the location and filename of the extraterrestrial solar source
                % ---------------------------------------------------------------------
                formatSpec = '%s %u %5s %s \n\n';
                fprintf(fileID, formatSpec,'day_of_year', inputs.RT.day_of_year, ' ', '# accounts for changing Earth-Sun distance');


                % Define the total precipitable water
                % -------------------------------------------------------------------------
                % Define override of total precipitable water. This will force the total
                % column of water vapor to be whatever value you define.
                % If you don't wish to change the default, define the variable with nan
                if inputs.RT.use_MODIS_aboveCloudWaterVapor == true

                    if isnan(data2compare.modis_column_water_vapor{dd}(nn))==false
                        formatSpec = '%s %s %f %s %5s %s \n';
                        fprintf(fileID, formatSpec,'mol_modify','H2O', data2compare.modis_column_water_vapor{dd}(nn),' MM', ' ', '# Total Precipitable Water');
                    end


                    % ----- TSETING COLUMN WATER VAPOR UNCERTAINTY -----
                    % Add +/- 75% uncertainty to the retrieved column water vapor amount
%                     rand_sign = 0;
%                     while rand_sign==0
%                         rand_sign = randi([-1, 1]);
%                     end
%                     col_water_vapor = data2compare.modis_column_water_vapor{dd}(nn)*(1 + rand_sign*0.75);
%                     
% 
%                     if isnan(col_water_vapor)==false
%                         formatSpec = '%s %s %f %s %5s %s \n';
%                         fprintf(fileID, formatSpec,'mol_modify','H2O', col_water_vapor,' MM', ' ', '# Total Precipitable Water');
%                     end
                    % ------------------------------------------------------------------

                end


                % Define the surface albedo
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'albedo', inputs.RT.surface_albedo, ' ', '# Surface albedo of the ocean');


                % Define the Water Cloud properties, if you want a cloud in your model
                % --------------------------------------------------------------------
                if inputs.RT.yesCloud==true

                    % Define the water cloud file
                    % ------------------------------------------------
                    formatSpec = '%s %s %5s %s \n';
                    fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{nn, dd}], ' ', '# Location of water cloud file');

                    % Define the percentage of horizontal cloud cover
                    % This is a number between 0 and 1
                    % ------------------------------------------------

                    formatSpec = '%s %f %5s %s \n';
                    fprintf(fileID, formatSpec,'cloudcover wc', inputs.RT.cloud_cover, ' ', '# Cloud cover percentage');


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
                fprintf(fileID, formatSpec,'wavelength', inputs.spec_response{ww}(1,1), inputs.spec_response{ww}(end,1), ' ', '# Wavelength range');




                if inputs.RT.use_coxMunk==true

                    % Define the wind speed for the Cox-Munk ocean surface bi-directional reflectance model
                    % be solve
                    % -------------------------------------------------------------------------
                    formatSpec = '%s %f %5s %s \n\n';
                    fprintf(fileID, formatSpec,'brdf_cam u10', inputs.RT.wind_speed, ' ', '# (m/s) Ocean Surface wind speed');

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
                    fprintf(fileID, formatSpec,'aerosol_modify tau set', inputs.RT.aerosol_opticalDepth, ' ', '# Optical Depth of aerosol layer');

                end




                % Define the sensor altitude
                % ------------------------------------------------
                formatSpec = '%s %s %5s %s \n';
                fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');

                % Define the solar zenith angle
                % ------------------------------------------------

                formatSpec = '%s %f %5s %s \n';
                fprintf(fileID, formatSpec,'sza', data2compare.modis_sza{dd}(nn), ' ', '# Solar zenith angle');

                % Define the solar azimuth angle
                % -------------------------------------------------------

                formatSpec = '%s %f %5s %s \n';
                fprintf(fileID, formatSpec,'phi0', data2compare.modis_phi0{dd}(nn), ' ', '# Solar azimuth angle');

                % Define the cosine of the zenith viewing angle
                % ------------------------------------------------

                formatSpec = '%s %f %5s %s \n';
                fprintf(fileID, formatSpec,'umu', round(cosd(data2compare.modis_vza{dd}(nn) ),4), ' ', '# Cosine of the zenith viewing angle');

                % Define the azimuth viewing angle
                % ------------------------------------------------

                formatSpec = '%s %f %5s %s \n\n';
                fprintf(fileID, formatSpec,'phi', data2compare.modis_phi{dd}(nn), ' ', '# Azimuthal viewing angle');



                %                 if inputs.RT.compute_reflectivity_uvSpec==true
                %                     % Set the output quantity to be reflectivity
                %                     % ------------------------------------------------
                %                     formatSpec = '%s %s %5s %s \n\n';
                %                     fprintf(fileID, formatSpec,'output_quantity', 'reflectivity', ' ', '# Output is reflectance');
                %                 end


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






            end



        else

            inputName{nn, dd, ww} = [num2str(floor((inputs.spec_response{ww}(end,1)-inputs.spec_response{ww}(1,1))/2 + inputs.spec_response{ww}(1,1))),...
                'nm_',num2str(dd*nn),'nn_',num2str(modis.solar.zenith(pixels2use{dd}.row(nn), pixels2use{dd}.col(nn))),...
                'sza_',num2str(round(double(modis.sensor.zenith(pixels2use{dd}.row(nn), pixels2use{dd}.col(nn))))),'vza_', inputs.RT.atm_file(1:end-4),'.INP'];

            outputName{nn, dd, ww} = ['OUTPUT_',inputName{nn, dd, ww}(1:end-4)];

        end








    end







end




%% Run all INP files

R_model = zeros(inputs.pixels.n_pixels, length(data2test), length(inputs.bands2run));


% do you want uv_spec to compute reflectivity?
compute_reflectivity_uvSpec = false;

% define the number of bbands to run outside the parfor loop
num_bands2run = length(inputs.bands2run);

% define the spectral response function outside the parfor loop
spec_response = inputs.spec_response;



tic

for dd = 1:length(data2test)

    parfor nn = 1:inputs.pixels.n_pixels

        for ww = 1:num_bands2run


            % compute INP file
            [inputSettings] = runUVSPEC(folder2save,inputName{nn, dd, ww}, outputName{nn, dd, ww});

            % read .OUT file
            [ds,~,~] = readUVSPEC(folder2save, outputName{nn, dd, ww},inputSettings(2,:), compute_reflectivity_uvSpec);

            if compute_reflectivity_uvSpec==false
                % save reflectance
                R_model(nn, dd, ww) = reflectanceFunction_4modis(inputSettings(2,:), ds, spec_response{ww}(:,2));

            else

                R_model(nn, dd, ww) = ds.reflectivity.value;
            end



        end

    end

end


toc


% Save all the calculations!
save(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/Fixing_Reflectance_Estimate/',...
    'modis_libRadTran_reflectance_comparison_', char(datetime('today')),'.mat'], "inputs", "data2test", "data2compare", "R_model")



%% Make a subplot of comparing LibRadTran reflectance to MODIS for each wavelength

% set the font label size
label_fontSize = 20;

% set the order of wavelengths from small to large
index_sort = [3, 4, 1, 2, 5, 6, 7];

% create 1 to 1 line
% find the minimum and maximum values to create a y=x line

one2one_line = linspace(0.1, 0.85, 150);

% define the subplot height and width
w = 0.35;
h = 0.35;

% define the x and y position of the first subplot
x0 = -0.01;
y0 = 0.6;

% define the x and y spacing between each plot
dx = 0.225;
dy = 0.45;



figure1 = figure;
set(figure1, 'Position', [0 0 2200 1200])


% compute the total average of all data points for each wavelength
avg_ratio_perWavelength = zeros(1, length(index_sort));

% compute the variance of all data points for each wavelength
var_ratio_perWavelength = zeros(1, length(index_sort));

for ww = 1:length(index_sort)
    
    %clear ax
    
    subplot(2,4,ww)
    ax(ww) = plot(one2one_line, one2one_line, 'k', 'LineWidth',1);


    hold on

    % Compute the average ratio for each day at each wavelength
    avg_ratio_perDay = zeros(1, length(data2compare));

    % Compute the varianve of the ratio for each day at each wavelength
    var_ratio_perDay = zeros(1, length(data2compare));



    % create an empty cell array for the legend string
    legend_str = cell(1, length(data2test)+1);
    % set the first legend entry to be the one-to-one line
    legend_str{1} = '$1:1$';

    
    % at each wavelength, we want to compute the avg ratio of the modis
    % reflectance and the libRadTran computed reflectacne
    modis_data2avg_total = [];
    model_data2avg_total = [];


    for dd = 1:length(data2test)



        % somtimes, the retrievials will meet all the requirements, but
        % some wavelengths will have high uncertainty. Create an index that
        % keeps that measurments which are less than the reflectance
        % uncertainty threshold defined above
        index_R_high_uncertainty = data2compare.modis_refl_uncert{dd}(:,index_sort(ww))<=(inputs.pixels.reflectance_uncertainty*0.01);


        % store all the data across different days for a single wavelength
        % in an array
        modis_data2avg_total = [modis_data2avg_total; data2compare.modis_refl{dd}(index_R_high_uncertainty, index_sort(ww))];
        model_data2avg_total = [model_data2avg_total; R_model(index_R_high_uncertainty,dd, index_sort(ww))];

        if sum(index_R_high_uncertainty)<inputs.pixels.n_pixels

            warning([newline, 'Data from ', data2test{dd}(6:7), '/', data2test{dd}(9:10), ' at wavelength: ',...
                num2str(round(median(inputs.spec_response{ww}(:,1)))), ' nm was found to have ',...
                num2str(inputs.pixels.n_pixels - sum(index_R_high_uncertainty)),...
                ' measurements with an uncertainty greater than ', num2str(inputs.pixels.reflectance_uncertainty), '%', newline])
        end


        % compute the average ratio for each day at each wavelength
        avg_ratio_perDay(dd) = mean(data2compare.modis_refl{dd}(index_R_high_uncertainty, index_sort(ww))./R_model(index_R_high_uncertainty,dd, index_sort(ww)));
        
        % comptue the variance of the ratio for each day and each
        % wavelength
        var_ratio_perDay(dd) = var(data2compare.modis_refl{dd}(index_R_high_uncertainty, index_sort(ww))./R_model(index_R_high_uncertainty,dd, index_sort(ww)));
        
        % --- Legend displays mean and Varaince (but it's very small) ---
        % define the legend string for each day
%         legend_str{dd+1} = [data2test{dd}(6:7), '/', data2test{dd}(9:10), ': $\mu$ = ', num2str(round(avg_ratio_perDay(dd), 2)), ...
%                             ' $\sigma^{2}$ = ', num2str(round(var_ratio_perDay(dd), 2))];

        % --- Legend displays mean  ---
        % define the legend string for each day
        legend_str{dd+1} = [data2test{dd}(6:7), '/', data2test{dd}(9:10), ': $\mu$ = ', num2str(round(avg_ratio_perDay(dd), 2))];
        
       errorbar(R_model(index_R_high_uncertainty,dd, index_sort(ww)), data2compare.modis_refl{dd}(index_R_high_uncertainty, index_sort(ww)),...
           data2compare.modis_refl_uncert{dd}(index_R_high_uncertainty, index_sort(ww)), 'vertical','.',...
           'Color', mySavedColors(dd, 'fixed'), 'MarkerSize',25)
    
       hold on
    
    end


    % Next, calculate the average over every day for each wavelength
    avg_ratio_perWavelength(ww) = mean(modis_data2avg_total./model_data2avg_total, 'all');
    
    % now compute the variance of all the data for each wavelength
    var_ratio_perWavelength(ww) = var(modis_data2avg_total./model_data2avg_total, 1, 'all');


    % we should print the legend with the corresponding coefficient value,
    % the average value for each data set
    legend(legend_str, 'Interpreter','latex', 'Location', 'northwest')

    hold on; grid on; grid minor
    xlim([one2one_line(1), one2one_line(end)])
    ylim([one2one_line(1), one2one_line(end)])

    if ww == 1 || ww == 5
        xlabel('LibRadTran Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
        ylabel('MODIS Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
    end
    title([num2str(round(median(inputs.spec_response{index_sort(ww)}(:,1)))), ' nm'])
    axis square





end



% ---- Change the size and position of each subplot! ----

for ww = 1:length(ax)

    % SET SUBPLOT POSITION
    if ww<5
        %         ax(ww).Parent.Position = [(x0 + (ww-1)*dx), y0, w, h];
        ax(ww).Parent.Position = [(x0 + (ww-1)*dx), y0, w, h];

        % insert a text box with the total average and the variance for all
        % days at each wavelength
        dim = [0.19+(ww-1)*dx, 0.635, 0.065, 0.046];
        str = ['$\mu_{tot}$ = ', num2str(round(avg_ratio_perWavelength(ww),2)), newline,...
            '$\sigma^{2}_{tot}$ = ', num2str(round(var_ratio_perWavelength(ww),5))];
        annotation('textbox',dim,'String',str,...
            'FontSize', 20, 'FontWeight', 'Bold', 'Interpreter', 'latex')

    else
        ax(ww).Parent.Position = [(x0 + (ww-5)*dx), y0 - dy, w, h];

        % insert a text box with the total average and the variance for all
        % days at each wavelength
        dim = [0.19+(ww-5)*dx, 0.635-dy, 0.065, 0.046];
        str = ['$\mu_{tot}$ = ', num2str(round(avg_ratio_perWavelength(ww),2)), newline,...
            '$\sigma^{2}_{tot}$ = ', num2str(round(var_ratio_perWavelength(ww),5))];
        annotation('textbox',dim,'String',str,...
            'FontSize', 20, 'FontWeight', 'Bold', 'Interpreter', 'latex')

    end

end



% ---------------------------------------------------------------------
% ---------------------------------------------------------------------




