%% This function will write .INP uvspec files for LibRadTran using EMIT data



%%

function inpNames = write_INP_file_4EMIT_Gauss_Newton(inputs, pixels2use, emit, wc_filename)



%%




% for each spectral bin, we have an image on the ground composed of 2030 *
% 1354 pixels. The swath on the ground is large enough that the solar
% zenith and solar azimuth change significantly. Ideally at each spectral bin, and
% each pixel, we could calculate the reflectance function and how it varies
% with changing tau and change re.



folder2save = [inputs.folder2save.libRadTran_INP_OUT]; % where the newly created .inp files will be saved
addpath(folder2save);

% If the folder doesn't exist, create it
if ~exist(folder2save, 'dir')
    mkdir(folder2save)
end




% -----------------------------------------------------------------------
% define the EMIT bands to run using the spectral response functions
wavelength = emit.spec_response.wavelength(inputs.bands2run, :);

% ------------------------------------------------------------------------


% we have 4 variables that can change for each INP file

%   1) pixel
%   2) modis band
%   3) re
%   4) tau_c


% step through each band, each effective raidus and each optical depth
inpNames = cell(length(pixels2use.row), length(inputs.bands2run));

for pp = 1:length(pixels2use.row)

    % ------------------------------------------------
    % ------ DEFINE SOLAR AND VIEWING GEOMETRY -------
    % ------------------------------------------------

    % Define the solar zenith angle
    % ------------------------------------------------
    sza = emit.obs.solar.zenith(pp);           % degree

    % Define the solar azimuth angle
    % -------------------------------------------------------
    % this is how we map EMIT-defined solar azimuth to the LibRadTran
    % definition.
    % LibRadTran defines the solar azimuth clockwise from South.
    % So at 0deg the Sun is due south, 90deg the Sun is due West,
    % and so on. EMIT defines the solar azimuth clockwise from due North.
    % So we need to add 180deg to the EMIT values, but modulo 360, since it
    % needs to wrap around.
    saz = mod(emit.obs.solar.azimuth(pp) + 180, 360);         % degree



    % ------ Define the Viewing Geometry ------
    % we need the cosine of the zenith viewing angle
    % positive values solve for upwelling radiance, where the sensor is
    % defined to be looking down towrads the Earth's surface. negative
    % values solve for downwelling radiance, where the sensor is looking
    % upwards towards the sky

    % Define the cosine of the zenith viewing angle
    % ------------------------------------------------
    % define the viewing zenith angle
    vza = double(emit.obs.sensor.zenith(pp)); % values are in degrees;                        % degree



    % Define the azimuth viewing angle
    % ------------------------------------------------
    % LibRadTran defines the viewing azimuth clockwise from North.
    % So at 0deg the sensor is in the north looking south, at 90deg
    % the sensor is in the east looking West, and so on.
    % EMIT defines the sensor azimuth clockwise from due North.
    % So we don't need to change the emit values
    vaz = emit.obs.sensor.azimuth(pp);



    % create the begining of the file name string
    fileBegin = ['pixel_',num2str(pixels2use.row(pp)),'r_',num2str(pixels2use.col(pp)),'c_sza_',num2str(round(sza)),'_vza_',num2str(round(vza)),'_'];

    for bb = 1:length(inputs.bands2run)


        % Define the file name
        inpNames{pp,bb} = [fileBegin,num2str(round(median(wavelength(bb,:)))),'nm.INP'];




        % ------------------------------------------------------------
        % --------------------- WRITE INP FILE -----------------------
        % ------------------------------------------------------------
        % fprintf writes lines in our text file from top to botom
        % wc.DAT files are written with the higher altitudes at the top, and the
        % surface at the bottom

        % to write column vectors in a text file, we have to store them as row
        % vectors

        % The first argument is the format specification
        % The designates what type of characters to print, such as floating point
        % arithmetic or string. The number before the character designates the
        % MINIMUM number of characters to print




        % ----------------- ******************** ---------------------
        % ------------------ Write the INP File --------------------
        % ----------------- ******************** ---------------------

        % Open the old file for writing
        fileID = fopen([folder2save,inpNames{pp,bb}], 'w');

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
            fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{1}], ' ', '# Location of water cloud file');
            %fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{1}(1)], ' ', '# Location of water cloud file');


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
        fprintf(fileID, formatSpec,'wavelength', wavelength(bb,1), wavelength(bb,end), ' ', '# Wavelength range');




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




        % Define the sensor altitude
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');

        % Define the solar zenith angle
        % ------------------------------------------------
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');

        % Define the solar azimuth angle
        % -------------------------------------------------------
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'phi0', saz, ' ', '# Solar azimuth angle');

        % Define the cosine of the zenith viewing angle
        % ------------------------------------------------
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'umu', round(cosd(vza),4), ' ', '# Cosine of the zenith viewing angle');

        % Define the azimuth viewing angle
        % ------------------------------------------------
        formatSpec = '%s %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'phi', vaz, ' ', '# Azimuthal viewing angle');




        % Set the error message to quiet of verbose
        % ------------------------------------------------
        formatSpec = '%s';
        fprintf(fileID, formatSpec, inputs.RT.err_msg);


        % Close the file!
        fclose(fileID);
        % ----------------------------------------------------
        % ----------------------------------------------------



































        % % ------------------------------------------------------
        % % ---------- Lets define the % of cloud cover ----------
        % % ------------------------------------------------------
        %
        % if strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1850/')==true
        %
        %     % For 11/11/2008 - 18:50 data
        %     cloud_cover = 0.795;
        %
        % elseif strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1430/')==true
        %
        %     % For 11/11/2008 - 14:30 data
        %     cloud_cover = 0.725;
        %
        % elseif strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_09/')==true
        %
        %     % For 11/09/2008 data
        %     cloud_cover = 0.70;
        %
        % else
        %
        %     cloud_cover = 0.775;
        %
        % end
        %
        %
        %
        %
        %
        %
        %
        %
        % % ------------------------------------------------------
        %
        %
        %
        % % write INP files for all 7 MODIS bands. This function writes files for a
        % % single pixel at a time, for now...
        % numPixels = length(pixels2use.row);
        % numBands = length(GN_inputs.bands2use);
        %
        % inpNames = cell(numPixels, numBands);
        %
        % for pp = 1:numPixels
        %
        %     sza = modis.solar.zenith(pixels2use.row(pp),pixels2use.col(pp));
        %     phi0 = modis.solar.azimuth(pixels2use.row(pp),pixels2use.col(pp));
        %
        %     % we need the cosine of the zenith viewing angle
        %     umu = round(cosd(double(modis.sensor.zenith(pixels2use.row(pp),pixels2use.col(pp)))),3); % values are in degrees
        %     phi = modis.sensor.azimuth(pixels2use.row(pp),pixels2use.col(pp));
        %
        %     % Modis defines the azimuth viewing angle as [0,180]
        %     % and [-180,0], whereas libradtran defines the azimuth
        %     % angle as [0,360]. So we need to make this adjustment
        %
        %     if phi<0
        %         phi = phi+360;
        %     end
        %
        %     % Modis defines the azimuth viewing angle as [0,180]
        %     % and [-180,0], whereas libradtran defines the azimuth
        %     % angle as [0,360]. So we need to make this adjustment
        %
        %     if phi0<0
        %         phi0 = phi0+360;
        %     end
        %
        %     % create the begining of the file name string
        %     fileBegin = ['GN_pixel_',num2str(pixels2use.row(pp)),'r_',num2str(pixels2use.col(pp)),'c_sza_',num2str(sza),'_saz_',num2str(phi0),'_band_'];
        %
        %     for bb = 1:numBands
        %
        %         band_num = GN_inputs.bands2use(bb);        % modis band number that defines the upper and lower wavelength boundaries used to compute the equation of radiative transfer
        %
        %
        %
        %
        %         % redefine the old file each time
        %         inpNames{pp,bb} = [fileBegin,num2str(band_num),'.INP'];
        %
        %         % ------------------------------------------------------------
        %         % --------------------- WRITE INP FILE -----------------------
        %         % ------------------------------------------------------------
        %
        %
        %
        %         % Create the water cloud file
        %         fileID = fopen([newFolder,inpNames{pp,bb}], 'w');
        %
        %         % fprintf writes lines in our text file from top to botom
        %         % wc.DAT files are written with the higher altitudes at the top, and the
        %         % surface at the bottom
        %
        %         % to write column vectors in a text file, we have to store them as row
        %         % vectors
        %
        %         % The first argument is the format specification
        %         % The designates what type of characters to print, such as floating point
        %         % arithmetic or string. The number before the character designates the
        %         % MINIMUM number of characters to print
        %
        %
        %         % Define which RTE solver to use
        %         formatSpec = '%s %s %5s %s \n';
        %         fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');
        %
        %
        %         % Define the number of streams to keep track of when solving the equation
        %         % of radiative transfer
        %         formatSpec = '%s %s %5s %s \n\n';
        %         fprintf(fileID, formatSpec,'number_of_streams','6',' ', '# Number of streams');
        %
        %
        %         % Define the location and filename of the atmopsheric profile to use
        %         formatSpec = '%s %s %5s %s \n';
        %         fprintf(fileID, formatSpec,'atmosphere_file','../data/atmmod/afglus.dat',' ', '# Location of atmospheric profile to use');
        %
        %         % Define the location and filename of the extraterrestrial solar source
        %         formatSpec = '%s %s %5s %s \n\n';
        %         fprintf(fileID, formatSpec,'source solar','../data/solar_flux/kurudz_1.0nm.dat', ' ', '# Bounds between 250 and 10000 nm');
        %
        %         % Define the ozone column
        %         formatSpec = '%s %s %s %s %5s %s \n';
        %         fprintf(fileID, formatSpec,'mol_modify','O3', '300.','DU', ' ', '# Set ozone column');
        %
        %         % Define the surface albedo
        %         formatSpec = '%s %s %5s %s \n\n';
        %         fprintf(fileID, formatSpec,'albedo','0.0600', ' ', '# Surface albedo of the ocean');
        %
        %
        %         % Define the water cloud file
        %         formatSpec = '%s %s %5s %s \n';
        %         fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{1}], ' ', '# Location of water cloud file');
        %
        %
        %         % Define the percentage of horizontal cloud cover
        %         % This is a number between 0 and 1
        %         formatSpec = '%s %s %5s %s \n';
        %         fprintf(fileID, formatSpec,'cloudcover wc', num2str(cloud_cover), ' ', '# Cloud cover percentage');
        %
        %         % Define the technique or parameterization used to convert liquid cloud
        %         % properties of r_eff and LWC to optical depth
        %         formatSpec = '%s %s %5s %s \n\n';
        %         fprintf(fileID, formatSpec,'wc_properties', wc_parameterization, ' ', '# optical properties parameterization technique');
        %
        %         % Define the wavelengths for which the equation of radiative transfer will
        %         % be solve
        %         formatSpec = '%s %f %f %5s %s \n\n';
        %         fprintf(fileID, formatSpec,'wavelength', wavelength(bb,2), wavelength(bb,3), ' ', '# Wavelength range');
        %
        %
        %         % Define the sensor altitude
        %         formatSpec = '%s %s %5s %s \n';
        %         fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');
        %
        %         % Define the solar zenith angle
        %         formatSpec = '%s %f %5s %s \n';
        %         fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');
        %
        %         % Define the solar azimuth angle
        %         formatSpec = '%s %f %5s %s \n';
        %         fprintf(fileID, formatSpec,'phi0', phi0, ' ', '# Solar azimuth angle');
        %
        %         % Define the cosine of the zenith viewing angle
        %         formatSpec = '%s %f %5s %s \n';
        %         fprintf(fileID, formatSpec,'umu', umu, ' ', '# Cosine of the zenith viewing angle');
        %
        %         % Define the azimuth viewing angle
        %         formatSpec = '%s %f %5s %s \n\n';
        %         fprintf(fileID, formatSpec,'phi', phi, ' ', '# Azimuthal viewing angle');
        %
        %
        %         % Set the error message to quiet of verbose
        %         formatSpec = '%s';
        %         fprintf(fileID, formatSpec,'quiet');
        %
        %
        %         fclose(fileID);
        %
        %




    end




end







end

