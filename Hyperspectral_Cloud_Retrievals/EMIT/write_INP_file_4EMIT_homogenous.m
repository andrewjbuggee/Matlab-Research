%% This function will write a .INP uvspec file for LibRadTran



%%

function [inpNames, inputs] = write_INP_file_4EMIT_homogenous(inputs, pixels2use, emit)

% ------------------------------------------------
% ---------- INPUTS AND FUNCTION SET UP ----------
% ------------------------------------------------


% what computer are we using?

% a template file has been set up to be edited and saved as a new file
% determine which computer is being used
userName = whatComputer;

if strcmp(userName,'anbu8374')

    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];

elseif strcmp(userName,'andrewbuggee')

    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4'];

else
    error('I dont recognize this computer user name')
end


addpath(libRadTran_path);

%% DEFINE UVSPEC INP FILE INPUTS


% -----------------------------------------------------------------------
% -------------- Define the grid of r_e and tau_c values ----------------
% -----------------------------------------------------------------------
% MODIS only considers homogenous plane parallel clouds. Lets construct the
% re matrix needed to create homogenous water clouds using write_wc_file
re = 3:3:24;     % microns - values of re that we wish to model

tau_c = [1,3, 5:5:50];      % values of tau that we wish to model

% re = 12.965;     % microns - values of re that we wish to model
% tau_c = 4.8949;      % values of tau that we wish to model




% save these values to the inputs structure
inputs.RT.re = re;
inputs.RT.tau_c = tau_c;
% -----------------------------------------------------------------------




% -----------------------------------------------------------------------
% -- Use the spectral response function to define the wavelength range --
% -----------------------------------------------------------------------

wavelengths = emit.spec_response.wavelength;

% ------------------------------------------------------------------------



% for each spectral bin, we have an image on the ground composed of 2030 *
% 1354 pixels. The swath on the ground is large enough that the solar
% zenith and solar azimuth change significantly.
% Ideally at each spectral bin, and
% each pixel, we could calculate the reflectance function and how it varies
% with changing tau and change re.


pixel_row = pixels2use.row; % for everything I need in this code, we use the 1km resolution pixel locations
pixel_col = pixels2use.col; %

% ---- DEFINE THE FOLDER NAME WHERE INP AND OUT FILES ARE STORED ---
% If the folder doesnt exist, create it!
folder2save_libRadTran_INP_OUT = inputs.folder2save.libRadTran_INP_OUT; % where the newly created .inp files will be saved

% If the folder doesn't exist, create it
if ~exist(folder2save_libRadTran_INP_OUT, 'dir')
    mkdir(folder2save_libRadTran_INP_OUT)

else

    % If the folder exists, delete all files
    delete([folder2save_libRadTran_INP_OUT,'*'])
end







% we repeat the re vector for every single optical depth
% Therefore we create a grid of values where there is a water cloud file
% for each re and each tau_c
re_vec = repmat(re, 1, length(tau_c));
% The tau vector is a column vector the length of the re vector
tau_vec = repmat(tau_c, length(re), 1);
tau_vec = tau_vec(:)';


% % the above matrices must be rerun for every modis band!
% re = repmat(re,1,length(inputs.bands2run));
% tau_c = repmat(tau_c, 1, length(inputs.bands2run));

% create the lambda matrix with all values needed to create INP files

% for the write_wc_file function we need to define the wavelength that is
% used to compute the optical depth. It should be bands 1 3 or 4.



% For now, lets only use the midpoint of the EMIT band in the visible
% spectrum used for this analysis
wavelength_forTau = emit.radiance.wavelength(inputs.bands2run(1));                % nm



% ---------------------------------------
% For now lets hard code the cloud height
% ---------------------------------------
% cloud top height and thickness is set to 9km and 1km respectively,
% the same values described in: "Overview of the MODIS Collection 6 
% Cloud Optical Property (MOD06) Retrieval Look-up Tables" Amarasinghe
% et. al 2017
%z_topBottom = [9, 8];                     % km - altitude above surface for the cloud top and cloud bottom
z_topBottom = [3, 2];
% ---------------------------------------
% ---------------------------------------


% Let's integrate over a modified gamma droplet distribution
dist_str = inputs.RT.drop_distribution_str;
% libRadTran stats in their manual that a distribution variance value of 7
% is adequate for liquid water clouds. The distribution variance needs to
% be the same size as the effective radius vector
dist_var = repmat(inputs.RT.drop_distribution_var,size(re_vec,1), size(re_vec,2));

% Since this is writing files to do the TBLUT retrieval,
% we assume clouds are vertical homogeneous
vert_homogeneous_str = 'vert-homogeneous';


% define how liquid water content will be computed
parameterization_str = inputs.RT.parameterization_str;


% -----------------------------------
% ---- Write a Water Cloud file! ----
% -----------------------------------

% ------------------------------------------------------
% --------------------VERY IMPORTANT ------------------
% ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
% ------------------------------------------------------
wc_filename = write_wc_file(re_vec,tau_vec,z_topBottom, wavelength_forTau, dist_str,...
    dist_var, vert_homogeneous_str, parameterization_str, 0);
% rehape so we can step through values for r and tau
wc_filename = reshape(wc_filename, length(re),length(tau_c));

% ------------------------------------------------------
% ------------------------------------------------------




% we have 4 variables that can change for each INP file

%   1) pixel
%   2) modis band
%   3) re
%   4) tau_c



% step through each pixel, and feed all re values, tau values, and modis
% bands of interest


% step through each band, each effective raidus and each optical depth
inpNames = cell(length(pixel_row), length(re),length(tau_c),length(inputs.bands2run));

for pp = 1:length(pixel_row)


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
    fileBegin = ['pixel_',num2str(pixel_row(pp)),'r_',num2str(pixel_col(pp)),'c_sza_',...
                 num2str(round(sza)),'_vza_',num2str(round(vza)),'_band_'];

    for bb = 1:length(inputs.bands2run)

        band_num = inputs.bands2run(bb);        % modis band number that defines the upper and lower wavelength boundaries used to compute the equation of radiative transfer


        for rr = 1:length(re)



            for tt = 1:length(tau_c)


                % redefine the old file each time
                inpNames{pp,rr,tt,bb} = [fileBegin,num2str(band_num),'_r_',num2str(re(rr)),'_T_',num2str(tau_c(tt)),'.INP'];





                % ------------------------------------------------------
                % --------------------VERY IMPORTANT ------------------
                % ----- TESTING IMPORTANCE OF CLOUD TOP HEIGHT --------
                % ------------------------------------------------------

%                 % define the geometric location of the cloud top and cloud bottom
%                 if inputs.RT.use_MODIS_cloudTopHeight==false
%                     z_topBottom = [9,8];          % km above surface
%             
%                 else
%             
%                     % if the cloud top height is below 1 km, make the lower altitude 0
%                     cloudTopHeight = modis.cloud.topHeight(pixel_row(pp),pixel_col(pp));
%             
%                     if cloudTopHeight>=1000
%                         z_topBottom = [cloudTopHeight, cloudTopHeight - 1000]./1000;      % km above surface
%                     elseif cloudTopHeight<1000
%                         z_topBottom = [cloudTopHeight, 0]./1000;      % km above surface
%                     end
%             
%                 end
% 
%                 wc_filename = write_wc_file(re(rr),tau_c(tt),z_topBottom, wavelength_forTau, dist_str,...
%                                 inputs.RT.drop_distribution_var, vert_homogeneous_str, parameterization_str, 0);
                % ------------------------------------------------------
                % ------------------------------------------------------








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
                fileID = fopen([folder2save_libRadTran_INP_OUT,inpNames{pp,rr,tt,bb}], 'w');

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
                    fprintf(fileID, formatSpec,'disort_intcor phase', ' ', '# Apply the Nakajima and Tanaka radiance correction');
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
               
                if inputs.RT.use_MODIS_aboveCloudWaterVapor==true

                    total_h2O_column = inputs.RT.effectiveWaterThickness.vapor(pp) * 10;    % mm of precipitable water

                    if isnan(total_h2O_column)==false
                        formatSpec = '%s %s %f %s %5s %s \n';
                        fprintf(fileID, formatSpec,'mol_modify','H2O', total_h2O_column,' MM', ' ', '# Total Precipitable Water');
                    end
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
                    fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{rr,tt}], ' ', '# Location of water cloud file');
                    %fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{1}(1)], ' ', '# Location of water cloud file');

                    % Define the percentage of horizontal cloud cover
                    % This is a number between 0 and 1
                    % ------------------------------------------------
                    if inputs.RT.linear_cloudFraction==true
                        %percent_cloud_cover(nn) = -0.1082 * modis.EV1km.reflectance(row(nn), col(nn),1) + 0.92164;      % custom linear fit 1 (seventh pass figure)
                        cloud_cover = (0.89 - 0.8459)/(0.2 - 0.7) * emit.reflectance(inputs.bands2use(1),pp) +...
                            (0.8459 - 0.7*(0.89 - 0.8459)/(0.2 - 0.7));      % custom linear fit 2 (Ninth pass figure)
                    else
                        cloud_cover = inputs.RT.cloud_cover;
                    end
                    formatSpec = '%s %f %5s %s \n';
                    fprintf(fileID, formatSpec,'cloudcover wc', cloud_cover, ' ', '# Cloud cover percentage');


                    % Define the technique or parameterization used to convert liquid cloud
                    % properties of r_eff and LWC to optical depth
                    % ----------------------------------------------------------------------
                    formatSpec = '%s %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'wc_properties', inputs.RT.wc_parameterization, ' ', '# optical properties parameterization technique');

                end



                % Define the wavelengths for which the equation of radiative transfer will
                % be solve. Should be in nm
                % -------------------------------------------------------------------------
                formatSpec = '%s %f %f %5s %s \n\n';
                fprintf(fileID, formatSpec,'wavelength', wavelengths{inputs.bands2run(bb)}(1),...
                    wavelengths{inputs.bands2run(bb)}(end), ' ', '# Wavelength range');




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























%                 % Create the water cloud file
%                 fileID = fopen([newFolder,inpNames{pp,rr,tt,bb}], 'w');
% 
%                 % fprintf writes lines in our text file from top to botom
%                 % wc.DAT files are written with the higher altitudes at the top, and the
%                 % surface at the bottom
% 
%                 % to write column vectors in a text file, we have to store them as row
%                 % vectors
% 
%                 % The first argument is the format specification
%                 % The designates what type of characters to print, such as floating point
%                 % arithmetic or string. The number before the character designates the
%                 % MINIMUM number of characters to print
% 
% 
%                 % Define which RTE solver to use
%                 formatSpec = '%s %s %5s %s \n';
%                 fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');
% 
% 
%                 % Define the number of streams to keep track of when solving the equation
%                 % of radiative transfer
%                 formatSpec = '%s %s %5s %s \n\n';
%                 fprintf(fileID, formatSpec,'number_of_streams','6',' ', '# Number of streams');
% 
% 
%                 % Define the location and filename of the US standard
%                 % atmosphere that will be used for this analysis
%                 formatSpec = '%s %s %5s %s \n';
%                 fprintf(fileID, formatSpec,'atmosphere_file','../data/atmmod/afglus.dat',' ', '# Location of atmospheric profile to use');
% 
%                 % Define the location and filename of the extraterrestrial solar source
%                 formatSpec = '%s %s %5s %s \n\n';
%                 fprintf(fileID, formatSpec,'source solar',['../data/solar_flux/',solarFlux_file], ' ', '# Bounds between 250 and 10000 nm');
% 
%                 % Define the location and filename of the extraterrestrial solar source
%                 formatSpec = '%s %f %5s %s \n\n';
%                 fprintf(fileID, formatSpec,'day_of_year', day_of_year, ' ', '# accounts for changing Earth-Sun distance');
% 
%                 % Define the ozone column
%                 %                 formatSpec = '%s %s %s %s %5s %s \n';
%                 %                 fprintf(fileID, formatSpec,'mol_modify','O3', '300.','DU', ' ', '# Set ozone column');
% 
%                 % Define the water vapor column
%                 formatSpec = '%s %s %s %s %5s %s \n';
%                 fprintf(fileID, formatSpec,'mol_modify','H2O', num2str(column_vapor), 'MM', ' ', '# Total Precipitable Water');
% 
%                 % Define the surface albedo
%                 formatSpec = '%s %s %5s %s \n\n';
%                 fprintf(fileID, formatSpec,'albedo','0.0600', ' ', '# Surface albedo of the ocean');
% 
% 
%                 % Define the water cloud file
%                 formatSpec = '%s %s %5s %s \n';
%                 fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{rr,tt}], ' ', '# Location of water cloud file');
% 
% 
%                 % Define the percentage of horizontal cloud cover
%                 % This is a number between 0 and 1
%                 formatSpec = '%s %s %5s %s \n';
%                 fprintf(fileID, formatSpec,'cloudcover wc', num2str(cloud_cover), ' ', '# Cloud cover percentage');
% 
%                 % Define the technique or parameterization used to convert liquid cloud
%                 % properties of r_eff and LWC to optical depth
%                 formatSpec = '%s %s %5s %s \n\n';
%                 fprintf(fileID, formatSpec,'wc_properties', wc_parameterization, ' ', '# optical properties parameterization technique');
% 
%                 % Define the wavelengths for which the equation of radiative transfer will
%                 % be solve
%                 formatSpec = '%s %f %f %5s %s \n\n';
%                 fprintf(fileID, formatSpec,'wavelength', wavelength(bb,1), wavelength(bb,2), ' ', '# Wavelength range');
% 
% 
%                 % Define the sensor altitude
%                 formatSpec = '%s %s %5s %s \n';
%                 fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');
% 
%                 % Define the solar zenith angle
%                 formatSpec = '%s %f %5s %s \n';
%                 fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');
% 
%                 % Define the solar azimuth angle
%                 formatSpec = '%s %f %5s %s \n';
%                 fprintf(fileID, formatSpec,'phi0', saz, ' ', '# Solar azimuth angle');
% 
%                 % Define the cosine of the zenith viewing angle
%                 formatSpec = '%s %f %5s %s \n';
%                 fprintf(fileID, formatSpec,'umu', umu, ' ', '# Cosine of the zenith viewing angle');
% 
%                 % Define the azimuth viewing angle
%                 formatSpec = '%s %f %5s %s \n\n';
%                 fprintf(fileID, formatSpec,'phi', phi, ' ', '# Azimuthal viewing angle');
% 
% 
%                 % Set the error message to quiet of verbose
%                 formatSpec = '%s';
%                 fprintf(fileID, formatSpec,'quiet');
% 
% 
%                 fclose(fileID);










            end




        end


    end







end



end

