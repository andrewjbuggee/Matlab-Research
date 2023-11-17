%% How does cloud top height affect reflectance as measured from space?

% By Andrew Buggee

%% Loop through different cloud top heights

clear variables



%% Define the parameters of the INP file




% Define the number of streams to use in your radiative transfer model
num_streams = 16;
% ------------------------------------------------------------------------


% Define the source file
%source_file = '../data/solar_flux/lasp_TSIS1_hybrid_solar_reference_p01nm_resolution.dat';
%source_file = '../data/solar_flux/kurudz_0.1nm.dat';
source_file = '../data/solar_flux/kurudz_1.0nm.dat';
source_file_resolution = 1;           % nm


% define the wavelength range. If monochromatic, enter the same number
% twice
% ------------------------------------------------------------------------

wavelength = [300, 2500];
%wavelength_step = 300;          % nm - number of wavelengths between each monochromatic calculation


% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
band_parameterization = 'reptran coarse';
%band_parameterization = 'reptran_channel modis_terra_b07';
% ------------------------------------------------------------------------




% define the atmospheric data file
atm_file = 'afglus.dat';

% define the surface albedo
albedo = 0.05;

% day of the year
day_of_year = 182;

% ------------------------------------------------------------------------
% ------------------- Define cloud top height to use----------------------
cloud_top_heights = 0.5:0.5:10;             % km
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------------------- Define an optical depth to use----------------------
tau_c = 5;             % km
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------------- Define droplet size at cloud top and bottom --------------
re_topBottom = [9, 5];             % microns

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
linear_cloudFraction = false;
% if false, define the cloud cover percentage
cloud_cover = 1;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
use_custom_mie_calcs = false;
% ------------------------------------------------------------------------


% define the type of droplet distribution
distribution_str = 'gamma';
% define the distribution varaince
distribution_var = 10;                  % for kokhanovsky, 1 is broad, 10 is narrow
% define whether this is a vertically homogenous cloud or not
vert_homogeneous_str = 'vert-non-homogeneous';
% define how liquid water content will be computed
parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
lambda_forTau = 650;            % nm



% Define the parameterization scheme used to comptue the optical quantities
if use_custom_mie_calcs==false
    wc_parameterization = 'mie interpolate';
else
    %wc_parameterization = '../data/wc/mie/wc.mie_test.cdf interpolate';
    wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
end

% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
use_coxMunk = true;
wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
yesAerosols = true;

aerosol_type = 4;               % 4 = maritime aerosols
aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
compute_reflectivity_uvSpec = true;
% --------------------------------------------------------------



% --------------------------------------------------------------
% --------------- Define the solar zenith angle ----------------
sza = 0;
% --------------------------------------------------------------



% --------------------------------------------------------------
% --------------- Define the viewing zenith angle --------------
vza = 0;
% --------------------------------------------------------------


% --------------------------------------------------------------
% --------------- Define the solar azimuth angle --------------
phi0 = 0;
% --------------------------------------------------------------


% --------------------------------------------------------------
% --------------- Define the solar azimuth angle --------------
vaz = 0;
% --------------------------------------------------------------








%% Write each INP file using various MODIS values

% Define the folder path where all .INP files will be saved
folder2save = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/testing_UVSPEC/';

inputName = cell(1, length(cloud_top_heights));
outputName = cell(1, length(cloud_top_heights));
wc_filename = cell(1, length(cloud_top_heights));




for nn = 1:length(cloud_top_heights)


    % -----------------------------------
    % ---- Write a Water Cloud file! ----
    % -----------------------------------
    % most uncertainties for the modis optical retrieval are between 2
    % and 10 percent. So lets round off all re values to the 1000th decimal
    % place

    % define the geometric location of the cloud top and cloud bottom
    z_topBottom = [cloud_top_heights(nn), cloud_top_heights(nn) - 0.5];     % km

    % define some altitude vector
    z = linspace(z_topBottom(2), z_topBottom(1), 20);


    % create a cloud droplet profile
    re = create_droplet_profile2(re_topBottom, z, 'altitude', 'adiabatic');


    % ------------------------------------------------------
    % --------------------VERY IMPORTANT ------------------
    % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
    % ------------------------------------------------------
    wc_filename{nn} = write_wc_file(re, tau_c, z_topBottom, lambda_forTau, distribution_str,...
        linspace(distribution_var, distribution_var, length(re)), vert_homogeneous_str, parameterization_str, nn);
    wc_filename{nn} = wc_filename{nn}{1};


    % ------------------------------------------------
    % ---- Define the input and output filenames! ----
    % ------------------------------------------------
    % input_names need a unique identifier. Let's give them the nn value so
    % they can be traced, and are writen over in memory
    inputName{nn} = ['Testing_cloud_top_height_of_',num2str(cloud_top_heights(nn)),'km_',atm_file(1:end-4),'.INP'];



    outputName{nn} = ['OUTPUT_',inputName{nn}(1:end-4)];



    % ----------------- ******************** ---------------------
    % ------------------ Write the INP File --------------------
    % ----------------- ******************** ---------------------

    % Open the old file for writing
    fileID = fopen([folder2save,inputName{nn}], 'w');

    % Define which RTE solver to use
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');


    % Define the number of streams to keep track of when solving the equation
    % of radiative transfer
    % ------------------------------------------------
    formatSpec = '%s %u %5s %s \n\n';
    fprintf(fileID, formatSpec,'number_of_streams', num_streams,' ', '# Number of streams');


    % Use phase function correction?
    % ------------------------------------------------
    if use_nakajima_phaseCorrection==true
        % define the pahse correction to be true
        % ------------------------------------------------
        formatSpec = '%s %5s %s \n\n';
        fprintf(fileID, formatSpec,'disort_intcor phase', ' ', '# Apply the Nakajima and Tanka radiance correction');
    end


    % Define the band model to use
    % of radiative transfer
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'mol_abs_param', band_parameterization,' ', '# Band model');


    % Define the location and filename of the atmopsheric profile to use
    % ------------------------------------------------
    formatSpec = '%s %5s %s \n';
    fprintf(fileID, formatSpec,['atmosphere_file ','../data/atmmod/',atm_file],' ', '# Location of atmospheric profile');

    % Define the location and filename of the extraterrestrial solar source
    % ---------------------------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'source solar', source_file, ' ', '# Bounds between 250 and 10000 nm');


    % Define the location and filename of the extraterrestrial solar source
    % ---------------------------------------------------------------------
    formatSpec = '%s %u %5s %s \n\n';
    fprintf(fileID, formatSpec,'day_of_year', day_of_year, ' ', '# accounts for changing Earth-Sun distance');


    % Define the total precipitable water
    % -------------------------------------------------------------------------
    % Define override of total precipitable water. This will force the total
    % column of water vapor to be whatever value you define.
    % If you don't wish to change the default, define the variable with nan
%     total_h2O_column = modis.vapor.col_nir(row(nn),col(nn))*10;        % mm of precipitable water
%     if isnan(total_h2O_column)==false
%         formatSpec = '%s %s %f %s %5s %s \n';
%         fprintf(fileID, formatSpec,'mol_modify','H2O', total_h2O_column,' MM', ' ', '# Total Precipitable Water');
%     end


    % Define the surface albedo
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'albedo', albedo, ' ', '# Surface albedo of the ocean');

    
    % Define the Water Cloud properties, if you want a cloud in your model
    % --------------------------------------------------------------------
    if yesCloud==true

        % Define the water cloud file
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{nn}], ' ', '# Location of water cloud file');

        % Define the percentage of horizontal cloud cover
        % This is a number between 0 and 1
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'cloudcover wc', cloud_cover, ' ', '# Cloud cover percentage');


        % Define the technique or parameterization used to convert liquid cloud
        % properties of r_eff and LWC to optical depth
        % ----------------------------------------------------------------------
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'wc_properties', wc_parameterization, ' ', '# optical properties parameterization technique');

    end



    % Define the wavelengths for which the equation of radiative transfer will
    % be solve
    % -------------------------------------------------------------------------
    formatSpec = '%s %f %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'wavelength', wavelength(1), wavelength(2), ' ', '# Wavelength range');

    % Define the wavelength step between each value to be computed
    % -------------------------------------------------------------------------
%     formatSpec = '%s %f %5s %s \n\n';
%     fprintf(fileID, formatSpec,'wavelength_step', wavelength_step, ' ', '# Wavelength step');




    if use_coxMunk==true

        % Define the wind speed for the Cox-Munk ocean surface bi-directional reflectance model
        % be solve
        % -------------------------------------------------------------------------
        formatSpec = '%s %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'brdf_cam u10', wind_speed, ' ', '# (m/s) Ocean Surface wind speed');

    end



    % Define the Aerosol Layer properties, if you want a cloud in your model
    % --------------------------------------------------------------------
    if yesAerosols==true

        % Turn on default aersol layer, which occupies lower 2km of model
        % --------------------------------------------------------------
        formatSpec = '%s %5s %s \n';
        fprintf(fileID, formatSpec,'aerosol_default', ' ', '# turn on Shettle (1989) boundary layer aerosols');


        % Specify the Aerosl type
        % 1=rural aersols,  4=maritime aersols,  5=Urban aerosols,
        % 6=Tropospheric aerosols
        % ------------------------------------------------
        formatSpec = '%s %u %5s %s \n';
        fprintf(fileID, formatSpec,'aerosol_haze', aerosol_type, ' ', '# Aerosol type');


        % Define aerosol layer optical depth
        % ----------------------------------------------------------------------
        formatSpec = '%s %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'aerosol_modify tau set', aerosol_opticalDepth, ' ', '# Optical Depth of aerosol layer');

    end




    % Define the sensor altitude
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');

    % Define the solar zenith angle
    % ------------------------------------------------
    % define the solar zenith angle
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');

    % Define the solar azimuth angle
    % -------------------------------------------------------
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'phi0', phi0, ' ', '# Solar azimuth angle');

    % Define the cosine of the zenith viewing angle
    % ------------------------------------------------
    % define the viewing zenith angle
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'umu', round(cosd(vza),4), ' ', '# Cosine of the zenith viewing angle');

    % Define the azimuth viewing angle
    % ------------------------------------------------
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'phi', vaz, ' ', '# Azimuthal viewing angle');


    
    if compute_reflectivity_uvSpec==true
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


end



%% Run INP files


L_model = zeros(length(wavelength(1):wavelength(2)), length(cloud_top_heights));
R_model = zeros(length(wavelength(1):wavelength(2)), length(cloud_top_heights));




tic
for nn = 1:length(cloud_top_heights)


    % compute INP file
    [inputSettings] = runUVSPEC(folder2save,inputName{nn},outputName{nn});

    % read .OUT file
    [ds,~,~] = readUVSPEC(folder2save,outputName{nn},inputSettings(2,:), compute_reflectivity_uvSpec);

    if compute_reflectivity_uvSpec==false
        % save reflectance
        R_model(:, nn) = reflectanceFunction(inputSettings(2,:), ds, spec_response);

        % save the radiance calculation

        L_model(:,nn) = ds.radiance.value;          % (mW/m^2/nm/sr) -

    else

        R_model(:,nn) = ds.reflectivity.value;
    end



end
toc



save([folder2save,'Reflectance across solar spectrum for varrying cloud top heights with optical depth of ',...
        num2str(tau_c), '_', char(datetime("today")),'.mat'],...
        'R_model', 'tau_c')


%% PLOT THE REFLECTANCE AS A FUNCTION OF WAVELENGTH AND CLOUD TOP HEIGHT

figure;
imagesc(wavelength(1):wavelength(2), cloud_top_heights, R_model')
ylabel('Cloud Top Height ($km$)', 'Interpreter','latex')
xlabel('Wavelength ($nm$)', 'Interpreter','latex')
c = colorbar;
ylabel(c, 'Reflectance ($1/sr$)', 'Interpreter','latex', 'FontSize', 28)
set(gcf, 'Position', [0 0 1200 650])
title('Nadir Viewing Reflectance of a Cloudy Scene', 'Interpreter','latex')

% flip the y-axis
set(gca, 'YDir', 'normal')

% set the colorbar limits
clim([0, 0.7])

% --- Create a Textbox with tau, solar zenith angle, viewing zenith angle
% ---
annotation('textbox',[0.9025 0.0559 0.081 0.114],...
    'String',{['$\tau_c$ = ', num2str(tau_c), newline,...
    'SZA = ', num2str(sza), newline,...
    'VZA = ', num2str(vza), newline]},...
    'LineWidth',2,...
    'Interpreter','latex',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');
