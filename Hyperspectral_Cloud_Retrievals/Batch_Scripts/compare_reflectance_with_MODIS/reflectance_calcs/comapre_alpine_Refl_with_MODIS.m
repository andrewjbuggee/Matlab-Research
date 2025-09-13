%% Compare reflectances computed using DISORT on CU's Alpine cluster with MODIS observations

% By Andrew John Buggee

%% Loop through 100 different MODIS Pixels

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ***** Define the MODIS Folder *****

    modisFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];


    % ***** Define the VOCALS-REx File *****

    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/SPS_1/'];



elseif strcmp(whatComputer,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % ----- Define the MODIS folder name -----

    modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];


    % ***** Define the VOCALS-REx Folder *****

    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/SPS_1/'];


elseif strcmp(whatComputer,'curc')==true



    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


    % Define the MODIS folder name

    modisFolder = '/projects/anbu8374/MODIS_data/';


    % ***** Define the VOCALS-REx Folder *****

    vocalsRexFolder = '/projects/anbu8374/VOCALS_REx_data/';




end




%% LOAD MODIS DATA

% Load modis data and create input structure


% ----- November 9th at decimal time 0.611 (14:40) -----
%modisData = '2008_11_09/';


% ----- November 11th at decimal time 0.604 (14:30) -----
%modisData = '2008_11_11_1430/';


% ----- November 11th at decimal time 0.784 (18:50) -----
modisData = '2008_11_11_1850/';



[modis,L1B_fileName] = retrieveMODIS_data([modisFolder, modisData]);

%%

folder_paths = define_folderPaths_for_HySICS(1);


%% Define the location of the data and the save reflectances

which_computer = folder_paths.which_computer;

% Load MODIS data set

if strcmp(which_computer, 'curc')

    modisFolder = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/2024_01_27/'];

elseif strcmp(which_computer, 'andrewbuggee')

    % modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2024_01_27/'];

    % No pixels meet all constraints with H<=0.3
    % modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/'];

    % No pixels meet all constraints with H<=0.3
    % modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
    %     'Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/'];

    % No pixels meet all constraints with H<=0.3
    modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/'];


    folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/compare_reflectance_with_MODIS/reflectance_calcs'];

end


%%  Delete old files?

% First, delete files in the HySICS INP folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% delete old water vapor profiles
delete([folder_paths.atm_folder_path, '*-aboveCloud.DAT'])

% delete old MIE files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])


%% Load MODIS data

[modis,L1B_1km_fileName] = retrieveMODIS_data(modisFolder);

%% Print err and status messages?

% Define the parameters of the INP file
print_libRadtran_err = true;


%% Grab n random pixels from the suitablePixels set

n_pixels = 300;

% find pixels above a certain optical depth
inputs.pixels.tau_min_threshold = 4;

% find pixels below a certain optical depth
inputs.pixels.tau_max_threshold = 30;

% find pixels above a certain effective radius
inputs.pixels.re_min_threshold = 3;

% find pixels below a certain effective radius
inputs.pixels.re_max_threshold = 23;

% set the H-index threshold
% Zhang and Platnick (2011) state that a value below 0.1 will have few 3D
% effects
inputs.pixels.H_index_max = 2;


[pixels] = findSuitablePixel(modis,inputs);


%% Start parallel pool

start_parallel_pool(folder_paths.which_computer)



%% Set up for loop here because each pixel will need a unique set of inputs?

%% Define the INP parameters


% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';

% we model two free parameters, r_top and r_bot
inputs.RT.num_re_parameters = 2;



% Define some needed folder and file names
saveCalculations_fileName = GN_inputs.save_calcs_fileName; % where to save the computed data
INP_folderName = modisInputs.INP_folderName; % Where to save the INP files

% --- compute the forward model at our current estimate ---
r_top = current_guess(1);
r_bottom = current_guess(2);
tau_c = current_guess(3);

profile_type = GN_inputs.model.profile.type; % type of water droplet profile

% Using the same wavelength MODIS write_INP_file_4MODIS_2 uses to compute
% the cloud properties
wavelength_tau_c = modisBands(1);    % nm - Wavelength used for cloud optical depth calculation
% ----------------------------------------------------------

% --------------------------------------------
% create water cloud file with droplet profile
% --------------------------------------------

% Set up a few constants for the water cloud
H = GN_inputs.RT.cloudDepth;                                % km - geometric thickness of cloud
n_layers = GN_inputs.RT.cloud_layers;                          % number of layers to model within cloud

% Cloud top
z_top = GN_inputs.RT.cloudTop_height(pp);        % km -  cloud top height

z = linspace(z_top-H, z_top,n_layers);        % km - altitude above ground vector

indVar = 'altitude';                    % string that tells the code which independent variable we used

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
constraint = profile_type;              % string that tells the code which physical constraint to use



re = create_droplet_profile2([r_top, r_bottom], z, indVar, constraint);     % microns - effective radius vector


% Set the droplet distribution type
dist_str = GN_inputs.RT.drop_distribution_str;                                 % droplet distribution

% define the droplet distribution variance
% This should be the same length as re
% A distribution variance must be defined for each re value

% -- For now, lets assume this is constant --
dist_var = linspace(GN_inputs.RT.drop_distribution_var,GN_inputs.RT.drop_distribution_var, GN_inputs.RT.cloud_layers);              % distribution variance

vert_homogeneous_str = GN_inputs.RT.vert_homogeneous_str;     % This tells the function to create a multi-layered cloud
% define the boundaries of the cloud in Z-space
z_topBottom = [z(end), z(1)];                    % km - boundaries of the altitude vector. 

% Tell the code to use a pre-computed mie table for the extinction
% efficiency, or to use the value of the extinction paradox -> Qe = 2
parameterization_str = GN_inputs.RT.parameterization_str;


% -----------------------------------
% ---- Write a Water Cloud file! ----
% -----------------------------------

% ------------------------------------------------------
% --------------------VERY IMPORTANT ------------------
% ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
% ------------------------------------------------------
loop_var = 0;

wc_filename = write_wc_file(re,tau_c,z_topBottom, wavelength_tau_c(1,1), dist_str,...
    dist_var, vert_homogeneous_str, parameterization_str, loop_var);





%%

%% Define the parameters of the INP file


% Define the MODIs spectral band you wish to run
% ------------------------------------------------------------------------
band_num = 1;
% ------------------------------------------------------------------------



% Define the number of streams to use in your radiative transfer model
num_streams = 16;
% ------------------------------------------------------------------------


% Define the source file
%source_file = '../data/solar_flux/lasp_TSIS1_hybrid_solar_reference_p01nm_resolution.dat';
%source_file = '../data/solar_flux/kurudz_0.1nm.dat';
source_file = '../data/solar_flux/kurudz_1.0nm.dat';
source_file_resolution = 1;           % nm



% Define the spectral response function
% ------------------------------------------------------------------------
spec_response = modis_terra_specResponse_func(band_num, source_file_resolution);


% define the wavelength range. If monochromatic, enter the same number
% twice
% ------------------------------------------------------------------------
% band7 = modisBands(band_num);
% wavelength = [band7(2), band7(3)];              % nm - monochromatic wavelength calcualtion
wavelength = [spec_response{1}(1,1), spec_response{1}(end,1)];


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
day_of_year = str2double(L1B_1km_fileName{1}(15:17));

% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
use_MODIS_cloudTopHeight = true;
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
distribution_var = 20;                  % for kokhanovsky, 1 is broad, 10 is narrow
% define whether this is a vertically homogenous cloud or not
vert_homogeneous_str = 'vert-homogeneous';
% define how liquid water content will be computed
parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
band1 = modisBands(1);
lambda_forTau = band1(1);            % nm



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
compute_reflectivity_uvSpec = false;
% --------------------------------------------------------------









%% Write each INP file using various MODIS values

% Define the folder path where all .INP files will be saved
folder2save = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/solving_modeling_discrepancy_2/';
% folder2save = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
%     'LibRadTran/libRadtran-2.0.4/Matching_MODIS_Reflectance/'];

inputName = cell(1, length(idx));
outputName = cell(1, length(idx));
wc_filename = cell(1, length(idx));

% save the modis droplet size and optical depth
re = zeros(1, length(idx));
tau_c = zeros(1, length(idx));

% save the solar zenith angle
sza = zeros(1, length(idx));

% save the cloud cover percentage
percent_cloud_cover = zeros(1, length(idx));


for nn = 1:length(idx)


    % -----------------------------------
    % ---- Write a Water Cloud file! ----
    % -----------------------------------
    % most uncertainties for the modis optical retrieval are between 2
    % and 10 percent. So lets round off all re values to the 1000th decimal
    % place

    % define the total cloud optical depth
    tau_c(nn) = round(modis.cloud.optThickness17(row(nn),col(nn)), 3);
    % define the droplet size

    re(nn) = round(modis.cloud.effRadius17(row(nn),col(nn)), 3);                      % microns

    % define the geometric location of the cloud top and cloud bottom
    if use_MODIS_cloudTopHeight==false
        z_topBottom(nn,:) = [9,8];          % km above surface

    else

        % if the cloud top height is below 1 km, make the lower altitude 0
        % ----- FIX THE 2% SYSTEMATIC BIAS OF MY REFLECTANCE ----
        % My reflectances, when using MODIS retrieved cloud top height,
        % effective radius and cloud optical depth, are consistantly 2%
        % less than the MODIS measurements. 
        %cloudTopHeight = 1.75 * modis.cloud.topHeight(row(nn), col(nn));
        cloudTopHeight = modis.cloud.topHeight(row(nn), col(nn));

        if cloudTopHeight>=1000
            z_topBottom(nn,:) = [cloudTopHeight, cloudTopHeight - 1000]./1000;      % km above surface
        elseif cloudTopHeight<1000
            z_topBottom(nn,:) = [cloudTopHeight, 0]./1000;      % km above surface
        end

    end

    % ------------------------------------------------------
    % --------------------VERY IMPORTANT ------------------
    % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
    % ------------------------------------------------------
    wc_filename{nn} = write_wc_file(re(nn),tau_c(nn),z_topBottom(nn,:), lambda_forTau, distribution_str,...
        distribution_var, vert_homogeneous_str, parameterization_str, nn);
    wc_filename{nn} = wc_filename{nn}{1};


    % ------------------------------------------------
    % ---- Define the input and output filenames! ----
    % ------------------------------------------------
    % input_names need a unique identifier. Let's give them the nn value so
    % they can be traced, and are writen over in memory
    if yesCloud==true
        inputName{nn} = [num2str(floor((wavelength(2)-wavelength(1))/2 + wavelength(1))),...
            'nm_withCloudLayer_',num2str(nn),'nn_',num2str(modis.solar.zenith(row(nn),col(nn))),'sza_',num2str(round(double(modis.sensor.zenith(row(nn),col(nn))))),...
            'vza_',atm_file(1:end-4),'.INP'];
    else
        inputName{nn} = [num2str((wavelength(2)-wavelength(1))/2 + wavelength(1)),...
            'nm_',num2str(nn),'nn_',num2str(modis.solar.zenith(row(nn),col(nn))),'sza_',num2str(round(double(modis.sensor.zenith(row(nn),col(nn))))),...
            'vza_',atm_file(1:end-4),'.INP'];
    end


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
    % The default value for the US standard ATM is 14.295 mm of water
    total_h2O_column = modis.vapor.col_nir(row(nn),col(nn))*10;        % mm of precipitable water
    if isnan(total_h2O_column)==false
        formatSpec = '%s %s %f %s %5s %s \n';
        fprintf(fileID, formatSpec,'mol_modify','H2O', total_h2O_column,' MM', ' ', '# Total Precipitable Water');
    end


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
        % ------------------------------------------------
        if linear_cloudFraction==true
            %percent_cloud_cover(nn) = -0.1082 * modis.EV1km.reflectance(row(nn), col(nn),1) + 0.92164;      % custom linear fit 1 (seventh pass figure)
            percent_cloud_cover(nn) = (0.89 - 0.8459)/(0.2 - 0.7) * modis.EV1km.reflectance(row(nn), col(nn),1) +...
                                        (0.8459 - 0.7*(0.89 - 0.8459)/(0.2 - 0.7));      % custom linear fit 2 (Ninth pass figure)
        else
            percent_cloud_cover(nn) = cloud_cover;
        end
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'cloudcover wc', percent_cloud_cover(nn), ' ', '# Cloud cover percentage');


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
    sza(nn) = modis.solar.zenith(row(nn),col(nn));           % degree
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'sza', sza(nn), ' ', '# Solar zenith angle');

    % Define the solar azimuth angle
    % -------------------------------------------------------
    % this is how we map MODIS azimuth of the sun to the LibRadTran measurement
    phi0 = modis.solar.azimuth(row(nn),col(nn)) + 180;         % degree
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'phi0', phi0, ' ', '# Solar azimuth angle');

    % Define the cosine of the zenith viewing angle
    % ------------------------------------------------
    % define the viewing zenith angle
    vza = double(modis.sensor.zenith(row(nn),col(nn))); % values are in degrees;                        % degree
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'umu', round(cosd(vza),4), ' ', '# Cosine of the zenith viewing angle');

    % Define the azimuth viewing angle
    % ------------------------------------------------
    % define the viewing azimuth angle
    % to properly map the azimuth angle onto the reference plane used by
    % libRadTran, we need an if statement
    if modis.sensor.azimuth(row(nn),col(nn))<0
        vaz = 360 + modis.sensor.azimuth(row(nn),col(nn));
    else
        vaz = modis.sensor.azimuth(row(nn),col(nn));
    end
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

L_model = zeros(numel(spec_response{1}(:,1)), length(idx));
R_model = zeros(1, length(idx));

R_modis = zeros(1, length(idx));
L_modis = zeros(1, length(idx));
R_uncertainty = zeros(1, length(idx));
%L_uncertainty = zeros(1, length(idx));

% store the modis reflectance and uncertainties in a standalone array to
% avoid unecessary communication during the parallel for loop

modis_reflectance = modis.EV1km.reflectance(:, :, band_num);
modis_reflectance_uncertainty = modis.EV1km.reflectanceUncert(:,:,band_num);

modis_radiance = modis.EV1km.radiance(:, :, band_num);
%modis_radiance_uncertainty = modis.EV1km.radiance(:,:,band_num);

tic
for nn = 1:length(idx)

    % Store the modis reflectance value
    % ---------------------------------------------------------------------
    % MODIS reflectance is not divided by cos(solar zenith angle)!!
    R_modis(nn) =  modis_reflectance(row(nn), col(nn))./cosd(double(sza(nn)));
    % ---------------------------------------------------------------------
    

     % store the MODIS reflectance uncertainty
     R_uncertainty(nn) = 0.01*modis_reflectance_uncertainty(row(nn), col(nn)) *  R_modis(nn);

    % Store the modis radiance value
     L_modis(nn) =  modis_radiance(row(nn), col(nn));

     % store the MODIS reflectance uncertainty
     %L_uncertainty(nn) = 0.01*modis_reflectance_uncertainty(row(nn), col(nn)) *  R_modis(nn);

    % compute INP file
    [inputSettings] = runUVSPEC(folder2save,inputName{nn},outputName{nn});

    % read .OUT file
    [ds,~,~] = readUVSPEC(folder2save,outputName{nn},inputSettings(2,:), compute_reflectivity_uvSpec);

    if compute_reflectivity_uvSpec==false
        % save reflectance
        R_model(nn) = reflectanceFunction(inputSettings(2,:), ds, spec_response{1}(:,2));

    else

        R_model(nn) = ds.reflectivity.value;
    end

    % save the radiance calculation 

    L_model(:,nn) = ds.radiance.value;          % (mW/m^2/nm/sr) - 

end
toc

%% Plot the ratio of modeled reflectance with MODIS reflectance to the optical depth


figure; 
subplot(1,3,1)
plot(tau_c, R_model./R_modis, '.', 'MarkerSize', 25)


hold on; grid on; grid minor

xlabel('$\tau_c$','Interpreter','latex', 'FontSize',35)
ylabel('$\frac{R_{est}}{R_{modis}}$','Interpreter','latex', 'FontSize',35)
title('Ratio doesnt trend well with optical depth', 'FontSize',25)
set(gcf, 'Position', [0 0 2000 800])
axis square

subplot(1,3,2)
plot(R_modis, R_model./R_modis, '.', 'MarkerSize', 25)


hold on; grid on; grid minor

xlabel('Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',35)
ylabel('$\frac{R_{est}}{R_{modis}}$','Interpreter','latex', 'FontSize',35)
title('Ratio doesnt trend well with measured reflectance', 'FontSize',25)
axis square

subplot(1,3,3)
plot(R_model, R_modis, '.', 'MarkerSize', 25)


hold on; grid on; grid minor

xlabel('My Estimate of Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',25)
ylabel('MODIS Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',25)
title('Deviation from 1 to 1', 'FontSize',25)
axis square




%% Sub Plot comparing modeled and measured REFLECTANCE, re, tau_c, and a histogram

label_fontSize = 20;
% ------------------------ SUBPLOT 1 ----------------------------------
% First make a scatter plot comparing the modeled and measured
% reflectance and have the color of each marker represent the retrieved
% droplet size
% ---------------------------------------------------------------------


% create 1 to 1 line
% find the minimum and maximum values to create a y=x line

min_R_est = min(R_model);
min_R_modis = min(R_modis);

max_R_est = max(R_model);
max_R_modis = max(R_modis);

min_R_global = min([min_R_est,min_R_modis]);

max_R_global = max([max_R_est,max_R_modis]);

x_r = linspace((0.9 * min_R_global),(1.1*max_R_global),150);

% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
C = colormap(parula(length(R_model)));
% sort the droplet size values
[re_sort, idx_sort] = sort(re, 'ascend');

figure; subplot(1,3,1)
plot(x_r, x_r, 'k', 'LineWidth',1)
hold on

for nn = 1:length(re_sort)
    
   errorbar(R_model(idx_sort(nn)), R_modis(idx_sort(nn)), R_uncertainty(idx_sort(nn)),'vertical','Marker','.','Color',C(nn,:),'MarkerSize',25)

   hold on

end

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
clim([min(re), max(re)]);
% set colorbar title
cb.Label.String = '$r_e$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

hold on; grid on; grid minor
xlim([x_r(1), x_r(end)])
ylim([x_r(1), x_r(end)])
xlabel('My Estimate of Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
ylabel('MODIS Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
set(gcf, 'Position', [0 0 1800 900])
axis square


% ------------------------ SUBPLOT 2 ----------------------------------
% Next make a scatter plot comparing the modeled and measured
% reflectance and have the color of each marker represent the retrieved
% optical depth
% ---------------------------------------------------------------------



% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
C = colormap(parula(length(tau_c)));
% sort the droplet size values
[tau_c_sort, idx_sort] = sort(tau_c, 'ascend');

subplot(1,3,2)
plot(x_r, x_r, 'k', 'LineWidth',1)
hold on

for nn = 1:length(tau_c_sort)
    
   errorbar(R_model(idx_sort(nn)), R_modis(idx_sort(nn)), R_uncertainty(idx_sort(nn)),'vertical','Marker','.','Color',C(nn,:),'MarkerSize',25)

   hold on

end

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
clim([min(tau_c), max(tau_c)]);
% set colorbar title
cb.Label.String = '$\tau_c$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

hold on; grid on; grid minor
xlim([x_r(1), x_r(end)])
ylim([x_r(1), x_r(end)])
xlabel('My Estimate of Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
ylabel('MODIS Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
axis square

% inset title
if linear_cloudFraction==false
    title(['Band ', num2str(band_num), ' - X pass: ', num2str(n_pixels), ' pixels, ', num2str(num_streams), ' streams, ',...
        distribution_str, ' droplet distribution, cloud cover = ',num2str(cloud_cover), newline,...
        'used ',band_parameterization, ' band model', ' dist-variance = ',num2str(distribution_var)]);

else
    title(['Band ', num2str(band_num), ' - X pass: ', num2str(n_pixels), ' pixels', num2str(num_streams), ' streams',...
        distribution_str, ' droplet distribution, linear function to define cloud cover ', newline,...
        'used ',band_parameterization, ' band model', ' dist-variance = ',num2str(distribution_var)]);

end




% ------------------------ SUBPLOT 3 ----------------------------------
% Make a histogram of the ratio between my reflectance estiamte and the
% true measured reflectance
% ---------------------------------------------------------------------


subplot(1,3,3)

h = histogram(R_model./R_modis, 'NumBins',50);
% find the bin with the most counts - the mode
[max_count, idx_max] = max(h.BinCounts);
mode_val = h.BinEdges(idx_max) + h.BinWidth/2;
% let's also compute the weighted average
avg_weighted = sum((h.BinCounts./sum(h.BinCounts)) .* (h.BinEdges(1:end-1) + h.BinWidth/2))...
                / sum((h.BinCounts./sum(h.BinCounts)));


hold on
% Plot the modal value
% xline(mode_val,'Label',['Mode = ',num2str(round(mode_val,2))], 'FontSize', 15,...
%     'LineWidth',1)
% Plot the weighted average
xline(mode_val,'Label',['W-Mean = ',num2str(round(avg_weighted,2))], 'FontSize', 17,...
    'LineWidth',1.5, 'FontWeight','bold', 'LabelVerticalAlignment','bottom')

hold on; grid on; grid minor

xlabel('$\frac{R_{est}}{R_{modis}}$','Interpreter','latex', 'FontSize',30)
ylabel('Counts','Interpreter','latex', 'FontSize',label_fontSize)
axis square






%% Sub Plot comparing modeled and measured radiance, re, tau_c, and a histogram

label_fontSize = 20;
% ------------------------ SUBPLOT 1 ----------------------------------
% First make a scatter plot comparing the modeled and measured
% reflectance and have the color of each marker represent the retrieved
% droplet size
% ---------------------------------------------------------------------



L_int = zeros(1, n_pixels);

for nn = 1:n_pixels
    % L_model has units of W/m^2/micron/sr

    % units of W/m^2/sr
%     L_int(nn) = trapz(spec_response{1}(:,1)./1e3, L_model(:,nn) .* spec_response{1}(:,2));

    % units of mW/m^2/nm/sr (same as W/m^2/micron/sr)
%     L_int(nn) = trapz(spec_response{1}(:,1), L_model(:,nn) .* spec_response{1}(:,2))/...
%                 ((spec_response{1}(end,1) - spec_response{1}(1,1)));

    
    % Perfect Spectral Response Function
    % units of mW/m^2/nm/sr (same as W/m^2/micron/sr)
    L_int(nn) = trapz(spec_response{1}(:,1), L_model(:,nn) .* ones(numel(spec_response{1}(:,2)),1))/...
                ((spec_response{1}(end,1) - spec_response{1}(1,1)));

    % Using just the bands listed on their website for band 7 (2105 -
    % 2155)nm
    % units of mW/m^2/nm/sr (same as W/m^2/micron/sr)
%     L_int(nn) = trapz(spec_response{1}(48:end-20,1), L_model(48:end-20,nn) .* spec_response{1}(48:end-20,2))/...
%                 ((spec_response{1}(end-20,1) - spec_response{1}(48,1)));


    % Correcting the channel radiance for band 1
%     L_int(nn) = trapz(spec_response{1}(:,1), L_model(:,nn) .* spec_response{1}(:,2))/...
%                 (50);


    

%     L_int(nn) = trapz(spec_response2.wl./1e3, L_model(:,nn) .* spec_response2.val)/...
%                 ((spec_response2.wl(1) - spec_response2.wl(end))/1e3);
end





% create 1 to 1 line
% find the minimum and maximum values to create a y=x line

min_L_est = min(L_int);
min_L_modis = min(L_modis);

max_L_est = max(L_int);
max_L_modis = max(L_modis);

min_L_global = min([min_L_est,min_L_modis]);

max_L_global = max([max_L_est,max_L_modis]);

x_r = linspace((0.9 * min_L_global),(1.1*max_L_global),150);

% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
C = colormap(parula(length(L_int)));
% sort the droplet size values
[re_sort, idx_sort] = sort(re, 'ascend');

figure; subplot(1,3,1)
plot(x_r, x_r, 'k', 'LineWidth',1)
hold on

for nn = 1:length(re_sort)
    
   plot(L_int(idx_sort(nn)), L_modis(idx_sort(nn)),'Marker','.','Color',C(nn,:),'MarkerSize',25)

   hold on

end

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
clim([min(re), max(re)]);
% set colorbar title
cb.Label.String = '$r_e$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

hold on; grid on; grid minor
xlim([x_r(1), x_r(end)])
ylim([x_r(1), x_r(end)])
xlabel('Radiance Estimate $(W/m^2/\mu m/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
ylabel('MODIS Measured Radiance $(W/m^2/\mu m/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
set(gcf, 'Position', [0 0 1800 900])
axis square


% ------------------------ SUBPLOT 2 ----------------------------------
% Next make a scatter plot comparing the modeled and measured
% reflectance and have the color of each marker represent the retrieved
% optical depth
% ---------------------------------------------------------------------



% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
C = colormap(parula(length(tau_c)));
% sort the droplet size values
[tau_c_sort, idx_sort] = sort(tau_c, 'ascend');

subplot(1,3,2)
plot(x_r, x_r, 'k', 'LineWidth',1)
hold on

for nn = 1:length(tau_c_sort)
    
   plot(L_int(idx_sort(nn)), L_modis(idx_sort(nn)),'Marker','.','Color',C(nn,:),'MarkerSize',25)

   hold on

end

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
clim([min(tau_c), max(tau_c)]);
% set colorbar title
cb.Label.String = '$\tau_c$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

hold on; grid on; grid minor
xlim([x_r(1), x_r(end)])
ylim([x_r(1), x_r(end)])
xlabel('Radiance Estimate $(W/m^2/\mu m/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
ylabel('MODIS Measured Radiance $(W/m^2/\mu m/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
axis square

% inset title
if linear_cloudFraction==false
    title(['Band ', num2str(band_num), ' - X pass: ', num2str(n_pixels), ' pixels, ', num2str(num_streams), ' streams, ',...
        distribution_str, ' droplet distribution, cloud cover = ',num2str(cloud_cover), newline,...
        'used ',band_parameterization, ' band model']);

else
    title(['Band ', num2str(band_num), ' - X pass: ', num2str(n_pixels), ' pixels', num2str(num_streams), ' streams',...
        distribution_str, ' droplet distribution, linear function to define cloud cover ', newline,...
        'used ',band_parameterization, ' band model']);

end




% ------------------------ SUBPLOT 3 ----------------------------------
% Make a histogram of the ratio between my reflectance estiamte and the
% true measured reflectance
% ---------------------------------------------------------------------


subplot(1,3,3)

h = histogram(L_int./L_modis, 'NumBins',50);
% find the bin with the most counts - the mode
[max_count, idx_max] = max(h.BinCounts);
mode_val = h.BinEdges(idx_max) + h.BinWidth/2;
% let's also compute the weighted average
avg_weighted = sum((h.BinCounts./sum(h.BinCounts)) .* (h.BinEdges(1:end-1) + h.BinWidth/2))...
                / sum((h.BinCounts./sum(h.BinCounts)));


hold on
% Plot the modal value
% xline(mode_val,'Label',['Mode = ',num2str(round(mode_val,2))], 'FontSize', 15,...
%     'LineWidth',1)
% Plot the weighted average
xline(mode_val,'Label',['W-Mean = ',num2str(round(avg_weighted,2))], 'FontSize', 17,...
    'LineWidth',1.5, 'FontWeight','bold', 'LabelVerticalAlignment','bottom')

hold on; grid on; grid minor

xlabel('$\frac{L_{est}}{L_{modis}}$','Interpreter','latex', 'FontSize',30)
ylabel('Counts','Interpreter','latex', 'FontSize',label_fontSize)
axis square
