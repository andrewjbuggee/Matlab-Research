%% Does the droplet distribution value matter? Can it meaningfully change the relfectance measurement?


clear variables
% By Andrew John Buggee

%% Define the cloud parameters that will be changing during each reflectance calculation


% define the distribution varaince
distribution_var = 1:5:100;

% define the droplet effective radius
re = 5;                      % microns


%% Define the parameters of the INP file


% Define the MODIs spectral band you wish to run
% ------------------------------------------------------------------------
band_num = 7;
% ------------------------------------------------------------------------



% Define the number of streams to use in your radiative transfer model
num_streams = 16;
% ------------------------------------------------------------------------


% Define the source file
%source_file = '../data/solar_flux/lasp_TSIS1_hybrid_solar_reference_p01nm_resolution.dat';
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
day_of_year = 180;

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
percent_cloud_cover = 1;
% ------------------------------------------------------------------------


% save the modis droplet size and optical depth
tau_c = 20;

% define the geometric location of the cloud top and cloud bottom
z_topBottom = [4,3];          % km above surface



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
use_custom_mie_calcs = false;
% ------------------------------------------------------------------------

% define the type of droplet distribution
distribution_str = 'gamma';
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


% define the solar zenith angle
sza = 20;           % degree
% Define the solar azimuth measurement between values 0 and 360
phi0 = 30;         % degree
% define the viewing zenith angle
vza = 10;           % degrees;
% define the viewing azimuth angle
vaz = 0;            % deg


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
folder2save = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/testing_droplet_distribution/'];

inputName = cell(numel(re), numel(distribution_var));
outputName = cell(numel(re), numel(distribution_var));
wc_filename = cell(numel(re), numel(distribution_var));





for rr = 1:length(re)


    for vv = 1:length(distribution_var)


        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------
        % most uncertainties for the modis optical retrieval are between 2
        % and 10 percent. So lets round off all re values to the 1000th decimal
        % place


        % ------------------------------------------------------
        % --------------------VERY IMPORTANT ------------------
        % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
        % ------------------------------------------------------
        wc_filename{rr,vv} = write_wc_file(re(rr), tau_c, z_topBottom, lambda_forTau, distribution_str,...
            distribution_var(vv), vert_homogeneous_str, parameterization_str, rr*vv);
        wc_filename{rr,vv} = wc_filename{rr,vv}{1};


        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory
        inputName{rr,vv} = [num2str(floor((wavelength(2)-wavelength(1))/2 + wavelength(1))),...
            'nm_withCloudLayer_',num2str(rr*vv),'nn_',num2str(re(rr)),'effRadius_',num2str(distribution_var(vv)),...
            'variance_',atm_file(1:end-4),'.INP'];



        outputName{rr,vv} = ['OUTPUT_',inputName{rr,vv}(1:end-4)];



        % ----------------- ******************** ---------------------
        % ------------------ Write the INP File --------------------
        % ----------------- ******************** ---------------------

        % Open the old file for writing
        fileID = fopen([folder2save,inputName{rr,vv}], 'w');

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
            fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{rr,vv}], ' ', '# Location of water cloud file');

            % Define the percentage of horizontal cloud cover
            % This is a number between 0 and 1
            % ------------------------------------------------
            formatSpec = '%s %f %5s %s \n';
            fprintf(fileID, formatSpec,'cloudcover wc', percent_cloud_cover, ' ', '# Cloud cover percentage');


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
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');

        % Define the solar azimuth angle
        % -------------------------------------------------------
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'phi0', phi0, ' ', '# Solar azimuth angle');

        % Define the cosine of the zenith viewing angle
        % ------------------------------------------------
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


end



%% Run INP files

%L_model = zeros(numel(spec_response{1}(:,1)), length(idx));
R_model = zeros(numel(re), numel(distribution_var));

legend_str = cel(1, numel(re));


tic
for rr = 1:numel(re)

    legend_str{rr} = ['$r_e = $', num2str(re(rr)), '$mu m$'];

    parfor vv = 1:numel(distribution_var)


        % compute INP file
        [inputSettings] = runUVSPEC(folder2save,inputName{rr,vv},outputName{rr,vv});

        % read .OUT file
        [ds,~,~] = readUVSPEC(folder2save,outputName{rr,vv},inputSettings(2,:), compute_reflectivity_uvSpec);

        if compute_reflectivity_uvSpec==false
            % save reflectance
            R_model(rr,vv) = reflectanceFunction(inputSettings(2,:), ds, spec_response{1}(:,2));

        else

            R_model(rr,vv) = ds.reflectivity.value;
        end

        % save the radiance calculation

        %L_model(:,rr) = ds.radiance.value;          % (mW/m^2/nm/sr) -
    end

end
toc

%% Plot the reflectance versus the droplet distribution for each effective radius


figure;
plot(distribution_var, R_model, '.-', 'MarkerSize', 25, 'LineWidth',1.5)

grid on; grid minor

xlabel('$\sigma$','Interpreter','latex', 'FontSize',35)
ylabel('$R$ ($1/sr$)','Interpreter','latex', 'FontSize',35)
title('Reflectance versus Distribution variance', 'Interpreter','latex', 'FontSize',25)

% write the legend, where each curve has a different effective raidus
legend(legend_str, 'Interpreter', 'latex', 'location','best', 'Fontsize',20)

% set the plot size
set(gcf, 'Position', [0 0 1200 600])


