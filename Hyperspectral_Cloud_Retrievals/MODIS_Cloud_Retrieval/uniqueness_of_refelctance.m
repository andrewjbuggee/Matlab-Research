%% Is there always a unique r_top, r_bot and tau_c using the first 7 MODIS spectral channels?




clear variables
% By Andrew John Buggee

%% Define the cloud parameters that will be changing during each reflectance calculation


r_top = 6:12;       % microns
r_bot = 4:10;        % microns
tau_c = 5:5:35;

% r_top = 6:2:12;       % microns
% r_bot = 4:2:10;        % microns
% tau_c = 5:5:35;


% r_top = 9.03;       % microns
% r_bot = 9.03;        % microns
% tau_c = 6.34;



%% Want to use real MODIS geometry inputs?

% Load modis data and create input structure

% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the MODIS data folder path

    modisPath = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/';



    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the MODIS data folder path

    modisPath = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];

    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];

end


% -------------------------------------
% ------- PICK MODIS DATA SET  --------
% -------------------------------------

% ----- November 9th at decimal time 0.611 (14:40) -----
%modisFolder = '2008_11_09/';

% ----- November 11th at decimal time 0.604 (14:30) -----
%modisFolder = '2008_11_11_1430/';

% ----- November 11th at decimal time 0.784 (18:50) -----
modisFolder = '2008_11_11_1850/';



[modis,L1B_fileName] = retrieveMODIS_data([modisPath, modisFolder]);


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008
%modis_idx = 348140;    % for 9 nov 2008 - pixel overlapping with VOCALS
modis_idx = 1278681;        % for 11 Nov 2008 @ 18:50 - pixel overlapping with VOCALS     
%modis_idx = 110293;        % for 11 Nove 2008 @ 1430 - pixel overlapping with VOCALS

%% Grab the MODIS reflectances for the pixel used
[r,c] = ind2sub(size(modis.EV1km.reflectance), modis_idx);
R_modis = zeros(1, size(modis.EV1km.reflectance,3));
R_uncert_modis = zeros(1, size(modis.EV1km.reflectance, 3));

for bb = 1:size(modis.EV1km.reflectance, 3)

    % ****** DID YOU USE REFLECTANCE_4MODIS? ******
    % If not you need to divide the MODIS reflectance by cos(sza)
    R_modis(bb) = modis.EV1km.reflectance(r,c,bb);
    R_uncert_modis(bb) = R_modis(bb) * 0.01*double(modis.EV1km.reflectanceUncert(r,c,bb)); % converted from percentage to reflectance
end


%% Define the parameters of the INP file


% Define the MODIS spectral band you wish to run
% ------------------------------------------------------------------------
band_num = 1:7;
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
for ww = 1:length(band_num)
    wavelength(ww,:) = [spec_response{ww}(1,1), spec_response{ww}(end,1)];
end

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
day_of_year = str2double(L1B_fileName{1}(15:17));

% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
use_MODIS_cloudTopHeight = true;
cloud_depth = 500;                % meters

if use_MODIS_cloudTopHeight==true
    z_topBottom = [modis.cloud.topHeight(modis_idx), modis.cloud.topHeight(modis_idx) - cloud_depth]./1e3; %km above surface

else
    % define the geometric location of the cloud top and cloud bottom
    z_topBottom = [2.5, 2];          % km above surface

end

% Water Cloud depth
H = z_topBottom(1) - z_topBottom(2);                                % km - geometric thickness of cloud

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
linear_cloudFraction = false;
% if false, define the cloud cover percentage
percent_cloud_cover = 1;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
use_custom_mie_calcs = false;
% ------------------------------------------------------------------------

% define the type of droplet distribution
distribution_str = 'gamma';
% define whether this is a vertically homogenous cloud or not
vert_homogeneous_str = 'vert-non-homogeneous';
% define how liquid water content will be computed
parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
band1 = modisBands(1);
lambda_forTau = band1(1);            % nm


% ------------------------------------------------------------------------
% ------------------- Radius Profile attributes --------------------------
% ------------------------------------------------------------------------

profile_type = 'adiabatic'; % type of water droplet profile

n_layers = 10;                          % number of layers to model within cloud

z = linspace(z_topBottom(1), z_topBottom(2), n_layers);        % km - altitude above ground vector

indVar = 'altitude';                    % string that tells the code which independent variable we used

dist_var = linspace(20,20,n_layers);              % distribution variance
% ------------------------------------------------------------------------



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
sza = double(modis.solar.zenith(modis_idx));           % degree

% Define the solar azimuth measurement between values 0 and 360
% this is how we map MODIS azimuth of the sun to the LibRadTran measurement
phi0 = double(modis.solar.azimuth(modis_idx) + 180);         % degree

% define the viewing zenith angle
vza = double(modis.sensor.zenith(modis_idx)); % values are in degrees;                        % degree

% define the viewing azimuth angle
% define the viewing azimuth angle
% to properly map the azimuth angle onto the reference plane used by
% libRadTran, we need an if statement
if modis.sensor.azimuth(modis_idx)<0
    vaz = 360 + modis.sensor.azimuth(modis_idx);
else
    vaz = modis.sensor.azimuth(modis_idx);
end


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








%% Write each INP file and Calculate Reflectance for MODIS

inputName = cell(length(r_top), length(r_bot), length(tau_c), length(band_num));
outputName = cell(length(r_top), length(r_bot), length(tau_c), length(band_num));
wc_filename = cell(length(r_top), length(r_bot), length(tau_c), length(band_num));


R_model = zeros(length(r_top), length(r_bot), length(tau_c), length(band_num));


tic
for rt = 1:length(r_top)


    for rb = 1:length(r_bot)


        for tc = 1:length(tau_c)


            disp(['Iteration: [rt, rb, tc] = [', [num2str(rt),', ', num2str(rb), ', ', num2str(tc)], ']...', newline])


            parfor ww = 1:length(band_num)



                % -----------------------------------
                % ---- Write a Water Cloud file! ----
                % -----------------------------------
                % most uncertainties for the modis optical retrieval are between 2
                % and 10 percent. So lets round off all re values to the 1000th decimal
                % place

                re = create_droplet_profile2([r_top(rt), r_bot(rb)], z, indVar, profile_type);     % microns - effective radius vector


                % ------------------------------------------------------
                % --------------------VERY IMPORTANT ------------------
                % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
                % ------------------------------------------------------
                wc_filename{rt,rb,tc,ww} = write_wc_file(re, tau_c(tc), z_topBottom, lambda_forTau, distribution_str,...
                    dist_var, vert_homogeneous_str, parameterization_str, ww);
                wc_filename{rt,rb,tc,ww} = wc_filename{rt,rb,tc,ww}{1};


                % ------------------------------------------------
                % ---- Define the input and output filenames! ----
                % ------------------------------------------------
                % input_names need a unique identifier. Let's give them the nn value so
                % they can be traced, and are writen over in memory
                %                 inputName{rt,rb,tc,ww} = [num2str(floor((wavelength(ww,2)-wavelength(ww,1))/2 + wavelength(ww,1))),...
                %                     'nm_withCloudLayer_',num2str(r_top(rt)),'rTop_',num2str(r_bot(rb)),...
                %                     'rBot_', num2str(tau_c(tc)), 'tauC_' ,atm_file(1:end-4),'.INP'];

                inputName{rt,rb,tc,ww} = [num2str(floor((wavelength(ww,2)-wavelength(ww,1))/2 + wavelength(ww,1))),...
                    'reflectance', atm_file(1:end-4),'.INP'];



                outputName{rt,rb,tc,ww} = ['OUTPUT_',inputName{rt,rb,tc,ww}(1:end-4)];



                % ----------------- ******************** ---------------------
                % ------------------ Write the INP File --------------------
                % ----------------- ******************** ---------------------

                % Open the old file for writing
                fileID = fopen([folder2save,inputName{rt,rb,tc,ww}], 'w');

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
                    fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{rt,rb,tc,ww}], ' ', '# Location of water cloud file');

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
                fprintf(fileID, formatSpec,'wavelength', wavelength(ww,1), wavelength(ww,2), ' ', '# Wavelength range');




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




                % ----------------------------------------------------
                % --------------- RUN RADIATIVE TRANSFER -------------
                % ----------------------------------------------------


                % compute INP file
                [inputSettings] = runUVSPEC(folder2save,inputName{rt,rb, tc, ww},outputName{rt,rb, tc, ww});

                % read .OUT file
                [ds,~,~] = readUVSPEC(folder2save,outputName{rt,rb, tc, ww},inputSettings(2,:), compute_reflectivity_uvSpec);

                if compute_reflectivity_uvSpec==false
                    % compute reflectance in the MODIS style (without
                    % dividing by cos(sza)
                    R_model(rt,rb, tc, ww) = reflectanceFunction_4modis(inputSettings(2,:), ds, spec_response{ww}(:,2));

                else

                    R_model(rt,rb, tc, ww) = ds.reflectivity.value;
                end




            end


        end


    end


end


% ----------------------------------------------
% ---------- SAVE REFLECTANCE OUTPUT! ----------
% ----------------------------------------------

rev = 1;

if strcmp(whatComputer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/Reflectance_Uniqueness/'];



elseif strcmp(whatComputer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/Reflectance_Uniqueness/'];


elseif strcmp(whatComputer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


filename = [folderpath,'reflectance_calcs_MODIS-data-from-',modisFolder(1:end-1),...
    '_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];

while isfile(filename)
    rev = rev+1;
    filename = [folderpath,'reflectance_calcs_MODIS-data-from-',modisFolder(1:end-1),...
    '_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];
end

save(filename,"r_top", "r_bot", "tau_c", "wavelength", "R_model", "modisFolder", 'modis_idx');

toc


%% Subplots of reflectance across different optical depths for a single wavelength

wave_len_idx = 1;

% find min and max values of reflectance for this wavelength
[minR, ~] = min(R_model(:,:,:,wave_len_idx), [], 'all');
[maxR, ~] = max(R_model(:,:,:,wave_len_idx), [], 'all');

% set colorbar limits for each subplot
color_lim = [minR, maxR];



figure;
for tc = 1:length(tau_c)

    subplot(2,4,tc); imagesc(r_bot, r_top, R_model(:,:,tc, wave_len_idx));

    if tc==1
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
    end

    title(['$\tau_c$ = ', num2str(tau_c(tc))], 'Interpreter', 'latex')

    if tc==4
        cb = colorbar;
        set(get(cb, 'label'), 'string', 'Reflectance $(1/sr)$','Interpreter','latex', 'Fontsize',28)
        %cb.Limits = color_lim;

    else
        colorbar
        %clim(color_lim)
    end

    %clim([minR, maxR])



end


% set the plot size
set(gcf, 'Position', [0 0 1200 600])

% Create textbox
annotation('textbox',...
    [0.322666666666667 0.00166666666666667 0.358166666666667 0.0916666666666668],...
    'VerticalAlignment','middle',...
    'String',['$\lambda$ = ',num2str(round(mean(wavelength(wave_len_idx, :)))), ' $nm$'],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',35,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');





%% Make lineplots of reflectance versus optical depth at a constant wavelength


wave_len_idx = 1;



figure;
for rt = 1:length(r_top)
    for rb = 1:length(r_bot)


        plot(tau_c, reshape(R_model(rt,rb,:, wave_len_idx), 1, []));

        hold on

    end
end

xlabel('$\tau_{c}$', 'Interpreter', 'latex')
ylabel('Reflectance $(1/sr)$', 'Interpreter', 'latex')

title(['Reflectance $\lambda$ = ', num2str(round(mean(wavelength(wave_len_idx, :)))), ' $nm$'], 'Interpreter', 'latex')



% set the plot size
set(gcf, 'Position', [0 0 1200 600])
grid on; grid minor




%% Subplots of reflectance across different wavelengths for a single optical depth

tau_idx = 1;

% find min and max values of reflectance for this wavelength
[minR, ~] = min(R_model(:,:,tau_idx,:), [], 'all');
[maxR, ~] = max(R_model(:,:,tau_idx,:), [], 'all');



figure;
for ww = 1:length(band_num)

    subplot(2,4,ww); imagesc(r_bot, r_top, R_model(:,:,tau_idx, ww));

    if ww==1
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
    end

    title(['$\lambda$ = ', num2str(round(mean(wavelength(ww, :)))), ' ($nm$)'],...
        'Interpreter', 'latex')

    if ww==4
        cb = colorbar;
        set(get(cb, 'label'), 'string', 'Reflectance $(1/sr)$','Interpreter','latex', 'Fontsize',28)


    else
        colorbar
    end

    %clim([minR, maxR])



end


% set the plot size
set(gcf, 'Position', [0 0 1200 600])

% Create textbox
annotation('textbox',...
    [0.322666666666667 0.00166666666666667 0.358166666666667 0.0916666666666668],...
    'VerticalAlignment','middle',...
    'String',['$\tau_c$ = ',num2str(tau_c(tau_idx))],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',35,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');







%% Subplots of reflectance across different optical depths for all wavelengths




for ww = 1:length(band_num)

    % find min and max values of reflectance for this wavelength
    [minR, ~] = min(R_model(:,:,:,ww), [], 'all');
    [maxR, ~] = max(R_model(:,:,:,ww), [], 'all');

    figure;
    for ww = 1:length(tau_c)

        subplot(2,4,ww); imagesc(r_bot, r_top, R_model(:,:,ww,ww));

        if ww==1
            xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
            ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
        end

        title(['$\tau_c$ = ', num2str(tau_c(ww))], 'Interpreter', 'latex')

        if ww==3
            cb = colorbar;
            set(get(cb, 'label'), 'string', 'Reflectance $(1/sr)$','Interpreter','latex', 'Fontsize',28)


        else
            colorbar
        end

        %clim([minR, maxR])



    end


    % set the plot size
    set(gcf, 'Position', [0 0 1200 600])

    % Create textbox
    annotation('textbox',...
        [0.322666666666667 0.00166666666666667 0.358166666666667 0.0916666666666668],...
        'VerticalAlignment','middle',...
        'String',['$\lambda$ = ',num2str(round(mean(wavelength(ww, :)))), ' $nm$'],...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',35,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor','none');

end



%% Plot a 3D surf figure with r_bot r_top and reflectance as the variables. Plot
% for a single optical depth at a single wavelength

wave_len_idx = 1;
tau_c_idx = 4;

% ADD A PLANE OF CONSTANT REFLECTANCE TO SHOW DIFFERENT SOLUTIONS
R_constant = 0.65;


r_top_mat = repmat(r_top', 1, length(r_bot));
r_bot_mat = repmat(r_bot, length(r_top), 1);
R_mat = reshape(R_model(:,:, tau_c_idx, wave_len_idx), size(r_top_mat));

% Set the color of the

% plot the visible band first
f = figure;

% plot the first band
surf(r_bot_mat, r_top_mat, R_mat);

% interpolate between points to smooth the surface
%shading interp

% set up plot stuff
grid on; grid minor
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
zlabel('Reflectance ($1/sr$)', 'Interpreter','latex')

title(['$\tau_c$ = ', num2str(tau_c(tau_c_idx)), '   $\lambda$ = ',...
    num2str(round(mean(wavelength(wave_len_idx, :)))), ' $nm$'],Interpreter='latex')


% ADD A PLANE OF CONSTANT REFLECTANCE TO SHOW DIFFERENT SOLUTIONS
hold on
surf(r_bot_mat, r_top_mat, repmat(R_constant, length(r_top), length(r_bot)));



% set the plot size
set(gcf, 'Position', [0 0 1200 600])


%% Plot a 3D surf figure with r_bot r_top and optical depth as the variables. Plot
% for multiple optical depths at a single wavelength

wave_len_idx = 7;

r_top_mat = repmat(r_top', 1, length(r_bot));
r_bot_mat = repmat(r_bot, length(r_top), 1);



% designate the color of each surface as the optical depth
% We span the colormap only by the number of unique data values
C = parula(length(tau_c));


figure;
for ww = 1:length(tau_c)

    R_mat = reshape(R_model(:,:, ww, wave_len_idx), size(r_top_mat));


    % plot the first band
    surf(r_bot_mat, r_top_mat, R_mat, repmat(tau_c(ww), length(r_top), length(r_bot)));

    hold on


end


% set colorbar label
cb = colorbar;
set(get(cb, 'label'), 'string', '$\tau_c$','Interpreter','latex', 'Fontsize',28)

% interpolate between points to smooth the surface
%shading interp

% set up plot stuff
grid on; grid minor
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
zlabel('Reflectance ($1/sr$)', 'Interpreter','latex')

title(['$\lambda$ = ',num2str(round(mean(wavelength(wave_len_idx, :)))), ' $nm$'],...
    Interpreter='latex')



% set the plot size
set(gcf, 'Position', [0 0 1200 600])



%% Plot a 3D surf figure with r_bot r_top and wavelength as the variables. Plot
% for multiple wavelengths as different planes for a single optical depth

tau_idx = 1;

r_top_mat = repmat(r_top', 1, length(r_bot));
r_bot_mat = repmat(r_bot, length(r_top), 1);



% designate the color of each surface as the optical depth
% We span the colormap only by the number of unique data values
C = parula(length(tau_c));


figure;
for ww = 1:length(band_num)

    R_mat = reshape(R_model(:,:, tau_idx, ww), size(r_top_mat));


    % plot the first band
    wl = round(mean(wavelength(ww, :)));

    surf(r_bot_mat, r_top_mat, R_mat, repmat(wl, length(r_top), length(r_bot)));

    hold on


end


% set colorbar label
cb = colorbar;
set(get(cb, 'label'), 'string', '$\lambda$ ($nm$)','Interpreter','latex', 'Fontsize',28)

% interpolate between points to smooth the surface
%shading interp

% set up plot stuff
grid on; grid minor
xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
zlabel('Reflectance ($1/sr$)', 'Interpreter','latex')

title(['$\tau_c$ = ',num2str(tau_c(tau_idx))],...
    Interpreter='latex')



% set the plot size
set(gcf, 'Position', [0 0 1200 600])



%% Can I determine uniqueness by rounding the reflectance computed to the measurement uncertainty of MODIS?
% Then, see how redundant certain states are for all 7 wavelenghts
% OR, should I instead ask, how many sets of measurements I computed are
% within the uncertainty of the MODIS measurements? That is maybe a better
% estimate of uniqueness, because MODIS claims a specific measurement but
% then include confidence intervals. They can say for certain that there
% measurement lies somewhere within that uncertainty interval.

% According to MODIS data, the reflectance uncertainty is betwen 1.5 and
% 2.5%

% Let's assume the reflectance uncertainty is 1%
% Let's truncate the reflectance data to hundreds decimal point

% For our 7 spectral measurements, how any different states of (r_top,
% r_bot, tau_c) lead to the same 7 measurements

% Let's reshape R_model_round to be 7 rows and N number of columns
% corresponding the N number of unique states

%R_model_round_states = zeros(length(band_num), length(r_top)*length(r_bot)*length(tau_c));
R_model_round_states = [];

% Grab the MODIS reflectances for the pixel used
[r,c] = ind2sub(size(modis.EV1km.reflectance(:,:,1)), modis_idx);
R_modis = zeros(1, length(band_num));
R_uncert_modis = zeros(1, length(band_num));

for bb = 1:length(band_num)

    % ****** DID YOU USE REFLECTANCE_4MODIS? ******
    % If not you need to divide the MODIS reflectance by cos(sza)
    R_modis(bb) = modis.EV1km.reflectance(r,c,bb)/cosd(double(modis.solar.zenith(r,c)));
    R_uncert_modis(bb) = R_modis(bb) * 0.01*double(modis.EV1km.reflectanceUncert(r,c,bb)); % converted from percentage to reflectance
end

redundant_states = [];


for rt = 1:length(r_top)


    for rb = 1:length(r_bot)


        for tc = 1:length(tau_c)

            % Round reflectance calculations to the nearest hundreth
            % decimal place
            R_model_round_states = [R_model_round_states; reshape(round(R_model(rt,rb,tc,:),2), 1, [])];

            % Check to see if the reflectance computed by the model is
            % within the listed uncertainty for MODIS
            %             redundant_states = [redundant_states, all(abs(R_modis - reshape(R_model(rt,rb,tc,:), 1, [])) <= R_uncert_modis)];
            redundant_states = [redundant_states; abs(R_modis - reshape(R_model(rt,rb,tc,:), 1, [])) <= R_uncert_modis];


        end
    end
end


% Find the number of unique measurements
[R_model_unique, idx_original, idx_unique] = unique(R_model_round_states, 'rows');

% print the percentage of redundant states
disp([newline, num2str(100*(size(R_model_round_states,1) - size(R_model_unique,1))/size(R_model_round_states,1)),...
    '% of retireved states are redundant', newline])

% find the logical array of unique values
idx_unique_logical = ismember((1:size(R_model_round_states,1)), idx_original);



%% Interpolate the reflectance calculations on a finer grid

tau_c_fine = tau_c(1):tau_c(end);

R_model_fine = zeros(length(r_top), length(r_bot), length(tau_c_fine), length(band_num));



for rt = 1:length(r_top)


    for rb = 1:length(r_bot)


        for wl = 1:length(band_num)

            new_reflectance = interp1(tau_c, reshape(R_model(rt, rb, :, wl), 1, []), tau_c_fine);
            R_model_fine(rt,rb,:,wl) = reshape(new_reflectance, 1,1,[],1);


        end

    end

end



%% 3D Interpolate the reflectance calculations on a finer grid and compare with MODIS measurement


% Meshgrid is defined on x,y,z space, not row, column, depth space
% In 3D space, z = row, x = column, y = depth
[R_bot, R_top, Tau_c] = meshgrid(r_bot, r_top, tau_c);

% Create the new fine grid to interpolate on
% define the discrete step length of each variable
d_r_top = 0.05;      % microns
d_r_bot = 0.05;      % microns
d_tau_c = 0.1;

r_top_fine = r_top(1):d_r_top:r_top(end);
r_bot_fine = r_bot(1):d_r_bot:r_bot(end);
tau_c_fine = tau_c(1):d_tau_c:tau_c(end);

[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);

R_model_fine = zeros(length(r_top_fine), length(r_bot_fine), length(tau_c_fine), size(R_model,4));


% ***** IF USING REFLECTANCE CALCS USING STANDARD REFLECTANCE DEFINITION *****
% MULTIPLY EACH VALUE BY THE AIRMASS: COS(SZA)
warning([newline, 'Are you using reflectance_calcs_standardReflectance_with_mu0_9-nov-2008-data_15-Nov-2023.mat?',...
    newline, 'Make sure you multiply all reflectances by cos(sza)!', newline])


for wl = 1:size(R_model,4)

    R_model_fine(:,:,:,wl) = interp3(R_bot, R_top, Tau_c, R_model(:, :, :, wl),...
        R_bot_fine, R_top_fine, Tau_c_fine);


end



% Using the new fine grid, calculate how many sets of measurements are
% within the MODIS measurement and it's uncertainty

% Grab the MODIS reflectances for the pixel used
[r,c] = ind2sub(size(modis.EV1km.reflectance(:,:,1)), modis_idx);
R_modis = zeros(1, length(band_num));
R_uncert_modis = zeros(1, length(band_num));

for bb = 1:length(band_num)

    % ****** DID YOU USE REFLECTANCE_4MODIS? ******
    % If not you need to divide the MODIS reflectance by cos(sza)
    R_modis(bb) = modis.EV1km.reflectance(r,c,bb)/cosd(double(modis.solar.zenith(r,c)));
    R_uncert_modis(bb) = R_modis(bb) * 0.01*double(modis.EV1km.reflectanceUncert(r,c,bb)); % converted from percentage to reflectance
end


redundant_states = zeros(size(R_model_fine,1), size(R_model_fine,2), size(R_model_fine,3));


for rt = 1:size(R_model_fine,1)


    for rb = 1:size(R_model_fine,2)


        for tc = 1:size(R_model_fine,3)

            % Check to see if the reflectance computed by the model is
            % within the listed uncertainty for MODIS
            redundant_states(rt,rb,tc) = all(abs(R_modis - reshape(R_model_fine(rt,rb,tc,:), 1, [])) <= R_uncert_modis);
            %redundant_states = [redundant_states; abs(R_modis - reshape(R_model_fine(rt,rb,tc,:), 1, [])) <= R_uncert_modis];


        end
    end
end

% print the percentage of redundant states
disp([newline, 'There are ', num2str(sum(redundant_states, 'all')), ' sets of measurements within',...
    ' the MODIS measurement and uncertainty.', newline])

[r_redun, c_redun, d_redun] = ind2sub(size(redundant_states), find(redundant_states));

% ----- Plot all the redundant states on a scatter plot -----

% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
r_top_redundant = zeros(length(r_redun), 1);
r_bot_redundant = zeros(length(r_redun), 1);
tau_c_redundant = zeros(length(r_redun), 1);

for nn = 1:length(r_redun)
    r_top_redundant(nn) = R_top_fine(r_redun(nn), c_redun(nn), d_redun(nn));
    r_bot_redundant(nn) = R_bot_fine(r_redun(nn), c_redun(nn), d_redun(nn));
    tau_c_redundant(nn) = Tau_c_fine(r_redun(nn), c_redun(nn), d_redun(nn));

end

C = colormap(parula(length(tau_c_redundant)));
% sort the droplet size values
[tau_c_redundant_sort, idx_sort] = sort(tau_c_redundant, 'ascend');

figure;

for nn = 1:length(tau_c_redundant_sort)

    plot(r_bot_redundant(idx_sort(nn)), r_top_redundant(idx_sort(nn)),'Marker','.','Color',C(nn,:),'MarkerSize',25)

    hold on

end

% Plot a one-to-one line to show the boundary for homogenous profiles
[min_radius_val, ~] = min([r_top_redundant; r_bot_redundant]);
[max_radius_val, ~] = max([r_top_redundant; r_bot_redundant]);
plot([min_radius_val, max_radius_val], [min_radius_val, max_radius_val], 'k-', ...
    'linewidth', 1)

xlim([min(r_bot_redundant), max(r_bot_redundant)])
ylim([min(r_top_redundant), max(r_top_redundant)])

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
clim([min(tau_c_redundant_sort), max(tau_c_redundant_sort)]);
% set colorbar title
cb.Label.String = '$\tau_c$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

grid on; grid minor
xlabel('$r_{bot}$ $(\mu m)$','Interpreter','latex')
ylabel('$r_{top}$ $(\mu m)$','Interpreter','latex')
set(gcf, 'Position', [0 0 1000 500])

% Set title as the resolution of each variable
title(['$\triangle r_{top} = $', num2str(d_r_top), ' $\mu m$',...
    '    $\triangle r_{bot} = $', num2str(d_r_bot), ' $\mu m$',...
    '    $\triangle \tau_{c} = $', num2str(d_tau_c)], ...
    'Fontsize', 25, 'Interpreter', 'latex')

% ---- plot the 2D space or r-top and r-bot showing area of redundancy ---

% for every r-top, what is the largest and smallest r_bot that results in a
% measurement within the MODIS measurement and uncertainty?
[r_top_unique, idx_original] = unique(r_top_redundant);

top_boundary = zeros(length(r_top_unique), 2);
bottom_boundary = zeros(length(r_top_unique), 2);
tau_c_points_top = cell(length(r_top_unique),1);
tau_c_points_top_minVal = zeros(length(r_top_unique), 1);
tau_c_points_bottom = cell(length(r_top_unique),1);
tau_c_points_bottom_minVal = zeros(length(r_top_unique), 1);

for nn = 1:length(r_top_unique)

    % find set of r_bottom values for each unique r_top
    r_bottom_set = r_bot_redundant(r_top_redundant==r_top_unique(nn));

    % If there is more than 1 value in the set, find the highest and lowest
    % value. These make up the upper and lower boundaries, respectively
    % The locations should be (r_bot,r_top)
    bottom_boundary(nn,:) = [min(r_bottom_set), r_top_unique(nn)];
    top_boundary(nn,:) = [max(r_bottom_set), r_top_unique(nn)];

    % grab the optical depth for each of these points
    tau_c_points_bottom{nn} = tau_c_redundant(r_top_redundant==r_top_unique(nn) & r_bot_redundant==min(r_bottom_set));
    tau_c_points_bottom_minVal(nn) = min(tau_c_points_bottom{nn});

    tau_c_points_top{nn} = tau_c_redundant(r_top_redundant==r_top_unique(nn) & r_bot_redundant==max(r_bottom_set));
    tau_c_points_top_minVal(nn) = min(tau_c_points_top{nn});

end

% Create a polyshape using the vertices above
figure;
p = patch([bottom_boundary(:,1); flipud(top_boundary(:,1))], ...
    [bottom_boundary(:,2); flipud(top_boundary(:,2))],...
    [tau_c_points_bottom_minVal; tau_c_points_top_minVal], 'FaceColor','interp');


xlim([r_bot(1), r_bot(end)])
ylim([r_top(1), r_top(end)])

% create a one-to-one line to delineate between profiles where r-top>r-bot
% and those where this isn't true
hold on;
plot([r_top(1), r_bot(end)], [r_top(1), r_bot(end)], 'k-', 'linewidth', 2)

p.EdgeAlpha = 0;

cb = colorbar;
% set colorbar title
cb.Label.String = '$\tau_c$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

% create legend
legend('Region of Redundant Solutions', 'Vertically homogenous droplet profile',...
    'Interpreter', 'latex', 'Fontsize', 18, 'Location', 'best')

title('State space where adiabatic profiles lead to reflectances within MODIS uncertainty', ...
    'Fontsize', 23)

grid on; grid minor

xlabel('$r_{bot}$ $(\mu m)$','Interpreter','latex')
ylabel('$r_{top}$ $(\mu m)$','Interpreter','latex')
set(gcf, 'Position', [0 0 1200 600])

%% Subplots of reflectance across different wavelengths for a single optical depth

tau_idx = tau_c_fine==8.5;

% find min and max values of reflectance for this wavelength
[minR, ~] = min(R_model_fine(:,:,tau_idx,:), [], 'all');
[maxR, ~] = max(R_model_fine(:,:,tau_idx,:), [], 'all');



figure;
for ww = 1:length(band_num)

    subplot(2,4,ww); imagesc(r_bot, r_top, R_model_fine(:,:,tau_idx, ww));

    if ww==1
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
    end

    title(['$\lambda$ = ', num2str(round(mean(wavelength(ww, :)))), ' ($nm$)'],...
        'Interpreter', 'latex')

    if ww==4
        cb = colorbar;
        set(get(cb, 'label'), 'string', 'Reflectance $(1/sr)$','Interpreter','latex', 'Fontsize',28)


    else
        colorbar
    end

    %clim([minR, maxR])



end


% set the plot size
set(gcf, 'Position', [0 0 1200 600])

% Create textbox
annotation('textbox',...
    [0.322666666666667 0.00166666666666667 0.358166666666667 0.0916666666666668],...
    'VerticalAlignment','middle',...
    'String',['$\tau_c$ = ',num2str(tau_c_fine(tau_idx))],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',35,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

