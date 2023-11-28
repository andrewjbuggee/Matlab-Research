%% Is there always a unique r_top, r_bot and tau_c using the first 7 MODIS spectral channels?




clear variables
% By Andrew John Buggee

%% Define the cloud parameters that will be changing during each reflectance calculation


r_top = 6:0.5:12;       % microns
r_bot = 4:0.5:10;        % microns
tau_c = 5:5:35;



%% Want to use real MODIS geometry inputs?

% Load modis data and create input structure

% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the MODIS data folder path

    modisPath = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/';


elseif strcmp(computer_name,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the MODIS data folder path

    modisPath = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];

end


% -------------------------------------
% ------- PICK MODIS DATA SET  --------
% -------------------------------------

% ----- November 9th at decimal time 0.611 (14:40) -----
modisFolder = '2008_11_09/';

% ----- November 11th at decimal time 0.604 (14:30) -----
%modisFolder = '2008_11_11_1430/';

% ----- November 11th at decimal time 0.784 (18:50) -----
%modisFolder = '2008_11_11_1850/';



[modis,L1B_fileName] = retrieveMODIS_data([modisPath, modisFolder]);


% Define an index to use
modis_idx = 110292;
%modis_idx = 348140;

%% Define the parameters of the INP file


% Define the MODIs spectral band you wish to run
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

if use_MODIS_cloudTopHeight==true
    z_topBottom = [modis.cloud.topHeight(modis_idx), modis.cloud.topHeight(modis_idx) - 500]./1e3; %km above surface

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

dist_var = linspace(7,7,n_layers);              % distribution variance
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
sza = modis.solar.zenith(modis_idx);           % degree

% Define the solar azimuth measurement between values 0 and 360
% this is how we map MODIS azimuth of the sun to the LibRadTran measurement
phi0 = modis.solar.azimuth(modis_idx) + 180;         % degree

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








%% Write each INP file using various MODIS values

% Define the folder path where all .INP files will be saved
folder2save = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];

inputName = cell(length(r_top), length(r_bot), length(tau_c), length(band_num));
outputName = cell(length(r_top), length(r_bot), length(tau_c), length(band_num));
wc_filename = cell(length(r_top), length(r_bot), length(tau_c), length(band_num));


R_model = zeros(length(r_top), length(r_bot), length(tau_c), length(band_num));


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
                    % save reflectance
                    R_model(rt,rb, tc, ww) = reflectanceFunction(inputSettings(2,:), ds, spec_response{ww}(:,2));

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

save(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/',...
    'Reflectance_Uniqueness/reflectance_calcs_',char(datetime("today")),'.mat'],...
    "r_top", "r_bot", "tau_c", "wavelength", "R_model");


%% Run INP files

%L_model = zeros(numel(spec_response{1}(:,1)), length(idx));
R_model = zeros(length(r_top), length(r_bot), length(tau_c), length(band_num));

%legend_str = cel(1, numel(re));


tic
for rt = 1:length(r_top)


    for rb = 1:length(r_bot)


        for tc = 1:length(tau_c)


            parfor ww = 1:length(band_num)

                %legend_str{rt} = ['$r_e = $', num2str(re(rt)), '$mu m$'];



                % compute INP file
                [inputSettings] = runUVSPEC(folder2save,inputName{rt,rb, tc, ww},outputName{rt,rb, tc, ww});

                % read .OUT file
                [ds,~,~] = readUVSPEC(folder2save,outputName{rt,rb, tc, ww},inputSettings(2,:), compute_reflectivity_uvSpec);

                if compute_reflectivity_uvSpec==false
                    % save reflectance
                    R_model(rt,rb, tc, ww) = reflectanceFunction(inputSettings(2,:), ds, spec_response{ww}(:,2));

                else

                    R_model(rt,rb, tc, ww) = ds.reflectivity.value;
                end

                % save the radiance calculation

                %L_model(:,rr) = ds.radiance.value;          % (mW/m^2/nm/sr) -

            end

        end
    end

end
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




%% Subplots of reflectance across different wavelengths for a single optical depth

tau_idx = 4;

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


for rt = 1:length(r_top)


    for rb = 1:length(r_bot)


        for tc = 1:length(tau_c)


            R_model_round_states = [R_model_round_states; reshape(round(R_model(rt,rb,tc,:),2), 1, [])];

            
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



