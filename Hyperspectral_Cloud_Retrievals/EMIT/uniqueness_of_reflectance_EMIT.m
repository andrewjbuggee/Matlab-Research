%% Is there always a unique r_top, r_bot and tau_c using the N spectral channels from EMIT?




clear variables
% By Andrew John Buggee

%% Define the cloud parameters that will be changing during each reflectance calculation


r_top = 5:13;       % microns
r_bot = 4:10;        % microns

tau_c = 5:5:35;

% r_top = 6:2:12;       % microns
% r_bot = 4:2:10;        % microns
% tau_c = 5:5:35;


% r_top = 9.03;       % microns
% r_bot = 9.03;        % microns
% tau_c = 6.34;




%% Want to use real EMIT geometry inputs?

% Load modis data and create input structure

% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the EMIT data folder path

    emitPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';



    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the EMIT data folder path

    emitPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];


end


% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitFolder = '17_Jan_2024_coast/';



[emit,L1B_fileName] = retrieveEMIT_data([emitPath, emitFolder]);


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% row = 1112;
% col = 974;

% 17_Jan_2024_coast - large optical depth
row = 912;
col = 929;

emit_idx = sub2ind(size(emit.radiance.measurements), row, col);    % for 9 nov 2008 - pixel overlapping with VOCALS


%% Grab the EMIT radiances for the pixel used
Rad_emit = reshape(emit.radiance.measurements(row, col, :), [],1);      % microW/cm^2/nm/sr



%% Define the parameters of the INP file


% Define the number of streams to use in your radiative transfer model
num_streams = 16;
% ------------------------------------------------------------------------


% ------- Define the source file ------
source_file_resolution = 0.1;           % nm

%source_file = '../data/solar_flux/lasp_TSIS1_hybrid_solar_reference_p01nm_resolution.dat';
if source_file_resolution==0.1
    
    source_file = '../data/solar_flux/kurudz_0.1nm.dat';

elseif source_file_resolution==1

    source_file = '../data/solar_flux/kurudz_1.0nm.dat';

end



% Define the spectral response function
% ------------------------------------------------------------------------

% ****************************************************************
% *-*-*-*- Only keep wavelengths that avoid water vapor -*-*-*-*-*

% The following indexes are for wavelengths that avoid water vapor
% absopriton according to figure 5 from King and Vaughan, which shows the
% information content for r_top, r_bot, tau_c, and water vapor across
% wavelengths from 500 to 2500 nm
wavelength_idx = [17, 24, 31, 40, 52, 65, 86, 92, 93, 115, 117, 118, 119, 121,...
    158, 159, 164, 165, 166, 167, 168, 174, 175, 221, 222, 226, 232, 234, 238,...
    248, 252, 259]';
Rad_emit = Rad_emit(wavelength_idx);
% -------------------------------------------------


% define the wavelength range. If monochromatic, enter the same number
% twice
% ------------------------------------------------------------------------

spec_response = cell(length(wavelength_idx), 1);
wavelength = zeros(length(wavelength_idx), 2);

for ww = 1:length(wavelength_idx)

    % first we will create and store the spectral response function from
    % the full-wdith-half-max provided for each spectral channel
    % the spectral response function is gaussian
    % the emit wavelength vector is the center wavelength
    
    % set the center wavelength as the mean of the distribution
    mu = emit.radiance.wavelength(wavelength_idx(ww));      % nm
    % compute the standard deviation from the FWHM

    sigma = emit.radiance.fwhm(wavelength_idx(ww))/(2*sqrt(2*log(2)));      % std

    % create a wavelength vector
    if source_file_resolution==0.1
        
        wl = round(mu-(1.5*emit.radiance.fwhm(wavelength_idx(ww))), 1):...
            source_file_resolution:round(mu+(1.5*emit.radiance.fwhm(wavelength_idx(ww))), 1);

    elseif source_file_resolution==1
        
        wl = round(mu-(1.5*emit.radiance.fwhm(wavelength_idx(ww)))):...
            round(mu+(1.5*emit.radiance.fwhm(wavelength_idx(ww))));

    end

    % compute the gaussian spectral response function
    spec_response{ww} = pdf('Normal', wl, mu, sigma)';

    % The wavelength vector for libRadTran is simply the lower and upper
    % bounds
    wavelength(ww,:) = [wl(1), wl(end)];

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
day_of_year = 17;

% ------------------------------------------------------------------------
cloud_depth = 500;                % meters

% define the geometric location of the cloud top and cloud bottom
z_topBottom = [2.5, 2];          % km above surface



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
sza = double(emit.obs.solar.zenith(emit_idx));           % degree

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
phi0 = mod(double(emit.obs.solar.azimuth(emit_idx) + 180), 360);         % degree

% define the viewing zenith angle
vza = double(emit.obs.sensor.zenith(emit_idx)); % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

vaz = emit.obs.sensor.azimuth(emit_idx);     % degree


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

% --------------------------------------------------------
% --------- What is column water vapor amount? -----------

% Using measurements from the AMSR2 instrument, a passive microwave
% radiometer for 17 Jan 2024
yesModify_waterVapor = true;

waterVapor_column = 16;              % mm - milimeters of water condensed in a column
% ------------------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
compute_reflectivity_uvSpec = false;
% --------------------------------------------------------------








%% Write each INP file and Calculate Reflectance for MODIS

inputName = cell(length(r_top), length(r_bot), length(tau_c), size(wavelength,1));
outputName = cell(length(r_top), length(r_bot), length(tau_c), size(wavelength,1));
wc_filename = cell(length(r_top), length(r_bot), length(tau_c), size(wavelength,1));


Rad_model = zeros(length(r_top), length(r_bot), length(tau_c), size(wavelength,1));
Refl_model = zeros(length(r_top), length(r_bot), length(tau_c), size(wavelength,1));

Refl_emit = zeros(size(Rad_emit));

tic
for rt = 1:length(r_top)


    for rb = 1:length(r_bot)


        for tc = 1:length(tau_c)


            disp(['Iteration: [rt, rb, tc] = [', [num2str(rt),', ', num2str(rb), ', ', num2str(tc)], ']...', newline])


            parfor ww = 1:size(wavelength,1)



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
                    'nm_radiance_', atm_file(1:end-4),'.INP'];



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



                % Define the column water vapor amount
                % --------------------------------------------------------------------
                if yesModify_waterVapor==true

                    % Turn on default aersol layer, which occupies lower 2km of model
                    % --------------------------------------------------------------
                    formatSpec = '%s %f %s %5s %s \n\n';
                    fprintf(fileID, formatSpec,'mol_modify H2O ', waterVapor_column, ' MM', ' ', '# Column water vapor amount');


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

                % compute the reflectance
                Refl_model(rt, rb, tc, ww) = reflectanceFunction(inputSettings(2,:), ds, spec_response{ww});

                % integrate the radiance over the wavelength channel and
                % convert the output to the same units as the EMIT data
                % LibRadTran reports radiance in units of mW/m^2/nm/sr
                % EMIT reports radiance in units of microW/cm^2/nm/sr
                % Divide the libRadTran values to get the output into the
                % same units as EMIT
                % Make sure to integrate with the spectral response
                % function!
                Rad_model(rt, rb, tc, ww) = trapz(ds.wavelength, spec_response{ww}.*ds.radiance.value)/10;        % microW/cm^2/sr

                % Now divide by the wavelength range of the channel to get
                % units of micro-watts/cm^2/sr/nm

                Rad_model(rt, rb, tc, ww) = Rad_model(rt, rb, tc, ww)/(ds.wavelength(end)-ds.wavelength(1));       % microW/cm^2/nm/sr

                % Convert the EMIT radiance into reflectance
                if rt==1 && rb==1 && tc==1
                    
                    % The source function has units of mW/m^2/nm
                    % Divide this by 10 to get units of microW/cm^2/nm
                    Refl_emit(ww) = pi*Rad_emit(ww)/(cosd(inputSettings{2,4}) * ...
                        trapz(inputSettings{2,7}(:,1), spec_response{ww} .* inputSettings{2,7}(:,2))/10);

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

    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/EMIT/Reflectance_Uniqueness/'];



elseif strcmp(whatComputer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'EMIT/Reflectance_Uniqueness/'];


elseif strcmp(whatComputer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    warning([newline, 'No folder to store things in!', newline])



end


filename = [folderpath,'reflectance_calcs_MODIS-data-from-',emitFolder(1:end-1),...
    '_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];

while isfile(filename)
    rev = rev+1;
    filename = [folderpath,'reflectance_calcs_MODIS-data-from-',emitFolder(1:end-1),...
        '_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];
end

save(filename,"r_top", "r_bot", "tau_c", "wavelength", "Rad_model", "emitFolder", 'emit_idx');

toc


%% Subplots of Radiance across different optical depths for a single wavelength

wave_len_idx = 1;

% find min and max values of reflectance for this wavelength
[minR, ~] = min(Rad_model(:,:,:,wave_len_idx), [], 'all');
[maxR, ~] = max(Rad_model(:,:,:,wave_len_idx), [], 'all');

% set colorbar limits for each subplot
color_lim = [minR, maxR];



figure;
for tc = 1:length(tau_c)

    subplot(2,4,tc); imagesc(r_bot, r_top, Rad_model(:,:,tc, wave_len_idx));

    if tc==1
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
    end

    title(['$\tau_c$ = ', num2str(tau_c(tc))], 'Interpreter', 'latex')

    if tc==4
        cb = colorbar;
        set(get(cb, 'label'), 'string', 'Radiance $(\mu W/cm^{2}/nm/sr)$','Interpreter','latex', 'Fontsize',28)
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





%% Make lineplots of Radiance versus optical depth at a constant wavelength


wave_len_idx = 1;



figure;
for rt = 1:length(r_top)
    for rb = 1:length(r_bot)


        plot(tau_c, reshape(Rad_model(rt,rb,:, wave_len_idx), 1, []));

        hold on

    end
end

xlabel('$\tau_{c}$', 'Interpreter', 'latex')
ylabel('Radiance $(\mu W/cm^{2}/nm/sr)$', 'Interpreter', 'latex')

title(['Reflectance $\lambda$ = ', num2str(round(mean(wavelength(wave_len_idx, :)))), ' $nm$'], 'Interpreter', 'latex')



% set the plot size
set(gcf, 'Position', [0 0 1200 600])
grid on; grid minor




%% Subplots of radiance across different wavelengths for a single optical depth

tau_idx = 4;

% find min and max values of reflectance for this wavelength
[minR, ~] = min(Rad_model(:,:,tau_idx,:), [], 'all');
[maxR, ~] = max(Rad_model(:,:,tau_idx,:), [], 'all');



figure;
for ww = 1:size(wavelength, 1)

    subplot(2,4,ww); imagesc(r_bot, r_top, Rad_model(:,:,tau_idx, ww));

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







%% Subplots of radiance across different optical depths for all wavelengths




for ww = 1:size(wavelength, 1)

    % find min and max values of reflectance for this wavelength
    [minR, ~] = min(Rad_model(:,:,:,ww), [], 'all');
    [maxR, ~] = max(Rad_model(:,:,:,ww), [], 'all');

    figure;
    for tc = 1:length(tau_c)

        subplot(2,4,tc); imagesc(r_bot, r_top, Rad_model(:,:,tc,ww));

        if tc==1
            xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter', 'latex')
            ylabel('$r_{top}$ ($\mu m$)', 'Interpreter', 'latex')
        end

        title(['$\tau_c$ = ', num2str(tau_c(tc))], 'Interpreter', 'latex')

        if tc==3
            cb = colorbar;
            set(get(cb, 'label'), 'string', 'Radiance $(\mu W/cm^{2}/nm/sr)$',...
                'Interpreter','latex', 'Fontsize',28)


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



%% Plot a 3D surf figure with r_bot r_top and radiance as the variables. Plot
% for a single optical depth at a single wavelength

wave_len_idx = 1;
tau_c_idx = 5;

% ADD A PLANE OF CONSTANT REFLECTANCE TO SHOW DIFFERENT SOLUTIONS
R_constant = Rad_emit(wave_len_idx);


r_top_mat = repmat(r_top', 1, length(r_bot));
r_bot_mat = repmat(r_bot, length(r_top), 1);
R_mat = reshape(Rad_model(:,:, tau_c_idx, wave_len_idx), size(r_top_mat));

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

    R_mat = reshape(Rad_model(:,:, tau_idx, ww), size(r_top_mat));


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



%% Can I determine uniqueness by rounding the reflectance computed to the measurement uncertainty of EMIT?
% Then, see how redundant certain states are for all EMIT wavelenghts
% OR, should I instead ask, how many sets of measurements I computed are
% within the uncertainty of the EMIT measurements? That is maybe a better
% estimate of uniqueness, because EMIT claims a specific measurement but
% then includes confidence intervals. They can say for certain that their
% measurement lies somewhere within that uncertainty interval.

% We have to compute the EMIT reflectance uncertainty using 3 spectrally
% indexed functions provided in the metadata

% ******** but I can't find these metadata ***********
% So for now, let's say it's 2% across all channels for now
emit_uncert = 0.02;

% Let's assume the reflectance uncertainty is 1%
% Let's truncate the reflectance data to hundreds decimal point

% For our spectral measurements, how many different states of (r_top,
% r_bot, tau_c) lead to the same set of spectral measurements within the
% uncertainty

% Let's reshape R_model_round to be W rows, where W is the number of wavelengths
% and N number of columns corresponding the number of unique states

% Grab the EMIT radiances for the pixel used
Refl_emit_uncert = Refl_emit * emit_uncert; % microW/cm^2/nm/sr - converted from percentage to radiance


redundant_states = [];


for rt = 1:length(r_top)


    for rb = 1:length(r_bot)


        for tc = 1:length(tau_c)

            % Check to see if the radaiance computed by the model is
            % within the listed uncertainty for EMIT
            redundant_states = [redundant_states, abs(Refl_emit - reshape(Refl_model(rt,rb,tc,:), [], 1)) <= Refl_emit_uncert];


        end
    end
end


% Find the number states that lead to modeled measurements within the EMIT
% measurement uncertainty
num_states = sum(all(redundant_states, 1));

% print the percentage of redundant states
disp([newline, num2str(100* (num_states/(numel(r_top)*numel(r_bot)*numel(tau_c)))),...
    '% of modeled states are redundant', newline])





%% 3D Interpolate the radiance calculations on a finer grid and compare with EMIT measurement


% Meshgrid is defined on x,y,z space, not row, column, depth space
% In 3D space, z = row, x = column, y = depth
[R_bot, R_top, Tau_c] = meshgrid(r_bot, r_top, tau_c);

% Create the new fine grid to interpolate on
% define the discrete step length of each variable
d_r_top = 0.25;      % microns
d_r_bot = 0.25;      % microns
d_tau_c = 0.1;

r_top_fine = r_top(1):d_r_top:r_top(end);
r_bot_fine = r_bot(1):d_r_bot:r_bot(end);
tau_c_fine = tau_c(1):d_tau_c:tau_c(end);

[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);

Refl_model_fine = zeros(length(r_top_fine), length(r_bot_fine), length(tau_c_fine), size(Refl_model,4));



for wl = 1:size(Refl_model,4)

    Refl_model_fine(:,:,:,wl) = interp3(R_bot, R_top, Tau_c, Refl_model(:, :, :, wl),...
        R_bot_fine, R_top_fine, Tau_c_fine);


end



% Using the new fine grid, calculate how many sets of measurements are
% within the EMIT measurement and it's uncertainty


redundant_states = [];
rms_residual = zeros(length(r_top_fine), length(r_bot_fine), length(tau_c_fine));



for rt = 1:size(Refl_model_fine,1)


    for rb = 1:size(Refl_model_fine,2)


        parfor tc = 1:size(Refl_model_fine,3)

            % Check to see if the radiance computed by the model is
            % within the listed uncertainty for EMIT
            %redundant_states(rt,rb,tc) = all(abs(R_emit - reshape(R_model_fine(rt,rb,tc,:), 1, [])) <= R_emit_uncert);
            redundant_states = [redundant_states, abs(Refl_emit - reshape(Refl_model_fine(rt,rb,tc,:), [], 1)) <= Refl_emit_uncert];
            rms_residual(rt, rb, tc) = sqrt(mean( (Refl_emit - reshape(Refl_model_fine(rt,rb,tc,:), [], 1)).^2) );


        end
    end
end

% Find the number states that lead to modeled measurements within the EMIT
% measurement uncertainty
num_states = sum(all(redundant_states, 1));

% print the percentage of redundant states
disp([newline, num2str(100* (num_states/(numel(r_top_fine)*numel(r_bot_fine)*numel(tau_c_fine)))),...
    '% of modeled states that produce measurements within the EMIT uncertainty', newline])


if num_states>1

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

end


%% Find the states with the lowest rms residul

% find n smallest rms states
n_states = 10;

% store the state values at each minimum
r_top_min = zeros(n_states, 1);
r_bot_min = zeros(n_states, 1);
tau_c_min = zeros(n_states, 1);

% store the rms value and the index
min_val = zeros(n_states, 1);
idx_min = zeros(n_states, 1);


for nn = 1:n_states
    
    % find the smallest rms residual value, omitting nans
    [min_val(nn), idx_min(nn)] = min(rms_residual, [], 'all', 'omitnan');

    r_top_min(nn) = R_top_fine(idx_min(nn));
    r_bot_min(nn) = R_bot_fine(idx_min(nn));
    tau_c_min(nn) = Tau_c_fine(idx_min(nn));

    % set the minimum value to nan and omit
    rms_residual(idx_min(nn)) = nan;


end



%% Create Contour plot of rms residual between true EMIT measurements and the libRadTran modeled measurements


% define the optical depth slice you'd like to plot
idx_tauC = tau_c_fine == tau_c_min(1);

% Create figure
figure;
colormap(hot);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create contour
[c1,h1] = contour(r_bot_fine, r_top_fine, rms_residual(:,:, idx_tauC),'LineWidth',3);
clabel(c1,h1,'FontSize',15,'FontWeight','bold');

% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex');

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex');

% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC))],'Interpreter','latex');

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% Create colorbar
colorbar(axes1);




% Also try a filled contour plot


% Create figure
figure;
colormap(hot);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create contour
[c1,h1] = contourf(r_bot_fine, r_top_fine, rms_residual(:,:, idx_tauC),'LineWidth',3);
clabel(c1,h1,'FontSize',15,'FontWeight','bold');

% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex');

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex');

% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC))],'Interpreter','latex');

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% Create colorbar
colorbar(axes1);


%% Create an imagesc plot of the rms residual


% define the optical depth slice you'd like to plot
idx_tauC = tau_c_fine == tau_c_min(1);

% Create figure
figure;
colormap(parula);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create contour
imagesc(r_bot_fine, r_top_fine, rms_residual(:,:, idx_tauC));

% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex');

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex');

% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC))],'Interpreter','latex');

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% Create colorbar
colorbar(axes1);

