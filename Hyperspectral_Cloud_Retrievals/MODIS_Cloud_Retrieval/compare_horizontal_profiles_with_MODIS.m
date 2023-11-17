%% Find MODIS pixels near same path as Horizontal Transects and compare microphysics

% By Andrew John Buggee

%% LOAD MODIS DATA

clear variables


% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer, 'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % Define the MODIS folder name

    % ----- November 9th at decimal time 0.611 (14:40) -----
    modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/';


    % ----- November 11th at decimal time 0.604 (14:30) -----
    %modisFolder = ['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/'];


    % ----- November 11th at decimal time 0.784 (18:50) -----
    %modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/';






elseif strcmp(whatComputer,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % Define the MODIS folder name

    % ----- October 18th at decimal time 0.6354 (15:15) -----
    %modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_10_18/';

    % ----- November 2nd at decimal time 0.607 (14:35) -----
    %modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_02/';


    % ----- November 9th at decimal time 0.611 (14:40) -----
    modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/';


    % ----- November 11th at decimal time 0.604 (14:30) -----
    %modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/'];


    % ----- November 11th at decimal time 0.784 (18:50) -----
    %modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/';




end




[modis,L1B_fileName] = retrieveMODIS_data(modisFolder);

% ----- Create a structure defining inputs of the problem -----

modisInputs = create_modis_inputs(modisFolder, L1B_fileName);

%% LOAD VOCALS REX DATA


% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    folder_path = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/';

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
        'VOCALS_REx/vocals_rex_data/SPS_1/'];


end

% Oct-15-2008 Data
%filename = 'RF01.20081015.164800_201200.PNI.nc';

% Oct-18-2008 Data
%filename = 'RF02.20081018.130300_213000.PNI.nc';

% Oct-21-2008 Data
%filename = 'RF03.20081021.060000_142400.PNI.nc';

% ----- October 23 data -----
%filename = 'RF04.20081023.055800_142200.PNI.nc';

% Oct-25-2008 Data
%filename = 'RF05.20081025.062900_152500.PNI.nc';

% ----- October 28 data -----
%filename = 'RF06.20081028.061800_151200.PNI.nc';

% ----- October 31 data -----
%filename = 'RF07.20081031.060000_150000.PNI.nc';

% ----- November 2 data -----
%filename = 'RF08.20081102.055700_152100.PNI.nc';

% ----- November 4 data -----
%filename = 'RF09.20081104.060000_145600.PNI.nc';

% ----- November 6 data -----
%filename = 'RF10.20081106.074800_142000.PNI.nc';

% ----- November 9 data -----
filename = 'RF11.20081109.125700_213600.PNI.nc';

% ------ November 11 data -----
%filename = 'RF12.20081111.125000_214500.PNI.nc';

% ----- November 13 data -----
%filename = 'RF13.20081113.125700_215700.PNI.nc';

% ----- November 15 data -----
%filename = 'RF14.20081115.125800_220400.PNI.nc';


vocalsRex = readVocalsRex([folder_path, filename]);



%% Find Horizontal Profile Closest to MODIS Overpass


% ---- set thresholds for the LWC and Nc ---
lwc_threshold = 0.03;       % g/m^3
Nc_threshold = 1;           % cm^{-3}
max_vertical_displacement = 20;     % meters


% Lets only keep the vocalsRex data we need
% Time is measured in seconds since the startTime

% ---- DO YOU WANT TO USE ADVECTION? -----
modisInputs.flags.useAdvection = true;

tic

vocalsRex = cropVocalsRex_horzProfs2MODIS(vocalsRex, lwc_threshold, Nc_threshold,...
    max_vertical_displacement, modis, modisInputs);
toc

% Let's save this data, since it can take a long time to calculate
if modisInputs.flags.useAdvection == true

    save([modisInputs.savedCalculations_folderName,'horz_prof_closest2MODIS_withAdvection_',...
        char(datetime("today")),'.mat'],...
        "vocalsRex", "lwc_threshold", "Nc_threshold", "max_vertical_displacement");

else

    save([modisInputs.savedCalculations_folderName,'horz_prof_closest2MODIS_withoutAdvection_',...
        char(datetime("today")),'.mat'],...
        "vocalsRex", "lwc_threshold", "Nc_threshold", "max_vertical_displacement");

end



%%  Plot Vocals-Rex Horizontal Profile along with MODIS estimates

% plot the effective radius - we want to understand horizontal variability

normalize_distance = false;

% define the colors of each curve
C = mySavedColors([8,9], 'fixed');

legend_str = cell(1, 2);

std_val = zeros(1, 1);

mean_val = zeros(1, 1);

figure;

% if normalize distance is true, all distance vectors will be
% normalized between 0 and 1

if normalize_distance==true

    norm_dist = (vocalsRex.horz_dist - min(vocalsRex.horz_dist))./...
        (max(vocalsRex.horz_dist) - min(vocalsRex.horz_dist));

    % plot the effective radius
    % if the 2DC data is compliant, plot the effective radius computed
    % using both instruments
    if vocalsRex.flag_2DC_data_is_conforming==true

        plot(norm_dist, vocalsRex.re, 'Color',C(nn,:));

        hold on;

        % plot the average value as a dashed line
        mean_val(nn) = mean(vocalsRex.re);
        constant_y_val = linspace(mean_val(nn), mean_val(nn), length(vocalsRex.horz_dist));
        plot(norm_dist, constant_y_val,'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

        % compute the standard deviation
        std_val(nn) = std(vocalsRex.re);


    else

        % if the 2DC data is non-conforming, use only the CDP data and
        % make a note of it
        plot(norm_dist, vocalsRex.re_CDP, 'Color',C(nn,:));

        % plot the average value as a dashed line
        mean_val(nn) = mean(vocalsRex.re_CDP);
        constant_y_val = linspace(mean_val, mean_val, length(vocalsRex.horz_dist));
        plot(norm_dist, constant_y_val,'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

        % compute the standard deviation
        std_val(nn) = std(vocalsRex.re_CDP);

    end
    hold on



else

    % --- DATA IS IN METERS - PLOT IN KILOMETERS ---

    % plot the effective radius
    % if the 2DC data is compliant, plot the effective radius computed
    % using both instruments

    if vocalsRex.flag_2DC_data_is_conforming==true

        plot(vocalsRex.horz_dist./1e3, vocalsRex.re,...
            'Color',C(1,:));

        hold on;

        % plot the average value as a dashed line
        mean_val = mean(vocalsRex.re);
        constant_y_val = linspace(mean_val, mean_val, length(vocalsRex.horz_dist));
        plot(vocalsRex.horz_dist./1e3, constant_y_val,...
            'LineStyle','--','LineWidth',1.5,'Color', C(1,:));

        % compute the standard deviation
        std_val = std(vocalsRex.re);

    else

        % if the 2DC data is non-conforming, use only the CDP data and
        % make a note of it
        plot(vocalsRex.horz_dist./1e3, vocalsRex.re_CDP,...
            'Color',C(1,:));

        hold on

        % plot the average value as a dashed line
        mean_val = mean(vocalsRex.re_CDP);
        constant_y_val = linspace(mean_val, mean_val, length(vocalsRex.horz_dist));
        plot(vocalsRex.horz_dist./1e3, constant_y_val,...
            'LineStyle','--','LineWidth',1.5,'Color', C(1,:));

        % compute the standard deviation
        std_val = std(vocalsRex.re_CDP);

    end
    hold on


end

legend_str{1} = ['$\sigma_{in-situ}$ = ', num2str(round(std_val, 2)), ' $\mu m$'];
% skip one because of the mean value
legend_str{2} = ['$\left<r_e^{in-situ} \right>$ = ', num2str(round(mean_val, 2)), ' $\mu m$'];

% ----------------------------------------
% ----- Now Lets Plot the MODIS Data -----
% ----------------------------------------

% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with a length unit of meters.
wgs84 = wgs84Ellipsoid("m");

modis_re = zeros(1, length(vocalsRex.modisIndex_minDist));
modis_tau_c = zeros(1, length(vocalsRex.modisIndex_minDist));
modis_lat = zeros(1, length(vocalsRex.modisIndex_minDist));
modis_long = zeros(1, length(vocalsRex.modisIndex_minDist));

for nn = 1:length(vocalsRex.modisIndex_minDist)

    modis_re(nn) = modis.cloud.effRadius17(vocalsRex.modisIndex_minDist(nn));
    modis_tau_c(nn) = modis.cloud.optThickness17(vocalsRex.modisIndex_minDist(nn));
    modis_lat(nn) = modis.geo.lat(vocalsRex.modisIndex_minDist(nn));
    modis_long(nn) = modis.geo.long(vocalsRex.modisIndex_minDist(nn));

end

modis_dist_traveled_fromVR = distance(vocalsRex.latitude(1), vocalsRex.longitude(1), modis_lat, modis_long, wgs84);

hold on
plot(modis_dist_traveled_fromVR./1e3, modis_re, '.', 'Color', C(2,:), 'MarkerSize', 17)

% plot the average value as a dashed line
mean_val = mean(modis_re);
constant_y_val = linspace(mean_val, mean_val, length(modis_dist_traveled_fromVR));
plot(modis_dist_traveled_fromVR./1e3, constant_y_val,...
    'LineStyle','--','LineWidth',1.5,'Color', C(2,:));

legend_str{3} = ['$\sigma_{modis}$ = ', num2str(round(std(modis_re), 2)), ' $\mu m$'];
legend_str{4} = ['$\left<r_e^{modis} \right>$ = ', num2str(round(mean(modis_re), 2)), ' $\mu m$'];



% ---- Make Plot Pretty! ----


grid on; grid minor;
% if the 2DC data is compliant, plot the effective radius computed
% using both instruments
if vocalsRex.flag_2DC_data_is_conforming==true
    ylabel('$r_e$ ($\mu m$)', 'Interpreter','latex')
else
    % if the 2DC data is non-conforming, use only the CDP data and
    % make a note of it
    ylabel('$r_e$ ($\mu m$) - (CDP only)', 'Interpreter','latex')
end

% include a title in the middle plot
if isfield(vocalsRex, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(vocalsRex.LWC_threshold),' $g/m^{3}$',...
        '     $N_c \geq$ ', num2str(vocalsRex.Nc_threshold), ' $cm^{-3}$',...
        '     Max vert displacement: ', num2str(vocalsRex.max_vert_displacement), ' $m$'], 'interpreter', 'latex')

elseif isfield(vocalsRex.inputs, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(vocalsRex.inputs.LWC_threshold),' $g/m^{3}$',...
        '     $N_c \geq$ ', num2str(vocalsRex.inputs.Nc_threshold), ' $cm^{-3}$',...
        '     Max vert displacement: ', num2str(vocalsRex.inputs.max_vert_displacement), ' $m$'], 'interpreter', 'latex')

end

% Include an x axis label on the middle plot
if normalize_distance==true

    xlabel('Normalized Horizontal Distance Travelled', 'Interpreter','latex');
else

    xlabel('Horizontal Distance Travelled ($km$)', 'Interpreter','latex');
end




% in the third subplot, define the indices_2_plot being plotted
legend(legend_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20)

% set plot size
set(gcf, 'Position', [0 0 1200 625])







%%  Plot Vocals-Rex Horizontal Profile along with MODIS estimates
% Set the color of the MODIS marker to be the value of the optical depth

% plot the effective radius - we want to understand horizontal variability

normalize_distance = false;

% define the colors of each curve
C = mySavedColors([5,9], 'fixed');

legend_str = cell(1, 2);

std_val = zeros(1, 1);

mean_val = zeros(1, 1);

figure;

% if normalize distance is true, all distance vectors will be
% normalized between 0 and 1

if normalize_distance==true

    norm_dist = (vocalsRex.horz_dist - min(vocalsRex.horz_dist))./...
        (max(vocalsRex.horz_dist) - min(vocalsRex.horz_dist));

    % plot the effective radius
    % if the 2DC data is compliant, plot the effective radius computed
    % using both instruments
    if vocalsRex.flag_2DC_data_is_conforming==true

        plot(norm_dist, vocalsRex.re, 'Color',C(nn,:));

        hold on;

        % plot the average value as a dashed line
        mean_val(nn) = mean(vocalsRex.re);
        constant_y_val = linspace(mean_val(nn), mean_val(nn), length(vocalsRex.horz_dist));
        plot(norm_dist, constant_y_val,'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

        % compute the standard deviation
        std_val(nn) = std(vocalsRex.re);


    else

        % if the 2DC data is non-conforming, use only the CDP data and
        % make a note of it
        plot(norm_dist, vocalsRex.re_CDP, 'Color',C(nn,:));

        % plot the average value as a dashed line
        mean_val(nn) = mean(vocalsRex.re_CDP);
        constant_y_val = linspace(mean_val, mean_val, length(vocalsRex.horz_dist));
        plot(norm_dist, constant_y_val,'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

        % compute the standard deviation
        std_val(nn) = std(vocalsRex.re_CDP);

    end
    hold on



else

    % --- DATA IS IN METERS - PLOT IN KILOMETERS ---

    % plot the effective radius
    % if the 2DC data is compliant, plot the effective radius computed
    % using both instruments

    if vocalsRex.flag_2DC_data_is_conforming==true

        plot(vocalsRex.horz_dist./1e3, vocalsRex.re,...
            'Color',C(1,:));

        hold on;

        % plot the average value as a dashed line
        mean_val = mean(vocalsRex.re);
        constant_y_val = linspace(mean_val, mean_val, length(vocalsRex.horz_dist));
        plot(vocalsRex.horz_dist./1e3, constant_y_val,...
            'LineStyle','--','LineWidth',1.5,'Color', C(1,:));

        % compute the standard deviation
        std_val = std(vocalsRex.re);

    else

        % if the 2DC data is non-conforming, use only the CDP data and
        % make a note of it
        plot(vocalsRex.horz_dist./1e3, vocalsRex.re_CDP,...
            'Color',C(1,:));

        hold on

        % plot the average value as a dashed line
        mean_val = mean(vocalsRex.re_CDP);
        constant_y_val = linspace(mean_val, mean_val, length(vocalsRex.horz_dist));
        plot(vocalsRex.horz_dist./1e3, constant_y_val,...
            'LineStyle','--','LineWidth',1.5,'Color', C(1,:));

        % compute the standard deviation
        std_val = std(vocalsRex.re_CDP);

    end
    hold on


end

legend_str{1} = ['$\sigma_{in-situ}$ = ', num2str(round(std_val, 2)), ' $\mu m$'];
% skip one because of the mean value
legend_str{2} = ['$\left<r_e^{in-situ} \right>$ = ', num2str(round(mean_val, 2)), ' $\mu m$'];

% ----------------------------------------
% ----- Now Lets Plot the MODIS Data -----
% ----------------------------------------

% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with a length unit of meters.
wgs84 = wgs84Ellipsoid("m");

modis_re = zeros(1, length(vocalsRex.modisIndex_minDist));
modis_tau_c = zeros(1, length(vocalsRex.modisIndex_minDist));
modis_lat = zeros(1, length(vocalsRex.modisIndex_minDist));
modis_long = zeros(1, length(vocalsRex.modisIndex_minDist));

for nn = 1:length(vocalsRex.modisIndex_minDist)

    modis_re(nn) = modis.cloud.effRadius17(vocalsRex.modisIndex_minDist(nn));
    modis_tau_c(nn) = modis.cloud.optThickness17(vocalsRex.modisIndex_minDist(nn));
    modis_lat(nn) = modis.geo.lat(vocalsRex.modisIndex_minDist(nn));
    modis_long(nn) = modis.geo.long(vocalsRex.modisIndex_minDist(nn));

end

% ----- Compute the projected distance along the VR profile -----
% Using the MODIS pixel position, we want to project this vector onto the
% VR in-situ vector asd compute the distance travelled along this
% trajectory.
% The origin of this coordinate system is the first VR data point



% REMEMBER - if a simple advection model was used, we need the new lat long
% position.
if modisInputs.flags.useAdvection==false

    % compute the distance (magnitude) from the first VR point to each MODIS
    % pixel location. (in km)
    modis_dist_fromVR = distance(vocalsRex.latitude(1), vocalsRex.longitude(1), modis_lat, modis_long, wgs84)./1e3;        % meters
    
    % Compute the angle between the first VR data point and each MODIS pixel
    % This angle is measured with respect to due north
    az_MODIS = azimuth(vocalsRex.latitude(1), vocalsRex.longitude(1), modis_lat, modis_long, wgs84);
    
    % Compute the angle with respect to noraml of the vocals-rex profile
    az_VR = azimuth(vocalsRex.latitude(1), vocalsRex.longitude(1), vocalsRex.latitude(end), vocalsRex.longitude(end), wgs84);


else

    % compute the distance (magnitude) from the first VR point to each MODIS
    % pixel location. (in km)
    modis_dist_fromVR = distance(vocalsRex.lat_withAdvection(1), vocalsRex.long_withAdvection(1), modis_lat, modis_long, wgs84)./1e3;        % meters
    

    % Compute the angle between the first VR data point and each MODIS pixel
    % This angle is measured with respect to due north
    az_MODIS = azimuth(vocalsRex.lat_withAdvection(1), vocalsRex.long_withAdvection(1), modis_lat, modis_long, wgs84);
    
    % Compute the angle with respect to noraml of the vocals-rex profile
    az_VR = azimuth(vocalsRex.lat_withAdvection(1), vocalsRex.long_withAdvection(1), vocalsRex.latitude(end), vocalsRex.longitude(end), wgs84);

end



% compute the angle between the vector connected the first VR data point
% and MODIS and the VR profile
az_between_MODIS_VR = abs(az_MODIS - az_VR);

% the distance of the MODIS pixel location along the VR flight path is the
% magnitude of the MODIS location multiplied by the cosine of the angle
% between the two 
dist_MODIS_along_VR = modis_dist_fromVR .* cosd(az_between_MODIS_VR);

% For now, we will assume that the distances between the first data point
% taken by vocals rex and the MODIS pixels is so small that we can ignore
% spherical trigonometry


% ------------------------------------------------------------------

% Define a Colormap for the optical depth values
[modis_tau_c_sorted, idx_sort] = sort(modis_tau_c);

% Check to see if all the data is unique
[modis_tau_c_unique, ~, idx_unique] = unique(modis_tau_c(idx_sort));

% We span the colormap only by the number of unique data values
C_tau = parula(length(modis_tau_c_unique));

% And now let's repeat the same colors for the redundant data points
C_tau = C_tau(idx_unique,:);


% sort the distance travelled and the modis effective radius by the optical
% depth
modis_dist_sort = dist_MODIS_along_VR(idx_sort);
modis_re_sort = modis_re(idx_sort);



hold on
for nn = 1:length(modis_re_sort)
    plot(modis_dist_sort(nn), modis_re_sort(nn), '.', 'Color', C_tau(nn,:), 'MarkerSize', 28)
end

% create colorbar and colorbar label
cb = colorbar;
clim([modis_tau_c_sorted(1), modis_tau_c_sorted(end)])

set(get(cb, 'label'), 'string', '$\tau_c$','Interpreter','latex', 'Fontsize',28)


% plot the average value as a dashed line
mean_val = mean(modis_re);
constant_y_val = linspace(mean_val, mean_val, length(dist_MODIS_along_VR));
plot(dist_MODIS_along_VR, constant_y_val,...
    'LineStyle','--','LineWidth',1.5,'Color', C_tau(1,:));

legend_str{3} = ['$\sigma_{modis}$ = ', num2str(round(std(modis_re), 2)), ' $\mu m$'];
legend_str{4} = ['$\left<r_e^{modis} \right>$ = ', num2str(round(mean(modis_re), 2)), ' $\mu m$'];



% ---- Make Plot Pretty! ----


grid on; grid minor;
% if the 2DC data is compliant, plot the effective radius computed
% using both instruments
if vocalsRex.flag_2DC_data_is_conforming==true
    ylabel('$r_e$ ($\mu m$)', 'Interpreter','latex')
else
    % if the 2DC data is non-conforming, use only the CDP data and
    % make a note of it
    ylabel('$r_e$ ($\mu m$) - (CDP only)', 'Interpreter','latex')
end

% include a title in the middle plot
if isfield(vocalsRex, 'LWC_threshold')==true

    if modisInputs.flags.useAdvection == false

        title(['Without Advection - $LWC \geq$ ', num2str(vocalsRex.LWC_threshold),' $g/m^{3}$',...
            '     $N_c \geq$ ', num2str(vocalsRex.Nc_threshold), ' $cm^{-3}$',...
            '     Max vert displacement: ', num2str(vocalsRex.max_vert_displacement), ' $m$'],...
            'FontSize', 30, 'interpreter', 'latex')

    else

        title(['With Advection - $LWC \geq$ ', num2str(vocalsRex.LWC_threshold),' $g/m^{3}$',...
            '     $N_c \geq$ ', num2str(vocalsRex.Nc_threshold), ' $cm^{-3}$',...
            '     Max vert displacement: ', num2str(vocalsRex.max_vert_displacement), ' $m$'],...
            'FontSize', 30, 'interpreter', 'latex')

    end

elseif isfield(vocalsRex.inputs, 'LWC_threshold')==true

    if modisInputs.flags.useAdvection == false

        title(['Without Advection - $LWC \geq$ ', num2str(vocalsRex.inputs.LWC_threshold),' $g/m^{3}$',...
            '     $N_c \geq$ ', num2str(vocalsRex.inputs.Nc_threshold), ' $cm^{-3}$',...
            '     Max vert displacement: ', num2str(vocalsRex.inputs.max_vert_displacement), ' $m$'],...
                'FontSize', 30, 'interpreter', 'latex')

    else

        title(['With Advection - $LWC \geq$ ', num2str(vocalsRex.inputs.LWC_threshold),' $g/m^{3}$',...
            '     $N_c \geq$ ', num2str(vocalsRex.inputs.Nc_threshold), ' $cm^{-3}$',...
            '     Max vert displacement: ', num2str(vocalsRex.inputs.max_vert_displacement), ' $m$'],...
                'FontSize', 30, 'interpreter', 'latex')

    end

end

% Include an x axis label on the middle plot
if normalize_distance==true

    xlabel('Normalized Horizontal Distance Travelled', 'Interpreter','latex');
else

    xlabel('Horizontal Distance Travelled ($km$)', 'Interpreter','latex');
end




% in the third subplot, define the indices_2_plot being plotted
legend(legend_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20)

% set plot size
set(gcf, 'Position', [0 0 1200 625])


% Create textbox
annotation('textbox',[0.781833333333335 0.008 0.0864999999999987 0.048],...
    'String',{[filename(6:9),'/',filename(10:11), '/', filename(12:13)]},...
    'LineWidth',2,...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');








%% Plot the VOCALS-REX flight path and the overlapping MODIS pixels

 

figure; 


% load the across and along pixel growth curve fits
load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/',...
    'MODIS_Cloud_Retrieval/along_across_pixel_growth.mat'])

% Assume these are the center coordinates of the pixel
% First, create a perfect square assuming, when MODIS is looking nadir,
% that the pixel is 1km by 1 km on the ground.
% Thus, the diagonal from the center to any corner is 1/sqrt(2) km long
%dist_to_corners = 1/sqrt(2);        % km

% Using the 'reckon' function, we con move to each corner from the center
% position. But we need to know the orientation of the MODIS instrument
% with respect to the meriodional lines.

terra_inclination = 98.2098;        % degrees
aqua_inclination = 	98.1987;        % degrees

azimuth_2_corners = [45, 135, 225, 315, 45];


% grab the indexes of the ordered set of effective radius for the data
% being used
data2plot = modis.cloud.effRadius17(vocalsRex.modisIndex_minDist);
% Set NaNs to 0
data2plot(isnan(data2plot))=0;

[~, index_sort] = sort(reshape(data2plot, [],1), 'ascend');
%[r,c] = ind2sub([n_rows, n_cols], index_sort);
C = parula(length(index_sort));


for ii = 1:length(data2plot)

        % the distance to the corners depends on the viewing zenith angle
        along_length = along_scan(modis.sensor.zenith(vocalsRex.modisIndex_minDist(ii)));
        across_length = across_scan(modis.sensor.zenith(vocalsRex.modisIndex_minDist(ii)));

        % compute the distance from the center to each corner
        dist_to_corner = sqrt((along_length/2)^2 + (across_length/2)^2) * 1e3;  % meters


        % The azimuth angle in the reckon function is with respect to the local
        % meridian (north). So we have to add the inclination angle to the azimuth
        % angles listed above.

        [lat_corner, long_corner] = reckon(modis.geo.lat(vocalsRex.modisIndex_minDist(ii)),...
            modis.geo.long(vocalsRex.modisIndex_minDist(ii)),linspace(dist_to_corner, dist_to_corner,5),...
            azimuth_2_corners+terra_inclination, wgs84Ellipsoid);

        % Create a geopolyshape
        modis_polyshape = geopolyshape(lat_corner, long_corner);

        gp = geoplot(modis_polyshape);
        gp.EdgeAlpha = 0;
        % The color corresponds to the linear index_sort
        gp.FaceColor = C(ii==index_sort,:);
        gpFaceAlpha = 0.9;

        

        hold on


    
end

% Plot the Vocals-Rex data
geoscatter(vocalsRex.latitude, vocalsRex.longitude, 10, "red",'*')


cb = colorbar;
%clim([min(modis.cloud.effRadius17(1:n_rows, 1:n_cols)), max(modis.cloud.effRadius17(1:n_rows, 1:n_cols))])
clim([min(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist), [], 'all'),...
    max(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist), [], 'all')])

set(get(cb, 'label'), 'string', '$r_e \; (\mu m)$','Interpreter','latex', 'Fontsize',22)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])

title(['MODIS Effective Radius - ', modisFolder(113:end-1)],'Interpreter','latex', 'FontSize', 35)



%% Compute the Cross Correlation between the MODIS retrieval of re and VOCALS-Rex

% To make the two signals the same length, let's find the Vocals-Rex data
% points closest to the MODIS pixel values used

for nn = 1:vocalsRex.modisIndex_minDist

    % For each MODIS Pixel, use the VOCALS-REx data point closest to it





end