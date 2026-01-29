%% Step through each VOCAL-REx vertical profile and find the closest ERA5 data set in space and time


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    % ***** Define the ERA5 data directory *****
    folderpath_era5 = [''];


    % ***** Define the ensemble profiles folder *****
    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

elseif strcmp(which_computer,'andrewbuggee')==true


    % ***** Define the ERA5 data directory *****
    folderpath_era5 = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/'];



    % ***** Define the ensemble profiles folder *****
    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


end



% ---------------------
% define the ensemble filename
% profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat';
profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_04-Dec-2025.mat';



%% Load the ensemble profiles

load([vocalsRexFolder, profiles]);





%% Step through each profile and find the radiosonde closest in space and time
% The radiosonde should be over the ocean, rather than over land


% the composite radiosonde data is saved with the day it was recorded
% on. So first find the radiosonde file for the day the current
% ensemble profile was recorded


total_column_pw = zeros(length(ensemble_profiles), 1);
above_cloud_pw = zeros(length(ensemble_profiles), 1);

cloud_topHeight = zeros(length(ensemble_profiles), 1);
cloud_baseHeight = zeros(length(ensemble_profiles), 1);

cloud_topPressure = zeros(length(ensemble_profiles), 1);
cloud_basePressure = zeros(length(ensemble_profiles), 1);

temp_prof = cell(length(ensemble_profiles), 1);
watVap_prof = cell(length(ensemble_profiles), 1);
pressure_prof = cell(length(ensemble_profiles), 1);
altitude = cell(length(ensemble_profiles), 1);





for nn = 1:length(ensemble_profiles)


    % extract the date of the nth profile
    date_profile = ensemble_profiles{nn}.dateOfFlight;

    % determine the month
    m = month(date_profile, 'name');
    era5_month_folderpath = [m{1}, '_2008/'];

    % define the filename using the day given by vocals rex
    d = day(date_profile, "dayofmonth");
    era5_filename = ['era5_vocalsrex_' m{1}, '_2008_day', num2str(d), '.nc'];




    % ----------------------------------------------------------
    % ----------------- read in ERA5 data ----------------------
    % ----------------------------------------------------------

    info = ncinfo([folderpath_era5, era5_month_folderpath, era5_filename]);

    % read the time
    time = double(ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'valid_time'));            % seconds since 1970

    % Convert to UTC datetime
    utcTime = datetime(time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

    % ----------------------------------------------------------
    % ----------------------------------------------------------




    % ----------------------------------------------------------
    % ----- First, find the closest time to the VR profile -----
    % ----------------------------------------------------------

    % Extract that time at the middle of the in-situ profile
    vr_time_mid = ensemble_profiles{nn}.time_utc( round( length(ensemble_profiles{nn}.time_utc)/2 ) );  % UTC decimal time
    date_profile.Hour = floor(vr_time_mid);  % The hour
    date_profile.Minute = floor(60*(vr_time_mid - floor(vr_time_mid)));  % The minute
    % convert this to UTC time
    date_profile.TimeZone = "UTC";

    % Now, find the ERA5 data set closest in time
    [min_timeDiff, idx_time] = min( abs( utcTime - date_profile ));



    % -----------------------------------------------------------------------------
    % ----- Next, find the ERA5 data point closest in space to the VR profile -----
    % -----------------------------------------------------------------------------

    % It's time to make a 2D mesh grid using the 4 independent variables:
    % lat, long, presssure, and time

    % read the ERA5 lat, long
    lat = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'latitude');                                                        % Meausred in degrees North
    long = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'longitude');

    % create the 2D meshgrid
    [Lat, Lon] = meshgrid(lat, long);

    % we will be computing the arclength between points on an ellipsoid
    % Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with units of meters.
    wgs84 = wgs84Ellipsoid("m");

    % Compute the distance between the mid point and all ERA5 data points
    % in meters
    dist_btwn_era5_and_VR = distance(Lat, Lon, ensemble_profiles{nn}.latitude(round(end/2)),...
        ensemble_profiles{nn}.longitude(round(end/2)), wgs84);

    [min_spatialDiff, idx_minDist] = min(dist_btwn_era5_and_VR, [], 'all');            % m - minimum distance
    [r,c] = ind2sub(size(dist_btwn_era5_and_VR), idx_minDist);



    % -----------------------------------------------------------------------------------
    % ----- With the minimum idx's save the ERA5 data in a cell array for each prof -----
    % -----------------------------------------------------------------------------------

    % extract pressure levels
    pressure = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'pressure_level');            % hPa

    % Extract temperature data from the netCDF file
    temperature = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 't');                      % K
    % grab only the data you need and delete the extra
    T = reshape( temperature(r,c, :, idx_time), [], 1);
    clear temperature

    % Extract cloud liquid water content data from the netCDF file
    clwc_specific = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'clwc');               % kg of water droplets / kg of total mass of moist air - The 'total mass of moist air' is the sum of the dry air, water vapour, cloud liquid, cloud ice, rain and falling snow.
    % grab only the data you need and delete the extra
    clwc_2keep = reshape( clwc_specific(r, c, :, idx_time), [], 1);
    clear clwc_specific

    % extract relative humidity
    RH = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'r');
    % grab only the data you need and delete the extra
    RH_2keep = reshape( RH(r, c, :, idx_time), [], 1);
    clear RH

    % extract specific humidity
    q = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'q');
    % grab only the data you need and delete the extra
    q_2keep = reshape( q(r, c, :, idx_time), [], 1);
    clear q


    % Compute cloud liquid water content from specific cloud LWC
    clwc_2keep = (clwc_2keep .* p) ./...
        (computeMoistAirGasConstant(q_2keep) .* T);        % kg of liquid water droplets/ m^3 of air
    
    % convert to grams per m^3
    clwc_2keep = clwc_2keep * 1000; % Convert to grams per m^3




end
















%%

% read in ERA5 data

info = ncinfo([foldername, filename]);

% read the time
time = double(ncread([foldername, filename], 'valid_time'));            % seconds since 1970
% Convert to UTC datetime
utcTime = datetime(time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% read the lat, long
lat = ncread([foldername, filename], 'latitude');                                                        % Meausred in degrees North
long = ncread([foldername, filename], 'longitude');

% Extract temperature data from the netCDF file
temperature = ncread([foldername, filename], 't');                      % K

% Extract cloud liquid water content data from the netCDF file
clwc_specific = ncread([foldername, filename], 'clwc');               % kg of water droplets / kg of total mass of moist air - The 'total mass of moist air' is the sum of the dry air, water vapour, cloud liquid, cloud ice, rain and falling snow.


% extract pressure level, the independent variable
pressure = ncread([foldername, filename], 'pressure_level');            % hPa


% extract relative humidity
RH = ncread([foldername, filename], 'r');

% extract specific humidity
q = ncread([foldername, filename], 'q');


%% Compute cloud liquid water content from specific LWC

clwc = clwc_specific .* repmat(reshape(pressure.*100, 1,1, [], 1), length(lat), length(long), 1, length(time)) ./...
    (computeMoistAirGasConstant(q) .* temperature);        % kg of liquid water droplets/ m^3 of air

% convert to grams per m^3
clwc = clwc * 1000; % Convert to grams per m^3

%% Convert temperature from Kelvin to Celcius

temperature = temperature - 273.15; % Convert temperature to Celsius

%% Plot just the ERA5 data

% Create a figure for the vertical profile
figure;
hold on;

time_idx = 10;
lat_idx = 10;
long_idx = 10;

fnt_sz = 20;



% Plot temperature
subplot(1,3,1)
semilogy(reshape(temperature(lat_idx, long_idx, :, time_idx), [], 1),...
    pressure, 'Color', mySavedColors(61,'fixed'));
xlabel('Temperature ($C$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
grid on; grid minor

% flip y axis
set(gca, 'YDir', 'reverse')




% Plot specific cloud liquid water content
subplot(1,3,2)
semilogy(reshape(clwc(lat_idx, long_idx, :, time_idx), [], 1),...
    pressure, 'Color', mySavedColors(61,'fixed'));
% plot(reshape(clwc_specific(lat_idx, long_idx, :, time_idx), [], 1),...
%     pressure, 'Color', mySavedColors(61,'fixed'));
xlabel('cloud LWC ($g/m^3$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)

grid on; grid minor
% flip y axis
set(gca, 'YDir', 'reverse')



% Plot relative humidity
subplot(1,3,3)
semilogy(reshape(RH(lat_idx, long_idx, :, time_idx), [], 1),...
    pressure, 'Color', mySavedColors(61,'fixed'));
xlabel('Relative Humidity ($\%$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
grid on; grid minor

% flip y axis
set(gca, 'YDir', 'reverse')




%% LOAD VOCALS-REx data in-situ aircraft data


% ----- November 9 data -----
%vocalsRexFile = 'RF11.20081109.125700_213600.PNI.nc';


% ----- November 11 data -----
vocalsRexFile = 'RF12.20081111.125000_214500.PNI.nc';




% ---------------------------------------------------------------------
% ------------ Do you want to load VOCALS-REx data? -------------------
% ---------------------------------------------------------------------


loadVOCALSdata = true;

% % ---------------------------------------------------------------------
% % ----------- Do you want to use VOCALS-REx pixels? -------------------
% % ---------------------------------------------------------------------
%
% % If flag below is true, only suitable pixels within the VOCALS-REx data set
% % will be used
% useVOCALS_pixelLocations = true;



if loadVOCALSdata==true
    vocalsRex = readVocalsRex([vocalsRexFolder,vocalsRexFile]);

end


%% Find vertical profiles

% ------------------------------------------------
% ----- We dont need all of this data ------------
% ------------------------------------------------

% ----- Find all vertical profiles within VOCALS-REx data ------
% find all vertical profiles in the vocals-REx data set. Vertical profiles
% qualify if the total number concentration is above 10^0 and relatively
% flat, the plane is climbing in altitude, and clears the cloud deck from
% bottom to top. If this is all true, save the vocals rex data
lwc_threshold = 0.03;           % g/m^3
stop_at_max_lwc = false;         % truncate profile after the maximum lwc value
Nc_threshold = 25;               % # droplets/cm^3

% ----- Find all vertical profiles within VOCALS-REx data ------
vert_prof = find_verticalProfiles_VOCALS_REx_ver2(vocalsRex, lwc_threshold, stop_at_max_lwc,...
    Nc_threshold, whatComputer);

clear vocalsRex


%% First find the reanalysis data point closest to the VOCALS-Rex measurement

% Let's step through each vertical profile and find the closest point to
% the reanalysis data



% we will be computing the arclength between points on an ellipsoid
% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with units of meters.
wgs84 = wgs84Ellipsoid("m");

% Set up an empty array for each vocals-rex profile
era5_minDist = zeros(1, length(vert_prof));
time_diff_era5_VR = zeros(1, length(vert_prof));

% Step through each vertical profile and find the MODIS pixel that overlaps
% with the mid point of the in-situ sample
for nn = 1:length(vert_prof)




    dist_btwn_era5_and_VR = distance(lat, long, vert_prof(nn).latitude(round(end/2)),...
        vert_prof(nn).longitude(round(end/2)), wgs84);

    [era5_minDist(nn), index_minDist] = min(dist_btwn_era5_and_VR, [], 'all');            % m - minimum distance

    %
    % % compute the time between the modis pixel closest to the VR sampled
    % % path and the time VOCALS was recorded
    % time_diff_era5_VR(nn) = abs(utcTime - vert_prof(nn).time_utc(round(end/2))) * 60;                         % minutes


end


%% Let's plot an example vertical profile from the same day as one of my retrievals from paper 1
% find the radiosonde closest in time

% step through each vertical profile and plot the radiosonde measurements
% closest in time

fnt_sz = 26;
lgnd_fnt = 20;

era5_clr = 'k';
C130_clr = 61;



for nn = 1:length(vert_prof)

    % Find the index of the closest launch time to a specific retrieval time
    targetTime = datetime(2008, 11, 11, floor(vert_prof(nn).time_utc(1)),...
        round(60*(vert_prof(nn).time_utc(1) - floor(vert_prof(nn).time_utc(1)))), 0,...
        'TimeZone', 'UTC'); % Example target time

    [~, closestIndex] = min(abs(utcTime - targetTime));


    % plot the temperature, pressure and relative humidity as a function of
    % height for the index determined above


    % Create a figure for the vertical profile
    figure;
    hold on;

    % Plot temperature
    subplot(2,3,1)
    semilogy(reshape(temperature(lat_idx, long_idx, :, closestIndex), [], 1),...
        pressure, 'Color', 'k');
    xlabel('Temperature ($C$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)

    grid on; grid minor
    % flip y axis
    set(gca, 'YDir', 'reverse')

    hold on
    plot(vert_prof(nn).temp, vert_prof(nn).pres, 'Color', mySavedColors(C130_clr, 'fixed'))
    legend('ERA-5', 'Aircraft','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
        'Color', 'white', 'TextColor', 'k')




    % Plot CLWC
    subplot(2,3,2)
    semilogy(reshape(clwc(lat_idx, long_idx, :, closestIndex), [], 1),...
        pressure, 'Color', era5_clr);
    xlabel('cloud LWC ($g/m^3$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)

    grid on; grid minor
    % flip y axis
    set(gca, 'YDir', 'reverse')

    hold on
    plot(vert_prof(nn).lwc, vert_prof(nn).pres, 'Color', mySavedColors(C130_clr, 'fixed'))
    legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
        'Color', 'white', 'TextColor', 'k')


    % title(['Radiosonde - (', num2str(lat(closestIndex)), ',', num2str(long(closestIndex)),...
    %     ') - ', num2str(month(closestIndex)), '/', num2str(day(closestIndex)), '/2008 - ',...
    %     num2str(hour(closestIndex)), ':', num2str(minute(closestIndex)), ' UTC'], ...
    %     'FontSize', 20, 'Interpreter', 'latex')




    % Plot relative humidity
    subplot(2,3,3)
    semilogy(reshape(RH(lat_idx, long_idx, :, closestIndex), [], 1),...
        pressure, 'Color', era5_clr);
    xlabel('Relative Humidity ($\%$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    grid on; grid minor

    % flip y axis
    set(gca, 'YDir', 'reverse')

    hold on
    plot(vert_prof(nn).relative_humidity, vert_prof(nn).pres, 'Color', mySavedColors(C130_clr, 'fixed'))
    legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
        'Color', 'white', 'TextColor', 'k')





    % Plot the effective radius
    subplot(2,3,4)
    plot(vert_prof(nn).re, vert_prof(nn).altitude, 'Color', mySavedColors(C130_clr, 'fixed'));
    xlabel('Effective Radius ($\mu m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    title('In-situ droplet size - aircraft', 'Interpreter','latex', 'FontSize', fnt_sz)
    grid on; grid minor
    % ylim([0, altitude(end)])


    % Show the position of the radiosonde launch and the airplane when
    % sampling the cloud
    subplot(2,3,5)
    % Plot the location of the radiosonde launch
    geoscatter(lat(closestIndex), long(closestIndex), 100, 'k', 'filled', 'DisplayName', 'ERA-5 Reanalysis');
    hold on;

    % Plot the location of the sampled cloud
    geoscatter(vert_prof(nn).latitude(1), vert_prof(nn).longitude(1), 100, 'b', 'filled', 'DisplayName', 'Sampled Cloud');
    legend('show', 'Location', 'best', 'FontSize', lgnd_fnt, 'Interpreter', 'latex');


    set(gcf, 'Position', [0,0, 2200, 1000])





    % ---- Second plot --------
    % Make geoscatter plot showing the location of the radiosonde launch and
    % the location of the sampled cloud (vertical profile)
    % figure;
    % geoscatter()


end





