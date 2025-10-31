%% Read in ERA5 Reanalysis data and compare with VOCALS-REx in-situ microphsyics measurements


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    foldername = [];

    % ***** Define the VOCALS-REx File *****

    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/ERA5_reanalysis/ERA5_data/'];

elseif strcmp(which_computer,'andrewbuggee')==true



end

filename = '2008_11_Nov_VOCALS_region.nc';


%%

% read in vocals rex radiosonde data - all data occured in 2008

info = ncinfo([foldername, filename]);

% read the time
time = double(ncread([foldername, filename], 'valid_time'));

% read the lat, long
lat = ncread([foldername, filename], 'latitude');                                                        % Meausred in degrees North
long = ncread([foldername, filename], 'longitude');

% Extract temperature data from the netCDF file
temperature = ncread([foldername, filename], 't');                               % K

% Extract cloud liquid water content data from the netCDF file
clwc_specific = ncread([foldername, filename], 'clwc');                           % kg of water droplets / kg of total mass of moist air - The 'total mass of moist air' is the sum of the dry air, water vapour, cloud liquid, cloud ice, rain and falling snow.

% Extract rain liquid water content data from the netCDF file
rwc_specific = ncread([foldername, filename], 'crwc');

% extract pressure level, the independent variable
pressure = ncread([foldername, filename], 'pressure_level');                       % hPa


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

%% Plot just the radiosonde data

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
parfor nn = 1:length(vert_prof)

  


    dist_btwn_era5_and_VR = distance(lat, long, vert_prof(nn).latitude(round(end/2)),...
        vert_prof(nn).longitude(round(end/2)), wgs84);

    [era5_minDist(nn), index_minDist] = min(dist_btwn_era5_and_VR, [], 'all');            % m - minimum distance


    % compute the time between the modis pixel closest to the VR sampled
    % path and the time VOCALS was recorded
    time_diff_era5_VR(nn) = abs(modis_pixel_time(index_minDist) - ...
        vert_prof(nn).time_utc(round(end/2))) * 60;                         % minutes

    % A single MODIS granule takes 5 minutes to be collected. First, lets
    % determine if any profiles were recorded within the the 5 minute window
    % when MODIS collected data.

    within_5min_window(nn) = any(vert_prof(nn).time_utc >= modis_pixel_time(1,1) & ...
        vert_prof(nn).time_utc <= modis_pixel_time(end,end));

end


%% Let's plot an example vertical profile from the same day as one of my retrievals from paper 1
% find the radiosonde closest in time

% step through each vertical profile and plot the radiosonde measurements
% closest in time

fnt_sz = 26;
lgnd_fnt = 20;

rdoSnde_clr = 'k';
C130_clr = 61;

% Determine the time of the radiosonde launches
launchTimes = datetime(2008, month, day, hour, minute, 0);


for nn = 1:length(vert_prof)

    % Find the index of the closest launch time to a specific retrieval time
    targetTime = datetime(2008, 11, 11, floor(vert_prof(nn).time_utc(1)),...
        round(60*(vert_prof(nn).time_utc(1) - floor(vert_prof(nn).time_utc(1)))), 0); % Example target time
    [~, closestIndex] = min(abs(launchTimes - targetTime));


    % plot the temperature, pressure and relative humidity as a function of
    % height for the index determined above


    % Create a figure for the vertical profile
    figure;
    hold on;

    % Plot temperature
    subplot(2,3,1)
    plot(temperature(:, closestIndex), altitude, 'Color', rdoSnde_clr)
    xlabel('Temperature ($C$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    grid on; grid minor

    hold on
    plot(vert_prof(nn).temp, vert_prof(nn).altitude, 'Color', mySavedColors(C130_clr, 'fixed'))
    legend('Radiosonde', 'Aircraft','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')




    % Plot pressure
    subplot(2,3,2)
    plot(pressure(:, closestIndex), altitude,'Color', rdoSnde_clr);
    xlabel('Pressure ($mb$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    grid on; grid minor

    hold on
    plot(vert_prof(nn).pres, vert_prof(nn).altitude, 'Color', mySavedColors(C130_clr, 'fixed'))
    legend('Radiosonde', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')


    title(['Radiosonde - (', num2str(lat(closestIndex)), ',', num2str(long(closestIndex)),...
        ') - ', num2str(month(closestIndex)), '/', num2str(day(closestIndex)), '/2008 - ',...
        num2str(hour(closestIndex)), ':', num2str(minute(closestIndex)), ' UTC'], ...
        'FontSize', 20, 'Interpreter', 'latex')




    % Plot relative humidity
    subplot(2,3,3)
    plot(RH(:, closestIndex), altitude, 'Color', rdoSnde_clr);
    xlabel('Relative Humidity', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    grid on; grid minor

    hold on
    plot(vert_prof(nn).relative_humidity, vert_prof(nn).altitude, 'Color', mySavedColors(C130_clr, 'fixed'))
    legend('Radiosonde', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')





    % Plot the effective radius
    subplot(2,3,4)
    plot(vert_prof(nn).re, vert_prof(nn).altitude, 'Color', 'k');
    xlabel('Effective Radius ($\mu m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
    grid on; grid minor
    % ylim([0, altitude(end)])


    % Show the position of the radiosonde launch and the airplane when
    % sampling the cloud
    subplot(2,3,5)
    % Plot the location of the radiosonde launch
    geoscatter(lat(closestIndex), long(closestIndex), 100, 'k', 'filled', 'DisplayName', 'Radiosonde Launch');
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





