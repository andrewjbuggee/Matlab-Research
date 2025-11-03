%% Read in Vocals Rex Radiosonde data


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    foldername = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/'];

    % ***** Define the VOCALS-REx File *****

    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/SPS_1/'];

elseif strcmp(which_computer,'andrewbuggee')==true


    foldername = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/'];


    % ***** Define the VOCALS-REx File *****

    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


end

filename = 'VOCALS2008_soundings_z_v4.2.nc';


%%

% read in vocals rex radiosonde data - all data occured in 2008

info = ncinfo([foldername, filename]);

% read the month, day, hour and minute of each launch
month = double(ncread(filename, 'month'));
day = double(ncread(filename, 'day'));
hour = double(ncread(filename, 'hour'));            % UTC hour of balloon launch
minute = double(ncread(filename, 'minute'));

% read the lat, long
lat = ncread(filename, 'lat');                                                        % Meausred in degrees North
long = ncread(filename, 'lon');

% Extract temperature and pressure data from the netCDF file
temperature = double(ncread(filename, 'T'));
pressure = double(ncread(filename, 'pres'));

% extract altitude from the netCDF file
altitude = double(ncread(filename, 'height'));

% extract relative humidity
RH = double(ncread(filename, 'RH'));



%% Plot just the radiosonde data

% Create a figure for the vertical profile
figure;
hold on;

plt_idx = 10;
fnt_sz = 20;

% Plot temperature
subplot(1,3,1)
plot(temperature(:, plt_idx), altitude, 'Color', mySavedColors(61,'fixed'));
xlabel('Temperature ($C$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
grid on; grid minor


% Plot pressure
subplot(1,3,2)
plot(pressure(:, plt_idx), altitude,'Color', mySavedColors(62,'fixed'));
xlabel('Pressure ($mb$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
grid on; grid minor
title(['Radiosonde - (', num2str(lat(plt_idx)), ',', num2str(long(plt_idx)),...
    ') - ', num2str(month(plt_idx)), '/', num2str(day(plt_idx)), '/2008 - ',...
    num2str(hour(plt_idx)), ':', num2str(minute(plt_idx)), ' UTC'], ...
    'FontSize', 20, 'Interpreter', 'latex')

% Plot relative humidity
subplot(1,3,3)
plot(RH(:, plt_idx), altitude, 'Color', mySavedColors(63,'fixed'));
xlabel('Relative Humidity', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
grid on; grid minor



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





