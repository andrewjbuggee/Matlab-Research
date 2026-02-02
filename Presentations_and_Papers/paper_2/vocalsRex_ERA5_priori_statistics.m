%% Step through each VOCAL-REx vertical profile and find the closest ERA5 data set in space and time


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    % ***** Define the ERA5 data directory *****
    folder_paths.era5_pres_lvl = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/'];


    % ***** Define the ensemble profiles folder *****
    folder_paths.vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];



elseif strcmp(which_computer,'andrewbuggee')==true


    % ***** Define the ERA5 pressure level data directory *****
    folder_paths.era5_pres_lvl = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/'];


    % ***** Define the ERA5 single level data directory *****
    folder_paths.era5_single_lvl = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/single_level_instantaneous/'];



    % ***** Define the ensemble profiles folder *****
    folder_paths.vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


end



% ---------------------
% define the ensemble filename
folder_paths.VR_profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat';
% profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_04-Dec-2025.mat';



%% Load the ensemble profiles

load([folder_paths.vocalsRexFolder, folder_paths.VR_profiles]);



%% Some notes on the ERA5 variables

% Specific humidity (q) is the mass of water vapor per unit mass of moist air -
% typically expressed in kg/kg or g/kg. It's one of several ways to express
% atmospheric moisture content, and it has the advantage of being conserved
%  during adiabatic processes (unlike relative humidity).



%% Step through each profile and find the radiosonde closest in space and time
% The radiosonde should be over the ocean, rather than over land


% the composite radiosonde data is saved with the day it was recorded
% on. So first find the radiosonde file for the day the current
% ensemble profile was recorded


total_column_pw = zeros(length(ensemble_profiles), 1);
above_cloud_pw = zeros(length(ensemble_profiles), 1);
total_column_pw_SL = zeros(length(ensemble_profiles), 1);
above_cloud_pw_usingVR = zeros(length(ensemble_profiles), 1);

cloud_topHeight = zeros(length(ensemble_profiles), 1);
cloud_baseHeight = zeros(length(ensemble_profiles), 1);
cloud_baseHeight_SL = zeros(length(ensemble_profiles), 1);
cloud_topHeight_usingVR = zeros(length(ensemble_profiles), 1);

cloud_topPressure = zeros(length(ensemble_profiles), 1);
cloud_basePressure = zeros(length(ensemble_profiles), 1);

temp_prof = cell(length(ensemble_profiles), 1);
RH_prof = cell(length(ensemble_profiles), 1);
pressure_prof = [];
altitude_prof = cell(length(ensemble_profiles), 1);
rho_v_prof = cell(length(ensemble_profiles), 1);
watVap_concentration_prof = cell(length(ensemble_profiles), 1);
clwc_prof = cell(length(ensemble_profiles), 1);

dist_2_VR = zeros(length(ensemble_profiles), 1);
time_2_VR = duration.empty(length(ensemble_profiles), 0);

dist_2_VR_SL = zeros(length(ensemble_profiles), 1);
time_2_VR_SL = duration.empty(length(ensemble_profiles), 0);


test = zeros(length(ensemble_profiles), 1);



for nn = 1:length(ensemble_profiles)

    disp([newline, 'Processing profile ', num2str(nn), '...', newline])


    % extract the date of the nth profile
    date_profile = ensemble_profiles{nn}.dateOfFlight;

    % determine the month
    m = month(date_profile, 'name');
    era5_month_folderpath = [m{1}, '_2008/'];

    % define the filename using the day given by vocals rex
    d = day(date_profile, "dayofmonth");
    era5_filename = ['era5_vocalsrex_' m{1}, '_2008_day', num2str(d), '.nc'];




   
    % ----------------------------------------------------------
    % ------------ read in ERA5 single level data --------------
    % ----------------------------------------------------------
    
    if strcmp(m{1}, 'October') == true

        era5_singleLevel_filename = 'Oct_2008_15-31.nc';

    elseif strcmp(m{1}, 'November') == true
            
        era5_singleLevel_filename = 'Nov_2008_2-15.nc';

    else

        error([newline, 'I only have single elvel ERA5 data for Oct. and Nov.', newline])

    end

    info_singleLevel = ncinfo([folder_paths.era5_single_lvl, era5_singleLevel_filename]);

    % read the time
    sl.time = double(ncread([folder_paths.era5_single_lvl, era5_singleLevel_filename], 'valid_time'));            % seconds since 1970

    % Convert to UTC datetime
    sl.utcTime = datetime(sl.time , 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

    % ----------------------------------------------------------
    % ----------------------------------------------------------





    % ----------------------------------------------------------
    % ---------- read in ERA5 pressure level data --------------
    % ----------------------------------------------------------

    info = ncinfo([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename]);

    % read the time
    time = double(ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'valid_time'));            % seconds since 1970

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




    % ----------------------------------------------------------
    % ----- Do the same but for the single level data -----
    % ----------------------------------------------------------

    % Now, find the ERA5 data set closest in time
    [sl.min_timeDiff, sl.idx_time] = min( abs( sl.utcTime - date_profile ));



    % -----------------------------------------------------------------------------
    % ----- Next, find the ERA5 data point closest in space to the VR profile -----
    % -----------------------------------------------------------------------------

    % It's time to make a 2D mesh grid using the 4 independent variables:
    % lat, long, presssure, and time

    % read the ERA5 lat, long
    lat = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'latitude');                                                        % Meausred in degrees North
    long = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'longitude');

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
    % -----------------------------------------------------------------------------


    % -----------------------------------------------------------------------------
    % ------------ Do this again but for the ERA5 single level data ---------------
    % -----------------------------------------------------------------------------
    % read the ERA5 lat, long
    sl.lat = ncread([folder_paths.era5_single_lvl, era5_singleLevel_filename], 'latitude');                                                        % Meausred in degrees North
    sl.long = ncread([folder_paths.era5_single_lvl, era5_singleLevel_filename], 'longitude');

    % create the 2D meshgrid
    [sl.Lat, sl.Lon] = meshgrid( sl.lat, sl.long);

    % Compute the distance between the mid point and all ERA5 data points
    % in meters
    sl.dist_btwn_era5_and_VR = distance(sl.Lat, sl.Lon, ensemble_profiles{nn}.latitude(round(end/2)),...
        ensemble_profiles{nn}.longitude(round(end/2)), wgs84);

    [sl.min_spatialDiff, sl.idx_minDist] = min(sl.dist_btwn_era5_and_VR, [], 'all');            % m - minimum distance
    [sl.r,sl.c] = ind2sub(size(sl.dist_btwn_era5_and_VR), sl.idx_minDist);
    % -----------------------------------------------------------------------------






    % -----------------------------------------------------------------------------------
    % ----- With the minimum idx's save the ERA5 data in a cell array for each prof -----
    % -----------------------------------------------------------------------------------

    % extract pressure levels
    p = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'pressure_level');            % hPa

    % Extract temperature data from the netCDF file
    temperature = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 't');                      % K
    % grab only the data you need and delete the extra
    T = reshape( temperature(r,c, :, idx_time), [], 1);
    clear temperature

    % Extract cloud liquid water content data from the netCDF file
    clwc_specific = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'clwc');               % kg of water droplets / kg of total mass of moist air - the mass of cloud liquid water droplets per kilogram of the total mass of moist air. The 'total mass of moist air' is the sum of the dry air, water vapour, cloud liquid, cloud ice, rain and falling snow.
    % grab only the data you need and delete the extra
    clwc_specific_2keep = reshape( clwc_specific(r, c, :, idx_time), [], 1);
    clear clwc_specific

    % extract relative humidity
    RH = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'r');
    % grab only the data you need and delete the extra
    RH_2keep = reshape( RH(r, c, :, idx_time), [], 1);
    clear RH

    % extract specific humidity
    q = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'q');    % [kg/kg]
    % grab only the data you need and delete the extra
    q_2keep = reshape( q(r, c, :, idx_time), [], 1);
    clear q


    % ----------------------------------------------------------------------
    % ----------------------------------------------------------------------
    % ***** Compute cloud liquid water content from specific cloud LWC *****
    % ----------------------------------------------------------------------
    % convert pressure to pascals

    clwc = specificCloudLiquidWater2LWC(clwc_specific_2keep, T, q_2keep, p.*100);      % [kg/m³]

    % convert to g/m^3
    clwc = clwc .* 1000; % Convert to grams per m^3

    % -----------------------------------------------------------------------
    % ----------------------------------------------------------------------


    % -----------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------
    % ---------------- compute the above cloud total precipitable water -----------------
    % -----------------------------------------------------------------------------------
    % Compute the mass of water vapor per unit volume at all altitudes
    % convert pressure from hPa to Pa
    use_virtual_temp = true;   % more accurate estimate, but the improvement is on the order of 1%
    % rho_v = [kg/m^3]
    % waterVapor_concentration_cm3 = [#/cm^3]
    [rho_v, waterVapor_concentration_cm3] = specificHumidity2waterVaporDensity(q_2keep, T, p.*100, use_virtual_temp);


    % --------------------------------------------------
    % % compute the geopotential height for each layer
    % --------------------------------------------------
    % ** my estiamtes vary from the ERA 5 estimates. Use theirs **
    % % we are only using data over ocean, so set Z_surface to 0
    % Z_sfc = 0;         % meters
    % Z = computeGeopotentialHeight(T, q_2keep, p, Z_sfc);   % meters

    % read in the geopotential and convert it to geopotential height
    Z = ncread([folder_paths.era5_pres_lvl, era5_month_folderpath, era5_filename], 'z');            % m^2/s^2
    Z = reshape( Z(r, c, :, idx_time), [], 1) ./ 9.80665;       % meters

    % read in the era5 single level surface pressure
    p0 = ncread([folder_paths.era5_single_lvl, era5_singleLevel_filename], 'sp');   % Pa
    p0 = p0(sl.r, sl.c, sl.idx_time);       % Pa



    % From ERA5's website, convert geopotential height to geometric height
    % (https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height)
    % This is likely to be very similar since geopotential height is close
    % to the geometric height in the lower atmosphere
    con = physical_constants;
    alt = con.R_earth * Z ./(con.R_earth - Z);       % meters - geometric height
    alt_with0 = [0; alt];

    % check total column amount
    total_column_pw(nn) = trapz(alt, rho_v);            % kg / m^2

    % total_column_pw(nn) = -1/con.g0 .* trapz(p.*100, q_2keep);     % kg / m^2 - basically gives me the same answer
    % -----------------------------------------------------------------------------------



    % ----------------------------------------------------------
    % Now extract the single level integrated column water vapor
    % ----------------------------------------------------------
    tot_col_pw_SL = ncread([folder_paths.era5_single_lvl, era5_singleLevel_filename], 'tcwv');   % kg/m^2
    total_column_pw_SL(nn) = tot_col_pw_SL( sl.r, sl.c, sl.idx_time);    % kg/m^2
    clear tot_col_pw_SL

    % ----------------------------------------------------------


    % ----------------------------------------------------------
    % ------ Extract the single level cloud base height --------
    % ----------------------------------------------------------
    cb_height = ncread([folder_paths.era5_single_lvl, era5_singleLevel_filename], 'cbh');   % meters
    cloud_baseHeight_SL(nn) = cb_height( sl.r, sl.c, sl.idx_time);    % m
    clear cb_height

    % ----------------------------------------------------------


    % Find the altitude level where CLWC is greater than 0
    is_cloud_idx = find(clwc > 0);

    if isempty(is_cloud_idx) == true

        error([newline, 'No cloud found for profile ', num2str(nn), newline])

    end


    % compute the above cloud precipitable water amount
    above_cloud_pw(nn) = trapz( alt( (is_cloud_idx(end) +1) : end ), rho_v( (is_cloud_idx(end) +1) : end ));            % kg / m^2

    % ** compute above cloud pw using the VR cloud top height **
    % ERA5 cloud top height matches VR cloud top heights for values below
    % 1500 meters. Above that, ERA 5 appears to overestimate cloud top
    % height
    [~, idx_cloudTop_VR] = min( abs( alt - max(ensemble_profiles{nn}.altitude) )); 
    above_cloud_pw_usingVR(nn) = trapz( alt( idx_cloudTop_VR : end ), rho_v( idx_cloudTop_VR : end ));            % kg / m^2
    % -----------------------------------------------------------------------------------


    % ----------------------------------------------------
    % ----------------------------------------------------
    % Define the cloud top and base height and pressure
    % ----------------------------------------------------

    % cloud_topHeight(nn) = alt(is_cloud_idx(end));       % meters - altitude above detected cloud without CLWC > 0
    % cloud_topPressure(nn) = p(is_cloud_idx(end));       % meters - pressure above detected cloud without CLWC > 0
    cloud_topHeight(nn) = alt(is_cloud_idx(end) +1);       % meters - altitude above detected cloud without CLWC > 0
    cloud_topPressure(nn) = p(is_cloud_idx(end) +1);       % meters - pressure above detected cloud without CLWC > 0

    cloud_topHeight_usingVR(nn) = alt( idx_cloudTop_VR);    % meters

    % check if cloud was detected at the lowest pressure level
    if is_cloud_idx(1)==1

        cloud_baseHeight(nn) = alt(is_cloud_idx(1));        % meters - altitude below detected cloud without CLWC > 0
        cloud_basePressure(nn) = p(is_cloud_idx(1));        % meters - pressure below detected cloud without CLWC > 0

    elseif is_cloud_idx(1)>1

        cloud_baseHeight(nn) = alt(is_cloud_idx(1)-1);        % meters - altitude below detected cloud without CLWC > 0
        cloud_basePressure(nn) = p(is_cloud_idx(1)-1);        % meters - pressure below detected cloud without CLWC > 0
    end



    % ----------------------------------------------------





    % store the temperature, pressure and RH profiles
    temp_prof{nn} = T;   % K
    pressure_prof = p;  % mb
    RH_prof{nn} = RH_2keep;  % percent
    rho_v_prof{nn} = rho_v;     % kg/m^3
    watVap_concentration_prof{nn} = waterVapor_concentration_cm3;   %  #/cm^3
    clwc_prof{nn} = clwc;       % g/m^3
    altitude_prof{nn} = alt;          % meters



    % ----------------------------------------------------
    % Keep the spatial distance and temporal differences
    dist_2_VR(nn) = min_spatialDiff;      % meters
    time_2_VR(nn) = min_timeDiff;       % duration object

    dist_2_VR_SL(nn) = sl.min_spatialDiff;      % meters
    time_2_VR_SL(nn) = sl.min_timeDiff;       % duration object



end








%% Save the outputs to a mat file

if strcmp(which_computer,'anbu8374')==true


    folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/'];



elseif strcmp(which_computer,'andrewbuggee')==true


    folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/'];


end



save([folderpath_2save,'ERA5_profiles_closest_to_', num2str(length(ensemble_profiles)),...
    '-VR_profiles_', char(datetime("today")),'.mat'],...
    "total_column_pw", "above_cloud_pw", "cloud_topHeight", "cloud_baseHeight", "cloud_basePressure",...
    "cloud_topPressure", "temp_prof", "RH_prof", "pressure_prof", "altitude_prof", "rho_v_prof",...
    "watVap_concentration_prof", "clwc_prof", "dist_2_VR", "time_2_VR", "sl", "total_column_pw_SL",...
    "cloud_baseHeight_SL", "cloud_topHeight_usingVR", "above_cloud_pw_usingVR", "folder_paths",...
    "dist_2_VR_SL", "time_2_VR_SL")






%% Make Plots



clear variables

if strcmp(whatComputer,'anbu8374')==true

    % ***** Define the EAR5 paper 2 data folder *****
    folderpath_era5 = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/'];

    % ***** Define the ensemble profiles folder *****
    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ***** Define the EAR5 paper 2 data folder *****
    folderpath_era5 = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/'];

    % ***** Define the ensemble profiles folder *****
    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

end

% era5 = load([folderpath_era5,...
%     'ERA5_profiles_closest_to_VR_profiles_30-Jan-2026.mat']);
% era5 = load([folderpath_era5,...
%     'ERA5_profiles_closest_to_73-VR_profiles_31-Jan-2026.mat']);
era5 = load([folderpath_era5,...
    'ERA5_profiles_closest_to_73-VR_profiles_01-Feb-2026.mat']);




% Load the ensemble profiles
load([era5.folder_paths.vocalsRexFolder, era5.folder_paths.VR_profiles]);




% Which profile would you like to plot?
idx_2plt = 3;

fnt_sz = 26;
lgnd_fnt = 20;

era5_clr = 'k';
C130_clr = 61;

% % Create a figure for the vertical profile
% fig1 = figure;
% hold on;
% 
% 
% 
% 
% % Plot temperature
% ax1 = subplot(1,4,1);
% semilogy(era5.temp_prof{idx_2plt} - 273.15, era5.pressure_prof, 'Color', 'k');
% xlabel('Temperature ($C$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% 
% grid on; grid minor
% % flip y axis
% set(gca, 'YDir', 'reverse')
% 
% hold on
% plot(ensemble_profiles{idx_2plt}.temp, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))
% 
% legend('ERA-5', 'Aircraft','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
%     'Color', 'white', 'TextColor', 'k')
% 
% 
% 
% 
% % Plot CLWC
% ax2 = subplot(1,4,2);
% semilogy(era5.clwc_prof{idx_2plt}, era5.pressure_prof, 'Color', 'k');
% xlabel('cloud LWC ($g/m^3$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% 
% grid on; grid minor
% % flip y axis
% set(gca, 'YDir', 'reverse')
% 
% hold on
% plot(ensemble_profiles{idx_2plt}.lwc, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))
% 
% legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
%     'Color', 'white', 'TextColor', 'k')
% 
% 
% % title(['Radiosonde - (', num2str(lat(closestIndex)), ',', num2str(long(closestIndex)),...
% %     ') - ', num2str(month(closestIndex)), '/', num2str(day(closestIndex)), '/2008 - ',...
% %     num2str(hour(closestIndex)), ':', num2str(minute(closestIndex)), ' UTC'], ...
% %     'FontSize', 20, 'Interpreter', 'latex')
% 
% 
% 
% 
% % Plot relative humidity
% ax3 = subplot(1,4,3);
% semilogy(era5.RH_prof{idx_2plt}, era5.pressure_prof, 'Color', 'k');
% xlabel('Relative Humidity ($\%$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% grid on; grid minor
% 
% % flip y axis
% set(gca, 'YDir', 'reverse')
% 
% hold on
% plot(ensemble_profiles{idx_2plt}.relative_humidity, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))
% 
% legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
%     'Color', 'white', 'TextColor', 'k')
% 
% 
% 
% 
% % Plot water vapor number concentration
% ax4 = subplot(1,4,4);
% semilogy(era5.watVap_concentration_prof{idx_2plt}, era5.pressure_prof, 'Color', 'k');
% xlabel('Vapor ($cm^{-3}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% grid on; grid minor
% 
% % flip y axis
% set(gca, 'YDir', 'reverse')
% 
% hold on
% plot(ensemble_profiles{idx_2plt}.Nc_vapor, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))
% 
% legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
%     'Color', 'white', 'TextColor', 'k')
% 
% 
% 
% % Link the y-axes of the specified axes handles
% linkaxes([ax1, ax2, ax3, ax4], 'y');
% 
% % set the yaxis limits to be a region close to the surface
% ylim([700, 1000])
% 
% set(gcf, 'Position', [0,0, 1400, 850])












% Create a figure for the vertical profile
fig1 = figure;



% Plot temperature
% Create axes
ax1 = axes('Parent',fig1,'Position',[0.0885714285714286 0.11 0.165607902735562 0.815]);
hold(ax1,'on');

semilogy(era5.temp_prof{idx_2plt} - 273.15, era5.pressure_prof, 'Color', 'k');
xlabel('Temperature ($C$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)

grid on; grid minor
% flip y axis
set(gca, 'YDir', 'reverse')

hold on
plot(ensemble_profiles{idx_2plt}.temp, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))

legend('ERA-5', 'Aircraft','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')




% Plot CLWC
% Create axes
ax2 = axes('Parent',fig1,'Position',[0.328571428571429 0.11 0.165607902735562 0.815]);
hold(ax2,'on');

semilogy(era5.clwc_prof{idx_2plt}, era5.pressure_prof, 'Color', 'k');
xlabel('cloud LWC ($g/m^3$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)

grid on; grid minor
% flip y axis
set(gca, 'YDir', 'reverse')

hold on
plot(ensemble_profiles{idx_2plt}.lwc, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))

legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')


% title(['Radiosonde - (', num2str(lat(closestIndex)), ',', num2str(long(closestIndex)),...
%     ') - ', num2str(month(closestIndex)), '/', num2str(day(closestIndex)), '/2008 - ',...
%     num2str(hour(closestIndex)), ':', num2str(minute(closestIndex)), ' UTC'], ...
%     'FontSize', 20, 'Interpreter', 'latex')




% Plot relative humidity
% Create axes
ax3 = axes('Parent',fig1,'Position',[0.563900247870048 0.11 0.165607902735562 0.815]);
hold(ax3,'on');

semilogy(era5.RH_prof{idx_2plt}, era5.pressure_prof, 'Color', 'k');
xlabel('Relative Humidity ($\%$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
grid on; grid minor

% flip y axis
set(gca, 'YDir', 'reverse')

hold on
plot(ensemble_profiles{idx_2plt}.relative_humidity, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))

legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')




% Plot water vapor number concentration
% Create axes
ax4 = axes('Parent',fig1,'Position',[0.800017269146644 0.11 0.165607902735562 0.815]);
hold(ax4,'on');

semilogy(era5.watVap_concentration_prof{idx_2plt}, era5.pressure_prof, 'Color', 'k');
xlabel('Vapor ($cm^{-3}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Pressure ($hPa$)', 'Interpreter','latex', 'FontSize', fnt_sz)
grid on; grid minor

% flip y axis
set(gca, 'YDir', 'reverse')

hold on
plot(ensemble_profiles{idx_2plt}.Nc_vapor, ensemble_profiles{idx_2plt}.pres, 'Color', mySavedColors(62, 'fixed'))

legend('ERA-5', 'Aircraft', 'Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')



% Link the y-axes of the specified axes handles
linkaxes([ax1, ax2, ax3, ax4], 'y');

% set the yaxis limits to be a region close to the surface
% ylim([700, 1000])

set(gcf, 'Position', [0,0, 1400, 850])












% % Plot the effective radius
% subplot(2,3,4)
% plot(vert_prof(nn).re, vert_prof(nn).altitude, 'Color', mySavedColors(C130_clr, 'fixed'));
% xlabel('Effective Radius ($\mu m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('Altitude ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
% title('In-situ droplet size - aircraft', 'Interpreter','latex', 'FontSize', fnt_sz)
% grid on; grid minor
% % ylim([0, altitude(end)])
% 
% 
% % Show the position of the radiosonde launch and the airplane when
% % sampling the cloud
% subplot(2,3,5)
% % Plot the location of the radiosonde launch
% geoscatter(lat(closestIndex), long(closestIndex), 100, 'k', 'filled', 'DisplayName', 'ERA-5 Reanalysis');
% hold on;
% 
% % Plot the location of the sampled cloud
% geoscatter(vert_prof(nn).latitude(1), vert_prof(nn).longitude(1), 100, 'b', 'filled', 'DisplayName', 'Sampled Cloud');
% legend('show', 'Location', 'best', 'FontSize', lgnd_fnt, 'Interpreter', 'latex');
% 
% 
% set(gcf, 'Position', [0,0, 2200, 1000])




%% Compare ERA5 cloud top height with VR cloud to height


VR_cloudTopHeight = zeros(length(ensemble_profiles), 1);

for nn = 1:length(ensemble_profiles)

    VR_cloudTopHeight(nn) = max(ensemble_profiles{nn}.altitude);


end

figure; 
subplot(1,2,1)
plot(era5.cloud_topHeight, VR_cloudTopHeight, '.')
hold on
grid on; grid minor

x = linspace(min([era5.cloud_topHeight; VR_cloudTopHeight]),...
    max([era5.cloud_topHeight; VR_cloudTopHeight]), 200);
hold on
plot(x,x, 'k-', "LineWidth", 1)

xlabel('ERA5 Cloud Top Height ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('VOCALS-REx Cloud Top Height ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)


% Compare VR cloud top height with ERA5 cloud top height using VR adjustment
subplot(1,2,2)
plot(era5.cloud_topHeight_usingVR, VR_cloudTopHeight, '.')
hold on
grid on; grid minor

x = linspace(min([era5.cloud_topHeight; VR_cloudTopHeight]),...
    max([era5.cloud_topHeight; VR_cloudTopHeight]), 200);
hold on
plot(x,x, 'k-', "LineWidth", 1)

xlabel('ERA5 CTH w/ VR adjustment ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('VOCALS-REx CTH ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)


set(gcf, 'Position', [0, 0, 1300, 600])





%% Plot a one to one line with the two estimate for above cloud precipitable water

figure; 


plot(era5.above_cloud_pw, era5.above_cloud_pw_usingVR, '.')
grid on; grid minor
hold on

x = linspace(min([era5.above_cloud_pw; era5.above_cloud_pw_usingVR]),...
    max([era5.above_cloud_pw; era5.above_cloud_pw_usingVR]), 200);

plot(x,x, 'k-', "LineWidth", 1)

xlabel('$pw_{ac}$ using ERA5 CTH ($kg/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('$pw_{ac}$ using ERA5 CTH w/ VR adjustment ($kg/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)

set(gcf, 'Position', [0, 0, 1000, 600])



%% Plot a one to one line with pressure level estimated total PW and the single level value

figure; 


plot(era5.total_column_pw, era5.total_column_pw_SL, '.')
grid on; grid minor
hold on

x = linspace(min([era5.total_column_pw; era5.total_column_pw_SL]),...
    max([era5.total_column_pw; era5.total_column_pw_SL]), 200);

plot(x,x, 'k-', "LineWidth", 1)

xlabel('Pressure lvl $pw$ estiamte ($kg/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Single lvl $pw$ value ($kg/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)



%% Plot a one to one line with pressure level estimated cloud base height and the single level value

figure; 


plot(era5.cloud_baseHeight, era5.cloud_baseHeight_SL, '.')
grid on; grid minor
hold on

x = linspace(min([era5.cloud_baseHeight; era5.cloud_baseHeight_SL]),...
    max([era5.cloud_baseHeight; era5.cloud_baseHeight_SL]), 200);

plot(x,x, 'k-', "LineWidth", 1)

xlabel('Cloud Base pressure lvls ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Cloud Base single lvl ($m$)', 'Interpreter','latex', 'FontSize', fnt_sz)

