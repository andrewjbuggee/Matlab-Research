%% Step through each VOCAL-REx vertical profile and find the closest ERA5 data set in space and time


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    % ***** Define the ERA5 data directory *****
    folderpath_era5 = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/'];


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

cloud_topHeight = zeros(length(ensemble_profiles), 1);
cloud_baseHeight = zeros(length(ensemble_profiles), 1);

cloud_topPressure = zeros(length(ensemble_profiles), 1);
cloud_basePressure = zeros(length(ensemble_profiles), 1);

temp_prof = cell(length(ensemble_profiles), 1);
watVap_prof = cell(length(ensemble_profiles), 1);
pressure_prof = cell(length(ensemble_profiles), 1);
altitude_prof = cell(length(ensemble_profiles), 1);
rho_v_prof = cell(length(ensemble_profiles), 1);
watVap_concentration_prof = cell(length(ensemble_profiles), 1);
clwc_prof = cell(length(ensemble_profiles), 1);


dist_2_VR = zeros(length(ensemble_profiles), 1);
time_2_VR = cell(length(ensemble_profiles), 1);


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
    p = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'pressure_level');            % hPa

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

    % -----------------------------------------------------------------------------------
    % Compute cloud liquid water content from specific cloud LWC
    clwc_2keep = (clwc_2keep .* p) ./...
        (computeMoistAirGasConstant(q_2keep) .* T);        % kg of liquid water droplets/ m^3 of air
    
    % convert to grams per m^3
    clwc_2keep = clwc_2keep * 1000; % Convert to grams per m^3
    % -----------------------------------------------------------------------------------


    % -----------------------------------------------------------------------------------
    % compute the above cloud total precipitable water
    % Compute the mass of water vapor per unit volume at all altitudes
    % convert pressure from hPa to Pa
    use_virtual_temp = true;   % more accurate estimate, but the improvement is on the order of 1%
    % rho_v = [kg/m^3]
    % waterVapor_concentration_cm3 = [#/cm^3]
    [rho_v, waterVapor_concentration_cm3] = specificHumidity2waterVaporDensity(q_2keep, T, p.*100, use_virtual_temp);


    % % compute the geopotential height for each layer;
    % ** my estiamtes vary from the ERA 5 estimates. Use theirs **
    % % we are only using data over ocean, so set Z_surface to 0
    % Z_sfc = 0;         % meters
    % Z = computeGeopotentialHeight(T, q_2keep, p, Z_sfc);   % meters

    % read in the geopotential and convert it to geopotential height 
    Z = ncread([folderpath_era5, era5_month_folderpath, era5_filename], 'z');            % m^2/s^2
    Z = reshape( Z(r, c, :, idx_time), [], 1) ./ 9.80665;       % meters

    % From ERA5's website, convert geopotential height to geometric height
    % (https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height)
    % This is likely to be very similar since geopotential height is close
    % to the geometric height in the lower atmosphere
    con = physical_constants;
    alt = con.R_earth * Z ./(con.R_earth - Z);       % meters - geometric height

    % check total column amount
    total_column_pw(nn) = trapz(alt, rho_v);            % kg / m^2

    % ----------------------------------------------------
    % Find the altitude level where CLWC is greater than 0
    is_cloud_idx = find(clwc_2keep > 0);

    if isempty(is_cloud_idx) == true

        error([newline, 'No cloud found for profile ', num2str(nn), newline])

    end


    % compute the above cloud precipitable water amount
    above_cloud_pw(nn) = trapz( alt( (is_cloud_idx(end) +1) : end ), rho_v( (is_cloud_idx(end) +1) : end ));            % kg / m^2

    % -----------------------------------------------------------------------------------


    % ----------------------------------------------------
    % Define the cloud top and base height and pressure

    cloud_topHeight = alt(is_cloud_idx(end)+1);       % meters - altitude above detected cloud without CLWC > 0
    cloud_topPressure = p(is_cloud_idx(end)+1);       % meters - pressure above detected cloud without CLWC > 0
    
    % check if cloud was detected at the lowest pressure level
    if is_cloud_idx(1)==1

        cloud_baseHeight = alt(is_cloud_idx(1));        % meters - altitude below detected cloud without CLWC > 0
        cloud_basePressure = p(is_cloud_idx(1));        % meters - pressure below detected cloud without CLWC > 0

    elseif is_cloud_idx(1)>1
        
        cloud_baseHeight = alt(is_cloud_idx(1)-1);        % meters - altitude below detected cloud without CLWC > 0
        cloud_basePressure = p(is_cloud_idx(1)-1);        % meters - pressure below detected cloud without CLWC > 0
    end

    
    
    % ----------------------------------------------------

    % store the temperature, pressure and RH profiles 
    temp_prof{nn} = T;   % K
    pressure_prof{nn} = p;  % mb
    watVap_prof{nn} = RH_2keep;  % percent
    rho_v_prof{nn} = rho_v;     % kg/m^3 
    watVap_concentration_prof{nn} = waterVapor_concentration_cm3;   %  #/cm^3
    clwc_prof{nn} = clwc_2keep;       % g/m^3
    altitude_prof{nn} = alt;          % meters


    
    % ----------------------------------------------------
    % Keep the spatial distance and temporal differences
    dist_2_VR(nn) = min_spatialDiff;      % meters
    time_2_VR{nn} = min_timeDiff;       % duration object

end








%% Save the outputs to a mat file

if strcmp(which_computer,'anbu8374')==true


    folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/'];



elseif strcmp(which_computer,'andrewbuggee')==true


    folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/'];


end



save([folderpath_2save,'ERA5_profiles_closest_to_VR_profiles', char(datetime("today")),'.mat'],...
    "total_column_pw", "above_cloud_pw", "cloud_topHeight", "cloud_baseHeight", "cloud_basePressure",...
    "cloud_topPressure", "temp_prof", "watVap_prof", "pressure_prof", "altitude_prof", "rho_v_prof",...
    "watVap_concentration_prof", "clwc_prof", "dist_2_VR", "time_2_VR")





%%








