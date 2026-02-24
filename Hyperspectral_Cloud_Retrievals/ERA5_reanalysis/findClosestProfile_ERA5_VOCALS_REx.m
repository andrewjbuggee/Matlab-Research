%% Find the ERA5 profile closest in space and time to VOCALS-REx in-situ measurement
%
% This function opens the ERA5 data directory for VOCALS-REx-Terra coincident
% observations, finds the file for the same day as the VOCALS-REx scene,
% identifies the time step closest to the VOCALS-REx acquisition time, and
% extracts ERA5 atmospheric profiles at the ERA5 grid points nearest
% to each overlap pixel.
%
% INPUTS:
% (1) vocalsRex - VOCALS-REx data structure (from retrieveVOCALS-REx_data /
%            findOverlap_pixels_VOCALS-REx_Terra_coincident_data)
%            Required fields:
%            .radiance.geo.lat    - Latitude of each VOCALS-REx pixel [degrees]
%            .radiance.geo.long   - Longitude of each VOCALS-REx pixel [degrees]
%            .time                - [hour, minute] UTC acquisition time
%            .time_decimal        - UTC acquisition time in decimal hours
%            .day_of_year         - Day of year of the VOCALS-REx scene
%
% (2) overlap_pixels - Structure with VOCALS-REx-MODIS overlap pixel indices
%            Required fields:
%            .vocalsRex.row            - Row indices of overlap pixels in VOCALS-REx array
%            .vocalsRex.col            - Column indices of overlap pixels in VOCALS-REx array
%            .vocalsRex.linear_idx     - Linear indices of overlap pixels in VOCALS-REx array
%
% (3) print_status_updates - (optional) logical flag to enable/disable
%            printing of status messages to the command window.
%            Default: false
%
% OUTPUTS:
% (1) era5 - Structure containing ERA5 reanalysis data at the relevant
%            spatial and temporal locations. Fields:
%            .temperature         - Temperature [K] (n_unique_era5_pix x n_levels)
%            .specificHumidity    - Specific humidity [kg/kg] (n_unique x n_levels)
%            .geopotential        - Geopotential [m^2/s^2] (n_unique x n_levels)
%            .pressure            - Pressure levels [hPa] (1 x n_levels)
%            .geo.Latitude        - Latitude of each unique ERA5 pixel [degrees]
%            .geo.Longitude       - Longitude of each unique ERA5 pixel [degrees]
%            .metadata            - Metadata structure:
%               .start_year       - Year of the ERA5 file (cell)
%               .start_month      - Month of the ERA5 file (cell)
%               .start_day        - Day of the ERA5 file (cell)
%               .era5_filename    - Filename of the ERA5 file used
%               .era5_utcTime     - UTC datetime of the ERA5 time step used
%               .min_timeDiff     - Time difference between ERA5 and VOCALS-REx [duration]
%               .dist_btwn_era5_and_vocalsRex_m - Distance from ERA5 grid point to
%                                            each VOCALS-REx overlap pixel [m]
%
% (2) overlap_pixels - Updated structure with new field:
%            .era5.linear_idx     - Linear index into ERA5 lat/lon meshgrid
%                                   for each overlap pixel (n_pixels x 1).
%                                   Used by write_ERA5_radiosonde_DAT_with_multiPixels
%                                   to map overlap pixels to unique ERA5 profiles.
%
% NOTES:
%   - The ERA5 data directory is hard-coded for the VOCALS-REx-Terra overlap dataset
%   - Day of year matching is used to find files from the same calendar day
%   - Among files on the same day, the one closest in time to the VOCALS-REx
%     acquisition is selected
%   - The WGS84 ellipsoid is used for all geodetic distance calculations
%   - Multiple VOCALS-REx pixels may share the same ERA5 grid point (ERA5 has
%     coarser spatial resolution than VOCALS-REx)
%   - The era5.temperature rows correspond to unique ERA5 grid points in the
%     order returned by unique(overlap_pixels.era5.linear_idx)
%
% SEE ALSO:
%   write_ERA5_radiosonde_DAT_with_multiPixels
%   vocalsRex_ERA5_priori_statistics
%   retrieveVOCALS-REx_data
%
% By Andrew John Buggee
%%

function era5 = findClosestProfile_ERA5_VOCALS_REx(vocals_rex, print_status_updates, which_computer)

% Default: do not print status updates
if isempty(print_status_updates)
    print_status_updates = false;
end


%% Define the ERA5 data directory recorded on the same days as the VOCALS-REx data



% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % define where the ERA5 data is stored
    era5_dir = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/'];



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define where the ERA5 data is stored
    era5_dir = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/'];


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------
    
    % define where the ERA5 data is stored
    era5_dir = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'ERA5_reanalysis/ERA5_data/VOCALS_REx_overlap/'];

end


%% Find all ERA5 NetCDF files in the directory

% grab both months
era5_files = [dir([era5_dir, 'October_2008/era5_*.nc']); dir([era5_dir, 'November_2008/era5_*.nc'])];

if isempty(era5_files)
    error([newline, 'No ERA5 files found in: ', era5_dir, newline])
end

n_files = length(era5_files);

if print_status_updates
    disp([newline, 'Found ', num2str(n_files), ' ERA5 files in directory.', newline])
end


%% Find the ERA5 file for the same day and closest in time

% Extract VOCALS-REx time information
vocalsRex_time_decimal = vocals_rex.time_utc(round(end/2));    % decimal hours UTC
vocalsRex_day_of_month  = vocals_rex.dateOfFlight.Day;     % integer day of month
vocalsRex_month  = vocals_rex.dateOfFlight.Month;     % integer day of month
vocalsRex_day_of_year  = day(datetime(2008, vocalsRex_month, vocalsRex_day_of_month), 'dayofyear');  % vocals rex day of year

% Parse date and time from each ERA5 filename
% Filename format: era5_vocalsrex_Month_YYYY_MM_dayDD.nc
era5_day_of_year  = zeros(n_files, 1);

for ff = 1:n_files

    fname = era5_files(ff).name;

    % Split by underscore: {'era5', 'YYYY', 'MM', 'DD', 'HHMMUTCz.nc'}
    parts = strsplit(fname, '_');

    yr = str2double(parts{4});
    if strcmp(parts{3}, 'october')==true
        mo = 10;
    elseif strcmp(parts{3}, 'November')==true
        mo = 11;
    end
    dy = str2double(extractBetween(parts{5}, 'day', '.nc'));

    % % parts{5} is like '1500UTC.nc' -- extract just the HHMM portion
    % time_str = parts{5};
    % hr = str2double(time_str(1:2));
    % mn = str2double(time_str(3:4));
    % 
    era5_day_of_year(ff)  = day(datetime(yr, mo, dy), 'dayofyear');
    % era5_time_decimal(ff) = hr + mn/60;

end


% Find ERA5 files on the same day as the VOCALS-REx scene
idx_same_day = find(era5_day_of_year == vocalsRex_day_of_year);

if isempty(idx_same_day)
    error([newline, 'No ERA5 file found for VOCALS-REx day of year: ', num2str(vocalsRex_day_of_year), newline])
end

% define the file path for the nc file from the same day
era5_filepath = [era5_files(idx_same_day).folder, '/', era5_files(idx_same_day).name];

% Among same-day files, select the one closest in time to the VOCALS-REx scene
% first read the time indpendent variable from the era5 data strcuture from
% that day
% Read time steps in the ERA5 file (seconds since 1970 epoch)
time    = double(ncread(era5_filepath, 'valid_time'));
utcTime = datetime(time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% Construct the VOCALS-REx acquisition datetime using the ERA5 file's calendar date
vocalsRex_datetime = datetime(2008, vocalsRex_month, vocalsRex_day_of_month, ...
    floor(vocalsRex_time_decimal), floor(60*(vocalsRex_time_decimal - floor(vocalsRex_time_decimal))),...
    0, 'TimeZone', 'UTC');

% Find the ERA5 time step closest to the VOCALS-REx acquisition time
[min_ERA5_timeDiff, idx_time] = min( abs( utcTime - vocalsRex_datetime ) );
idx_time = idx_time(1);     % take first in case of tie



if print_status_updates
    disp([newline, 'Using ERA5 file: ', era5_files(idx_same_day).name])
    disp([newline, 'ERA5 time step selected: ', char(utcTime(idx_time)), ' UTC'])
    disp([newline, 'Time difference between VOCALS-REx acquisition and ERA5 file: ', ...
        num2str(minutes(min_ERA5_timeDiff), '%.1f'), ' minutes', newline])
end



%% Read the ERA5 NetCDF file

% Parse the date from the best-matching ERA5 filename
fname  = era5_files(idx_same_day).name;
parts = strsplit(fname, '_');

era5_year = str2double(parts{4});
if strcmp(parts{3}, 'october')==true
    era5_month = 10;
elseif strcmp(parts{3}, 'November')==true
    era5_month = 11;
end
era5_day_num = str2double(extractBetween(parts{5}, 'day', '.nc'));





%% Read the ERA5 spatial grid and pressure levels

lat      = ncread(era5_filepath, 'latitude');        % degrees north
long     = ncread(era5_filepath, 'longitude');       % degrees east
pressure = ncread(era5_filepath, 'pressure_level');  % hPa

% Create 2D meshgrid (size: n_lon x n_lat)
[Lat, Lon] = meshgrid(lat, long);

% WGS84 ellipsoid for geodetic distance calculations (distance in meters)
wgs84 = wgs84Ellipsoid('m');


%% For each overlap pixel, find the closest ERA5 grid point


% Get the lat/lon of this VOCALS-REx overlap pixel
% vocalsRex_row = overlap_pixels.vocalsRex.row(pp);
% vocalsRex_col = overlap_pixels.vocalsRex.col(pp);
VR_lat  = double( vocals_rex.latitude(round(end/2)) );
VR_lon  = double( vocals_rex.longitude(round(end/2)) );

% Compute geodetic distance from this VOCALS-REx pixel to all ERA5 grid points [m]
dist = distance(Lat, Lon, VR_lat, VR_lon, wgs84);

% Find the closest ERA5 grid point and its linear index in the meshgrid
[dist_era5_vocalsRex_m, idx_minDist] = min(dist, [], 'all');



%% Read ERA5 atmospheric profile variables
% Read the full 4D arrays (lon x lat x pressure_level x time) at once

T_4D = ncread(era5_filepath, 't');     % Temperature [K]
q_4D = ncread(era5_filepath, 'q');     % Specific humidity [kg/kg]
z_4D = ncread(era5_filepath, 'z');     % Geopotential [m^2/s^2]


%% Extract profiles closest to VR in situ measurement

n_levels        = length(pressure);
n_hrs           = length(utcTime);

% Convert linear meshgrid index to (lon, lat) subscripts
[r, c] = ind2sub(size(Lat), idx_minDist);

% Pre-allocate output arrays (n_unique_era5_pixels x n_pressure_levels)
era5.temperature      = reshape(T_4D(r, c, :, idx_time), 1, []);
era5.specificHumidity = reshape(q_4D(r, c, :, idx_time), 1, []);
era5.geopotential     = reshape(z_4D(r, c, :, idx_time), 1, []);
era5.geo.Latitude     = Lat(r, c);
era5.geo.Longitude    = Lon(r, c);



% Pressure levels are shared across all pixels [hPa]
era5.pressure = pressure(:)';   % row vector: 1 x n_levels


%% Store metadata

era5.metadata.start_year{1}  = era5_year;
era5.metadata.start_month{1} = era5_month;
era5.metadata.start_day{1}   = era5_day_num;

era5.metadata.era5_filename              = era5_files(idx_same_day).name;
era5.metadata.era5_utcTime               = utcTime(idx_time);
era5.metadata.min_timeDiff               = min_ERA5_timeDiff;
era5.metadata.dist_btwn_era5_and_vocalsRex_m  = dist_era5_vocalsRex_m;

if print_status_updates
    disp([newline, 'ERA5 profile closest to VOCALS-REx in-situ extracted.', newline])
    disp(['Distance between ERA5 grid point and VOCALS-REx measurement: ', ...
        num2str(dist_era5_vocalsRex_m / 1000, '%.2f'), ' km', newline])
end


end
