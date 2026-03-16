%% Find the ERA5 profile closest in space and time to an ORACLES in-situ profile
%
% Analogous to findClosestProfile_ERA5_VOCALS_REx.m, adapted for the ORACLES
% P-3 dataset.  Opens the ORACLES ERA5 data directory, finds the ERA5 file
% for the same calendar day as the ORACLES profile, identifies the hourly
% time step closest to the mid-profile UTC time, and extracts the ERA5
% atmospheric profile at the ERA5 grid point nearest to the profile location.
%
% INPUTS:
% -------
%   oracles_profile        - One element of the ORACLES ensemble_profiles
%                            cell array (output of find_verticalProfiles_ORACLES).
%                            Required fields:
%                              .latitude     - [degrees N] (1 x N_time)
%                              .longitude    - [degrees E] (1 x N_time)
%                              .time         - [s since midnight UTC] (1 x N_time)
%                              .dateOfFlight - datetime scalar
%
%   print_status_updates   - logical; if true, print progress messages
%
%   which_computer         - string identifying the machine
%                            ('anbu8374', 'andrewbuggee', or 'curc')
%
% OUTPUTS:
% --------
%   era5 - Structure with the following fields (same format as
%          findClosestProfile_ERA5_VOCALS_REx output):
%            .temperature         - [K]       (1 x n_levels)
%            .specificHumidity    - [kg/kg]   (1 x n_levels)
%            .geopotential        - [m^2/s^2] (1 x n_levels)
%            .pressure            - [hPa]     (1 x n_levels)
%            .geo.Latitude        - ERA5 grid point latitude [degrees]
%            .geo.Longitude       - ERA5 grid point longitude [degrees]
%            .metadata            - struct with date/file/distance info
%
% NOTES:
%   - ERA5 files follow the naming convention: era5_YYYY_MM_DD_HHMMUTC.nc
%   - The WGS84 ellipsoid is used for geodetic distance calculations
%
% SEE ALSO:
%   findClosestProfile_ERA5_VOCALS_REx, write_ERA5_radiosonde_DAT_with_multiPixels
%
% By Andrew John Buggee

%%

function era5 = findClosestProfile_ERA5_ORACLES(oracles_profile, print_status_updates, which_computer)


% Default: do not print status updates
if isempty(print_status_updates)
    print_status_updates = false;
end


%% Define the ERA5 data directory for ORACLES

if strcmp(which_computer, 'anbu8374')

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    era5_dir = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
        'ERA5_reanalysis/ERA5_data/ORACLES/'];

elseif strcmp(which_computer, 'andrewbuggee')

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    era5_dir = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
        'ERA5_reanalysis/ERA5_data/ORACLES/'];

elseif strcmp(which_computer, 'curc')

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    era5_dir = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
        'ERA5_reanalysis/ERA5_data/ORACLES/'];

end


%% Find all ERA5 NetCDF files in the ORACLES directory

era5_files = dir([era5_dir, 'era5_*.nc']);

if isempty(era5_files)
    error([newline, 'No ERA5 files found in: ', era5_dir, newline])
end

n_files = length(era5_files);

if print_status_updates
    disp([newline, 'Found ', num2str(n_files), ' ERA5 files in directory.', newline])
end


%% Parse date and time from each ERA5 filename
%
% ORACLES ERA5 filename format: era5_YYYY_MM_DD_HHMMUTC.nc
% Example: era5_2016_09_06_0700UTC.nc
% Split by '_' gives: {'era5', 'YYYY', 'MM', 'DD', 'HHMMUTC.nc'}

era5_yr           = zeros(n_files, 1);
era5_mo           = zeros(n_files, 1);
era5_dy           = zeros(n_files, 1);
era5_datetimes    = NaT(n_files, 1, 'TimeZone', 'UTC');

for ff = 1:n_files

    parts = strsplit(era5_files(ff).name, '_');

    yr = str2double(parts{2});
    mo = str2double(parts{3});
    dy = str2double(parts{4});

    % parts{5} is like '0700UTC.nc' -- extract HHMM portion
    time_str = parts{5};
    hr = str2double(time_str(1:2));
    mn = str2double(time_str(3:4));

    era5_yr(ff) = yr;
    era5_mo(ff) = mo;
    era5_dy(ff) = dy;
    era5_datetimes(ff) = datetime(yr, mo, dy, hr, mn, 0, 'TimeZone', 'UTC');

end


%% Identify the ORACLES profile location and time

% Use the mid-profile point for the spatial and temporal match
mid_idx = round(length(oracles_profile.latitude) / 2);

oracles_lat = double(oracles_profile.latitude(mid_idx));
oracles_lon = double(oracles_profile.longitude(mid_idx));

% Convert seconds since midnight UTC to a full datetime
oracles_date  = oracles_profile.dateOfFlight;
oracles_time_s = oracles_profile.time(mid_idx);       % s since midnight UTC

oracles_datetime = datetime(oracles_date.Year, oracles_date.Month, oracles_date.Day, ...
    0, 0, oracles_time_s, 'TimeZone', 'UTC');


%% Find ERA5 files on the same calendar day as the ORACLES profile

idx_same_day = find( ...
    era5_yr == oracles_date.Year  & ...
    era5_mo == oracles_date.Month & ...
    era5_dy == oracles_date.Day );

if isempty(idx_same_day)
    error([newline, 'No ERA5 file found for ORACLES date: ', ...
        char(oracles_date, 'yyyy-MM-dd'), newline])
end


%% Among same-day files, select the one closest in time

[min_ERA5_timeDiff, idx_best_rel] = min( abs(era5_datetimes(idx_same_day) - oracles_datetime) );
idx_best = idx_same_day(idx_best_rel(1));

era5_filepath = [era5_files(idx_best).folder, '/', era5_files(idx_best).name];

if print_status_updates
    disp([newline, 'Using ERA5 file: ', era5_files(idx_best).name])
    disp([newline, 'ERA5 time step selected: ', char(era5_datetimes(idx_best)), ' UTC'])
    disp([newline, 'Time difference between ORACLES profile and ERA5 file: ', ...
        num2str(minutes(min_ERA5_timeDiff), '%.1f'), ' minutes', newline])
end


%% Read the ERA5 spatial grid and pressure levels

lat      = ncread(era5_filepath, 'latitude');        % degrees north
lon      = ncread(era5_filepath, 'longitude');       % degrees east
pressure = ncread(era5_filepath, 'pressure_level');  % hPa

% Create 2D meshgrid (size: n_lon x n_lat)
[Lat, Lon] = meshgrid(lat, lon);

% WGS84 ellipsoid for geodetic distance calculations (meters)
wgs84 = wgs84Ellipsoid('m');


%% Find the ERA5 grid point closest to the ORACLES profile location

dist = distance(Lat, Lon, oracles_lat, oracles_lon, wgs84);
[dist_era5_oracles_m, idx_minDist] = min(dist, [], 'all');

% Convert linear meshgrid index to (lon, lat) subscripts
[r, c] = ind2sub(size(Lat), idx_minDist);


%% Read ERA5 atmospheric profile variables

% Read the ERA5 file time coordinate (seconds since 1970 epoch)
era5_time_raw  = double(ncread(era5_filepath, 'valid_time'));
era5_utcTimes  = datetime(era5_time_raw, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% Find the time index closest to the ORACLES mid-profile time
[min_timeDiff_step, idx_time] = min( abs(era5_utcTimes - oracles_datetime) );
idx_time = idx_time(1);

T_4D = ncread(era5_filepath, 't');     % Temperature [K]       (n_lon x n_lat x n_lev x n_time)
q_4D = ncread(era5_filepath, 'q');     % Specific humidity [kg/kg]
z_4D = ncread(era5_filepath, 'z');     % Geopotential [m^2/s^2]


%% Extract atmospheric profile at the closest ERA5 grid point

era5.temperature      = reshape(T_4D(r, c, :, idx_time), 1, []);    % 1 x n_levels
era5.specificHumidity = reshape(q_4D(r, c, :, idx_time), 1, []);    % 1 x n_levels
era5.geopotential     = reshape(z_4D(r, c, :, idx_time), 1, []);    % 1 x n_levels
era5.pressure         = pressure(:)';                                % 1 x n_levels  [hPa]

era5.geo.Latitude  = Lat(r, c);
era5.geo.Longitude = Lon(r, c);


%% Store metadata

era5.metadata.start_year{1}  = era5_yr(idx_best);
era5.metadata.start_month{1} = era5_mo(idx_best);
era5.metadata.start_day{1}   = era5_dy(idx_best);

era5.metadata.era5_filename               = era5_files(idx_best).name;
era5.metadata.era5_utcTime                = era5_datetimes(idx_best);
era5.metadata.min_timeDiff                = min_ERA5_timeDiff;
era5.metadata.dist_btwn_era5_and_oracles_m = dist_era5_oracles_m;


if print_status_updates
    disp([newline, 'ERA5 profile closest to ORACLES in-situ profile extracted.'])
    disp(['Distance between ERA5 grid point and ORACLES profile: ', ...
        num2str(dist_era5_oracles_m / 1000, '%.2f'), ' km', newline])
end


end
