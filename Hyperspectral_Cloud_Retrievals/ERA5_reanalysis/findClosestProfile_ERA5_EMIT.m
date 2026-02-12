%% Find the ERA5 profile closest in space and time to the EMIT scene
%
% This function opens the ERA5 data directory for EMIT-Terra coincident
% observations, finds the file for the same day as the EMIT scene,
% identifies the time step closest to the EMIT acquisition time, and
% extracts ERA5 atmospheric profiles at the ERA5 grid points nearest
% to each overlap pixel.
%
% INPUTS:
% (1) emit - EMIT data structure (from retrieveEMIT_data /
%            findOverlap_pixels_EMIT_Terra_coincident_data)
%            Required fields:
%            .radiance.geo.lat    - Latitude of each EMIT pixel [degrees]
%            .radiance.geo.long   - Longitude of each EMIT pixel [degrees]
%            .time                - [hour, minute] UTC acquisition time
%            .time_decimal        - UTC acquisition time in decimal hours
%            .day_of_year         - Day of year of the EMIT scene
%
% (2) overlap_pixels - Structure with EMIT-MODIS overlap pixel indices
%            Required fields:
%            .emit.row            - Row indices of overlap pixels in EMIT array
%            .emit.col            - Column indices of overlap pixels in EMIT array
%            .emit.linear_idx     - Linear indices of overlap pixels in EMIT array
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
%               .min_timeDiff     - Time difference between ERA5 and EMIT [duration]
%               .dist_btwn_era5_and_emit_m - Distance from ERA5 grid point to
%                                            each EMIT overlap pixel [m]
%
% (2) overlap_pixels - Updated structure with new field:
%            .era5.linear_idx     - Linear index into ERA5 lat/lon meshgrid
%                                   for each overlap pixel (n_pixels x 1).
%                                   Used by write_ERA5_radiosonde_DAT_with_multiPixels
%                                   to map overlap pixels to unique ERA5 profiles.
%
% NOTES:
%   - The ERA5 data directory is hard-coded for the EMIT-Terra overlap dataset
%   - Day of year matching is used to find files from the same calendar day
%   - Among files on the same day, the one closest in time to the EMIT
%     acquisition is selected
%   - The WGS84 ellipsoid is used for all geodetic distance calculations
%   - Multiple EMIT pixels may share the same ERA5 grid point (ERA5 has
%     coarser spatial resolution than EMIT)
%   - The era5.temperature rows correspond to unique ERA5 grid points in the
%     order returned by unique(overlap_pixels.era5.linear_idx)
%
% SEE ALSO:
%   write_ERA5_radiosonde_DAT_with_multiPixels
%   vocalsRex_ERA5_priori_statistics
%   retrieveEMIT_data
%
% By Andrew John Buggee
%%

function [era5, overlap_pixels] = findClosestProfile_ERA5_EMIT(emit, overlap_pixels)


%% Define the ERA5 data directory for EMIT-Terra overlap

era5_dir = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/ERA5_reanalysis/ERA5_data/EMIT_Terra_overlap/';


%% Find all ERA5 NetCDF files in the directory

era5_files = dir([era5_dir, 'era5_*.nc']);

if isempty(era5_files)
    error([newline, 'No ERA5 files found in: ', era5_dir, newline])
end

n_files = length(era5_files);

disp([newline, 'Found ', num2str(n_files), ' ERA5 files in directory.', newline])


%% Find the ERA5 file for the same day and closest in time

% Extract EMIT time information
emit_time_decimal = emit.time_decimal;    % decimal hours UTC
emit_day_of_year  = emit.day_of_year;     % integer day of year

% Parse date and time from each ERA5 filename
% Filename format: era5_YYYY_MM_DD_HHMMUTCz.nc
era5_day_of_year  = zeros(n_files, 1);
era5_time_decimal = zeros(n_files, 1);

for ff = 1:n_files

    fname = era5_files(ff).name;

    % Split by underscore: {'era5', 'YYYY', 'MM', 'DD', 'HHMMUTCz.nc'}
    parts = strsplit(fname, '_');

    yr = str2double(parts{2});
    mo = str2double(parts{3});
    dy = str2double(parts{4});

    % parts{5} is like '1500UTC.nc' -- extract just the HHMM portion
    time_str = parts{5};
    hr = str2double(time_str(1:2));
    mn = str2double(time_str(3:4));

    era5_day_of_year(ff)  = day(datetime(yr, mo, dy), 'dayofyear');
    era5_time_decimal(ff) = hr + mn/60;

end


% Find ERA5 files on the same day as the EMIT scene
idx_same_day = find(era5_day_of_year == emit_day_of_year);

if isempty(idx_same_day)
    error([newline, 'No ERA5 file found for EMIT day of year: ', num2str(emit_day_of_year), newline])
end


% Among same-day files, select the one closest in time to the EMIT scene
[min_timeDiff_hrs, idx_min_time] = min( abs( era5_time_decimal(idx_same_day) - emit_time_decimal ) );
idx_best = idx_same_day(idx_min_time);

era5_filepath = [era5_dir, era5_files(idx_best).name];

disp([newline, 'Using ERA5 file: ', era5_files(idx_best).name])
disp(['Time difference between EMIT acquisition and ERA5 file: ', ...
    num2str(min_timeDiff_hrs * 60, '%.1f'), ' minutes', newline])


%% Read the ERA5 NetCDF file

% Parse the date from the best-matching ERA5 filename
fname  = era5_files(idx_best).name;
parts  = strsplit(fname, '_');
era5_year    = str2double(parts{2});
era5_month   = str2double(parts{3});
era5_day_num = str2double(parts{4});

% Read time steps in the ERA5 file (seconds since 1970 epoch)
time    = double(ncread(era5_filepath, 'valid_time'));
utcTime = datetime(time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% Construct the EMIT acquisition datetime using the ERA5 file's calendar date
emit_datetime = datetime(era5_year, era5_month, era5_day_num, ...
    emit.time(1), emit.time(2), 0, 'TimeZone', 'UTC');

% Find the ERA5 time step closest to the EMIT acquisition time
[min_ERA5_timeDiff, idx_time] = min( abs( utcTime - emit_datetime ) );
idx_time = idx_time(1);     % take first in case of tie

disp(['ERA5 time step selected: ', char(utcTime(idx_time)), ' UTC'])
disp(['Time difference to EMIT acquisition: ', char(min_ERA5_timeDiff), newline])


%% Read the ERA5 spatial grid and pressure levels

lat      = ncread(era5_filepath, 'latitude');        % degrees north
long     = ncread(era5_filepath, 'longitude');       % degrees east
pressure = ncread(era5_filepath, 'pressure_level');  % hPa

% Create 2D meshgrid (size: n_lon x n_lat)
[Lat, Lon] = meshgrid(lat, long);

% WGS84 ellipsoid for geodetic distance calculations (distance in meters)
wgs84 = wgs84Ellipsoid('m');


%% For each overlap pixel, find the closest ERA5 grid point

n_overlap = length(overlap_pixels.emit.linear_idx);

overlap_pixels.era5.linear_idx = zeros(n_overlap, 1);
dist_era5_emit_m                = zeros(n_overlap, 1);

for pp = 1:n_overlap

    % Get the lat/lon of this EMIT overlap pixel
    emit_row = overlap_pixels.emit.row(pp);
    emit_col = overlap_pixels.emit.col(pp);
    pix_lat  = double( emit.radiance.geo.lat(emit_row, emit_col) );
    pix_lon  = double( emit.radiance.geo.long(emit_row, emit_col) );

    % Compute geodetic distance from this EMIT pixel to all ERA5 grid points [m]
    dist = distance(Lat, Lon, pix_lat, pix_lon, wgs84);

    % Find the closest ERA5 grid point and its linear index in the meshgrid
    [dist_era5_emit_m(pp), idx_minDist] = min(dist, [], 'all');

    % Store as linear index into the ERA5 meshgrid
    overlap_pixels.era5.linear_idx(pp) = idx_minDist;

end


%% Read ERA5 atmospheric profile variables
% Read the full 4D arrays (lon x lat x pressure_level x time) at once

T_4D = ncread(era5_filepath, 't');     % Temperature [K]
q_4D = ncread(era5_filepath, 'q');     % Specific humidity [kg/kg]
z_4D = ncread(era5_filepath, 'z');     % Geopotential [m^2/s^2]


%% Extract profiles for each unique ERA5 spatial pixel

unique_era5_pix = unique(overlap_pixels.era5.linear_idx);
n_unique        = length(unique_era5_pix);
n_levels        = length(pressure);

% Pre-allocate output arrays (n_unique_era5_pixels x n_pressure_levels)
era5.temperature      = zeros(n_unique, n_levels);
era5.specificHumidity = zeros(n_unique, n_levels);
era5.geopotential     = zeros(n_unique, n_levels);
era5.geo.Latitude     = zeros(n_unique, 1);
era5.geo.Longitude    = zeros(n_unique, 1);

for kk = 1:n_unique

    % Convert linear meshgrid index to (lon, lat) subscripts
    [r, c] = ind2sub(size(Lat), unique_era5_pix(kk));

    era5.temperature(kk, :)      = reshape(T_4D(r, c, :, idx_time), 1, []);
    era5.specificHumidity(kk, :) = reshape(q_4D(r, c, :, idx_time), 1, []);
    era5.geopotential(kk, :)     = reshape(z_4D(r, c, :, idx_time), 1, []);

    era5.geo.Latitude(kk)  = Lat(r, c);
    era5.geo.Longitude(kk) = Lon(r, c);

end

% Pressure levels are shared across all pixels [hPa]
era5.pressure = pressure(:)';   % row vector: 1 x n_levels


%% Store metadata

era5.metadata.start_year{1}  = era5_year;
era5.metadata.start_month{1} = era5_month;
era5.metadata.start_day{1}   = era5_day_num;

era5.metadata.era5_filename              = era5_files(idx_best).name;
era5.metadata.era5_utcTime               = utcTime(idx_time);
era5.metadata.min_timeDiff               = min_ERA5_timeDiff;
era5.metadata.dist_btwn_era5_and_emit_m  = dist_era5_emit_m;

disp([newline, 'ERA5 profiles extracted for ', num2str(n_unique), ' unique spatial pixels.'])
disp(['Average distance between ERA5 grid points and EMIT pixels: ', ...
    num2str(mean(dist_era5_emit_m) / 1000, '%.2f'), ' km', newline])


end
