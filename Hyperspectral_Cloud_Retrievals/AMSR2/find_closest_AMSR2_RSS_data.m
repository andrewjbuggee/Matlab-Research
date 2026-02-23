%% Find the closest AMSR2 RSS L3 data in time and space for a given EMIT pixel
%
%   This function searches the RSS AMSR2 L3 daily NetCDF files to find the
%   observation closest in both time and space to a specific EMIT pixel.
%   It reads the appropriate daily file based on the EMIT acquisition date,
%   finds the nearest grid cell, selects the pass closest in time, and
%   returns a structure containing only the data at that single location.
%
% INPUTS:
%   emit        - EMIT data structure (from retrieveEMIT_data.m) containing:
%                   .radiance.geo.lat  - latitude array
%                   .radiance.geo.long - longitude array
%                   .time_decimal      - UTC time in decimal hours
%
%   overlap_pixels - Structure with overlap pixel indices containing:
%                       .emit.row    - row indices into EMIT arrays
%                       .emit.col    - column indices into EMIT arrays
%
%   pixel_num   - Index into overlap_pixels arrays for the pixel of interest
%
%   amsr2_folder - (Optional) Path to folder containing RSS AMSR2 L3 .nc files
%                  Default: '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/
%                            Hyperspectral_Cloud_Retrievals/AMSR2/RSS_AMSR2_LWP/'
%
%   L1B_fileName - EMIT L1B filename string used to extract the acquisition date
%                  (e.g., 'EMIT_L1B_RAD_20240113T183456_...')
%
% OUTPUTS:
%   amsr2_closest - Structure containing the closest AMSR2 data point with fields:
%                     .lat                - latitude of closest grid cell (deg N)
%                     .lon                - longitude of closest grid cell (deg E, 0-360)
%                     .time               - observation time (fractional hours UTC)
%                     .pass               - pass number (1=ascending ~1:30PM, 2=descending ~1:30AM)
%                     .SST                - sea surface temperature (deg C)
%                     .wind_speed_LF      - low frequency wind speed (m/s)
%                     .wind_speed_MF      - medium frequency wind speed (m/s)
%                     .wind_speed_AW      - all weather wind speed (m/s)
%                     .water_vapor        - columnar water vapor (kg/m^2)
%                     .cloud_liquid_water - columnar cloud liquid water (kg/m^2)
%                     .rain_rate          - surface rain rate (mm/hr)
%                     .land_mask          - land mask flag
%                     .sea_ice_mask       - sea ice mask flag
%                     .coast_mask         - coast mask flag
%                     .noobs_mask         - no observation mask flag
%                     .distance_km        - distance to EMIT pixel (km)
%                     .time_diff_hrs      - time difference from EMIT (hours)
%                     .nc_filename        - name of the NetCDF file used
%                     .nc_global_attrs    - global attributes from the NetCDF file
%
%
% by Andrew J. Buggee
%%

function amsr2_closest = find_closest_AMSR2_RSS_data(emit, overlap_pixels, pixel_num, L1B_fileName, amsr2_folder)


%% Set default AMSR2 folder if not provided

if nargin < 5

    % Determine which computer you're using
    which_computer = whatComputer();

    if strcmp(which_computer, 'anbu8374') == true

        amsr2_folder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
            'Hyperspectral_Cloud_Retrievals/AMSR2/RSS_AMSR2_LWP/'];

    elseif strcmp(which_computer, 'andrewbuggee') == true

        amsr2_folder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
            'Hyperspectral_Cloud_Retrievals/AMSR2/RSS_AMSR2_LWP/'];

    end

end


%% Extract the EMIT pixel location and time

% Check if EMIT data has been filtered by remove_unwanted_emit_data()
% If so, lat/lon are 1D vectors indexed directly by pixel_num
% If not, lat/lon are 2D arrays indexed by (row, col)
if isvector(emit.radiance.geo.lat)

    % Data has been filtered - 1D indexing
    emit_lat = emit.radiance.geo.lat(pixel_num);
    emit_lon = emit.radiance.geo.long(pixel_num);

    % Use per-pixel UTC time if available (more accurate than scene-level)
    if isfield(emit, 'obs') && isfield(emit.obs, 'utc_time') && ...
            length(emit.obs.utc_time) >= pixel_num
        emit_time_utc = emit.obs.utc_time(pixel_num);
    else
        emit_time_utc = emit.time_decimal;
    end

else

    % Data has NOT been filtered - 2D indexing using overlap_pixels
    emit_row = overlap_pixels.emit.row(pixel_num);
    emit_col = overlap_pixels.emit.col(pixel_num);

    emit_lat = emit.radiance.geo.lat(emit_row, emit_col);
    emit_lon = emit.radiance.geo.long(emit_row, emit_col);

    % Use per-pixel UTC time if available
    if isfield(emit, 'obs') && isfield(emit.obs, 'utc_time')
        emit_time_utc = emit.obs.utc_time(emit_row, emit_col);
    else
        emit_time_utc = emit.time_decimal;
    end

end


%% Extract the EMIT acquisition date from the L1B filename

% EMIT L1B filename format: EMIT_L1B_RAD_YYYYMMDDTHHmmss_...
emit_year  = str2double(L1B_fileName(18:21));
emit_month = str2double(L1B_fileName(22:23));
emit_day   = str2double(L1B_fileName(24:25));


%% Find the matching AMSR2 RSS L3 daily file

% RSS AMSR2 filename format: RSS_AMSR2_ocean_L3_daily_YYYY-MM-DD_v08.2.nc
date_string = sprintf('%04d-%02d-%02d', emit_year, emit_month, emit_day);
nc_filename = ['RSS_AMSR2_ocean_L3_daily_', date_string, '_v08.2.nc'];
nc_filepath = fullfile(amsr2_folder, nc_filename);

% Check if the file exists
if ~isfile(nc_filepath)
    error([newline, 'AMSR2 RSS L3 file not found: ', nc_filepath, newline,...
        'Make sure the file for ', date_string, ' is in the folder.', newline])
end


%% Read the AMSR2 grid coordinates

amsr2_lat = ncread(nc_filepath, 'lat');       % (720 x 1) - degrees north
amsr2_lon = ncread(nc_filepath, 'lon');       % (1440 x 1) - degrees east (0-360)


%% Convert EMIT longitude to 0-360 range to match AMSR2

% EMIT uses -180 to 180, AMSR2 RSS uses 0 to 360
if emit_lon < 0
    emit_lon_360 = emit_lon + 360;
else
    emit_lon_360 = emit_lon;
end


%% Find the closest grid cell in space

% The AMSR2 RSS L3 grid is regular 0.25 deg, so we can find the closest
% indices directly
[~, lon_idx] = min(abs(amsr2_lon - emit_lon_360));
[~, lat_idx] = min(abs(amsr2_lat - emit_lat));

% Compute the actual distance using WGS84 ellipsoid
wgs84 = wgs84Ellipsoid('kilometer');
distance_km = distance(emit_lat, emit_lon, amsr2_lat(lat_idx), ...
    amsr2_lon(lon_idx) - 360*(amsr2_lon(lon_idx) > 180), wgs84);


%% Read time data for both passes at the closest grid cell

% NetCDF dimensions are (lon, lat, pass) in MATLAB due to column-major reading
% In the file: (pass, lat, lon) -> MATLAB reads as (lon, lat, pass)
amsr2_time = ncread(nc_filepath, 'time');     % (1440 x 720 x 2) in MATLAB

fill_value = -999.0;

% Extract time for both passes at the closest spatial location
time_pass1 = amsr2_time(lon_idx, lat_idx, 1);    % ascending (~1:30 PM local)
time_pass2 = amsr2_time(lon_idx, lat_idx, 2);    % descending (~1:30 AM local)


%% Find the pass closest in time to the EMIT observation

% Compute time differences, handling fill values
time_diff = nan(2, 1);

if time_pass1 ~= fill_value
    time_diff(1) = abs(emit_time_utc - time_pass1);
end

if time_pass2 ~= fill_value
    time_diff(2) = abs(emit_time_utc - time_pass2);
end

% Check that at least one pass has valid data
if all(isnan(time_diff))
    error([newline, 'No valid AMSR2 observations found at the closest grid cell ',...
        '(lat=', num2str(amsr2_lat(lat_idx)), ', lon=', num2str(amsr2_lon(lon_idx)), ...
        ') for date ', date_string, '.', newline,...
        'Both ascending and descending passes have no observations at this location.', newline])
end

% Select the pass with the smallest time difference
[min_time_diff, best_pass] = min(time_diff);


%% Read all variables at the closest location and best pass

% List of physical variables to read
phys_vars = {'SST', 'wind_speed_LF', 'wind_speed_MF', 'wind_speed_AW',...
    'water_vapor', 'cloud_liquid_water', 'rain_rate'};

% List of mask/flag variables to read
mask_vars = {'land_mask', 'sea_ice_mask', 'coast_mask', 'noobs_mask'};


%% Build the output structure

amsr2_closest = struct();

% Spatial and temporal info
amsr2_closest.lat = double(amsr2_lat(lat_idx));
amsr2_closest.lon = double(amsr2_lon(lon_idx));
amsr2_closest.time = double(amsr2_time(lon_idx, lat_idx, best_pass));
amsr2_closest.pass = best_pass;         % 1 = ascending, 2 = descending

% Read physical variables
for vv = 1:length(phys_vars)

    data = ncread(nc_filepath, phys_vars{vv});
    val = double(data(lon_idx, lat_idx, best_pass));

    % Replace fill values with NaN
    if val == fill_value
        val = NaN;
    end

    amsr2_closest.(phys_vars{vv}) = val;

end

% Read mask variables
for vv = 1:length(mask_vars)

    data = ncread(nc_filepath, mask_vars{vv});
    amsr2_closest.(mask_vars{vv}) = double(data(lon_idx, lat_idx, best_pass));

end

% Metadata
amsr2_closest.distance_km = distance_km;
amsr2_closest.time_diff_hrs = min_time_diff;
amsr2_closest.nc_filename = nc_filename;

% Read and store global attributes
nc_info = ncinfo(nc_filepath);
amsr2_closest.nc_global_attrs = struct();
for aa = 1:length(nc_info.Attributes)
    % Sanitize attribute name for use as a struct field
    attr_name = matlab.lang.makeValidName(nc_info.Attributes(aa).Name);
    amsr2_closest.nc_global_attrs.(attr_name) = nc_info.Attributes(aa).Value;
end


%% Print summary

fprintf('\n--- AMSR2 RSS L3 Closest Match ---\n');
fprintf('  EMIT pixel:  lat = %.3f, lon = %.3f, time = %.2f UTC\n', ...
    emit_lat, emit_lon, emit_time_utc);
fprintf('  AMSR2 match: lat = %.3f, lon = %.3f, time = %.2f UTC\n', ...
    amsr2_closest.lat, amsr2_closest.lon, amsr2_closest.time);
fprintf('  Pass: %d (%s)\n', best_pass, ...
    ternary(best_pass == 1, 'ascending ~1:30PM local', 'descending ~1:30AM local'));
fprintf('  Spatial distance: %.2f km\n', distance_km);
fprintf('  Time difference:  %.2f hours\n', min_time_diff);
fprintf('  Cloud LWP: %.4f kg/m^2\n', amsr2_closest.cloud_liquid_water);
fprintf('  File: %s\n', nc_filename);
fprintf('----------------------------------\n\n');


end



%% Helper function
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
