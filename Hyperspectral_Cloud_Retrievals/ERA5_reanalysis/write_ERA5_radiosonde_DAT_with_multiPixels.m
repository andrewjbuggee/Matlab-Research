%% Write ERA5 atmospheric profile to libRadtran .DAT format
%
% This function writes ERA5 reanalysis atmospheric profiles to a
% libRadtran-compatible .DAT file with pressure, temperature, and
% water vapor number density (molecules/cm³)
%
% Unlike AIRS data, ERA5 is model-based reanalysis data so no quality
% control flags are needed. ERA5 provides geopotential height directly,
% and specific humidity which we convert to water vapor number density.
%
% INPUTS:
% (1) era5 - structure containing ERA5 reanalysis data
%            Expected fields:
%            .temperature   - Temperature profile [K] (n_pixels x n_levels)
%            .pressure      - Pressure levels [hPa] (1 x n_levels) or (n_pixels x n_levels)
%            .specificHumidity - Specific humidity [kg/kg] (n_pixels x n_levels)
%            .geopotential  - Geopotential [m²/s²] (n_pixels x n_levels)
%            .geo.Latitude  - Latitude [degrees] (n_pixels x 1)
%            .geo.Longitude - Longitude [degrees] (n_pixels x 1)
%
% (2) folder_paths - structure containing the path to save the file
%            Must include 'atm_folder_path' field
%            Must include 'libRadtran_data' field (path to libRadtran data folder
%            containing atmmod/ with US standard atmosphere files)
%
% (3) idx - linear index of the pixel to extract profile (1 to num_pixels)
%
% (4) filename (optional) - custom filename. If not provided, generates automatic name
%
% (5) num_vars - if 2, use just temperature and pressure. if 3,
%                use temperature, pressure and water vapor
%
% (6) overlap_pixels - structure containing overlap pixel information
%                      Used to map between unique ERA5 pixels and EMIT pixels
%
% (7) us_std_atm_file (optional) - filename of US standard atmosphere file
%                                   Default: 'afglus.dat'
%
% (8) print_status_updates (optional) - flag to print status messages
%                                       Default: true
%
% OUTPUTS:
% (1) filename_fullPath - full path to the saved .DAT file
%
% (2) era5 - Updated structure with added field:
%       .datProfiles(idx) - Structure containing the profiles written to file:
%           .GP_height  - Geopotential height [m]
%           .T          - Temperature [K]
%           .p          - Pressure [hPa]
%           .q          - Specific humidity [kg/kg]
%           .vapor_concentration - Water vapor number density [molecules/cm³]
%           .vapor_massDensity   - Water vapor mass density [kg/m³]
%
% NOTES:
%   - ERA5 data is reanalysis/model data, so no quality control is needed
%   - Geopotential height is computed from ERA5's geopotential field
%   - Water vapor is written as number density (molecules/cm³) for libRadtran
%   - Pressure levels should be ordered from TOA to surface (low to high)
%     but the function will flip if needed for libRadtran format
%
% EXAMPLE USAGE:
%   [filename, era5] = write_ERA5_radiosonde_DAT_with_multiPixels(...
%       era5, folder_paths, 1, [], 3, overlap_pixels, 'afglus.dat', true);
%
% SEE ALSO:
%   write_AIRS_radiosonde_DAT_with_multiPixels
%   specificHumidity2waterVaporDensity
%   computeGeopotentialHeight
%
% By Andrew John Buggee
%%

function [filename_fullPath, era5] = write_ERA5_radiosonde_DAT_with_multiPixels(era5, folder_paths, idx, filename,...
    num_vars, overlap_pixels, us_std_atm_file, print_status_updates)

% Set default for print_status_updates if not provided
if nargin < 8
    print_status_updates = true;
end

% Check if folder_paths structure has the atmosphere folder field
if ~isfield(folder_paths, 'atm_folder_path')
    error([newline, 'folder_paths structure must contain "atm_folder_path" field with path to atmosphere folder', newline])
end

atm_folder_path = folder_paths.atm_folder_path;

% Make sure the path ends with a forward slash
if atm_folder_path(end) ~= '/'
    atm_folder_path = [atm_folder_path, '/'];
end

% Set default US standard atmosphere file if not provided
if isempty(us_std_atm_file)
    us_std_atm_file = 'afglus.dat';
end

%% Map overlap pixels to unique ERA5 pixels
% This handles the case where multiple EMIT pixels map to the same ERA5 pixel
% (ERA5 has coarser spatial resolution than EMIT)

% Check if overlap_pixels structure exists and has era5 field
if exist('overlap_pixels', 'var') && isfield(overlap_pixels, 'era5')
    unique_era5_pix = unique(overlap_pixels.era5.linear_idx);
    unique_pix_idx = zeros(1, length(overlap_pixels.era5.linear_idx));
    for xx = 1:length(unique_pix_idx)
        unique_pix_idx(xx) = find(unique_era5_pix == overlap_pixels.era5.linear_idx(xx));
    end
    % Get the actual ERA5 pixel index to use
    era5_pixel_idx = unique_pix_idx(idx);
else
    % No overlap structure provided, use idx directly
    era5_pixel_idx = idx;
end

%% Extract the profile data for the specified pixel

% Extract pressure levels
% Pressure is in hPa (mb)
if isfield(era5, 'pressure')
    if isvector(era5.pressure)
        pressure = double(era5.pressure(:));  % Make column vector
    else
        % Pressure varies spatially
        pressure = double(era5.pressure(era5_pixel_idx, :))';
    end
else
    error([newline, 'era5 structure must contain "pressure" field', newline])
end

% Extract temperature profile at the specified pixel
% Temperature is in Kelvin
if isfield(era5, 'temperature')
    temperature = era5.temperature(era5_pixel_idx, :)';  % (n_levels x 1)
else
    error([newline, 'era5 structure must contain "temperature" field', newline])
end

% Check for NaN values in temperature
if all(isnan(temperature))
    error([newline, 'All values of ERA5 temperature are NaN', newline])
end

% Extract specific humidity profile
% Specific humidity is in kg/kg
if isfield(era5, 'specificHumidity')
    q = era5.specificHumidity(era5_pixel_idx, :)';  % (n_levels x 1)
elseif isfield(era5, 'q')
    q = era5.q(era5_pixel_idx, :)';  % (n_levels x 1)
else
    error([newline, 'era5 structure must contain "specificHumidity" or "q" field', newline])
end

% Check for NaN values in specific humidity
if all(isnan(q))
    error([newline, 'All values of ERA5 specific humidity are NaN', newline])
end

%% Extract or compute geopotential height

if isfield(era5, 'geopotential')
    % ERA5 provides geopotential in m²/s²
    % Convert to geopotential height by dividing by standard gravity
    geopotential = era5.geopotential(era5_pixel_idx, :)';  % (n_levels x 1)
    Z_geopotential = geopotential / 9.80665;  % meters

    % Convert geopotential height to geometric height
    % From ERA5 documentation: h = R*Z/(R-Z)
    con = physical_constants;
    altitude = con.R_earth * Z_geopotential ./ (con.R_earth - Z_geopotential);  % meters

elseif isfield(era5, 'z')
    % Alternative field name for geopotential
    geopotential = era5.z(era5_pixel_idx, :)';  % (n_levels x 1)
    Z_geopotential = geopotential / 9.80665;  % meters

    % Convert to geometric height
    con = physical_constants;
    altitude = con.R_earth * Z_geopotential ./ (con.R_earth - Z_geopotential);  % meters

else
    % Compute geopotential height from temperature and specific humidity
    if print_status_updates
        warning('Geopotential field not found. Computing from hypsometric equation.');
    end
    Z_sfc = 0;  % Assume sea level (ERA5 data is typically over ocean for cloud studies)
    altitude = computeGeopotentialHeight(temperature, q, pressure .* 100, Z_sfc);  % meters
end

%% Convert specific humidity to water vapor mass density and number density

% Convert specific humidity to water vapor mass density using virtual temperature
% rho_v = [kg/m³], waterVapor_concentration_cm3 = [molecules/cm³]
use_virtual_temp = true;  % More accurate (~1% improvement)
[rho_v, waterVapor_concentration_cm3] = specificHumidity2waterVaporDensity(q, temperature, pressure .* 100, use_virtual_temp);

%% Verify total precipitable water (TPW)
% Compute TPW by integrating water vapor mass density over altitude

% For integration, we need data sorted from surface to TOA (increasing altitude)
% First, determine the current ordering
if altitude(1) > altitude(end)
    % Currently TOA to surface, need to flip
    altitude_sorted = flipud(altitude);
    rho_v_sorted = flipud(rho_v);
    pressure_sorted = flipud(pressure);
    temperature_sorted = flipud(temperature);
    q_sorted = flipud(q);
    waterVapor_numDensity_sorted = flipud(waterVapor_concentration_cm3);
else
    % Already surface to TOA
    altitude_sorted = altitude;
    rho_v_sorted = rho_v;
    pressure_sorted = pressure;
    temperature_sorted = temperature;
    q_sorted = q;
    waterVapor_numDensity_sorted = waterVapor_concentration_cm3;
end

% Compute total column water vapor
TPW_computed = trapz(altitude_sorted, rho_v_sorted);  % kg/m²

if print_status_updates
    fprintf('\n--- Total Column Water Vapor ---\n');
    fprintf('  Integrated from ERA5 profile: %.4f kg/m² (mm)\n', TPW_computed);
    fprintf('-------------------------------\n\n');
end

%% Handle missing data (NaN values)

% Check for NaN values in the profiles
if any(isnan(temperature)) || any(isnan(waterVapor_numDensity_sorted))
    warning(['Profile at pixel index ', num2str(idx), ...
        ' contains NaN values. These will be written to file but may cause issues in libRadtran.'])
end

%% Create filename

% Get latitude and longitude for this pixel
if isfield(era5, 'geo') && isfield(era5.geo, 'Latitude')
    lat = era5.geo.Latitude(era5_pixel_idx);
    lon = era5.geo.Longitude(era5_pixel_idx);
elseif isfield(era5, 'latitude')
    lat = era5.latitude(era5_pixel_idx);
    lon = era5.longitude(era5_pixel_idx);
else
    lat = NaN;
    lon = NaN;
end

if isempty(filename)
    % Generate automatic filename based on location and time
    if isfield(era5, 'metadata') && isfield(era5.metadata, 'start_year')
        year = era5.metadata.start_year{1};
        month = era5.metadata.start_month{1};
        day = era5.metadata.start_day{1};

        if num_vars == 3
            filename_radiosonde = sprintf('ERA5_profile_T-P-WV_lat%.2f_lon%.2f_%04d-%02d-%02d.dat', lat, lon, year, month, day);
        elseif num_vars == 2
            filename_radiosonde = sprintf('ERA5_profile_T-P_lat%.2f_lon%.2f_%04d-%02d-%02d.dat', lat, lon, year, month, day);
        end

    else
        if num_vars == 3
            filename_radiosonde = sprintf('ERA5_profile_T-P-WV_lat%.2f_lon%.2f.dat', lat, lon);
        elseif num_vars == 2
            filename_radiosonde = sprintf('ERA5_profile_T-P_lat%.2f_lon%.2f.dat', lat, lon);
        end
    end
else
    % User provided a filename — use it
    filename_radiosonde = filename;
    % Make sure it ends with .dat
    if ~contains(filename_radiosonde, '.dat') && ~contains(filename_radiosonde, '.DAT')
        filename_radiosonde = [filename_radiosonde, '.dat'];
    end
end

%% Write the radiosonde file

if num_vars == 2
    % libRadtran radiosonde files must have pressure increasing
    % (altitude decreasing) from top to bottom of file
    % We need TOA to surface (decreasing altitude, increasing pressure)
    if pressure_sorted(1) > pressure_sorted(end)
        % Currently surface to TOA, need to flip for libRadtran
        pressure_2write = flipud(pressure_sorted);
        temperature_2write = flipud(temperature_sorted);
    else
        % Already TOA to surface
        pressure_2write = pressure_sorted;
        temperature_2write = temperature_sorted;
    end

    fileID = fopen([atm_folder_path, filename_radiosonde], 'w');

    if fileID == -1
        error([newline, 'Could not open file for writing: ', atm_folder_path, filename_radiosonde, newline])
    end

    fprintf(fileID, '# extracted from ERA5 reanalysis data, Lat=%.4f, Lon=%.4f\n', lat, lon);
    fprintf(fileID, '#   p(hPa)  T(K)\n');

    for ii = 1:length(pressure_2write)
        fprintf(fileID, '%12.5f %5.1f\n', pressure_2write(ii), temperature_2write(ii));
    end

    fclose(fileID);

elseif num_vars == 3
    % libRadtran radiosonde files must have pressure increasing
    % (altitude decreasing) from top to bottom of file
    if pressure_sorted(1) > pressure_sorted(end)
        % Currently surface to TOA, need to flip
        pressure_2write = flipud(pressure_sorted);
        temperature_2write = flipud(temperature_sorted);
        vapor_density_2write = flipud(waterVapor_numDensity_sorted);
    else
        % Already TOA to surface
        pressure_2write = pressure_sorted;
        temperature_2write = temperature_sorted;
        vapor_density_2write = waterVapor_numDensity_sorted;
    end

    fileID = fopen([atm_folder_path, filename_radiosonde], 'w');

    if fileID == -1
        error([newline, 'Could not open file for writing: ', atm_folder_path, filename_radiosonde, newline])
    end

    fprintf(fileID, '# extracted from ERA5 reanalysis data, Lat=%.4f, Lon=%.4f\n', lat, lon);
    fprintf(fileID, '# Total column water vapor = %f kg/m2 (mm)\n', TPW_computed);
    fprintf(fileID, '#   p(hPa)  T(K)  H2O(cm-3)\n');

    for ii = 1:length(pressure_2write)
        fprintf(fileID, '%f %f %e\n', pressure_2write(ii), temperature_2write(ii), vapor_density_2write(ii));
    end

    fclose(fileID);
end

%% Define output filename with full path

filename_fullPath = [atm_folder_path, filename_radiosonde];

% Store the profiles in era5 structure
if num_vars == 3
    era5.datProfiles(idx).GP_height = altitude_sorted;                          % meters - geopotential height
    era5.datProfiles(idx).T = temperature_sorted;                               % K - temperature
    era5.datProfiles(idx).p = pressure_sorted;                                  % hPa - pressure
    era5.datProfiles(idx).q = q_sorted;                                         % kg/kg - specific humidity
    era5.datProfiles(idx).vapor_concentration = waterVapor_numDensity_sorted;   % molecules/cm³ - water vapor number density
    era5.datProfiles(idx).vapor_massDensity = rho_v_sorted;                     % kg/m³ - water vapor mass density
end

% Display success message
if print_status_updates
    disp(newline)
    fprintf('ERA5 radiosonde file written to: %s\n', filename_fullPath);
    fprintf('Pressure range: %.2f - %.2f hPa\n', min(pressure_2write), max(pressure_2write));
    fprintf('Temperature range: %.2f - %.2f K\n', min(temperature_2write), max(temperature_2write));
    if num_vars == 3
        fprintf('Water vapor density range: %.4e - %.4e molecules/cm³\n', ...
            min(vapor_density_2write), max(vapor_density_2write));
    end
    disp(newline)
end

end
