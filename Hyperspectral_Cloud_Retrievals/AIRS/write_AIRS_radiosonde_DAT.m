%% Write AIRS atmospheric profile to libRadtran .DAT format

% This function writes AIRS L2 retrieved atmospheric profiles to a 
% libRadtran-compatible .DAT file with pressure, temperature, and 
% relative humidity

% INPUTS:
% (1) AIRS - structure containing AIRS L2 data from readAIRS_L2_data()
%            Note: This should be the reduced AIRS structure after using
%            remove_unwanted_airs_data(), so all spatial data is 1D
% (2) folder_paths - structure containing the path to save the file
% (3) idx - linear index of the pixel to extract profile (1 to num_pixels)
% (4) filename (optional) - custom filename. If not provided, generates automatic name
% (5) num_vars - if 2, use just temperature and pressure from AIRS. if 3,
%                use temperature, pressure and relative humidity

% OUTPUTS:
% (1) filename_fullPath - full path to the saved .DAT file

% By Andrew John Buggee

%%

function [filename_fullPath] = write_AIRS_radiosonde_DAT(airs, folder_paths, idx, filename, num_vars)

% Check if folder_paths structure has the atmosphere folder field
if ~isfield(folder_paths, 'atm_folder_path')
    error([newline, 'folder_paths structure must contain "atm_folder_path" field with path to atmosphere folder', newline])
end

atm_folder_path = folder_paths.atm_folder_path;

% Make sure the path ends with a forward slash
if atm_folder_path(end) ~= '/'
    atm_folder_path = [atm_folder_path, '/'];
end

%% Extract the profile data for the specified pixel

% Extract pressure levels (standard pressure for temperature)
% Pressure is in hPa (mb)
pressure = double(airs.pressStd{1});  % (StdPressureLev = 28)

% Extract temperature profile at the specified pixel
% After using remove_unwanted_airs_data, temp.prof_std is (num_pixels x StdPressureLev)
% Temperature is in Kelvin
temperature = airs.temp.prof_std(idx, :)';  % (StdPressureLev x 1)

% Check to see if there are nan values. If so, spline fit the nan values
if sum(isnan(temperature))==1
    
    % remove the nan value
    idx_nan = isnan(temperature);

    temperature = exp(interp1(log(pressure(~idx_nan)), log(temperature(~idx_nan)),...
        log(pressure), 'linear', 'extrap')); 

elseif sum(isnan(relHum_H2O_levels))>1

    error([newline, 'There are a lot of nans in the Relative humidity data!', newline])

end

% Extract relative humidity profile
% After using remove_unwanted_airs_data, H2O.RelHum is (num_pixels x H2OPressureLev)
% RelHum is on H2O pressure levels (15 levels), need to map to pressStd (28 levels)
relHum_H2O_levels = airs.H2O.RelHum(idx, :)';  % (H2OPressureLev x 1)
pressH2O = double(airs.pressH2O{1});  % (H2OPressureLev = 15)

% Check to see if there are nan values. If so, spline fit the nan values
if sum(isnan(relHum_H2O_levels))==1
    
    % remove the nan value
    idx_nan = isnan(relHum_H2O_levels);

    relHum_H2O_levels = exp(interp1(log(pressH2O(~idx_nan)), log(relHum_H2O_levels(~idx_nan)),...
        log(pressH2O), 'linear', 'extrap')); 

elseif sum(isnan(relHum_H2O_levels))>1

    error([newline, 'There are a lot of nans in the Relative humidity data!', newline])

end


% Interpolate relative humidity to standard pressure levels
% Use linear interpolation in log-pressure space for atmospheric profiles
relHum = exp(interp1(log(pressH2O), log(relHum_H2O_levels), log(pressure), 'linear', 'extrap'));  % (StdPressureLev x 1)

% Ensure relative humidity is between 0 and 100%
relHum(relHum < 0) = 0;
relHum(relHum > 100) = 100;

% Convert relative humidity from percent to fraction (0-1 range)
% libRadtran uses fractional relative humidity, not percentage
relHum = relHum / 100;  % Convert from % to fraction

%% Handle missing data (NaN values)

% Check for NaN values in the profile
if any(isnan(temperature)) || any(isnan(relHum))
    warning(['Profile at pixel index ', num2str(idx), ...
        ' contains NaN values. These will be written to file but may cause issues in libRadtran.'])
end

%% Create filename if not provided

if nargin < 4 || isempty(filename)
    % Generate automatic filename based on location and time
    lat = airs.geo.Latitude(idx);
    lon = airs.geo.Longitude(idx);
    
    % Get time information
    if isfield(airs, 'metadata') && isfield(airs.metadata, 'start_year')
        year = airs.metadata.start_year{1};
        month = airs.metadata.start_month{1};
        day = airs.metadata.start_day{1};

        if num_vars==2

            filename = sprintf('AIRS_profile_T-P_lat%.2f_lon%.2f_%04d%02d%02d.dat', lat, lon, year, month, day);
        elseif num_vars==3

            filename = sprintf('AIRS_profile_T-P-RH_lat%.2f_lon%.2f_%04d%02d%02d.dat', lat, lon, year, month, day);
        end

    else
        
        if num_vars==2

            filename = sprintf('AIRS_profile_T-P_lat%.2f_lon%.2f.dat', lat, lon);

        elseif num_vars==3

            filename = sprintf('AIRS_profile_T-P-RH_lat%.2f_lon%.2f.dat', lat, lon);

        end

    end
end

% Make sure filename ends with .dat (lowercase to match example)
if ~contains(filename, '.dat') && ~contains(filename, '.DAT')
    filename = [filename, '.dat'];
end

%% Prepare data for writing

% libRadtran radiosonde files must have pressure increasing 
% (altitude decreasing) from top to bottom of file
% Check if pressure is in ascending or descending order
if pressure(1) > pressure(end)
    % Pressure is ascending (surface to top), need to flip
    pressure_2write = fliplr(pressure);
    temperature_2write = fliplr(temperature);
    relHum_2write = fliplr(relHum);
else
    % Pressure is already descending (top to surface)
    pressure_2write = pressure;
    temperature_2write = temperature;
    relHum_2write = relHum;
end

%% Write the file

if num_vars==2

    % Open file for writing
    fileID = fopen([atm_folder_path, filename], 'w');

    if fileID == -1
        error([newline, 'Could not open file for writing: ', atm_folder_path, filename, newline])
    end

    % Write header comment - single line format like the example
    fprintf(fileID, '# extracted from AIRS L2 data, Lat=%.4f, Lon=%.4f\n', ...
        airs.geo.Latitude(idx), airs.geo.Longitude(idx));

    % Write column headers - exactly matching the example format
    fprintf(fileID, '#   p(hPa)  T(K)\n');

    % Write the data - three columns with specific formatting to match example
    % Format: pressure with 5 decimal places, temperature with 1 decimal, RH in scientific notation
    for ii = 1:length(pressure_2write)
        fprintf(fileID, '%12.5f %5.1f\n', pressure_2write(ii), temperature_2write(ii));
    end

    % Close the file
    fclose(fileID);


elseif num_vars==3


    % Open file for writing
    fileID = fopen([atm_folder_path, filename], 'w');

    if fileID == -1
        error([newline, 'Could not open file for writing: ', atm_folder_path, filename, newline])
    end

    % Write header comment - single line format like the example
    fprintf(fileID, '# extracted from AIRS L2 data, Lat=%.4f, Lon=%.4f\n', ...
        airs.geo.Latitude(idx), airs.geo.Longitude(idx));

    % Write column headers - exactly matching the example format
    fprintf(fileID, '#   p(hPa)  T(K)  h2o(relative humidity)\n');

    % Write the data - three columns with specific formatting to match example
    % Format: pressure with 5 decimal places, temperature with 1 decimal, RH in scientific notation
    for ii = 1:length(pressure_2write)
        fprintf(fileID, '%12.5f %5.1f %e\n', pressure_2write(ii), temperature_2write(ii), relHum_2write(ii));
    end

    % Close the file
    fclose(fileID);



end

%% Define output filename with full path

filename_fullPath = [atm_folder_path, filename];

% Display success message
fprintf('AIRS profile written to: %s\n', filename_fullPath);
fprintf('Pressure range: %.2f - %.2f hPa\n', min(pressure_2write), max(pressure_2write));
fprintf('Temperature range: %.2f - %.2f K\n', min(temperature_2write), max(temperature_2write));
fprintf('Relative humidity range: %.4f - %.4f (fraction)\n', min(relHum_2write), max(relHum_2write));

end