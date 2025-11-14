function data = read_cls_file(filename)
% READ_CLS_FILE Reads a .cls sounding data file with multiple soundings
%
% Input:
%   filename - path to the .cls file
%
% Output:
%   data - structure array where each element contains one sounding:
%       data(i).data_type      - type of sounding data (string)
%       data(i).project_id     - project identifier (string)
%       data(i).release_site   - name of release site (string)
%       data(i).longitude      - release longitude in decimal degrees
%       data(i).latitude       - release latitude in decimal degrees
%       data(i).altitude       - release altitude in meters
%       data(i).release_time   - UTC release time as datetime object
%       data(i).sonde_info     - sonde/radiosonde information (string, or 'non existent')
%       data(i).quality        - quality information (string, or 'non existent')
%       data(i).variable_names - cell array of variable names
%       data(i).units          - cell array of units for each variable
%       data(i).sounding_data  - matrix of sounding data (rows x 21 columns)
%                                with erroneous data rows removed

    % Open the file and read all lines
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Read all lines into a cell array
    all_lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            all_lines{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);
    
    % Find all lines that start with "Data Type:" to identify soundings
    sounding_starts = [];
    for i = 1:length(all_lines)
        if startsWith(strtrim(all_lines{i}), 'Data Type:')
            sounding_starts(end+1) = i; %#ok<AGROW>
        end
    end
    
    num_soundings = length(sounding_starts);
    fprintf('Found %d sounding(s) in file: %s\n', num_soundings, filename);
    
    % Initialize structure array
    data = struct('data_type', {}, 'project_id', {}, 'release_site', {}, ...
                  'longitude', {}, 'latitude', {}, 'altitude', {}, ...
                  'release_time', {}, 'sonde_info', {}, 'quality', {}, ...
                  'variable_names', {}, 'units', {}, 'sounding_data', {});
    
    % Process each sounding
    soundings_processed = 0;
    for s = 1:num_soundings
        fprintf('\n--- Processing Sounding %d/%d ---\n', s, num_soundings);
        
        % Determine the range of lines for this sounding
        start_line = sounding_starts(s);
        if s < num_soundings
            end_line = sounding_starts(s+1) - 1;
        else
            end_line = length(all_lines);
        end
        
        % Extract lines for this sounding
        sounding_lines = all_lines(start_line:end_line);
        
        % Parse this sounding
        sounding_data = parse_single_sounding(sounding_lines, s);
        
        % Check if sounding was skipped due to excessive bad data
        if isempty(sounding_data)
            fprintf('  Sounding %d skipped due to excessive erroneous data\n', s);
            continue;  % Skip this sounding
        end
        
        % Add to data array
        soundings_processed = soundings_processed + 1;
        data(soundings_processed) = sounding_data;
        
        % Display summary
        fprintf('Data Type: %s\n', data(soundings_processed).data_type);
        fprintf('Project: %s\n', data(soundings_processed).project_id);
        fprintf('Release Site: %s\n', data(soundings_processed).release_site);
        fprintf('Location: %.3f°, %.3f° (%.1f m)\n', ...
                data(soundings_processed).longitude, data(soundings_processed).latitude, data(soundings_processed).altitude);
        fprintf('Release Time: %s UTC\n', datestr(data(soundings_processed).release_time));
        fprintf('Sonde Info: %s\n', data(soundings_processed).sonde_info);
        fprintf('Quality: %s\n', data(soundings_processed).quality);
        fprintf('Number of data points: %d\n', size(data(soundings_processed).sounding_data, 1));
    end
    
    fprintf('\n=== Successfully read %d out of %d sounding(s) ===\n', soundings_processed, num_soundings);
end

function sounding = parse_single_sounding(lines, sounding_num)
% PARSE_SINGLE_SOUNDING Parse a single sounding from a cell array of lines
%
% Input:
%   lines        - cell array of strings containing one sounding
%   sounding_num - sounding number (for display purposes)
%
% Output:
%   sounding - structure containing the parsed sounding data

    sounding = struct();
    
    % Find key header lines by searching for their identifiers
    line_indices = struct();
    
    % Lines 1-5 are always consistent
    line_indices.data_type = 1;      % "Data Type:"
    line_indices.project_id = 2;     % "Project ID:"
    line_indices.release_site = 3;   % "Release Site Type/Site ID:"
    line_indices.location = 4;       % "Release Location (lon,lat,alt):"
    line_indices.release_time = 5;   % "UTC Release Time (y,m,d,h,m,s):"
    
    % Search for optional header lines
    line_indices.sonde_info = [];    % "Sonde Id/Sonde Type:" or "Radiosonde Serial Number:"
    line_indices.quality = [];       % "System Operator/Comments:"
    line_indices.variable_names = [];
    line_indices.units = [];
    line_indices.data_start = [];
    
    for i = 1:length(lines)
        line = strtrim(lines{i});
        
        % Look for sonde information (either format)
        if startsWith(line, 'Sonde Id/Sonde Type:') || ...
           startsWith(line, 'Radiosonde Serial Number:')
            line_indices.sonde_info = i;
        end
        
        % Look for System Operator/Comments line
        if startsWith(line, 'System Operator/Comments:')
            line_indices.quality = i;
        end
        
        % Look for variable names (contains "Time" "Press" "Temp" etc.)
        if contains(line, 'Time') && contains(line, 'Press') && contains(line, 'Temp') && ...
           contains(line, 'RH') && ~contains(line, 'UTC') && ~contains(line, 'Nominal')
            line_indices.variable_names = i;
        end
        
        % Look for units line (contains "sec" "mb" "C" etc.)
        if contains(line, 'sec') && contains(line, 'mb') && contains(line, 'deg')
            line_indices.units = i;
        end
        
        % Look for separator line (dashes)
        if startsWith(line, '---')
            line_indices.data_start = i + 1;
        end
    end
    
    % Parse Line 1: Data Type
    line1 = lines{line_indices.data_type};
    tokens = strsplit(line1, ':');
    sounding.data_type = strtrim(tokens{2});
    
    % Parse Line 2: Project ID
    line2 = lines{line_indices.project_id};
    tokens = strsplit(line2, ':');
    sounding.project_id = strtrim(tokens{2});
    
    % Parse Line 3: Release Site
    line3 = lines{line_indices.release_site};
    tokens = strsplit(line3, ':');
    sounding.release_site = strtrim(tokens{2});
    
    % Parse Line 4: Location (lon, lat, alt)
    line4 = lines{line_indices.location};
    tokens = strsplit(line4, ':');
    location_str = strtrim(tokens{2});
    % Extract the decimal degree values (last three comma-separated values)
    parts = strsplit(location_str, ',');
    sounding.longitude = str2double(strtrim(parts{end-2}));
    sounding.latitude = str2double(strtrim(parts{end-1}));
    sounding.altitude = str2double(strtrim(parts{end}));
    
    % Parse Line 5: UTC Release Time
    line5 = lines{line_indices.release_time};
    tokens = strsplit(line5, ':');
    time_str = strtrim(strjoin(tokens(2:end), ':'));
    % Extract year, month, day, hour, minute, second
    time_parts = strsplit(time_str, ',');
    year = str2double(strtrim(time_parts{1}));
    month = str2double(strtrim(time_parts{2}));
    day = str2double(strtrim(time_parts{3}));
    % Handle the time part which contains ':'
    time_component = strtrim(time_parts{4});
    time_vals = sscanf(time_component, '%d:%d:%d');
    hour = time_vals(1);
    minute = time_vals(2);
    second = time_vals(3);
    
    % Create datetime object
    sounding.release_time = datetime(year, month, day, hour, minute, second, 'TimeZone','UTC');
    
    % Parse Sonde Information (optional - line 6)
    if ~isempty(line_indices.sonde_info)
        line_sonde = lines{line_indices.sonde_info};
        tokens = strsplit(line_sonde, ':');
        sounding.sonde_info = strtrim(tokens{2});
    else
        sounding.sonde_info = 'non existent';
    end
    
    % Parse Quality (System Operator/Comments) - optional
    if ~isempty(line_indices.quality)
        line_quality = lines{line_indices.quality};
        tokens = strsplit(line_quality, ':');
        sounding.quality = strtrim(tokens{2});
    else
        sounding.quality = 'non existent';
        fprintf('  Note: Quality information not found for sounding %d\n', sounding_num);
    end
    
    % Parse Variable Names
    if ~isempty(line_indices.variable_names)
        line_vars = lines{line_indices.variable_names};
        var_names_str = strtrim(line_vars);
        sounding.variable_names = strsplit(var_names_str);
    else
        error('Could not find variable names line in sounding %d', sounding_num);
    end
    
    % Parse Units
    if ~isempty(line_indices.units)
        line_units = lines{line_indices.units};
        units_str = strtrim(line_units);
        sounding.units = strsplit(units_str);
    else
        error('Could not find units line in sounding %d', sounding_num);
    end
    
    % Parse Data
    if isempty(line_indices.data_start)
        error('Could not find data start line in sounding %d', sounding_num);
    end
    
    % Read all data lines
    sounding.sounding_data = [];
    for i = line_indices.data_start:length(lines)
        line = lines{i};
        if ~isempty(strtrim(line))
            % Try to parse the line as numeric data
            values = str2num(line); %#ok<ST2NM>
            if ~isempty(values) && length(values) == 21
                sounding.sounding_data = [sounding.sounding_data; values];
            end
        end
    end
    
    % Remove erroneous data rows
    % Check for rows where pressure >= 9999 or temperature >= 999 or <= -999 or dewpoint >= 999 or <= -999
    if ~isempty(sounding.sounding_data)
        num_rows_before = size(sounding.sounding_data, 1);
        
        % Find column indices
        press_idx = find(strcmp(sounding.variable_names, 'Press'), 1);
        temp_idx = find(strcmp(sounding.variable_names, 'Temp'), 1);
        dewpt_idx = find(strcmp(sounding.variable_names, 'Dewpt'), 1);
        
        if isempty(press_idx) || isempty(temp_idx) || isempty(dewpt_idx)
            % If column names don't match exactly, assume standard positions
            press_idx = 2;  % Pressure is column 2
            temp_idx = 3;   % Temperature is column 3
            dewpt_idx = 4;  % Dewpoint is column 4
        end
        
        % Identify bad rows (pressure >= 9999 OR temperature >= 999 OR temperature <= -999 
        %                     OR dewpoint >= 999 OR dewpoint <= -999)
        bad_rows = (sounding.sounding_data(:, press_idx) >= 9999) | ...
                   (sounding.sounding_data(:, temp_idx) >= 999) | ...
                   (sounding.sounding_data(:, temp_idx) <= -999) | ...
                   (sounding.sounding_data(:, dewpt_idx) >= 999) | ...
                   (sounding.sounding_data(:, dewpt_idx) <= -999);
        
        num_bad = sum(bad_rows);
        percent_bad = (num_bad / num_rows_before) * 100;
        
        % Check if more than 50% of data is bad
        if percent_bad > 50
            fprintf('  WARNING: %.1f%% of data is erroneous (>50%%) - SKIPPING this sounding\n', percent_bad);
            sounding = [];  % Return empty to signal this sounding should be skipped
            return;
        end
        
        % Remove bad rows
        sounding.sounding_data(bad_rows, :) = [];
        
        num_rows_after = size(sounding.sounding_data, 1);
        num_removed = num_rows_before - num_rows_after;
        
        if num_removed > 0
            fprintf('  Removed %d erroneous data row(s) from sounding %d (%.1f%% of data)\n', ...
                    num_removed, sounding_num, percent_bad);
        end
    end
end
