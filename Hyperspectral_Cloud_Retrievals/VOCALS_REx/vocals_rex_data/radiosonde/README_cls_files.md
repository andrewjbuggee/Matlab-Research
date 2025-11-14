# CLS File Reader for MATLAB

## Overview
This MATLAB function reads `.cls` sounding data files that may contain multiple soundings. The function automatically detects how many soundings are in the file and returns them in a structure array.

## Files
- `read_cls_file.m` - Main function to read .cls files
- `example_usage.m` - Example script demonstrating usage

## Data Structure

The function returns a **structure array** where each element represents one sounding:

```matlab
data = read_cls_file('filename.cls');

% data is a structure array with the following fields:
% data(i).data_type      - Type of sounding (e.g., "GAUS SOUNDING DATA/Ascending")
% data(i).project_id     - Project identifier (e.g., "VOCALS_2008")
% data(i).release_site   - Release site name (string)
% data(i).longitude      - Longitude in decimal degrees
% data(i).latitude       - Latitude in decimal degrees
% data(i).altitude       - Altitude in meters
% data(i).release_time   - UTC release time (datetime object)
% data(i).sonde_info     - Sonde/radiosonde information (string, or 'non existent')
% data(i).quality        - Quality information (string, or 'non existent')
% data(i).variable_names - Cell array of variable names (21 elements)
% data(i).units          - Cell array of units (21 elements)
% data(i).sounding_data  - Numeric data matrix (N rows x 21 columns)
```

## Basic Usage

### Read a file with multiple soundings:
```matlab
data = read_cls_file('VOCALS_2008_5mb_20081005.cls');
num_soundings = length(data);  % Get number of soundings
```

### Access individual soundings:
```matlab
% Access first sounding
first_sounding = data(1);

% Access second sounding
second_sounding = data(2);

% Get release time of first sounding
release_time = data(1).release_time;
```

### Loop through all soundings:
```matlab
for i = 1:length(data)
    fprintf('Sounding %d: %s at %s\n', i, ...
            data(i).release_site, datestr(data(i).release_time));
end
```

### Extract specific variables:
```matlab
% Find the column index for a variable
temp_idx = find(strcmp(data(1).variable_names, 'Temp'));
alt_idx = find(strcmp(data(1).variable_names, 'Alt'));

% Extract the data
temperature = data(1).sounding_data(:, temp_idx);
altitude = data(1).sounding_data(:, alt_idx);

% Plot
plot(temperature, altitude);
```

### Access metadata fields:
```matlab
% Check data type
if contains(data(1).data_type, 'GAUS')
    fprintf('This is a GAUS sounding\n');
end

% Check if quality information exists
if ~strcmp(data(1).quality, 'non existent')
    fprintf('Quality: %s\n', data(1).quality);
end

% Filter soundings by project
vocals_soundings = data(strcmp({data.project_id}, 'VOCALS_2008'));
```

## Advantages of Structure Array vs Cell Array

**Structure Array** (current implementation):
```matlab
data(1).release_time    % Easy to access
data(2).temperature     % Natural indexing
length(data)            % Get number of soundings
```

**Cell Array** (alternative):
```matlab
data{1}.release_time    % Requires curly braces
data{2}.temperature     % Less intuitive
length(data)            % Still works
```

The structure array is preferred because:
1. More intuitive syntax with parentheses
2. Better integration with MATLAB tools (tables, etc.)
3. Can still use array operations: `[data.longitude]` gets all longitudes
4. Standard convention for multiple similar objects in MATLAB

## Converting to Cell Array (if needed)

If you prefer a cell array, you can easily convert:
```matlab
data_struct = read_cls_file('filename.cls');
data_cell = num2cell(data_struct);
```

## Variable Names in the Data

The 21 columns in `sounding_data` correspond to:
1. Time (sec)
2. Press (mb)
3. Temp (C)
4. Dewpt (C)
5. RH (%)
6. Ucmp (m/s)
7. Vcmp (m/s)
8. spd (m/s)
9. dir (deg)
10. Wcmp (m/s)
11. Lon (deg)
12. Lat (deg)
13. Ele (deg)
14. Azi (deg)
15. Alt (m)
16. Qp (code)
17. Qt (code)
18. Qrh (code)
19. Qu (code)
20. Qv (code)
21. QdZ (code)

## Tips

1. **Check number of soundings first:**
   ```matlab
   data = read_cls_file('file.cls');
   fprintf('This file contains %d soundings\n', length(data));
   ```

2. **Compare soundings:**
   ```matlab
   if length(data) > 1
       time_diff = data(2).release_time - data(1).release_time;
       fprintf('Time between soundings: %.2f hours\n', hours(time_diff));
   end
   ```

3. **Extract all metadata at once:**
   ```matlab
   all_times = [data.release_time];
   all_lons = [data.longitude];
   all_lats = [data.latitude];
   ```

## Data Quality Control

The function automatically cleans the data by:

1. **Removing erroneous rows**: Rows where pressure ≥ 9999 OR temperature ≥ 999 OR temperature ≤ -999 are automatically removed
2. **Handling missing quality information**: If the "System Operator/Comments" header is not present, quality is set to "non existent"

Example:
```matlab
data = read_cls_file('file.cls');
% If erroneous rows were found, you'll see:
% "Removed 3 erroneous data row(s) from sounding 2"
```

## Flexible Header Parsing

The function handles variable header formats:
- **Required headers** (lines 1-5): Data Type, Project ID, Release Site, Location, Release Time
- **Optional headers**: 
  - Sonde Id/Sonde Type OR Radiosonde Serial Number
  - System Operator/Comments (quality information)
  - Reference Launch Data Source/Time
  - Post Processing Comments
- **Automatic detection**: Finds variable names, units, and data start by searching for keywords

This means the function works even if:
- Some soundings have different numbers of header lines
- Some soundings use "Sonde Id/Sonde Type" while others use "Radiosonde Serial Number"
- Header lines appear in slightly different orders
- Some optional headers are missing
- Some lines contain only forward slashes `/` as placeholders

**Example header variations handled:**

*Full header (GAUS soundings):*
```
Data Type:                         GAUS SOUNDING DATA/Ascending
Project ID:                        VOCALS_2008
Release Site Type/Site ID:         R/V Jose Olaya
Release Location (lon,lat,alt):    ...
UTC Release Time (y,m,d,h,m,s):    ...
Sonde Id/Sonde Type:               082033635/Vaisala RS92-SGP (ccGPS)
Reference Launch Data Source/Time: Manual Entry/01:55:09.76
System Operator/Comments:          pc/none, Good Sounding
Post Processing Comments:          Aspen Version
```

*Minimal header (Ron Brown soundings):*
```
Data Type:                         Ron Brown Soundings/Ascending
Project ID:                        VOCALS_2008
Release Site Type/Site ID:         R/V Ron Brown/WTEC
Release Location (lon,lat,alt):    ...
UTC Release Time (y,m,d,h,m,s):    ...
Radiosonde Serial Number:          Z4127134
/
/
/
```

## Error Handling

The function will display an error if:
- File cannot be opened
- Variable names or units lines cannot be found
- Data start line cannot be found

The function automatically handles:
- Any number of soundings (1, 2, 3, or more)
- Different numbers of data points in each sounding
- Different header formats between soundings in the same file
- Missing quality information
- Erroneous data values
