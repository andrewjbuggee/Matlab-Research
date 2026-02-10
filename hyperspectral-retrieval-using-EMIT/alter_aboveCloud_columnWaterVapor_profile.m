% Function to alter the above cloud column amount of the water vapor
% density profile

% This function supports two types of atmosphere files:
%   1. US Standard Atmosphere .DAT files (9 columns: z, p, T, air, O3, O2, H2O, CO2, NO2)
%   2. AIRS radiosonde .DAT files (2-3 columns: p [mbar], T [K], WV [CM-3])
%      Created by write_AIRS_radiosonde_DAT_with_multiPixels()

% INPUTS:

% (1) inputs - the input structure for creating the input files

% (2) aboveCloudTotal - This is the total column water vapor amount above
% cloud used to alter the profile - kg/m^2

% (3) atm_folder_path - the location of libRadtrans atmosphere files and
% where to save the water vapor file

% By Andrew John Buggee

%%

function [filename_fullPath] = alter_aboveCloud_columnWaterVapor_profile(inputs, aboveCloudTotal,...
    atm_folder_path, airs_datProfiles, pixel_num)


if isfield(inputs.RT, 'use_radiosonde_file') == true && inputs.RT.use_radiosonde_file == true

    % extract the AIRS filename
    str = extractBetween(inputs.RT.radiosonde_file_T_P_WV, '/AIRS', '.dat');

    % create the modified atm file name
    filename = ['AIRS', str{1}, '_H2O_MODIFIED_', num2str(aboveCloudTotal), 'mm-aboveCloud',...
        '_cloudTop_', num2str(inputs.RT.z_topBottom(1)), 'km.DAT'];

else

    % create the modified atm file name
    filename = [inputs.RT.atm_file(1:end-4), '_H2O_MODIFIED_', num2str(aboveCloudTotal), 'mm-aboveCloud',...
        '_cloudTop_', num2str(inputs.RT.z_topBottom(1)), 'km.DAT'];

end


%% Read in the atmosphere profile and determine its type

% Check if this is an AIRS file by looking at the filename
useAIRS_file = isfield(inputs.RT, 'use_radiosonde_file') == true && inputs.RT.use_radiosonde_file == true && ...
    contains(inputs.RT.radiosonde_file_T_P_WV, 'AIRS', 'IgnoreCase', true);

if useAIRS_file
    %% ---- AIRS radiosonde file ----
    % These files have 2-3 columns: pressure (hPa), temperature (K), [relative humidity (fraction)]
    % Created by write_AIRS_radiosonde_DAT_with_multiPixels()

    % [z, waterVapor_column, airs_data] = read_AIRS_profile_for_H2O_scaling(inputs, atm_folder_path);
    % check to make sure the airs_datProfiles exits
    if exist("airs_datProfiles", "var") == false

        error([newline, 'Need the AIRS dat file profiles.', newline])

    else

        z = airs_datProfiles(pixel_num).GP_new;  % meters - geopotential height
        waterVapor_column = airs_datProfiles(pixel_num).vapor_concentration;    % #/cm^3 - water vapor number concentration

        % convert water vapor densities to m^(-3)
        waterVapor_column = waterVapor_column * 1e6;  % molecules/m^3
        

    end

else
    %% ---- US Standard Atmosphere file ----
    % These files have 9 columns and are stored in the BASE atmmod folder

    % these files are always stored in the BASE folder
    position = strfind(atm_folder_path, 'atmmod');
    if isnumeric(position)==true
        % then the substring exists. Good! It should
        % stop at the end of atmmod, the base folder
        atm = read_libRadtran_atm_dat_profiles_ver2([atm_folder_path(1:position+5),'/', inputs.RT.atm_file]);

    elseif isempty(position)

        error([newline, 'I cant find the base atmmod folder where the unaltered atmosphere .DAT files are stored.', newline])

    end

    % altitude values are listed in the 1st column (km)
    z = atm(:,1);
    % convert to meters
    z = z*1e3; % (m)

    % Water vapor densities are listed in the 7th column (cm^(-3))
    waterVapor_column = atm(:, 7);
    % convert water vapor densities to m^(-3)
    waterVapor_column = waterVapor_column * 1e6;  % molecules/m^3

end


%% Solve for the scalar value that alters the above cloud column water vapor amount

% First, interpolate the profile so that the cloud top height and the
% sensor altitude are included
if ischar(inputs.RT.sensor_altitude) && strcmp(inputs.RT.sensor_altitude, 'toa')==true

    % In addition, add another data z value 1 meter above cloud top. This
    % is the height where the scaling will begin so that, if libRadtran
    % interpolates within the cloud to get the value of water vapor density
    % within cloud, there isn't a difference between the original profile
    % and the new one.
    new_z = sort([z; inputs.RT.z_topBottom(1)*1e3; (inputs.RT.z_topBottom(1)*1e3)+1], 'descend');  % m

    % water vapor displays a more linear behavior in log space. Interpolate
    % the log of water vapor density
    waterVapor_column_interp = exp(interp1(z, log(waterVapor_column), new_z, "linear"));    % moleules/m^3

    % set the index for the sensor location
    idx_sensor = find(new_z==max(new_z));

elseif isdouble(inputs.RT.sensor_altitude)

    % In addition, add another data z value 1 meter above cloud top. This
    % is the height where the scaling will begin so that, if libRadtran
    % interpolates within the cloud to get the value of water vapor density
    % within cloud, there isn't a difference between the original profile
    % and the new one.
    new_z = sort([z; inputs.RT.z_topBottom(1)*1e3; (inputs.RT.z_topBottom(1)*1e3)+1;...
        inputs.RT.sensor_altitude], 'descending');

    waterVapor_column_interp = interp1(z, waterVapor_column, new_z, "linear");    % moleules/m^3

    % set the index for the sensor location
    idx_sensor = find(new_z==inputs.RT.sensor_altitude);

end


% Solve for the scalar value by integrating column water vapor from JUST above cloud
% top to sensor location

idx_justAboveCloudTop = find(new_z==((inputs.RT.z_topBottom(1)*1e3)+1));

% integrate from cloud top to sensor location and convert the density
% profile from molecules/m^3 to kg/m^2
con = physical_constants;
aboveCloud_columnAmount = -(con.Mol_mass_h2o_vap/con.N_A) *...
    trapz(new_z(idx_sensor:idx_justAboveCloudTop), waterVapor_column_interp(idx_sensor:idx_justAboveCloudTop));  % kg/m^2

% solve for the scalar constant
a = aboveCloudTotal / aboveCloud_columnAmount;


%% rescale the water vapor profile using the new scalar constant

waterVapor_column_2Write = waterVapor_column_interp;      % moleules/m^3

waterVapor_column_2Write(idx_sensor:idx_justAboveCloudTop) = a.*waterVapor_column_interp(idx_sensor:idx_justAboveCloudTop);  % moleules/m^3

% convert this back to molcules per cubic centimeter
waterVapor_column_2Write = waterVapor_column_2Write ./ 1e6;  % moleules/cm^3

% Convert the z vector to km
z_2write = new_z./1e3;  % (km)


%% Write the file

if useAIRS_file == false

    % Then write a .dat file that will update the H2O concentration profile
    % from some US standard atmosphere

    % Create the denisty file
    fileID = fopen([atm_folder_path,filename], 'w');

    % fprintf writes lines in our text file from top to botom
    % wc.DAT files are written with the higher altitudes at the top, and the
    % surface at the bottom

    % to write column vectors in a text file, we have to store them as row
    % vectors
    toWrite = [z_2write'; waterVapor_column_2Write'];

    % Create the opening comment lines of the water vapor DAT file

    fprintf(fileID, '%s %s %s \n','#','z', 'water-vapor-density');
    fprintf(fileID, '%s %s %s \n','#','(km)','(molecules/cm^3)');

    % Write in the data
    fprintf(fileID,'%f %f \n', toWrite);
    fclose(fileID);


    %% If using AIRS file, also write a new radiosonde .dat file with scaled RH
    %  The new radiosonde file uses the INTERPOLATED grid (new_z) which includes
    %  cloud top height, so that pressure, temperature, and RH are all on the
    %  same grid as the scaled water vapor density profile.

else


    % % Get original AIRS data for interpolation
    % pressure_orig = airs_data.pressure;      % hPa
    % temperature_orig = airs_data.temperature; % K
    % z_orig = airs_data.z;                     % m (sorted descending, high alt at top)
    % 
    % % Interpolate pressure and temperature to the new_z grid (which includes cloud top height)
    % % new_z is in meters, sorted descending (high altitude at top)
    % pressure_interp = exp(interp1(z_orig, log(pressure_orig), new_z, 'linear', 'extrap'));  % hPa (log-linear interp for pressure)
    % temperature_interp = interp1(z_orig, temperature_orig, new_z, 'linear', 'extrap');      % K
    % 
    % % Convert scaled water vapor density to relative humidity on the new_z grid
    % % waterVapor_column_2Write is in molecules/cm^3, need to convert to molecules/m^3
    % waterVapor_scaled = waterVapor_column_2Write .* 1e6;  % molecules/m^3
    % 
    % % Convert water vapor number density back to vapor pressure
    % % n = (e * N_A) / (R * T)  =>  e = (n * R * T) / N_A
    % R_universal = 8.314462;  % J/(mol·K)
    % e_Pa_scaled = (waterVapor_scaled .* R_universal .* temperature_interp) ./ con.N_A;  % Pa
    % e_hPa_scaled = e_Pa_scaled / 100;  % hPa
    % 
    % % Compute saturation vapor pressure (Bolton 1980)
    % T_celsius = temperature_interp - 273.15;
    % e_sat = 6.112 * exp((17.67 * T_celsius) ./ (T_celsius + 243.5));  % hPa
    % 
    % % Compute scaled relative humidity
    % relHum_scaled = e_hPa_scaled ./ e_sat;  % fraction (0-1)
    % 
    % % Ensure RH is between 0 and 1
    % relHum_scaled(relHum_scaled < 0) = 0;
    % relHum_scaled(relHum_scaled > 1) = 1;
    % 
    % % Prepare data for writing - need pressure increasing (low to high, i.e., TOA to surface)
    % % new_z is sorted descending (high alt at top), so pressure_interp is increasing (low to high)
    % % which is already in the correct order for radiosonde format (TOA to surface)
    % if pressure_interp(1) > pressure_interp(end)
    %     % Surface to TOA, need to flip for radiosonde format
    %     pressure_2write_rs = flipud(pressure_interp);
    %     temperature_2write_rs = flipud(temperature_interp);
    %     relHum_2write_rs = flipud(relHum_scaled);
    % else
    %     % Already TOA to surface (low pressure to high pressure)
    %     pressure_2write_rs = pressure_interp;
    %     temperature_2write_rs = temperature_interp;
    %     relHum_2write_rs = relHum_scaled;
    % end




    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % ** If using the radiosonde file, interpolate temperature and pressure
    % to the new altitude grid! **
    % ---------------------------------------------------------------------
    % The new altitude grid (new_z) includes two new levels:
    %   - cloud top height: inputs.RT.z_topBottom(1)*1e3 (meters)
    %   - 1 meter above cloud top: inputs.RT.z_topBottom(1)*1e3 + 1 (meters)
    % We need to interpolate pressure and temperature to these new levels.

    % Extract original profiles from airs_datProfiles
    z_orig = airs_datProfiles(pixel_num).GP_new;    % meters (sorted surface-to-TOA, ascending)
    p_orig = airs_datProfiles(pixel_num).p';        % hPa (sorted surface-to-TOA, descending pressure)
    T_orig = airs_datProfiles(pixel_num).T;         % K (sorted surface-to-TOA)

    new_p = exp(interp1(z_orig, log(p_orig), new_z, 'linear', 'extrap'));       % hPa
    new_T = interp1(z_orig, T_orig, new_z, 'linear', 'extrap');       % K

   

    % Prepare data for writing - need pressure increasing (low to high, i.e., TOA to surface)
    % new_z is sorted descending (high alt at top), so pressure_new should already be
    % increasing from low (TOA) to high (surface)
    if new_p(1) > new_p(end)
        % Surface to TOA order, need to flip for radiosonde format
        pressure_2write_rs = flipud(new_p);
        temperature_2write_rs = flipud(new_T);
    else
        % Already TOA to surface (low pressure to high pressure) - correct order
        pressure_2write_rs = new_p;
        temperature_2write_rs = new_T;
        % waterVapor_column_2Write is already in the correct order
    end
    % ---------------------------------------------------------------------

    % Write the radiosonde file (same format as write_AIRS_radiosonde_DAT_with_multiPixels)
    fileID_rs = fopen([atm_folder_path,filename], 'w');

    if fileID_rs == -1
        error([newline, 'Could not open file for writing: ', [atm_folder_path,filename], newline])
    end

    % Write header comment
    fprintf(fileID_rs, '# Modified AIRS profile with scaled above-cloud water vapor (%.2f mm)\n', aboveCloudTotal);
    fprintf(fileID_rs, '# Interpolated to include cloud top height at %.3f km\n', inputs.RT.z_topBottom(1));

    % Write column headers
    fprintf(fileID_rs, '#   p(hPa)  T(K)  H2O(CM-3)\n');

    % Write the data
    for ii = 1:length(pressure_2write_rs)
        fprintf(fileID_rs, '%f %f %e\n', pressure_2write_rs(ii), temperature_2write_rs(ii), waterVapor_column_2Write(ii));
    end

    fclose(fileID_rs);

end



filename_fullPath = [atm_folder_path, filename];




end




%% ========================================================================
%  Helper function to read AIRS radiosonde files and convert to z/H2O format
%  ========================================================================

function [z, vapor_concentration, rho_v] = read_AIRS_profile_for_H2O_scaling(inputs, atm_folder_path)
% READ_AIRS_PROFILE_FOR_H2O_SCALING reads an AIRS radiosonde .DAT file and
% converts it to altitude (m) and water vapor density (molecules/m^3)
%
% AIRS files created by write_AIRS_radiosonde_DAT_with_multiPixels() have:
%   - Column 1: pressure (hPa)
%   - Column 2: temperature (K)
%   - Column 3 (optional): relative humidity (fraction, 0-1)
%
% This function:
%   1. Reads the AIRS file
%   2. Converts pressure to altitude using the hypsometric equation
%   3. Converts relative humidity to water vapor number density
%
% OUTPUTS:
%   z               - altitude profile (m), sorted descending (high alt at top)
%   waterVapor_column - water vapor number density (molecules/m^3)
%   airs_data       - structure containing original AIRS data:
%                     .pressure (hPa), .temperature (K), .relHum (fraction), .z (m)



% Read the AIRS file
fid = fopen(inputs.RT.radiosonde_file_T_P_WV, 'r');
if fid == -1
    error(['Could not open AIRS file: ', inputs.RT.radiosonde_file_T_P_WV]);
end

% Read all lines, skipping comments
allLines = {};
tline = fgetl(fid);
while ischar(tline)
    allLines{end+1, 1} = tline; %#ok<AGROW>
    tline = fgetl(fid);
end
fclose(fid);

% Filter out empty lines and comment lines
validLines = {};
for ii = 1:length(allLines)
    trimmedLine = strtrim(allLines{ii});
    if ~isempty(trimmedLine) && trimmedLine(1) ~= '#'
        validLines{end+1, 1} = trimmedLine; %#ok<AGROW>
    end
end

% Determine number of columns by reading first data line
firstLineData = sscanf(validLines{1}, '%f');
numCols = length(firstLineData);

% Parse data
data = zeros(length(validLines), numCols);
for ii = 1:length(validLines)
    data(ii, :) = sscanf(validLines{ii}, '%f', numCols);
end

% Extract columns
% first row starts at TOA, last row is at the surface
pressure = data(:, 1);      % hPa
temperature = data(:, 2);   % K

if numCols >= 3
    vapor_concentration = data(:, 3);    % #/cm^3
else
    % If no RH column, use a default profile or error
    error([newline, 'AIRS file does not contain relative humidity column. ', ...
        'Cannot compute water vapor profile.', newline]);
end

% Convert pressure from hPa to Pa for calculations
pressure_Pa = pressure * 100;  % Pa

%% Convert pressure to altitude using hypsometric equation
% Uses mean virtual temperature between layers
% Z2 - Z1 = (R_d / g) * T_v_mean * ln(p1/p2)

% Physical constants
R_d = 286.933175;   % Specific gas constant for dry air - same value used by libRadtran v 2.0.6 [J/(kg·K)]
g = 9.799999;         % Value for gravity at z=0 km used by libRadtran v 2.0.6 [m/s²]
epsilon = 0.622;     % Ratio R_d/R_v

% Convert relative humidity to specific humidity for virtual temperature
% First, compute saturation vapor pressure using Clausius-Clapeyron (Bolton 1980)
T_celsius = temperature - 273.15;
e_sat = 6.112 * exp((17.67 * T_celsius) ./ (T_celsius + 243.5));  % hPa

% Actual vapor pressure
e = vapor_concentration .* e_sat;  % hPa

% Specific humidity: q = epsilon * e / (p - (1-epsilon)*e)
q = epsilon * e ./ (pressure - (1 - epsilon) * e);  % kg/kg

% Virtual temperature
T_v = temperature ./ (1 - (1 - epsilon) .* q);  % K

% Integrate from surface (highest pressure) upward
% AIRS files are ordered from low pressure (TOA) to high pressure (surface)
% We need to flip for integration from surface upward
if pressure(1) < pressure(end)
    % Data is TOA to surface, flip it
    pressure_Pa_sorted = flipud(pressure_Pa);
    T_v_sorted = flipud(T_v);
    needsFlip = true;
else
    % Data is surface to TOA
    pressure_Pa_sorted = pressure_Pa;
    T_v_sorted = T_v;
    needsFlip = false;
end

% Compute altitude profile
nlevels = length(pressure_Pa_sorted);
z_sorted = zeros(nlevels, 1);
z_sorted(1) = 0;  % Surface altitude = 0 m

for kk = 1:(nlevels-1)
    T_v_mean = (T_v_sorted(kk) + T_v_sorted(kk+1)) / 2;
    thickness = (R_d / g) * T_v_mean * log(pressure_Pa_sorted(kk) / pressure_Pa_sorted(kk+1));
    z_sorted(kk+1) = z_sorted(kk) + thickness;
end

% Flip back to original order if needed
if needsFlip
    z = flipud(z_sorted);  % m
else
    z = z_sorted;  % m
end



% ----- Compute altitude using hypsometric equation -----
% Z2 - Z1 = (R_d / g) * T_v_mean * ln(p1/p2)
R_d = 286.933175;   % Specific gas constant for dry air - same value used by libRadtran v 2.0.6 [J/(kg*K)]
g = 9.799999;         % Value for gravity at z=0 km used by libRadtran v 2.0.6 [m/s^2]

% Virtual temperature
T_v_28 = temperature ./ (1 - (1 - epsilon) .* q_28);  % K

% Sort pressure from high (surface) to low (TOA) for upward integration
[pressure_sorted, sort_idx_hyps] = sort(pressure, 'descend');
T_v_sorted = T_v_28(sort_idx_hyps);

% Integrate upward from surface
nlevels = length(pressure_sorted);
z_hyps_sorted = zeros(nlevels, 1);
z_hyps_sorted(1) = 0;  % Surface altitude = 0 m

for kk = 1:(nlevels-1)
    T_v_mean = (T_v_sorted(kk) + T_v_sorted(kk+1)) / 2;
    thickness = (R_d / g) * T_v_mean * log(pressure_sorted(kk) / pressure_sorted(kk+1));
    z_hyps_sorted(kk+1) = z_hyps_sorted(kk) + thickness;  % meters
end

%% Convert relative humidity to water vapor number density (molecules/m^3)
% Using the relationship: n = (e * N_A) / (R * T)
% where e is vapor pressure in Pa, R is universal gas constant

% ---------------------------------------------------------------------------------
% ** DONT NEED AT THE MOMENT BECAUSE RADIOSONDE FILES REPORT NUMBER CONCENTRATION **
% ---------------------------------------------------------------------------------

% con = physical_constants;
% R_universal = 8.314462;  % J/(mol·K) - universal gas constant
% 
% % Vapor pressure in Pa
% e_Pa = e * 100;  % Convert hPa to Pa
% 
% % Number density: n = (e * N_A) / (R * T)
% % molecules/m^3
% waterVapor_column = (e_Pa .* con.N_A) ./ (R_universal .* temperature);  % molecules/m^3
% 
% % Ensure z is sorted in descending order (high altitude at top) to match
% % the expected format for the rest of the function
% if z(1) < z(end)
%     z = flipud(z);
%     waterVapor_column = flipud(waterVapor_column);
%     % Also flip the airs_data to match
%     pressure = flipud(pressure);
%     temperature = flipud(temperature);
%     vapor_concentration = flipud(vapor_concentration);
% end

% ---------------------------------------------------------------------------------
% convert to mass per unit volume
con = physical_constants;

rho_v = (vapor_concentration * 1e6) .* (con.Mol_mass_h2o_vap/con.N_A);  % kg/m^3
% ---------------------------------------------------------------------------------

% Store original AIRS data in output structure (in same order as z output)
% airs_data.pressure = pressure;       % hPa
% airs_data.temperature = temperature; % K
% airs_data.vapCon = vapor_concentration;           % fraction (0-1)
% airs_data.z = z;                     % m

end
