%% Write AIRS atmospheric profile to libRadtran .DAT format

% This function writes AIRS L2 retrieved atmospheric profiles to a
% libRadtran-compatible .DAT file with pressure, temperature, and
% relative humidity

% Quality control is performed using AIRS quality flags:
%   - airs.temp.prof_std_QC: Temperature quality flag
%   - airs.H2O.RelHum_QC: Relative humidity quality flag
%   Quality flag values: 0 = Highest Quality, 1 = Good Quality, 2 = Do Not Use
%   Pressure levels with QC = 2 are replaced with US standard atmosphere values

% INPUTS:
% (1) AIRS - structure containing AIRS L2 data from readAIRS_L2_data()
%            Note: This should be the reduced AIRS structure after using
%            remove_unwanted_airs_data(), so all spatial data is 1D
% (2) folder_paths - structure containing the path to save the file
%            Must include 'atm_folder_path' field
%            Must include 'libRadtran_data' field (path to libRadtran data folder
%            containing atmmod/ with US standard atmosphere files)
% (3) idx - linear index of the pixel to extract profile (1 to num_pixels)
% (4) filename (optional) - custom filename. If not provided, generates automatic name
% (5) num_vars - if 2, use just temperature and pressure from AIRS. if 3,
%                use temperature, pressure and relative humidity
% (6) overlap_pixels - structure containing overlap pixel information
% (7) us_std_atm_file (optional) - filename of US standard atmosphere file
%                                   Default: 'afglus.dat'

% OUTPUTS:
% (1) filename_fullPath - full path to the saved .DAT file

% By Andrew John Buggee

%%

function [filename_fullPath, airs] = write_AIRS_radiosonde_DAT_with_multiPixels(airs, folder_paths, idx, filename,...
    num_vars, overlap_pixels, us_std_atm_file, print_status_updates)

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

%% Read US Standard Atmosphere for quality control replacement
% This will be used to replace bad quality AIRS measurements

% Check if libRadtran_data path exists
if ~isfield(folder_paths, 'libRadtran_data')
    error([newline, 'folder_paths structure must contain "libRadtran_data" field for US standard atmosphere', newline])
end

% Construct path to US standard atmosphere file
us_std_atm_path = [folder_paths.libRadtran_data, 'atmmod/', us_std_atm_file];

% Read the US standard atmosphere
% Columns: z(km), p(hPa), T(K), air, O3, O2, H2O(molecules/cm^3), CO2, NO2
us_std_atm = read_libRadtran_atm_dat_profiles_ver2(us_std_atm_path);

% Extract relevant columns from US standard atmosphere
us_std_z = us_std_atm(:, 1);        % altitude (km)
us_std_p = us_std_atm(:, 2);        % pressure (hPa)
us_std_T = us_std_atm(:, 3);        % temperature (K)
us_std_H2O = us_std_atm(:, 7);      % water vapor density (molecules/cm^3)

%% Extract the profile data for the specified pixel

% Extract pressure levels (standard pressure for temperature)
% Pressure is in hPa (mb)
pressure = double(airs.pressStd{1});  % (StdPressureLev = 28)

% Extract temperature profile at the specified pixel
% After using remove_unwanted_airs_data, temp.prof_std is (num_pixels x StdPressureLev)
% Temperature is in Kelvin
% ** use the AIRS measurement closest to EMIT **
unique_airs_pix = unique(overlap_pixels.airs.linear_idx);
unique_pix_idx = zeros(1, length(overlap_pixels.airs.linear_idx));
for xx = 1:length(unique_pix_idx)

    unique_pix_idx(xx) = find(unique_airs_pix==overlap_pixels.airs.linear_idx(xx));

end

temperature = airs.temp.prof_std(unique_pix_idx(idx), :)';  % (StdPressureLev x 1)

%% Check temperature quality flags and replace bad values with US standard atmosphere
% Quality flag: 0 = Highest Quality, 1 = Good Quality, 2 = Do Not Use

if isfield(airs.temp, 'prof_std_QC')
    temp_QC = airs.temp.prof_std_QC(unique_pix_idx(idx), :)';  % (StdPressureLev x 1)

    % Find pressure levels with QC = 2 (Do Not Use)
    idx_bad_temp = (temp_QC == 2);

    if any(idx_bad_temp)
        % Interpolate US standard atmosphere temperature to AIRS pressure levels
        % Use log-pressure interpolation for atmospheric profiles
        us_std_T_interp = interp1(log(us_std_p), us_std_T, log(pressure), 'linear', 'extrap');

        % Replace bad quality temperature values with US standard atmosphere
        temperature(idx_bad_temp) = us_std_T_interp(idx_bad_temp);

        if print_status_updates == true
            fprintf('Replaced %d temperature values (QC=2) with US standard atmosphere\n', sum(idx_bad_temp));
        end
    end
end

% if all values are nan, throw an error
if all(isnan(temperature))

    error([newline, 'All values of AIRS temperature are NaN', newline])

end

% Check to see if there are nan values. If so, remove them and do a spline
% fit to replace them
if sum(isnan(temperature))>=1

    % remove the nan value
    idx_nan = isnan(temperature);
    idx_nan_num = find(idx_nan);

    new_temps = zeros(1, length(idx_nan_num));

    for ii = 1:length(idx_nan_num)

        % Check if we need to extrapolate
        if idx_nan_num(ii)==1 || idx_nan_num(ii)==length(temperature)

            % we could linearly extrapolate
            new_temps(ii) = exp(interp1(log(pressure(~idx_nan)), log(temperature(~idx_nan)),...
                log(pressure(idx_nan_num(ii))), 'linear', 'extrap'));

            % Or we could just use the previous value it is is extrapolation
            % temperature = [];

        else

            % We can interpolate
            new_temps(ii) = exp(interp1(log(pressure(~idx_nan)), log(temperature(~idx_nan)),...
                log(pressure(idx_nan_num(ii))), 'linear'));

            % sometimes this value is right next to the first or last value
            % and also needs extrapolation
            if isnan(new_temps(ii))==true

                new_temps(ii) = exp(interp1(log(pressure(~idx_nan)), log(temperature(~idx_nan)),...
                    log(pressure(idx_nan_num(ii))), 'linear', 'extrap'));

            end



        end

    end

    % replace NaNs with the new temperatures
    temperature(idx_nan) = new_temps;


end



% Extract relative humidity profile
% After using remove_unwanted_airs_data, H2O.RelHum is (num_pixels x H2OPressureLev)
% RelHum is on H2O pressure levels (15 levels), need to map to pressStd (28 levels)
relHum_H2O_levels = airs.H2O.RelHum(unique_pix_idx(idx), :)';  % (H2OPressureLev x 1)
pressH2O = double(airs.pressH2O{1});  % (H2OPressureLev = 15)

%% Check relative humidity quality flags and replace bad values with US standard atmosphere
% Quality flag: 0 = Highest Quality, 1 = Good Quality, 2 = Do Not Use

if isfield(airs.H2O, 'RelHum_QC')
    relHum_QC = airs.H2O.RelHum_QC(unique_pix_idx(idx), :)';  % (H2OPressureLev x 1)

    % Find pressure levels with QC = 2 (Do Not Use)
    idx_bad_RH = (relHum_QC == 2);

    if any(idx_bad_RH)
        % Convert US standard atmosphere water vapor density to relative humidity
        % First, interpolate US std atm T and H2O to the H2O pressure levels
        us_std_T_H2O = interp1(log(us_std_p), us_std_T, log(pressH2O), 'linear', 'extrap');
        us_std_H2O_H2O = interp1(log(us_std_p), us_std_H2O, log(pressH2O), 'linear', 'extrap');

        % Convert water vapor number density (molecules/cm^3) to vapor pressure
        % n = (e * N_A) / (R * T)  =>  e = (n * R * T) / N_A
        % where n is in molecules/m^3, so convert from cm^3 to m^3
        con = physical_constants;
        R_universal = 8.314462;  % J/(mol·K)
        us_std_H2O_m3 = us_std_H2O_H2O * 1e6;  % molecules/m^3
        e_Pa = (us_std_H2O_m3 .* R_universal .* us_std_T_H2O) ./ con.N_A;  % Pa
        e_hPa = e_Pa / 100;  % hPa

        % Compute saturation vapor pressure (Bolton 1980)
        T_celsius = us_std_T_H2O - 273.15;
        e_sat = 6.112 * exp((17.67 * T_celsius) ./ (T_celsius + 243.5));  % hPa

        % Compute relative humidity from US standard atmosphere (as percentage)
        us_std_RH_H2O = (e_hPa ./ e_sat) * 100;  % percentage

        % Ensure RH is between 0 and 100%
        us_std_RH_H2O(us_std_RH_H2O < 0) = 0;
        us_std_RH_H2O(us_std_RH_H2O > 100) = 100;

        % Replace bad quality RH values with US standard atmosphere
        relHum_H2O_levels(idx_bad_RH) = us_std_RH_H2O(idx_bad_RH);

        if print_status_updates == true
            fprintf('Replaced %d relative humidity values (QC=2) with US standard atmosphere\n', sum(idx_bad_RH));
        end
    end
end

if all(isnan(relHum_H2O_levels))

    error([newline, 'All values of AIRS relative humidity are NaN', newline])

end

% ----- Track which H2O pressure levels originally had NaN values -----
% These are typically near-surface levels below cloud where AIRS cannot
% retrieve humidity. We'll use this later to identify which levels on the
% 28-level grid were filled by extrapolation, and scale them to match
% the AIRS total column water vapor.
idx_nan_H2O_original = (relHum_QC == 2);

% Determine the valid pressure range of the original (non-NaN) RH data
% on the 15-level H2O grid. Any 28-level pressures outside this range
% were extrapolated rather than interpolated.
if any(idx_nan_H2O_original)
    valid_pressH2O_range = [min(pressH2O(~idx_nan_H2O_original)), max(pressH2O(~idx_nan_H2O_original))];
else
    valid_pressH2O_range = [min(pressH2O), max(pressH2O)];
end

% Check to see if there are nan values. If so, remove them and do a spline
% fit to replace them
if sum(isnan(relHum_H2O_levels))>=1

    % remove the nan value
    idx_nan = isnan(relHum_H2O_levels);
    idx_nan_num = find(idx_nan);

    new_RH = zeros(1, length(idx_nan_num));

    for ii = 1:length(idx_nan_num)

        % Check if we need to extrapolate
        if idx_nan_num(ii)==1 || idx_nan_num(ii)==length(relHum_H2O_levels)

            % we could linearly extrapolate
            new_RH(ii) = exp(interp1(log(pressH2O(~idx_nan)), log(relHum_H2O_levels(~idx_nan)),...
                log(pressH2O(idx_nan_num(ii))), 'linear', 'extrap'));

            % Or we could just use the previous value it is is extrapolation
            % relHum_H2O_levels = [];

        else

            % We can interpolate
            new_RH(ii) = exp(interp1(log(pressH2O(~idx_nan)), log(relHum_H2O_levels(~idx_nan)),...
                log(pressH2O(idx_nan_num(ii))), 'linear'));

            % sometimes this value is right next to the first or last value
            % and also needs extrapolation
            if isnan(new_RH(ii))==true

                new_RH(ii) = exp(interp1(log(pressH2O(~idx_nan)), log(relHum_H2O_levels(~idx_nan)),...
                    log(pressH2O(idx_nan_num(ii))), 'linear', 'extrap'));

            end


        end

    end

    % relative humidity cannot be above 100%
    new_RH(new_RH>100) = 100;

    % replace NaNs with the new temperatures
    relHum_H2O_levels(idx_nan) = new_RH;



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

% ----- Identify which of the 28 standard pressure levels were extrapolated -----
% These are levels where pressure falls outside the valid range of the
% original (non-NaN) H2O pressure levels. These are the levels we'll scale
% to force the integrated TPW to match airs.H2O.totCol_Std.
idx_extrapolated_28 = (pressure < valid_pressH2O_range(1)) | (pressure > valid_pressH2O_range(2));

if print_status_updates == true
    fprintf('  %d of %d standard pressure levels were extrapolated from original RH data\n', ...
        sum(idx_extrapolated_28), length(pressure));
end


%% Convert RH profile to water vapor mass density and compute altitude
% These are needed for TPW verification, scaling, and writing the output file

% ----- Convert RH profile to water vapor mass density (kg/m^3) -----
% Saturation vapor pressure (Bolton 1980)
T_celsius_28 = temperature - 273.15;
e_sat_28 = 6.112 * exp((17.67 * T_celsius_28) ./ (T_celsius_28 + 243.5));  % hPa

% Actual vapor pressure from relative humidity (relHum is now a fraction)
e_hPa_28 = relHum' .* e_sat_28;  % hPa

% Water vapor mass density using ideal gas law for water vapor
% rho_v = e / (R_v * T), where e is in Pa and R_v = 461.52 J/(kg*K)
R_v = 461.52;  % Specific gas constant for water vapor [J/(kg*K)]
rho_v = (e_hPa_28 * 100) ./ (R_v .* temperature);  % kg/m^3

% Also compute specific humidity for the hypsometric method
epsilon = 0.622;  % R_d / R_v
q_28 = epsilon * e_hPa_28 ./ (pressure' - (1 - epsilon) * e_hPa_28);  % kg/kg


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

% Map rho_v to the sorted (surface-to-TOA) ordering for integration
rho_v_sorted = rho_v(sort_idx_hyps);

% Map the extrapolated index to the sorted ordering
idx_extrapolated_sorted = idx_extrapolated_28(sort_idx_hyps);

% Compute initial (unscaled) TPW
valid_hyps = ~isnan(z_hyps_sorted) & ~isnan(rho_v_sorted);
TPW_unscaled = trapz(z_hyps_sorted(valid_hyps), rho_v_sorted(valid_hyps));  % kg/m^2


%% Scale rho_v at extrapolated levels to match AIRS total column water vapor
% The extrapolated levels (near-surface, below cloud) have uncertain RH
% values. We scale rho_v only at these levels so the total integrated
% column water vapor matches airs.H2O.totCol_Std.

TPW_airs = airs.H2O.totCol_Std(unique_pix_idx(idx));  % kg/m^2

if print_status_updates == true
    fprintf('\n--- Total Column Water Vapor Verification (before scaling) ---\n');
    fprintf('  AIRS retrieval (totCol_Std):              %.4f kg/m^2\n', TPW_airs);
    fprintf('  Integrated (hypsometric equation):        %.4f kg/m^2  (%.1f%% difference)\n', ...
        TPW_unscaled, 100*(TPW_unscaled - TPW_airs)/TPW_airs);
end

if any(idx_extrapolated_sorted) && TPW_airs > 0

    % Compute partial integrals: contribution from non-extrapolated and
    % extrapolated levels separately
    % We need to integrate rho_v over the sub-regions defined by the
    % extrapolated levels. The scaling factor 'a' satisfies:
    %   TPW_airs = TPW_non_extrap + a * TPW_extrap

    % Create the scaled rho_v profile
    rho_v_scaled_sorted = rho_v_sorted;

    % Compute the contribution from non-extrapolated levels
    rho_v_non_extrap = rho_v_sorted;
    rho_v_non_extrap(idx_extrapolated_sorted) = 0;
    TPW_non_extrap = trapz(z_hyps_sorted(valid_hyps), rho_v_non_extrap(valid_hyps));

    % Compute the contribution from extrapolated levels
    rho_v_extrap = rho_v_sorted;
    rho_v_extrap(~idx_extrapolated_sorted) = 0;
    TPW_extrap = trapz(z_hyps_sorted(valid_hyps), rho_v_extrap(valid_hyps));

    if TPW_extrap > 0
        % Solve for the scalar: TPW_airs = TPW_non_extrap + a * TPW_extrap
        a_scale = (TPW_airs - TPW_non_extrap) / TPW_extrap;

        % Apply scaling only to the extrapolated levels
        rho_v_scaled_sorted(idx_extrapolated_sorted) = a_scale .* rho_v_sorted(idx_extrapolated_sorted);

        if print_status_updates == true
            fprintf('  Scaling factor for extrapolated levels:   %.4f\n', a_scale);
        end
    else
        if print_status_updates == true
            warning('Extrapolated levels contribute zero TPW. Cannot scale to match AIRS total.');
        end
    end

    % Verify the scaled TPW
    TPW_scaled = trapz(z_hyps_sorted(valid_hyps), rho_v_scaled_sorted(valid_hyps));

    if print_status_updates == true
        fprintf('  Integrated (after scaling):               %.4f kg/m^2  (%.1f%% difference)\n', ...
            TPW_scaled, 100*(TPW_scaled - TPW_airs)/TPW_airs);
    end

else

    % No extrapolated levels or invalid AIRS total — use unscaled profile
    rho_v_scaled_sorted = rho_v_sorted;
    if print_status_updates == true
        fprintf('  No extrapolated levels to scale.\n');
    end

end

if print_status_updates == true
    fprintf('--------------------------------------------------------------\n\n');
end


%% Convert scaled rho_v to water vapor number density (molecules/cm^3)
% n = rho_v * N_A / M_w
% where M_w = molar mass of water vapor (kg/mol), N_A = Avogadro's number
con = physical_constants;
% rho_v is in kg/m^3, convert to number density in molecules/cm^3
% n (molecules/m^3) = rho_v * N_A / M_w
% n (molecules/cm^3) = n (molecules/m^3) / 1e6
waterVapor_numDensity_sorted = (rho_v_scaled_sorted .* con.N_A ./ con.Mol_mass_h2o_vap) ./ 1e6;  % molecules/cm^3

% Convert altitude from meters to km for the output file
% z_2write_km = z_hyps_sorted ./ 1e3;  % km


%% Compute relative humidity for reference (not written to file)
% Back-convert the scaled rho_v to RH for diagnostic purposes

% Scaled vapor pressure
% e_Pa_scaled = rho_v_scaled_sorted .* R_v .* T_v_sorted;  % Pa (using T at sorted levels)
% Actually use temperature, not virtual temperature, for the RH calculation
temperature_sorted = temperature(sort_idx_hyps);
e_Pa_scaled_ref = rho_v_scaled_sorted .* R_v .* temperature_sorted;  % Pa
e_hPa_scaled_ref = e_Pa_scaled_ref / 100;  % hPa

% Saturation vapor pressure at each level
T_celsius_sorted = temperature_sorted - 273.15;
e_sat_sorted = 6.112 * exp((17.67 * T_celsius_sorted) ./ (T_celsius_sorted + 243.5));  % hPa

% Relative humidity (fraction)
relHum_scaled_ref = e_hPa_scaled_ref ./ e_sat_sorted;  % fraction (0-1)
relHum_scaled_ref(relHum_scaled_ref < 0) = 0;
relHum_scaled_ref(relHum_scaled_ref > 1) = 1;


%% Handle missing data (NaN values)

% Check for NaN values in the profile
if any(isnan(temperature)) || any(isnan(waterVapor_numDensity_sorted))
    warning(['Profile at pixel index ', num2str(idx), ...
        ' contains NaN values. These will be written to file but may cause issues in libRadtran.'])
end


%% Create filenames

lat = airs.geo.Latitude(unique_pix_idx(idx));
lon = airs.geo.Longitude(unique_pix_idx(idx));

if isempty(filename)
    % Generate automatic filename based on location and time
    if isfield(airs, 'metadata') && isfield(airs.metadata, 'start_year')
        year = airs.metadata.start_year{1};
        month = airs.metadata.start_month{1};
        day = airs.metadata.start_day{1};

        if num_vars == 3

            filename_radiosonde = sprintf('AIRS_profile_T-P-WV_lat%.2f_lon%.2f_%04d-%02d-%02d.dat', lat, lon, year, month, day);

        elseif num_vars == 2
            filename_radiosonde = sprintf('AIRS_profile_T-P_lat%.2f_lon%.2f_%04d-%02d-%02d.dat', lat, lon, year, month, day);
        end

    else

        if num_vars ==3 
            filename_radiosonde = sprintf('AIRS_profile_T-P-WV_lat%.2f_lon%.2f.dat', lat, lon);
        elseif num_vars == 2
            filename_radiosonde = sprintf('AIRS_profile_T-P_lat%.2f_lon%.2f.dat', lat, lon);
        end


    end
else
    % User provided a filename — use it for T-P, derive H2O filename
    filename_radiosonde = filename;
    % Make sure it ends with .dat
    if ~contains(filename_radiosonde, '.dat') && ~contains(filename_radiosonde, '.DAT')
        filename_radiosonde = [filename_radiosonde, '.dat'];
    end
    % Create H2O filename by replacing T-P-RH or T-P with H2O
    filename_H2O = strrep(filename_radiosonde, 'T-P-RH', 'H2O');
    filename_H2O = strrep(filename_H2O, 'T-P', 'H2O');
    % If no replacement was made, prepend H2O_
    if strcmp(filename_H2O, filename_radiosonde)
        [~, name_only, ext] = fileparts(filename_radiosonde);
        filename_H2O = [name_only, '_H2O', ext];
    end
end


%% Write the radiosonde file (T and P only, 2 variables)
% This defines the temperature and pressure profile for libRadtran

if num_vars == 2

    % libRadtran radiosonde files must have pressure increasing
    % (altitude decreasing) from top to bottom of file
    if pressure(1) > pressure(end)
        pressure_2write = fliplr(pressure);
        temperature_2write = fliplr(temperature');
    else
        pressure_2write = pressure;
        temperature_2write = temperature';
    end

    fileID = fopen([atm_folder_path, filename_radiosonde], 'w');

    if fileID == -1
        error([newline, 'Could not open file for writing: ', atm_folder_path, filename_radiosonde, newline])
    end

    fprintf(fileID, '# extracted from AIRS L2 data, Lat=%.4f, Lon=%.4f\n', lat, lon);
    fprintf(fileID, '#   p(hPa)  T(K)\n');

    for ii = 1:length(pressure_2write)
        fprintf(fileID, '%12.5f %5.1f\n', pressure_2write(ii), temperature_2write(ii));
    end

    fclose(fileID);


elseif num_vars == 3



     % libRadtran radiosonde files must have pressure increasing
    % (altitude decreasing) from top to bottom of file
    if pressure(1) > pressure(end)

        pressure_2write = fliplr(pressure);
        temperature_2write = fliplr(temperature');
        vapor_density_2write = fliplr(waterVapor_numDensity_sorted');

    else
        pressure_2write = pressure;
        temperature_2write = temperature';
        vapor_density_2write = waterVapor_numDensity_sorted';
    end

    fileID = fopen([atm_folder_path, filename_radiosonde], 'w');

    if fileID == -1
        error([newline, 'Could not open file for writing: ', atm_folder_path, filename_radiosonde, newline])
    end

    fprintf(fileID, '# extracted from AIRS L2 data, Lat=%.4f, Lon=%.4f\n', lat, lon);
    fprintf(fileID, '# AIRS total column water vapor = %f kg/m2 (mm)\n', TPW_airs);
    fprintf(fileID, '#   p(hPa)  T(K)  H2O(cm-3)\n');

    for ii = 1:length(pressure_2write)
        fprintf(fileID, '%f %f %e\n', pressure_2write(ii), temperature_2write(ii), vapor_density_2write(ii));
    end

    fclose(fileID);



end

%% Write the water vapor density file (z and molecules/cm^3)
% This defines the water vapor profile for libRadtran via mol_file H2O
% z_2write_km and waterVapor_numDensity_sorted are in surface-to-TOA order
% (ascending altitude). libRadtran mol_file expects highest altitude first.

% fileID_h2o = fopen([atm_folder_path, filename_H2O], 'w');
% 
% if fileID_h2o == -1
%     error([newline, 'Could not open file for writing: ', atm_folder_path, filename_H2O, newline])
% end
% 
% fprintf(fileID_h2o, '# Water vapor density from AIRS L2 data, Lat=%.4f, Lon=%.4f\n', lat, lon);
% fprintf(fileID_h2o, '# Extrapolated levels scaled to match AIRS totCol_Std = %.4f kg/m^2\n', TPW_airs);
% fprintf(fileID_h2o, '# z(km)  H2O(molecules/cm^3)\n');
% 
% % Write from highest altitude to lowest (flip the surface-to-TOA ordering)
% for ii = length(z_2write_km):-1:1
%     fprintf(fileID_h2o, '%f %e\n', z_2write_km(ii), waterVapor_numDensity_sorted(ii));
% end
% 
% fclose(fileID_h2o);


%% Define output filenames with full paths

filename_fullPath = [atm_folder_path, filename_radiosonde];
% filename_fullPath.waterVapor_H2O = [atm_folder_path, filename_H2O];

% store the new and scaled profiles in airs
if num_vars == 3
    airs.datProfiles(idx).GP_new = z_hyps_sorted;    % meters - geopotential height for new scaled vapor concentration
    airs.datProfiles(idx).Tv = T_v_sorted;           % K - virtual temperature
    airs.datProfiles(idx).T = temperature_sorted;           % K - temperature
    airs.datProfiles(idx).p = pressure_sorted;       % mb - 28 pressure levels
    airs.datProfiles(idx).vapor_concentration = waterVapor_numDensity_sorted;   % #/cm^3 - water vapor number concentration
    airs.datProfiles(idx).vapor_massDensity = rho_v_scaled_sorted;             % kg/m^3 - water vapor mass density
end



% Display success message
if print_status_updates == true
    disp(newline)
    fprintf('AIRS T-P radiosonde file written to: %s\n', filename_fullPath);
    % fprintf('AIRS H2O density file written to:    %s\n', filename_fullPath.waterVapor_H2O);
    fprintf('Pressure range: %.2f - %.2f hPa\n', min(pressure_2write), max(pressure_2write));
    fprintf('Temperature range: %.2f - %.2f K\n', min(temperature_2write), max(temperature_2write));
    fprintf('Water vapor density range: %.4e - %.4e molecules/cm^3\n', ...
        min(waterVapor_numDensity_sorted), max(waterVapor_numDensity_sorted));
    fprintf('Relative humidity range (reference): %.4f - %.4f (fraction)\n', ...
        min(relHum_scaled_ref), max(relHum_scaled_ref));
    disp(newline)
end

end