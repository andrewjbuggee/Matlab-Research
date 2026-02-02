%% Write atmospheric profile to libRadtran .DAT format using VOCALS-REx Radiosonde measurements

% This function writes radiosonde measured atmospheric profiles to a
% libRadtran-compatible .DAT file with pres_prof, temp_prof, and
% relative humidity

% INPUTS:
% (1) temp_prof - radiosonde measured temp_prof profile. Must be in units
% of Kelvin
% (2) press_prof - radiosonde measured pres_prof profile. Must be in units
% of mb
% (3) rh_prof - radiosonde measured relative humidity profile. Must be
% defined as a fraction
% (4) ancillary_data - contains lat, long, date and time of radiosonde
% release
% (5) folder_paths - structure containing the path to save the file
% (6) idx - linear index of the pixel to extract profile (1 to num_pixels)
% (7) num_vars - if 2, use just temp_prof and pres_prof from AIRS. if 3,
%                use temp_prof, pres_prof and relative humidity

% OUTPUTS:
% (1) filename_fullPath - full path to the saved .DAT file

% By Andrew John Buggee

%%

function [filename_fullPath] = write_VR_radiosonde_DAT_for_libRadtran(temp_prof, pres_prof, rh_prof,...
    ancillary_data, folder_paths, VR_profile_idx,num_vars)



% Check if folder_paths structure has the atmosphere folder field
if ~isfield(folder_paths, 'atm_folder_path')
    error([newline, 'folder_paths structure must contain "atm_folder_path" field with path to atmosphere folder', newline])
end

atm_folder_path = folder_paths.atm_folder_path;

% Make sure the path ends with a forward slash
if atm_folder_path(end) ~= '/'
    atm_folder_path = [atm_folder_path, '/'];
end



%%


% Check to see if there are nan values. If so, remove them and do a spline
% fit to replace them

if sum(isnan(temp_prof))>=1

    error ([newline, 'Radiosonde temperature profile has NaN values', newline])

    % remove the nan value
    idx_nan = isnan(temp_prof);
    idx_nan_num = find(idx_nan);

    new_temps = zeros(1, length(idx_nan_num));

    for ii = 1:length(idx_nan_num)

        % Check if we need to extrapolate
        if idx_nan_num(ii)==1 || idx_nan_num(ii)==length(temp_prof)

            % we could linearly extrapolate
            new_temps(ii) = exp(interp1(log(pres_prof(~idx_nan)), log(temp_prof(~idx_nan)),...
                log(pres_prof(idx_nan_num(ii))), 'linear', 'extrap'));

            % Or we could just use the previous value it is is extrapolation
            % temp_prof = [];

        else

            % We can interpolate
            new_temps(ii) = exp(interp1(log(pres_prof(~idx_nan)), log(temp_prof(~idx_nan)),...
                log(pres_prof(idx_nan_num(ii))), 'linear'));

            % sometimes this value is right next to the first or last value
            % and also needs extrapolation
            if isnan(new_temps(ii))==true

                new_temps(ii) = exp(interp1(log(pres_prof(~idx_nan)), log(temp_prof(~idx_nan)),...
                    log(pres_prof(idx_nan_num(ii))), 'linear', 'extrap'));

            end



        end

    end

    % replace NaNs with the new temp_profs
    temp_prof(idx_nan) = new_temps;


end







if all(isnan(rh_prof))

    error([newline, 'All values of radiosonde relative humidity are NaN', newline])

end




% Check to see if there are nan values. If so, remove them and do a spline
% fit to replace them
if sum(isnan(rh_prof))>=1

    error ([newline, 'Radiosonde Relative humidity profile has NaN values', newline])

    % remove the nan value
    idx_nan = isnan(rh_prof);
    idx_nan_num = find(idx_nan);

    new_RH = zeros(1, length(idx_nan_num));

    for ii = 1:length(idx_nan_num)

        % Check if we need to extrapolate
        if idx_nan_num(ii)==1 || idx_nan_num(ii)==length(rh_prof)

            % we could linearly extrapolate
            new_RH(ii) = exp(interp1(log(pressH2O(~idx_nan)), log(rh_prof(~idx_nan)),...
                log(pressH2O(idx_nan_num(ii))), 'linear', 'extrap'));

            % Or we could just use the previous value it is is extrapolation
            % rh_prof = [];

        else

            % We can interpolate
            new_RH(ii) = exp(interp1(log(pressH2O(~idx_nan)), log(rh_prof(~idx_nan)),...
                log(pressH2O(idx_nan_num(ii))), 'linear'));

            % sometimes this value is right next to the first or last value
            % and also needs extrapolation
            if isnan(new_RH(ii))==true

                new_RH(ii) = exp(interp1(log(pressH2O(~idx_nan)), log(rh_prof(~idx_nan)),...
                    log(pressH2O(idx_nan_num(ii))), 'linear', 'extrap'));

            end


        end

    end

    % relative humidity cannot be above 100%
    new_RH(new_RH>100) = 100;

    % replace NaNs with the new temp_profs
    rh_prof(idx_nan) = new_RH;



end



% Ensure relative humidity is between 0 and 100%
rh_prof(rh_prof < 0) = 0;
rh_prof(rh_prof > 100) = 100;

% Convert relative humidity from percent to fraction (0-1 range)
% libRadtran uses fractional relative humidity, not percentage
rh_prof = rh_prof ./ 100;  % Convert from % to fraction

%% Handle missing data (NaN values)

% Check for NaN values in the profile
if any(isnan(temp_prof)) || any(isnan(rh_prof))
    warning(['Profile at pixel index ', num2str(VR_profile_idx), ...
        ' contains NaN values. These will be written to file but may cause issues in libRadtran.'])
end

%% Create filename if not provided

formatSpec = 'dd-MMM-yyyy_HH-mm-ss';
dtstr_formatted = string(ancillary_data.dateTime_release{VR_profile_idx}, formatSpec);

if num_vars==2

    filename = join(['VR_radiosonde_profiles_T-P_releaseDateTime_', dtstr_formatted,...
        '_VR_profileNum-',num2str(VR_profile_idx), '.dat'], "");

elseif num_vars==3

    filename = join(['VR_radiosonde_profiles_T-P-RH_releaseDateTime_', dtstr_formatted,...
        '_VR_profileNum-',num2str(VR_profile_idx), '.dat'], "");
end




% Make sure filename ends with .dat (lowercase to match example)
if ~contains(filename, '.dat') && ~contains(filename, '.DAT')
    filename = [filename, '.dat'];
end

%% Prepare data for writing

% libRadtran radiosonde files must have pres_prof increasing
% (altitude decreasing) from top to bottom of file
% Check if pres_prof is in ascending or descending order
if pres_prof(1) > pres_prof(end)
    % pres_prof is ascending (surface to top), need to flip
    pres_prof_2write = flipud(pres_prof);   % mb
    temp_prof_2write = flipud(temp_prof);   % K
    rh_prof_2write = flipud(rh_prof);       % fraction
else
    % pres_prof is already descending (top to surface)
    pres_prof_2write = pres_prof;
    temp_prof_2write = temp_prof;
    rh_prof_2write = rh_prof;
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
        airs.geo.Latitude(unique_pix_idx(VR_profile_idx)), airs.geo.Longitude(unique_pix_idx(VR_profile_idx)));

    % Write column headers - exactly matching the example format
    fprintf(fileID, '#   p(hPa)  T(K)\n');

    % Write the data - three columns with specific formatting to match example
    % Format: pres_prof with 5 decimal places, temp_prof with 1 decimal, RH in scientific notation
    for ii = 1:length(pres_prof_2write)
        fprintf(fileID, '%12.5f %5.1f\n', pres_prof_2write(ii), temp_prof_2write(ii));
    end

    % Close the file
    fclose(fileID);


elseif num_vars==3


    % Open file for writing
    fileID = fopen(strcat(atm_folder_path, filename), 'w');

    if fileID == -1
        error(strcat(newline, 'Could not open file for writing: ', atm_folder_path, filename, newline))
    end

    % Write header comment - single line format like the example
    fprintf(fileID, '# extracted from VOCALS-REx radiosonde data, Lat=%.4f, Lon=%.4f\n', ...
        ancillary_data.lat_release(VR_profile_idx), ancillary_data.long_release(VR_profile_idx));

    % Write column headers - exactly matching the example format
    fprintf(fileID, '#   p(hPa)  T(K)  h2o(relative humidity)\n');

    % Write the data - three columns with specific formatting to match example
    % Format: pres_prof with 5 decimal places, temp_prof with 1 decimal, RH in scientific notation
    for ii = 1:length(pres_prof_2write)
        fprintf(fileID, '%12.5f %5.1f %e\n', pres_prof_2write(ii), temp_prof_2write(ii), rh_prof_2write(ii));
    end

    % Close the file
    fclose(fileID);



end

%% Define output filename with full path

filename_fullPath = strcat(atm_folder_path, filename);



end