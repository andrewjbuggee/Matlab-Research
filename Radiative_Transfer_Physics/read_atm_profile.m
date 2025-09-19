%% This function will read an atm .txt profile

% INPUTS:
%   (1) altitude - vector containing the altitude range of interest
%   (km) - user will provide an altitude range in a vector with two
%   quantities [min, max]. The function will find all values within this
%   range

%   (2) file_name - name of uvspec atm txt file (string) - there are
%   two options to choose from:
%       (a) 'atm_profile_no_cloud.txt' - These data consist of the 1976 US
%       standard atmosphere under clear sky conditions. There data file
%       includes temperature, pressure, the number density of air, ozone,
%       O2, H2O, CO2, NO2, and O4

%       (b) 'atm_profile_with_homogenous_cloud.txt' - These data consist of
%       the US 1976 standard atmosphere with a homogenous cloud between 2
%       and 3 km and a constant droplet radius of 10 microns.

%   (3) 'interp_option' - true or flase flag that, if true, tells the code
%   to interpolate the US standard atmosphere at the values defined by the
%   altitude input. If true, altitude values can be a vector, or a single
%   value. If false, the altitude input has to be the boundaries of the
%   range of interest.


% OUTPUTS:
%   (1) altitude (km) - a vector containing the altitude levels
%   corresponding to each layer

%   (2) pressure (hPa) - a vector containing pressure at each level

%   (3) temperature (K) - a vector containing temperature at each level

%   (4) air_Nc (cm^(-3)) - a vector containing the number density of air at
%   each level

%   (5) O3_Nc (cm^(-3)) - a vector containing the number density of ozone at
%   each level

%   (6) O2_Nc (cm^(-3)) - a vector containing the number density of
%   molecular oxygen at each level

%   (7) H2O_Nc (cm^(-3)) - a vector containing the number density of water
%   vapor at each level

%   (8) CO2_Nc (cm^(-3)) - a vector containing the number density of carbon
%   dioxide at each level

%   (9) NO2_Nc (cm^(-3)) - a vector containing the number density of
%   nitrogen dioxide at each level

%   (10) O4_Nc (cm^(-3)) - a vector containing the number density of
%   tetraoxygen at each level




% NOTES:
%   LibRadTran comes supplied with a folder of atmospheric constituents for
%   the US 1976 standard atmosphere. This text file was created by using
%   the table that is shown in the 'verbose' error message, which describes
%   in detail the atmosphere under simulation

% By Andrew John Buggee
%%

function atm_prof = read_atm_profile(altitude_boundaries, file_name, interp_option)


% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, altitude bounds and a file_name
% to read


if nargin~=3
    error([newline,'Not enough inputs. Need 3: altitude bounds and a filename', newline])
end

% Check to make sure altitude has a length of 2

if interp_option==false && length(altitude_boundaries)~=2

    error([newline,'The altitude input is a vector consisting of: [altitude_min, altitude_max]', newline])

end


% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    atm_profile_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/atm_profiles/';

elseif strcmp(computer_name,'anbu8374')==true

    atm_profile_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/atm_profiles/';

end


%%

% ------------------------------------------------------------
% ---------------------- Read Source File --------------------
% ------------------------------------------------------------



if strcmp(file_name,'atm_profile_no_cloud.txt')

    % Lets check to make sure the altitude input is within bounds of the
    % file selected

    altitude_regime = [0, 120];            % km - altitude boundaries

    if any(altitude_boundaries<altitude_regime(1)) || any(altitude_boundaries>altitude_regime(2))
        error([newline, 'Altitude is out of the range of atm_profile_no_cloud.txt. Must be between [0, 120] km.', newline])
    end


    % ---- Open the File ----

    file_id = fopen([atm_profile_folder,file_name], 'r');   % 'r' tells the function to open the file for reading

    % ------------------------------------------------------
    % -------- Reading .txt file using fscanf --------------
    % ------------------------------------------------------

    % define how the data is written
    % '%' starts a new character
    % '%*s' tells the code to skip string characters
    % 'f' tells us to look for floating point numbers
    % 'e' tells the code to look for exponential numbers
    % '5.1' tells us to look for a floating point number with 5 entities, 1
    % of which comes after the decimal point

    %     format_spec = '  %f %e';
    %     shape_output = [2 Inf];
    %
    %     % lets skip the first 11 lines, since they are commented
    %     for ii = 1:11
    %         fgets(file_id);
    %     end
    %
    %     % now the file pointer will be at the data
    %     A = fscanf(file_id, format_spec, shape_output)'; % sxtract data!

    % ------------------------------------------------------
    % -------- Reading .txt file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore

    format_spec = '%f %f %f %f %f %f %f %f %f %f %f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'Delimiter',' ',...
        'MultipleDelimsAsOne',1, 'CommentStyle','#');

    % read the values between the user defined boundaries if true
    if interp_option==false
        index_altitude = B{2}>=altitude_boundaries(1) & B{1}<=altitude_boundaries(2);

        atm_prof.altitude = B{2}(index_altitude);                % altitudes within the user specified range

        atm_prof.pressure = B{3}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.temperature = B{4}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.air_Nc = B{5}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.O3_Nc = B{6}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.O2_Nc = B{7}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.H2O_Nc = B{8}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.CO2_Nc = B{9}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.NO2_Nc = B{10}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.O4_Nc = B{11}(index_altitude);                % hPa - pressure at the corresponding altitude values

    else

        % If true, interpolate at the point of interest 

        atm_prof.altitude = altitude_boundaries;        % km

        atm_prof.pressure = interp1(B{2}, B{3}, altitude_boundaries, 'linear');         % hPa - pressure at the corresponding altitude values

        atm_prof.temperature = interp1(B{2}, B{4}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.air_Nc = interp1(B{2}, B{5}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.O3_Nc = interp1(B{2}, B{6}, altitude_boundaries, 'linear');               % hPa - pressure at the corresponding altitude values

        atm_prof.O2_Nc = interp1(B{2}, B{7}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.H2O_Nc = interp1(B{2}, B{8}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.CO2_Nc = interp1(B{2}, B{9}, altitude_boundaries, 'linear');               % hPa - pressure at the corresponding altitude values

        atm_prof.NO2_Nc = interp1(B{2}, B{10}, altitude_boundaries, 'linear');             % hPa - pressure at the corresponding altitude values

        atm_prof.O4_Nc = interp1(B{2}, B{11}, altitude_boundaries, 'linear');              % hPa - pressure at the corresponding altitude values

    end

elseif strcmp(file_name,'atm_profile_with_homogenous_cloud.txt')

    % Lets check to make sure the altitude input is within bounds of the
    % file selected

    altitude_regime = [0, 120];            % km - altitude boundaries

    if any(altitude_boundaries<altitude_regime(1)) || any(altitude_boundaries>altitude_regime(2))
        error([newline, 'Altitude is out of the range of atm_profile_no_cloud.txt. Must be between [0, 120] km.', newline])
    end


    % ---- Open the File ----

    file_id = fopen([atm_profile_folder,file_name], 'r');   % 'r' tells the function to open the file for reading


    % ------------------------------------------------------
    % -------- Reading .txt file using textscan ------------
    % ------------------------------------------------------

    format_spec = '%f %f %f %f %f %f %f %f %f %f %f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'Delimiter',' ',...
        'MultipleDelimsAsOne',1, 'CommentStyle','#');

        % read the values between the user defined boundaries if true
    if interp_option==false
        index_altitude = B{2}>=altitude_boundaries(1) & B{1}<=altitude_boundaries(2);

        atm_prof.altitude = B{2}(index_altitude);                % altitudes within the user specified range

        atm_prof.pressure = B{3}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.temperature = B{4}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.air_Nc = B{5}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.O3_Nc = B{6}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.O2_Nc = B{7}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.H2O_Nc = B{8}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.CO2_Nc = B{9}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.NO2_Nc = B{10}(index_altitude);                % hPa - pressure at the corresponding altitude values

        atm_prof.O4_Nc = B{11}(index_altitude);                % hPa - pressure at the corresponding altitude values

    else

        % If true, interpolate at the point of interest 

        atm_prof.altitude = altitude_boundaries;        % km

        atm_prof.pressure = interp1(B{2}, B{3}, altitude_boundaries, 'linear');         % hPa - pressure at the corresponding altitude values

        atm_prof.temperature = interp1(B{2}, B{4}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.air_Nc = interp1(B{2}, B{5}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.O3_Nc = interp1(B{2}, B{6}, altitude_boundaries, 'linear');               % hPa - pressure at the corresponding altitude values

        atm_prof.O2_Nc = interp1(B{2}, B{7}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.H2O_Nc = interp1(B{2}, B{8}, altitude_boundaries, 'linear');                % hPa - pressure at the corresponding altitude values

        atm_prof.CO2_Nc = interp1(B{2}, B{9}, altitude_boundaries, 'linear');               % hPa - pressure at the corresponding altitude values

        atm_prof.NO2_Nc = interp1(B{2}, B{10}, altitude_boundaries, 'linear');             % hPa - pressure at the corresponding altitude values

        atm_prof.O4_Nc = interp1(B{2}, B{11}, altitude_boundaries, 'linear');              % hPa - pressure at the corresponding altitude values

    end



else

    error([newline,'I dont recognize the atm profile you entered!',newline])

end




end
