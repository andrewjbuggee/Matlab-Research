%% Analyze composite soundings from VOCALS-REx


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/allSoundings_composite_5mb/'];

    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

elseif strcmp(which_computer,'andrewbuggee')==true


    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/radiosonde/allSoundings_composite_5mb/'];


    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


end



% ---------------------
% define the ensemble filename
profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat';



%% Load the ensemble profiles

load([vocalsRexFolder, profiles]);


%% What do I need to preform this analysis?

% Ideally I would have a radiosonde measurement close to each vertical
% droplet profile, thus I could estimate the water vapor concentration
% profile for each droplet profile measured and therefore estimate
% statistics of the above-cloud precipitable water amount.

% To do this, let's step through each ensemble profile, find the radiosonde
% closest to the sampled cloud, and store the radiosonde sampled vertical
% profiles

%% extract the filenames of each sounding composite file

sounding_files = dir(folderpath);

% extract the names
% sounding_file_names = {sounding_files.name};

%% Step through each profile and find the radiosonde closest in space and time
% The radiosonde should be over the ocean, rather than over land


% the composite radiosonde data is saved with the day it was recorded
% on. So first find the radiosonde file for the day the current
% ensemble profile was recorded


for nn = 1:length(profiles)

    % extract the date of the nth profile
    date_profile = ensemble_profiles{nn}.dateOfFlight;

    % Find the sounding date for the same day
    for mm = 1:length(sounding_files)

        % step through each filename
        sounding_file_name = sounding_files(mm).name;

        % check that the filename is longer than 2 characters
        if length(sounding_file_name)>2 && strcmp(sounding_file_name(1), 'V')==true

            % determine the date of the sounding
            date_sounding = datetime(str2double(sounding_file_name(17:20)), str2double(sounding_file_name(21:22)),...
                str2double(sounding_file_name(23:24)));

            % check to see if this is the same day as the profile
            date_difference = date_profile-date_sounding;

            % if the difference is 0, break the loop.
            if date_difference==0

                idx_sounding = mm;

                break

            end


        end


    end


    % Load the radiosonde data
    radiosonde = read_cls_file([folderpath, sounding_files(mm).name]);

    % % each radiosonde cls file contains multiple radiosondes each day. Find
    % % which profile is closest in space to the vertical profile
    % % store the MODIS latitude and longitude
    % radiosonde_lat = [radiosonde.latitude];
    % radiosonde_long = [radiosonde.longitude];
    % 
    % 
    % % we will be computing the arclength between points on an ellipsoid
    % % Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with units of meters.
    % wgs84 = wgs84Ellipsoid("m");
    % 
    % 
    % % Step through each vertical profile and find the MODIS pixel that overlaps
    % % with the mid point of the in-situ sample
    % 
    % dist_btwn_C130_radiosonde = distance(radiosonde_lat, radiosonde_long, ensemble_profiles{nn}.latitude(round(end/2)),...
    %     ensemble_profiles{nn}.longitude(round(end/2)), wgs84);                                  % m - minimum distance
    % 
    % [radiosonde_minDist, index_minDist] = min(dist_btwn_C130_radiosonde, [], 'all');            % m - minimum distance


 


    % % Let's weight time and space equally by adding the spatial and
    % % temporal differences in quadrature and finding the minium value
    % % NOPE, you can't do that. Different units. Think of another way to
    % % weight these two equally when making a decision to pick a single
    % % radiosonde.
    % % Idea: instead of taking the difference, let's normalize the time of
    % % release for each baloon with the time the profile was sampled. 
    % % Then, use these normalized times to weight the distance values. Hmm,
    % % but even for a large distance, it will have a value of 0 if weighted
    % % with the exact time of the profile. 
    % [test, idx_sort_time] = sort(abs(time_diff_C130_radiosonde), 2, 'ascend');
    % [test2, idx_sort_dist] = sort(abs(dist_btwn_C130_radiosonde), 2, 'ascend');
    % % create vectors of the same length as the time_diff and dist vectors
    % % that use integers between 1 and length(vector), where the entry is
    % % the order sorted integer where 1 is the smallest value 
    % time_diff_sort = zeros(1, length(time_diff_C130_radiosonde));
    % dist_btwn_sort = zeros(1, length(dist_btwn_C130_radiosonde));
    % for ii = 1:length(time_diff_C130_radiosonde)
    % 
    %     time_diff_sort(ii) = find(idx_sort_time==ii);
    %     dist_btwn_sort(ii) = find(idx_sort_dist==ii);
    % 
    % end
    % 
    % % multiply these two together and find the smallest value. 
    % [~, idx_min] = min(time_diff_sort.*dist_btwn_sort);


    %% Find the radiosonde closest in time to the profile

    % compute the time between the radiosonde launch and the C130 measured
    % profile
    time_diff_C130_radiosonde = datetime(ensemble_profiles{nn}.dateOfFlight.Year, ensemble_profiles{nn}.dateOfFlight.Month,...
        ensemble_profiles{nn}.dateOfFlight.Day, floor(ensemble_profiles{nn}.time_utc(round(end/2))),...
        floor(60*(ensemble_profiles{nn}.time_utc(round(end/2)) - floor(ensemble_profiles{nn}.time_utc(round(end/2))))),...
        round(60*(60*(ensemble_profiles{nn}.time_utc(round(end/2)) - floor(ensemble_profiles{nn}.time_utc(round(end/2)))) - ...
        floor(60*(ensemble_profiles{nn}.time_utc(round(end/2)) - floor(ensemble_profiles{nn}.time_utc(round(end/2))))))),...
        'TimeZone','UTC') - [radiosonde.release_time];

    % Find the minimum 
    [radiosonde_minTime, idx_minTime] = min(abs(time_diff_C130_radiosonde));            % m - minimum distance

    % Check that the radiosonde profile selected is over ocean
    while strcmp(radiosonde(idx_minTime).release_site(1:3), 'R/V') ~= true

        % this radiosonde is over land, find the radiosonde closest in time
        % to the vertical profile over ocean
        time_diff_C130_radiosonde(idx_minTime) = duration(24,0,0);                     % set to 24 hours
        % Find the minimum 
        [radiosonde_minTime, idx_minTime] = min(abs(time_diff_C130_radiosonde));            % m - minimum distance

    end
    % check that the radiosonde selected has a cloud

    % Define the upper boundary of the cloud. Radiosonde data includes
    % relative humidity, dew point temperature, environmental temperature
    % and pressure. From these values, how do I define whether or not the
    % balloon is sampling without a cloud or not? 








end

