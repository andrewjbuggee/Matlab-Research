

function vocalsRex = cropVocalsRex_vertProfs2MODIS(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold, modis, modisInputs)


% ----- Find all vertical profiles within VOCALS-REx data ------
vert_prof = find_verticalProfiles_VOCALS_REx(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold);


%% Let's step through each vertical profile and find the MODIS pixel that overlaps


% --------------------------------------------------------------------
% ------ Find the MODIS pixels to use for the vertical retrieval -----
% --------------------------------------------------------------------

% store the MODIS latitude and longitude
modis_lat = modis.geo.lat;
modis_long = modis.geo.long;

% store modis pixel time
modis_pixel_time = modis.EV1km.pixel_time_decimal;


% we will be computing the arclength between points on an ellipsoid
% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with units of meters.
wgs84 = wgs84Ellipsoid("m");

% Set up an empty array for each vocals-rex data point
modis_minDist = zeros(1, length(vert_prof.latitude));
time_diff_MODIS_VR = zeros(1, length(vert_prof.latitude));
within_5min_window = zeros(1, length(vert_prof.time_utc));

% Step through each vertical profile and find the MODIS pixel that overlaps
% with the mid point of the in-situ sample
parfor nn = 1:length(vert_prof.latitude)


    dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, vert_prof.latitude{nn}(round(end/2)),...
        vert_prof.longitude{nn}(round(end/2)), wgs84);

    [modis_minDist(nn), index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');            % m - minimum distance


    % compute the time between the modis pixel closest to the VR sampled
    % path and the time VOCALS was recorded
    time_diff_MODIS_VR(nn) = abs(modis_pixel_time(index_minDist) - ...
        vert_prof.time_utc{nn}(round(end/2))) * 60;                         % minutes

    % A single MODIS granule takes 5 minutes to be collected. First, lets
    % determine if any profiles were recorded within the the 5 minute window
    % when MODIS collected data.

    within_5min_window(nn) = any(vert_prof.time_utc{nn} >= modis_pixel_time(1,1) & ...
        vert_prof.time_utc{nn} <= modis_pixel_time(end,end));

end

%% Check the time difference and the distance between the closest pixel and VR


% First, keep all profiles recorded within the MODIS granule time window
% and have a minimum distance with the closest MODIS pixel of less than 1km
if any(within_5min_window==true & modis_minDist<=1000)
    
    index_vertProfs_2keep = find(within_5min_window==true & modis_minDist<=1000);

% Next, keep profiles with the smallest time difference as long as the
% distance between the closest MODIS pixel and the VR profile is less than
% 1km
elseif any(time_diff_MODIS_VR<5 & modis_minDist<1000)
    
%     % First, check to see if this is the MODIS file from 11 Nov. 2008 recorded at 1850 UTC.
%     % This file is tricky because two vertical profiles recorded right
%     % before and right after the MODIS granule are very different. The
%     % first has an adiabatic droplet profile, the second is nearly
%     % homogeneous.
%     if strcmp(modisInputs.L1B_filename, 'MYD021KM.A2008316.1850.061.2018039033053.hdf')==true
%         % pick the second closest profile
%         [~, idx_min] = min(time_to_halfwayPoint);
%         % get rid of this value
%         time_to_halfwayPoint(idx_min) = inf;
%         % grab the next vertical profile closest in time
%         [~, index_vertProfs_2keep] = min(time_to_halfwayPoint);
%     else
% 
%         [~, index_vertProfs_2keep] = min(time_to_halfwayPoint);
% 
%     end


    % save vertical profiles that were recorded within 5 minutes of the
    % MODIS recording and the distance between VR and the closest pixel is
    % less than 1km
    index_vertProfs_2keep = find(time_diff_MODIS_VR<5 & modis_minDist<1000);

end

% I dont' know what to do if there is more than 1 profile that satisfies
% the above criteria
if length(index_vertProfs_2keep)>1
    error([newline, "Code isn't set up for more than 1 profile.", newline])
end

% Let's keep the vertical profile that is closest to MODIS
clear vocalsRex;

% to get the data we want, we need to convert the structure to a cell array
fields = fieldnames(vert_prof);
vert_prof_cell = struct2cell(vert_prof);

data2keep = cell(1, numel(vert_prof_cell));


% step through each field. If its a cell, only keep the index found above
for ii = 1:length(vert_prof_cell)
    if iscell(vert_prof_cell{ii})==true

        data2keep{ii} = vert_prof_cell{ii}{index_vertProfs_2keep};

    else

        data2keep{ii} = vert_prof_cell{ii};

    end
end

% Convert back to a structure and keep the data closest in time to MODIS
vocalsRex = cell2struct(data2keep, fields, 2);




%%
% --------------------------------------------------------------------
% ------ Find the MODIS pixels to use for the vertical retrieval -----
% --------------------------------------------------------------------


% Set up an empty array for each vocals-rex data point
modisIndex_minDist = zeros(1, length(vocalsRex.latitude));
modis_minDist = zeros(1, length(vocalsRex.latitude));

% Store Vocals-Rex lat and long
vr_lat = vocalsRex.latitude;
vr_long = vocalsRex.longitude;

% store length of data points for VOCALS-REx
n_data_VR = length(vr_lat);


% First, let's find the MODIS pixel closest to ALL VOCALS-REx locations
% throughout the sampled vertical profile
parfor nn = 1:n_data_VR


    dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, vr_lat(nn), vr_long(nn), wgs84);

    [modis_minDist(nn), modisIndex_minDist(nn)] = min(dist_btwn_MODIS_and_VR, [], 'all');       % m - minimum distance between VR and closest MODIS pixel


end

% Did the VOCALS-REx measurement take place before or after the MODIS
% sample? Use the pixel closest to VOCALS-REx
[global_min_dist, min_idx] = min(modis_minDist);
VR_sampled_before_MODIS = vocalsRex.time_utc(min_idx) < modis.EV1km.pixel_time_decimal(modisIndex_minDist(min_idx));


if modisInputs.flags.useAdvection==false

    % ------------------------ NO ADVECTION -------------------------
    % if advection flag is false, simply find the MODIS pixels closest in
    % space to the vocals-rex in-situ measurements




else


    % ----------------------- USE ADVECTION -----------------------
    % if advection is true, use the measured windspeed and direction to
    % project where the cloud was at the time of the MODIS overpass



    % project EACH vocalsRex data point using the measured wind speed and
    % direction
    horz_wind_speed = reshape(vocalsRex.horz_wind_speed,[], 1);    % m/s

    % use EACH wind direction
    % This is the direction the wind is coming from, NOT the direction the
    % wind is blowing torwards
    wind_from_direction = vocalsRex.horz_wind_direction;         % degrees from north (0 deg)

    % compute the direction the wind is blowing towards using modulo
    % arithmetic
    wind_direction = mod(wind_from_direction + 180, 360);               % degrees from north (0 deg)

    % If VOCALS-REx sampled after the MODIS pixel closest to the in-situ
    % measured path was recorded, then in order to line up the same cloudy
    % region for comparison, we need to move the VOCALS-REx flight backwards
    % in time along the direction of where the wind is comming from. That's
    % because the we want to use the cloud that moved into the VOCALS-REx
    % position once VOCALS-REx made it's measurement

    % compute the time in seconds between the MODIS overpass and the
    % vocalsRex in-situ measurement
    d_time_sec = abs(vocalsRex.time_utc(min_idx) - modis.EV1km.pixel_time_decimal(modisIndex_minDist(min_idx)))*3600;           % sec


    % compute the horizontal distance travelled by the cloud during this
    % time
    dist_m = horz_wind_speed.*d_time_sec;                % meters travelled


    % First, find whether or not MODIS passed overhead before or after the
    % VOCALS-REX made in-situ measurements.
    % if true, then MODIS passed overhead before VOCALS sampled the cloud
    % set this to be a value of 1
    modis_before_vocals = modis.EV1km.pixel_time_decimal(modisIndex_minDist(min_idx))<vocalsRex.time_utc;


    % if position change is a logical 1, we have to move the VOCALS-REx
    % in-situ cloud position backward in time into the direction the
    % wind is coming from
    azimuth_angle = zeros(n_data_VR,1);
    azimuth_angle(modis_before_vocals) = wind_from_direction(modis_before_vocals);

    % if position change is logical 0, we have to move the VOCALS-REx
    % in-situ cloud position forwards in time. This time we move the data
    % along the direction the wind is moving into
    azimuth_angle(~modis_before_vocals) = wind_direction(~modis_before_vocals);


    % Use the reckon function to compute the new lat long position.
    % Inputs are the original lat/long, the distance travelled and the
    % azimuth, which is the angle between the direction of travel and
    % true north. Or, as MATLAB puts it, it is the angle between the
    % local meridian line and the direction of travel, where the angle
    % is swept out in the clockwise direction.
    % (0 - due north, 90 - due east, 180 - due west etc.)

    [lat_withAdvection, long_withAdvection] = reckon(vr_lat', vr_long', dist_m, azimuth_angle, wgs84);

    
    % Set up an empty array for each vocals-rex data point
    modisIndex_minDist = zeros(1, length(vocalsRex.latitude));
    modis_minDist = zeros(1, length(vocalsRex.latitude));
    parfor nn = 1:n_data_VR

        dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, lat_withAdvection(nn), long_withAdvection(nn), wgs84);

        [min_dist, index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');
        % save this index so we know what MODIs pixel to use in comparison
        modisIndex_minDist(nn) = index_minDist;
        modis_minDist(nn) = min_dist;            % meters

    end

    % Save the new lat and long
    vocalsRex.lat_withAdvection = lat_withAdvection;
    vocalsRex.long_withAdvection = long_withAdvection;




end




% Store the modis index values and the distance from VOCALS to the pixel
vocalsRex.modisIndex_minDist = modisIndex_minDist;
vocalsRex.modis_minDist = modis_minDist;            % meters




% ---- Check to see if all indices are unique. If not, delete the
% redundancy

[~, idx_unique] = unique(vocalsRex.modisIndex_minDist);
idx_unique_logic = ismember(1:length(vocalsRex.modisIndex_minDist), idx_unique);

% delete the indices that are not found above

vocalsRex.modisIndex_minDist = vocalsRex.modisIndex_minDist(idx_unique_logic);





% ---------------------------------------------------------------------
% % If VOCALS-REx sampled some location after MODIS sampled the same
% % location, we need to move the VOCALS-REx profile in the direction of the
% % wind. Move the profile using the average horizontal wind speed and the
% % time between when MODIS sampled and VOCALS-Rex sampled.
% 
% % Let's test this by moving the VOCALS-REx path by +/- 5 minutes
% dist_m = mean(horz_wind_speed).*(-300:30:300);
% 
% azimuth_angle = zeros(size(dist_m));
% 
% % for moving backwards in time, move along the direction the wind is
% % comming from
% azimuth_angle(dist_m<0) = mean(wind_from_direction);
% 
% % for moving forwards in time, move along the direction the wind is blowing
% azimuth_angle(dist_m>0) = mean(wind_direction);
% 
% % Use the reckon function to compute the new lat long position.
% % Inputs are the original lat/long, the distance travelled and the
% % azimuth, which is the angle between the direction of travel and
% % true north. Or, as MATLAB puts it, it is the angle between the
% % local meridian line and the direction of travel, where the angle
% % is swept out in the clockwise direction.
% % (0 - due north, 90 - due east, 180 - due west etc.)
% 
% lat_withAdvection = zeros(length(vr_lat), length(dist_m));
% long_withAdvection = zeros(length(vr_long), length(dist_m));
% 
% for nn = 1:length(dist_m)
% 
%     if dist_m(nn)==0
%        
%         lat_withAdvection(:,nn) = vr_lat';
%         long_withAdvection(:,nn) = vr_long';
%     
%     else
% 
%         [lat_withAdvection(:,nn), long_withAdvection(:,nn)] = reckon(vr_lat', vr_long', abs(dist_m(nn)), azimuth_angle(nn), wgs84);
% 
%     end
% 
% end
% ---------------------------------------------------------------------



% ---------------------------------------------------------------------
% -------------------------- OLD CODE --------------------------------
% ---------------------------------------------------------------------







% ---------------------------------------------------------------------
% ---- Below is the times associated with a MANUAL profile search -----
% ---------------------------------------------------------------------

% define the time within the vocals-rex profile where a nice cloud profile
% exists

% % convert the MODIS time to decimal time (fraction of 1 day)
% sec_per_hour = 3600;                                                                    % sec
% sec_per_min = 60;                                                                       % sec
% sec_per_day = 86400;                                                                    % sec

% % Convert start time of VOCALS-REx flight to same decimal time (fraction of
% % 1 day)
% startTime_dec = (vocalsRex.startTime(1)*sec_per_hour + vocalsRex.startTime(2)*sec_per_min)/sec_per_day;      % fraction of a day
%
% % Convert the rest of the Vocals-REx data record time into fractions of a
% % single day
% vocalsRex_time = startTime_dec + double(vocalsRex.time)./sec_per_day;



% if strcmp(modisFolder(end-16:end),'/2008_11_11_1430/')==true && strcmp(vocalsRexFolder(end-11:end),'/2008_11_11/')==true
%
%     cloud_profile_time = 0.6069;
%
%     index_delay = 4;
%
% elseif strcmp(modisFolder(end-11:end),'/2008_11_09/')==true && strcmp(vocalsRexFolder(end-11:end),'/2008_11_09/')==true
%
%     cloud_profile_time = 0.6120;
%
%     index_delay = 0;
%
% elseif strcmp(modisFolder(end-16:end),'/2008_11_11_1850/')==true && strcmp(vocalsRexFolder(end-11:end),'/2008_11_11/')==true
%
%     cloud_profile_time = 0.7816;
%
%     index_delay = -8;
%
% else
%
%     error([newline, 'What is the decimal time of the cloud profile you wish to look at?', newline])
%
% end
%
%
% cloudProfile_secSinceStart = floor((cloud_profile_time - startTime_dec)*sec_per_day);                               % seconds since flight started up to the time indicated by painemal and zudema
%
%
% % Define the number of discrete points in time you wish to look at. The
% % center point will be at the start time calculated above
% windowLength = 100;
%
% % Using the time defined above, find the location in space closest to
% % the airplanes location
%
% % Using lat and long with can minimize the euclidean norm
% % ***** There are some special cases *****
%
% % The values are closer on Novemeber 11th 2008 if we compare one time step
% % ahead. The re values are spatially similar, but the optical thickness
% % change quite a bit!
% dist_btwn_PZ_startTime_and_MODIS = sqrt((modis.geo.lat - vocalsRex.latitude(cloudProfile_secSinceStart+index_delay)).^2 + (modis.geo.long - vocalsRex.longitude(cloudProfile_secSinceStart+index_delay)).^2);
% [~, index_minDist] = min(dist_btwn_PZ_startTime_and_MODIS, [], 'all');
%
%
%
%
% % ****** There are two ways we could define Cloud Top *****
% % Cloud top is tricky. The vocals rex data shows droplets being recorded
% % for several data points past the peak LWC value. But MODIS estimates for
% % optical depth better align with the cloud top being defined as the peak
% % LWC value. After this maximum, LWC tends to fall off dramatically
%
% % Cloud top = location where re goes to 0 after LWC> 0.03 g/m^3
% % Cloud top = maximum value of LWC after LWC> 0.03 g/m^3
%
% % First find the index cloud bottom using the definition above
%
% % ---------------------------------------------------------------------
% % Cloud bottom = location where LWC = 0.03 g/m^3 - Painemal and Zuidema
% % (2011) defintion
% % ---------------------------------------------------------------------
%
% lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base
%
% LWC_window_index = cloudProfile_secSinceStart - windowLength/2 : cloudProfile_secSinceStart + windowLength/2;
% LWC_window = vocalsRex.lwc(LWC_window_index);
% index_minVal = find(LWC_window >= lwc_lim);
% index_cloudBottom = LWC_window_index(1) + index_minVal(1) -1;
%
%
%
% % Now lets find the index for cloud top by finding the value of effective
% % radius goes to 0
%
% window_data = vocalsRex.re(index_cloudBottom : index_cloudBottom + 100);
% index_nan = find(isnan(window_data));
%
% index_cloudTop = index_cloudBottom + index_nan(1) - 2;
%
% % Lets also find the index where LWC is a maxium after the index found for
% % cloud bottom
%
% window_data = vocalsRex.lwc(index_cloudBottom : index_cloudBottom + 100);
% [~,index_maxLWC] = max(window_data);
%
% index_cloudTop2 = index_cloudBottom + index_maxLWC - 1;
%
%
%
% vocalsRex.total_Nc = vocalsRex.total_Nc(index_cloudBottom:index_cloudTop2);
% vocalsRex.time = vocalsRex.time(index_cloudBottom:index_cloudTop2);
% vocalsRex.latitude = vocalsRex.latitude(index_cloudBottom:index_cloudTop2);
% vocalsRex.longitude = vocalsRex.longitude(index_cloudBottom:index_cloudTop2);
% vocalsRex.altitude = vocalsRex.altitude(index_cloudBottom:index_cloudTop2);
% vocalsRex.re = vocalsRex.re(index_cloudBottom:index_cloudTop2);
% vocalsRex.lwc = vocalsRex.lwc(index_cloudBottom:index_cloudTop2);
% vocalsRex.startTime = vocalsRex.startTime;                               % We have to assume that this is in UTC time as well
% vocalsRex.modisIndex_minDist = index_minDist;




end