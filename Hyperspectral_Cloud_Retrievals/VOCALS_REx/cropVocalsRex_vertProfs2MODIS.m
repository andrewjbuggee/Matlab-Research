

function vocalsRex = cropVocalsRex_vertProfs2MODIS(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold, modis, modisInputs)


% ----- Find all vertical profiles within VOCALS-REx data ------
vert_prof = find_verticalProfiles_VOCALS_REx(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold);


min_per_hour = 60;                                                                      % min

modis_timeUTC = modis.time(1) + modis.time(2)/min_per_hour;



% Now combine all time vectors together so we can find which profile is
% closest in time to the MODIS data
time2modis = zeros(1, length(vert_prof.time_utc));
for nn = 1:length(vert_prof.time_utc)
    time2modis(nn) = min(abs(vert_prof.time_utc{nn} - modis_timeUTC));    % times are in hours
end

% Determine which profile is closest to the time MODIS data was recorded
[~, index_vertProf_closest2modis] = min(time2modis);


% Let's keep the vertical profile that is closest to MODIS
clear vocalsRex;

% to get the data we want, we need to convert the structure to a cell array
fields = fieldnames(vert_prof);
vert_prof_cell = struct2cell(vert_prof);

data2keep = cell(1, numel(vert_prof_cell));


% step through each field. If its a cell, only keep the index found above
for ii = 1:length(vert_prof_cell)
    if iscell(vert_prof_cell{ii})==true

        data2keep{ii} = vert_prof_cell{ii}{index_vertProf_closest2modis};

    else

        data2keep{ii} = vert_prof_cell{ii};

    end
end

% Convert back to a structure and keep the data closest in time to MODIS
vocalsRex = cell2struct(data2keep, fields, 2);





% --------------------------------------------------------------------
% ------ Find the MODIS pixels to use for the vertical retrieval -----
% --------------------------------------------------------------------

% store the MODIS latitude and longitude
modis_lat = modis.geo.lat;
modis_long = modis.geo.long;


% we will be computing the arclength between points on an ellipsoid
% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with a length unit of meters.
wgs84 = wgs84Ellipsoid("m");

% Set up an empty array for each vocals-rex data point
modisIndex_minDist = zeros(1, length(vocalsRex.latitude));
modis_minDist = zeros(1, length(vocalsRex.latitude));

% Store Vocals-Rex lat and long
vr_lat = vocalsRex.latitude;
vr_long = vocalsRex.longitude;

% store length of data points for VOCALS-REx
n_data_VR = length(vr_lat);



if modisInputs.flags.useAdvection==false

    % ------------------------ NO ADVECTION -------------------------
    % if advection flag is false, simply find the MODIS pixels closest in
    % space to the vocals-rex in-situ measurements



    % Find the MODIS pixel closest to every VOCALS-Rex location in the
    % horizontal profile
    parfor nn = 1:n_data_VR


        dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, vr_lat(nn), vr_long(nn), wgs84);

        [min_dist, index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');
        % save this index so we know what MODIs pixel to use in comparison
        modisIndex_minDist(nn) = index_minDist;
        modis_minDist(nn) = min_dist;            % meters

    end



else


    % ----------------------- USE ADVECTION -----------------------
    % if advection is true, use the measured windspeed and direction to
    % project where the cloud was at the time of the MODIS overpass


    % ----------------------------- OLD CODE -------------------------------
    % ----------------------------------------------------------------------

    %     % use the mean windspeed
    %     mean_wind_speed = mean(vocalsRex.horz_wind_speed);      % m/s
    %
    %     % use the median wind direction
    %     % This is the direction the wind is coming from, NOT the direction the
    %     % wind is blowing torwards
    %     median_wind_from_direction = median(vocalsRex.horz_wind_direction);         % degrees from north (0 deg)
    %
    %     % compute the direction the wind is blowing towards using modulo
    %     % arithmetic
    %     median_wind_direction = mod(median_wind_from_direction + 180, 360);               % degrees from north (0 deg)
    %
    %     % compute the time in seconds between the MODIS overpass and the
    %     % vocalsRex in-situ measurement
    %     d_time_sec = time2modis(index_vertProf_closest2modis)*3600;           % sec
    %
    %     % compute the horizontal distance travelled by the cloud during this
    %     % time
    %     dist_m = mean_wind_speed*d_time_sec;                % meters travelled

    %     % Now we need to incorporate the distance within the vocals-rex
    %     % position. We assume that the cloud that vocals-rex sampled moves in
    %     % time according to the mean wind speed and median direction. So we
    %     % project the VOCALS location to a new location using this heading.
    %
    %
    %     % First, find whether or not MODIS passed overhead before or after the
    %     % VOCALS-REX made in-situ measurements
    %
    %     % if true, then MODIS passed overhead after VOCALS sampled the cloud
    %     % set this to be a value of 1
    %     position_change = modis_timeUTC>vocalsRex.time_utc;
    %
    %
    %     % if position change is a logical 1, we have to move the VOCALS-REx
    %     % in-situ cloud position forward in time. Use the median_wind_direction
    %     azimuth_angle = zeros(1, n_data_VR);
    %     azimuth_angle(position_change) = median_wind_direction;
    %
    %     % if position change if logical 0, we have to move the VOCALS-REx
    %     % in-situ cloud position backwards in time. Use the
    %     % median_wind_from_direction
    %     azimuth_angle(~position_change) = median_wind_from_direction;
    %
    %
    %     % Use the reckon function to compute the new lat long position.
    %     % Inputs are the original lat/long, the distance travelled and the
    %     % azimuth, which is the angle between the direction of travel and
    %     % true north. Or, as MATLAB puts it, it is the angle between the
    %     % local meridian line and the direction of travel, where the angle
    %     % is swept out in the clockwise direction.
    %     % (0 - due north, 90 - due east, 180 - due west etc.)
    %     [new_lat, new_long] = reckon(vr_lat, vr_long, dist_m, azimuth_angle, wgs84);


    %     % Find the MODIS pixel closest to the vocalsRex data by using the median of
    %     % the VOCALS-REx latitude and longitude
    %
    %     %dist_btwn_MODIS_and_VR = sqrt((double(modis_lat) - new_lat_median).^2 + (double(modis_long) - new_long_median).^2);
    %     dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, median(new_lat), median(new_long), wgs84);
    %
    %     [min_dist, index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');
    %     % save this index so we know what MODIs pixel to use in comparison
    %     vocalsRex.modisIndex_minDist_median = index_minDist;
    %     vocalsRex.modis_minDist_median = min_dist;            % meters
    %
    %
    %
    %
    %
    %
    %     % ---- Using the first Position of VOCALS-REx during sampling ----
    %
    %
    %     % Find the MODIS pixel closest to the vocalsRex data by using the
    %     % position of VOCALS-REx at the first sample of the vert profile
    %
    %     %dist_btwn_MODIS_and_VR = sqrt((double(modis_lat) - new_lat_first).^2 + (double(modis_long) - new_long_first).^2);
    %     dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, new_lat(1), new_long(1), wgs84);
    %
    %     [min_dist, index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');
    %     % save this index so we know what MODIs pixel to use in comparison
    %     vocalsRex.modisIndex_minDist_first = index_minDist;
    %     vocalsRex.modis_minDist_first = min_dist;            % meters
    %
    %
    %
    %
    %
    %     % ---- Using the last Position of VOCALS-REx during sampling ----
    %
    %
    %     % Find the MODIS pixel closest to the vocalsRex data by using the
    %     % position of VOCALS-REx at the first sample of the vert profile
    %
    %     %dist_btwn_MODIS_and_VR = sqrt((double(modis_lat) - new_lat_last).^2 + (double(modis_long) - new_long_last).^2);
    %     dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, new_lat(end), new_long(end), wgs84);
    %
    %     [min_dist, index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');
    %     % save this index so we know what MODIs pixel to use in comparison
    %     vocalsRex.modisIndex_minDist_last = index_minDist;
    %     vocalsRex.modis_minDist_last = min_dist;            % meters

    % ----------------------------------------------------------------------
    % ----------------------------------------------------------------------

    % ---- Using the Median Position of VOCALS-REx during sampling ----


    % So I start with some lat and long position. To move along a sphere, I
    % need to invoke spherical trigonometry. But if I incorrectly use circular
    % trig, I can solve this right now.
    %
    % My justification for ignoring spherical trigonometry is that the
    % circumference of the Earth is much much larger than the arc lengths
    % travelled by a cloud over 20 minutes. The circumference of the Earth is
    % 4e7 meters long. The clouds measured by VOCALS-REx travel an average of
    % 7.5 m/s. Over 30 minutes, this would travel a distance of 13.5 km or
    % 1.3e4 meters, roughly 1/1000 th of Earth's circumference


    % ------------------------ OLD CODE --------------------------------
    %     % First let's find the median location of the original vocals rex data
    %     vr_lat_median = median(vocalsRex.latitude);
    %     vr_long_median = median(vocalsRex.longitude);
    %
    %     [~, idx] = min(abs(vocalsRex.latitude - vr_lat_median));
    %
    %     % did MODIS pass overhead before or after vocals sampled the cloud?
    %     if modis_timeUTC<vocalsRex.time_utc(idx)
    %         % then MODIS passed overhead before VOCALS sampled the cloud
    %         position_change = -1;
    %
    %     elseif modis_timeUTC>vocalsRex.time_utc(idx)
    %         % then MODIS passed overhead after VOCALS sampled the cloud
    %         position_change = 1;
    %
    %     else
    %
    %         error([newline, 'MODIS passed overhead during the VOCALS-REx measurement. I dont know what to do.', newline])
    %
    %     end
    %

    %     % nudge the latitude by a small amount proportional to the wind heading
    %     % direction
    %     dLat = 0.00001 * cosd(median_wind_direction);
    %     Lat2Add = dLat;
    %
    %     % nudge the longitude by a small amount proportional to the wind heading
    %     % direction
    %     dLong = 0.00001 * sind(median_wind_direction);
    %     Long2Add = dLong;
    %
    %     % compute the linear distance between two points on an ellipsoid
    %     d = distance(vr_lat_median, vr_long_median, vr_lat_median + position_change*Lat2Add, vr_long_median + position_change*Long2Add, wgs84);
    %
    %     % step d until it is equal to the distance travelled by the cloud
    %     % during the time inbetween the VOCALS-REx sampling time and the MODIS
    %     % overpass
    %     while d<dist_m
    %
    %         Lat2Add = Lat2Add + dLat;
    %         Long2Add = Long2Add + dLong;
    %         d = distance(vr_lat_median, vr_long_median, vr_lat_median + position_change*Lat2Add, vr_long_median + position_change*Long2Add, wgs84);
    %
    %     end
    %
    %     % record the new lat/long position of the cloud
    %     new_lat_median = vr_lat_median + position_change*Lat2Add;
    %     new_long_median = vr_long_median + position_change*Long2Add;
    % ----------------------------------------------------------------






    % project EACH vocalsRex data point using the measured wind speed
    horz_wind_speed = reshape(vocalsRex.horz_wind_speed,[], 1);    % m/s

    % use EACH wind direction
    % This is the direction the wind is coming from, NOT the direction the
    % wind is blowing torwards
    wind_from_direction = vocalsRex.horz_wind_direction;         % degrees from north (0 deg)

    % compute the direction the wind is blowing towards using modulo
    % arithmetic
    wind_direction = mod(wind_from_direction + 180, 360);               % degrees from north (0 deg)

    % compute the time in seconds between the MODIS overpass and the
    % vocalsRex in-situ measurement
    d_time_sec = abs(vocalsRex.time_utc - modis_timeUTC)*3600;           % sec

    % compute the horizontal distance travelled by the cloud during this
    % time
    dist_m = horz_wind_speed.*d_time_sec;                % meters travelled

    % First, find whether or not MODIS passed overhead before or after the
    % VOCALS-REX made in-situ measurements

    % if true, then MODIS passed overhead after VOCALS sampled the cloud
    % set this to be a value of 1
    position_change = modis_timeUTC>vocalsRex.time_utc;


    % if position change is a logical 1, we have to move the VOCALS-REx
    % in-situ cloud position forward in time. Use the median_wind_direction
    azimuth_angle = zeros(n_data_VR,1);
    azimuth_angle(position_change) = wind_direction(position_change);

    % if position change if logical 0, we have to move the VOCALS-REx
    % in-situ cloud position backwards in time. Use the
    % median_wind_from_direction
    azimuth_angle(~position_change) = wind_from_direction(~position_change);


    % Use the reckon function to compute the new lat long position.
    % Inputs are the original lat/long, the distance travelled and the
    % azimuth, which is the angle between the direction of travel and
    % true north. Or, as MATLAB puts it, it is the angle between the
    % local meridian line and the direction of travel, where the angle
    % is swept out in the clockwise direction.
    % (0 - due north, 90 - due east, 180 - due west etc.)

    [lat_withAdvection, long_withAdvection] = reckon(vr_lat', vr_long', dist_m, azimuth_angle, wgs84);


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




% ---- Check to see if all 3 indices are unique. If not, delete the
% redundancy

[~, idx_unique] = unique(vocalsRex.modisIndex_minDist);
idx_unique_logic = ismember(1:length(vocalsRex.modisIndex_minDist), idx_unique);

% delete the indices that are not found above

vocalsRex.modisIndex_minDist = vocalsRex.modisIndex_minDist(idx_unique_logic);







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