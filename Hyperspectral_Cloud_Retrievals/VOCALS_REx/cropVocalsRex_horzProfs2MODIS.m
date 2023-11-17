

function vocalsRex = cropVocalsRex_horzProfs2MODIS(vocalsRex, lwc_threshold, Nc_threshold, max_vertical_displacement, modis, modisInputs)


% ---- Find all Horizontal Profiles ---
horz_prof = find_horizontalProfiles_VOCALS_REx(vocalsRex, lwc_threshold, Nc_threshold, max_vertical_displacement);


min_per_hour = 60;                                                                      % min

modis_timeUTC = modis.time(1) + modis.time(2)/min_per_hour;



% Now combine all time vectors together so we can find which profile is
% closest in time to the MODIS data
time2modis = zeros(1, length(horz_prof.time_utc));
for nn = 1:length(horz_prof.time_utc)
    time2modis(nn) = min(abs(horz_prof.time_utc{nn} - modis_timeUTC));    % times are in hours
end

% Determine which profile is closest to the time MODIS data was recorded
[~, index_vertProf_closest2modis] = min(time2modis);


% Let's keep the vertical profile that is closest to MODIS
clear vocalsRex;

% to get the data we want, we need to convert the structure to a cell array
fields = fieldnames(horz_prof);
horz_prof_cell = struct2cell(horz_prof);

data2keep = cell(1, numel(horz_prof_cell));


% step through each field. If its a cell, only keep the index found above
for ii = 1:length(horz_prof_cell)
    if iscell(horz_prof_cell{ii})==true

        data2keep{ii} = horz_prof_cell{ii}{index_vertProf_closest2modis};

    else

        data2keep{ii} = horz_prof_cell{ii};

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



    % Now we need to incorporate the distance within the vocals-rex
    % position. We assume that the cloud that vocals-rex sampled moves in
    % time according to the mean wind speed and median direction. So we
    % project the VOCALS location to a new location using this heading.


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


    % ----------------------- OLD CODE --------------------------------



    %     % nudge the latitude by a small amount proportional to the wind heading
    %     % direction
    %     dLat = 0.00001 * cosd(median_wind_direction);
    %
    %     % for a parfor loop, the Lat2Add has to be a preallocated vector
    %     Lat2Add = ones(1, n_data_VR) * dLat;
    %
    %     % nudge the longitude by a small amount proportional to the wind heading
    %     % direction
    %     dLong = 0.00001 * sind(median_wind_direction);
    %
    %     % for a parfor loop, the Long2Add has to be a preallocated vector
    %     Long2Add = ones(1, n_data_VR) * dLong;
    %
    %
    %     % Find the MODIS pixel closest to every VOCALS-Rex location in the
    %     % horizontal profile
    %     for nn = 1:n_data_VR
    %
    %         % First, project the VOCALS-DATA
    %
    %         % compute the linear distance between two points on an ellipsoid
    %         d = distance(vr_lat, vr_long, vr_lat + position_change(nn)*Lat2Add(nn),...
    %                         vr_long + position_change(nn)*Long2Add(nn), wgs84);
    %
    %         % step d until it is equal to the distance travelled by the cloud
    %         % during the time inbetween the VOCALS-REx sampling time and the MODIS
    %         % overpass
    %         while d<dist_m
    %
    %             Lat2Add(nn) = Lat2Add(nn) + dLat;
    %             Long2Add(nn) = Long2Add(nn) + dLong;
    %             d = distance(vr_lat(nn), vr_long(nn), vr_lat(nn) + position_change(nn)*Lat2Add(nn),...
    %                                             vr_long(nn) + position_change(nn)*Long2Add(nn), wgs84);
    %
    %         end
    % ------------------------------------------------------------------------


    
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






end