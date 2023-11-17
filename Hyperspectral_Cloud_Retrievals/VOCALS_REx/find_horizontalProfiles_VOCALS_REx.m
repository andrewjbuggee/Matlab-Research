%% Finding VOCALS-REx data where the plane flies horizontally within a cloud layer
% Save these horizontal profiles


% INPUTS:
% -------
%       (1) vocals-rex data set

%       (2) LWC_threshold - (g/m^3) this is a threshold that helps clean
%       the data. If set to 0, the vertical profiles will have values on
%       either end of the true vertical droplet profile that represent
%       measurements outside of what we want. the LWC_threshold can help
%       clean the data by limiting the data to have a certain liquid water
%       content threshold. For example. Painemal and Zuidema (2011) defined
%       the cloud top boundary as the level where the LWC was 0.03 g/m^3.
%       Beyond this, the excluded the data. This value will be used to
%       truncate data before and after the vertical profile.


%       (3) Nc_threshold - (droplets/cm^3) this is a threshold that helps
%       the function find profiles with some tangible physical meaning. At
%       times there are confounding measurements where the LWC is greater
%       than the defined threshold but the total number concentration is
%       less than 1. Typically this coincides with erroneous droplet size
%       estimates. This value will be used to ensure only contiguous data
%       above this threshold will satisfy the profile search.


%       (4) max_vertical_displacement - (meters) the maximum vertical
%       displacement allowed throughout a horizontal profile. If set to x,
%       than profiles where the plane ascended or descended by more than x
%       will be thrown out.



function [horz_profs] = find_horizontalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, Nc_threshold, max_vertical_displacement)

% Define the length of consecutive values found above the liquid water
% content threshold that is required to deem it a vertical profile.
consecutive_length_threshold = 50;




% Store whether or not the 2DC was conforming
horz_profs.flag_2DC_data_is_conforming = vocalsRex.flag_2DC_data_is_conforming;

% Let's also store the time vector in UTC format

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

% IF the CDP data was sampled at 10Hz, make a time vector to use for
% indexing
if length(vocalsRex.re_CDP)>length(vocalsRex.time)
    time_sps10 = linspace(vocalsRex.time(1), vocalsRex.time(end), length(vocalsRex.re_CDP));
end


% ---- Vertical Profile Requirements ----
% dz/dt must be close to 0
% Total Nc has to start at a value below 1
% Total Nc has to end at a value below 1
% ----------------------------------------

% First lets crop the data only to those portions where the plane is
% ascending or descending

dz_dt = diff(vocalsRex.altitude)./diff(vocalsRex.time);            % m/s

% Compute the mean with a sliding window for every n data points
% This will smooth out the data and make the horizontal flight segments
% easier to find. It's important we ignore the horizontal segments and find
% only the vertical profiles

n_window = 20;

dz_dt_mean = movmean(dz_dt,n_window);


% The plane's vertical velocity when flying horizontally is typically less
% than 0.5 m/s, on average

index_horizontal = find(abs(dz_dt_mean)<= 0.5);


% Find consecutive logical true values, which represent stand alone
% profiles where the plane is flying horizontally
index_consec = find(diff(index_horizontal)~=1);

% include a 1 to start
index_consec = [0, index_consec];


% store the total number of droplets in each bin to calculate the LWC for
% each instrument
Nc_per_bin = cell(1, length(index_consec)-1);

% -----------------------------------------------------------------------
% For each break in the data, create a profile
for ii = 1:length(index_consec)-1

    % read in the total number of droplets per unit volume
    horz_profs.Nc{ii} = vocalsRex.total_Nc(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    horz_profs.Nc_CDP{ii} = vocalsRex.total_Nc_CDP(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    horz_profs.Nc_2DC{ii} = vocalsRex.total_Nc_2DC(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';

    % read in the total number of droplets for each size bin, but there is
    % no need to store this variable
    Nc_per_bin{ii} = vocalsRex.Nc(:,index_horizontal(index_consec(ii)+1:(index_consec(ii+1))));

    % read in the number of droplets for each size bin, and divide by the
    % width of the size bin to estimate the droplet distribution
    horz_profs.nr{ii} = vocalsRex.Nc(:,index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))./...
        repmat(diff(double(vocalsRex.drop_radius_bin_edges))', 1, length(index_horizontal(index_consec(ii)+1:(index_consec(ii+1)))));

    horz_profs.time{ii} = vocalsRex.time(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    horz_profs.time_utc{ii} = UTC_starttime + double(horz_profs.time{ii})./3600;
    horz_profs.latitude{ii} = vocalsRex.latitude(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    horz_profs.longitude{ii} = vocalsRex.longitude(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    horz_profs.altitude{ii} = vocalsRex.altitude(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';

    if vocalsRex.flag_2DC_data_is_conforming==true
        horz_profs.re{ii} = vocalsRex.re(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
        horz_profs.re_2DC{ii} = vocalsRex.re_2DC(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    else
        horz_profs.mean_r_2DC{ii} = vocalsRex.mean_r_2DC(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    end

    % check to see if the CDP data was sampled at 1Hz or 10 Hz
    if length(vocalsRex.re_CDP)>length(vocalsRex.time)
        % CDP data sampled at 10Hz
        indices_sps10 = time_sps10>=horz_profs.time{ii}(1) & time_sps10<=horz_profs.time{ii}(end);
        horz_profs.re_CDP{ii} = vocalsRex.re_CDP(indices_sps10)';
        horz_profs.lwc_CDP{ii} = vocalsRex.lwc_CDP(indices_sps10)';
    else
        % CDP data sampled at 1Hz
        horz_profs.re_CDP{ii} = vocalsRex.re_CDP(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
        horz_profs.lwc_CDP{ii} = vocalsRex.lwc_CDP(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    end

    horz_profs.lwc{ii} = vocalsRex.lwc(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    horz_profs.lwc_2DC{ii} = vocalsRex.lwc_2DC(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';

    horz_profs.horz_wind_speed{ii} = vocalsRex.horz_wind_speed(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';
    horz_profs.horz_wind_direction{ii} = vocalsRex.horz_wind_direction(index_horizontal(index_consec(ii)+1:(index_consec(ii+1))))';


    % ----------- These variables are not being used -----------------
    %     horz_profs.SWT{ii} = vocalsRex.SWT(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    %     horz_profs.SWB{ii} = vocalsRex.SWB(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    %     horz_profs.LWT{ii} = vocalsRex.LWT(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    %     horz_profs.LWB{ii} = vocalsRex.LWB(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';

end

% Grab the non-time-stamped data as well
horz_profs.drop_radius_bin_edges = double(vocalsRex.drop_radius_bin_edges);
horz_profs.drop_radius_bin_center = double(vocalsRex.drop_radius_bin_center);
horz_profs.startTime = vocalsRex.startTime;                                                                    % We have to assume that this is in UTC time as well




% -----------------------------------------------------------------------
% Now find which of these profiles start and end below the Nc threshold AND have some
% values above 10^7

index2delete = [];

for ii = 1:length(horz_profs.Nc)

    %[horz_profs.Nc{ii}(1), horz_profs.Nc{ii}(end), any(horz_profs.Nc{ii}>10^7), all(horz_profs.Nc{ii}<1)]

    % Check to see if any of these statements are met. If so, delete the
    % profile
    if any(horz_profs.Nc{ii}>10^7) || all(horz_profs.Nc{ii}<Nc_threshold)

        % if this is true, mark the index for deletion
        index2delete = [index2delete, ii];

    end




end




% Let's delete all cells that met the above conditions
horz_profs.Nc(index2delete) = [];
Nc_per_bin(index2delete) = [];
horz_profs.nr(index2delete) = [];
horz_profs.time(index2delete) = [];
horz_profs.time_utc(index2delete) = [];
horz_profs.lwc(index2delete) = [];
horz_profs.latitude(index2delete) = [];
horz_profs.longitude(index2delete) = [];
horz_profs.altitude(index2delete) = [];

if vocalsRex.flag_2DC_data_is_conforming==true
    horz_profs.re(index2delete) = [];
    horz_profs.re_2DC(index2delete) = [];
else
    horz_profs.mean_r_2DC(index2delete) = [];
end

horz_profs.re_CDP(index2delete) = [];
horz_profs.lwc_CDP(index2delete) = [];
horz_profs.lwc_2DC(index2delete) = [];
horz_profs.Nc_CDP(index2delete) = [];
horz_profs.Nc_2DC(index2delete) = [];
horz_profs.horz_wind_speed(index2delete) = [];
horz_profs.horz_wind_direction(index2delete) = [];

% ----------- These variables are not being used -----------------
% horz_profs.SWT(index2delete) = [];
% horz_profs.SWB(index2delete) = [];
% horz_profs.LWT(index2delete) = [];
% horz_profs.LWB(index2delete) = [];




% -----------------------------------------------------------------------
% Usually at the end of each profile there are several zeros. Let's delete
% all 0 values at the end of the vector ONLY when all remaining values are
% 0

% for ii = 1:length(horz_profs.Nc)
%
%
%
%     % -----------------------------------------------------------------------
%     % find the point where all remaining values are zero and truncate
%     % -----------------------------------------------------------------------
%
%
%
%
% end





% -----------------------------------------------------------------------
% Now sort through the vertical profiles and get rid of data points where
% the LWC and the total Nc is below some threshold
% -----------------------------------------------------------------------


% zero vector incase it is needed below
max_lwc = zeros(1, length(horz_profs.lwc));


% create a zero vector that defines which profiles don't meet requirements
% and need to be delted
%profile_idx_2delete = zeros(1, length(horz_profs.lwc));


for ii = 1:length(horz_profs.lwc)

    % Find the first data point and the last data point where the LWC
    % threshold is exceeded. That is, before the first index, several
    % values should be below the threshold. And several data points
    % after the last index should also be below the threshold

    indexes_above_LWC_Nc_threshold = horz_profs.lwc{ii}>=LWC_threshold & horz_profs.Nc{ii}>=Nc_threshold;


    % find the first 1, the first data point that exceeds the LWC
    % threshold. The first value in the index below is where the
    % profile will start. The last index below will be the end of the
    % profile
    indexes2keep = find(indexes_above_LWC_Nc_threshold);

    % We also need to check every profile to ensure the liquid water
    % content stays above our threshold for the entire profile. It should
    % only drop below our threshold before and after the profile found.
    % the below while loop searches for the index to end on. Everything
    % inbetween the first and last index should be greater than or equal to
    % the LWC threshold
    consecutive_length = zeros(length(indexes2keep), length(indexes2keep));


    % This is getting redundant, but oh well

    if isempty(indexes2keep)==true

        % do nothing, this profile will be thrown out

    else

        for firstIndex = indexes2keep(1):indexes2keep(end)-1
            for lastIndex = firstIndex+1:indexes2keep(end)

                % if all values between the first index and last index are
                % true, this means they are all above the minimum liquid water
                % content and number concentration threshold threshold
                if all(horz_profs.lwc{ii}(firstIndex:lastIndex)>=LWC_threshold & horz_profs.Nc{ii}(firstIndex:lastIndex)>=Nc_threshold)==true

                    % compute the length of the vector
                    consecutive_length(firstIndex, lastIndex) = length(horz_profs.lwc{ii}(firstIndex:lastIndex));

                end

            end
        end

        % Now we find the first and last index that corresponds to the longest
        % consecutive string of values above the LWC and Nc threshold
        [~, max_idx] = max(consecutive_length, [], 'all');

        % the row and column correspond the the first and last index,
        % respectively
        [firstIndex, lastIndex] = ind2sub(size(consecutive_length), max_idx);



        % find the max LWC and the index between the first and last index.
        % The max value wil be used to remove profiles that don't have a max
        % LWC above the minimum threshold.
        [max_lwc(ii), index_max_lwc] = max(horz_profs.lwc{ii}(firstIndex:lastIndex));

        [~, index_absolute_max_lwc] = max(horz_profs.lwc{ii});


    end






    % if none of the values within the vertical profile are above the LWC
    % threshold, delete this profile. But for now, just set it to be a
    % vector of 0's
    if isempty(indexes2keep)==true

        horz_profs.Nc{ii} = 0;
        Nc_per_bin{ii} = 0;
        horz_profs.nr{ii} = 0;
        horz_profs.time{ii} = 0;
        horz_profs.time_utc{ii} = 0;
        horz_profs.lwc{ii} = 0;
        horz_profs.latitude{ii} = 0;
        horz_profs.longitude{ii} = 0;
        horz_profs.altitude{ii} = 0;

        if vocalsRex.flag_2DC_data_is_conforming==true
            horz_profs.re{ii} = 0;
            horz_profs.re_2DC{ii} = 0;
        else
            horz_profs.mean_r_2DC{ii} = 0;
        end


        horz_profs.re_CDP{ii} = 0;
        horz_profs.lwc_CDP{ii} = 0;
        horz_profs.lwc_2DC{ii} = 0;
        horz_profs.Nc_CDP{ii} = 0;
        horz_profs.Nc_2DC{ii} = 0;

        horz_profs.horz_wind_speed{ii} = 0;
        horz_profs.horz_wind_direction{ii} = 0;

        % ----------- These variables are not being used -----------------
        %         horz_profs.SWT{ii} = 0;
        %         horz_profs.SWB{ii} = 0;
        %         horz_profs.LWT{ii} = 0;
        %         horz_profs.LWB{ii} = 0;


        % if the longest consecutive vector of measurements above the
        % defined liquid water content threshold is less than the
        % consecutive length threshold, we won't keep this vertical
        % profile.

    elseif length(firstIndex:lastIndex)<consecutive_length_threshold

        horz_profs.Nc{ii} = 0;
        Nc_per_bin{ii} = 0;
        horz_profs.nr{ii} = 0;
        horz_profs.time{ii} = 0;
        horz_profs.time_utc{ii} = 0;
        horz_profs.lwc{ii} = 0;
        horz_profs.latitude{ii} = 0;
        horz_profs.longitude{ii} = 0;
        horz_profs.altitude{ii} = 0;

        if vocalsRex.flag_2DC_data_is_conforming==true
            horz_profs.re{ii} = 0;
            horz_profs.re_2DC{ii} = 0;
        else
            horz_profs.mean_r_2DC{ii} = 0;
        end

        horz_profs.re_CDP{ii} = 0;
        horz_profs.lwc_CDP{ii} = 0;
        horz_profs.lwc_2DC{ii} = 0;
        horz_profs.Nc_CDP{ii} = 0;
        horz_profs.Nc_2DC{ii} = 0;

        horz_profs.horz_wind_speed{ii} = 0;
        horz_profs.horz_wind_direction{ii} = 0;

        % ----------- These variables are not being used -----------------
        %         horz_profs.SWT{ii} = 0;
        %         horz_profs.SWB{ii} = 0;
        %         horz_profs.LWT{ii} = 0;
        %         horz_profs.LWB{ii} = 0;




    else

        % delete all time-stamped data that meet the above logical statement

        horz_profs.Nc{ii} = horz_profs.Nc{ii}(firstIndex:lastIndex);
        Nc_per_bin{ii} = Nc_per_bin{ii}(:,firstIndex:lastIndex);
        horz_profs.nr{ii} = horz_profs.nr{ii}(:, firstIndex:lastIndex);
        horz_profs.time{ii} = horz_profs.time{ii}(firstIndex:lastIndex);
        horz_profs.time_utc{ii} = horz_profs.time_utc{ii}(firstIndex:lastIndex);
        horz_profs.lwc{ii} = horz_profs.lwc{ii}(firstIndex:lastIndex);
        horz_profs.latitude{ii} = horz_profs.latitude{ii}(firstIndex:lastIndex);
        horz_profs.longitude{ii} = horz_profs.longitude{ii}(firstIndex:lastIndex);
        horz_profs.altitude{ii} = horz_profs.altitude{ii}(firstIndex:lastIndex);

        if vocalsRex.flag_2DC_data_is_conforming==true
            horz_profs.re{ii} = horz_profs.re{ii}(firstIndex:lastIndex);
            horz_profs.re_2DC{ii} = horz_profs.re_2DC{ii}(firstIndex:lastIndex);
        else
            horz_profs.mean_r_2DC{ii} = horz_profs.mean_r_2DC{ii}(firstIndex:lastIndex);
        end

        % check to see if were dealing with SPS1 of SPS10 data
        if length(vocalsRex.re_CDP)>length(vocalsRex.time)
            error([newline, 'I dont know how to handle SPS10 data.', newline])
        else
            horz_profs.re_CDP{ii} = horz_profs.re_CDP{ii}(firstIndex:lastIndex);
            horz_profs.lwc_CDP{ii} = horz_profs.lwc_CDP{ii}(firstIndex:lastIndex);
        end

        horz_profs.lwc_2DC{ii} = horz_profs.lwc_2DC{ii}(firstIndex:lastIndex);
        horz_profs.Nc_CDP{ii} = horz_profs.Nc_CDP{ii}(firstIndex:lastIndex);
        horz_profs.Nc_2DC{ii} = horz_profs.Nc_2DC{ii}(firstIndex:lastIndex);

        horz_profs.horz_wind_speed{ii} = horz_profs.horz_wind_speed{ii}(firstIndex:lastIndex);
        horz_profs.horz_wind_direction{ii} = horz_profs.horz_wind_direction{ii}(firstIndex:lastIndex);

        % ----------- These variables are not being used -----------------
        %         horz_profs.SWT{ii} = horz_profs.SWT{ii}(firstIndex:lastIndex);
        %         horz_profs.SWB{ii} = horz_profs.SWB{ii}(firstIndex:lastIndex);
        %         horz_profs.LWT{ii} = horz_profs.LWT{ii}(firstIndex:lastIndex);
        %         horz_profs.LWB{ii} = horz_profs.LWB{ii}(firstIndex:lastIndex);


    end


    % Flag profiles where the LWC has been set to a single value comprised of a 0
    % These need to be deleted.

    if length(horz_profs.lwc{ii})==1 && horz_profs.lwc{ii}==0

        profile_idx_2delete(ii) = true;


        % Also check to see if the total profile found doesn't travel
        % vertically by more than the defined max_vertical_displacement input

    elseif (max(horz_profs.altitude{ii}) - min(horz_profs.altitude{ii}))>max_vertical_displacement

        profile_idx_2delete(ii) = true;

    else


        profile_idx_2delete(ii) = false;

    end



end






% the is_profile_zero logical array defines which profiles to delete


% if this is true, delete all entries for the index
if sum(profile_idx_2delete)>0

    % delete all time-stamped data that meet the above logical statement

    horz_profs.Nc(profile_idx_2delete) = [];
    Nc_per_bin(profile_idx_2delete) = [];
    horz_profs.nr(profile_idx_2delete) = [];
    horz_profs.time(profile_idx_2delete) = [];
    horz_profs.time_utc(profile_idx_2delete) = [];
    horz_profs.lwc(profile_idx_2delete) = [];
    horz_profs.latitude(profile_idx_2delete) = [];
    horz_profs.longitude(profile_idx_2delete) = [];
    horz_profs.altitude(profile_idx_2delete) = [];

    if vocalsRex.flag_2DC_data_is_conforming==true
        horz_profs.re(profile_idx_2delete) = [];
        horz_profs.re_2DC(profile_idx_2delete) = [];
    else
        horz_profs.mean_r_2DC(profile_idx_2delete) = [];
    end


    horz_profs.re_CDP(profile_idx_2delete) = [];
    horz_profs.lwc_CDP(profile_idx_2delete) = [];
    horz_profs.lwc_2DC(profile_idx_2delete) = [];
    horz_profs.Nc_CDP(profile_idx_2delete) = [];
    horz_profs.Nc_2DC(profile_idx_2delete) = [];

    horz_profs.horz_wind_speed(profile_idx_2delete) = [];
    horz_profs.horz_wind_direction(profile_idx_2delete) = [];

    % ----------- These variables are not being used -----------------
    %     horz_profs.SWT(profile_idx_2delete) = [];
    %     horz_profs.SWB(profile_idx_2delete) = [];
    %     horz_profs.LWT(profile_idx_2delete) = [];
    %     horz_profs.LWB(profile_idx_2delete) = [];

end


% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Store the LWC threshold used
horz_profs.LWC_threshold = LWC_threshold;            % g/m^3

% Store the Nc threshold used
horz_profs.Nc_threshold = Nc_threshold;            % g/m^3

% Save the maximum vertical displacement allowed
horz_profs.max_vert_displacement = max_vertical_displacement;       % meters




% ----------------------------------------------------------------------
% ----------- COMPUTE THE HORIZONTAL DISTANCE TRAVELLED ----------------
% ----------------------------------------------------------------------

% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with a length unit of meters.
wgs84 = wgs84Ellipsoid("m");        % units of meters


for nn = 1:length(horz_profs.Nc)

    % Step through each point to get the linear distance travelled as a
    % vector
    horz_distance_travelled = zeros(1, length(horz_profs.latitude{nn}));

    for xx = 2:length(horz_profs.latitude{nn})

        % Find the linear distance between the start and end of the horizontal profile.
        % When you specify a reference ellipsoid as input to the distance function,
        % the function returns linear distances in the units of the semimajor axis of the ellipsoid.
        horz_distance_travelled(xx) = distance(horz_profs.latitude{nn}(1), horz_profs.longitude{nn}(1),...
            horz_profs.latitude{nn}(xx), horz_profs.longitude{nn}(xx),wgs84);

    end

    horz_profs.horz_dist{nn} = horz_distance_travelled;


end





% ----------------------------------------------------------------------
% ------------------ Compute Liquid Water Path -------------------------
% ----------------------------------------------------------------------



% we want to compute the total LWP, and the LWP for the two instruments
% used to measure droplets.
% The CDP probe measures radii between 0.875 and 25.055 microns.
% The 2DC probe measures radii between 31.25 and 793 microns

if iscell(horz_profs.horz_dist)==true

    num_profiles = length(horz_profs.horz_dist);


    % step through each profile
    for nn = 1:num_profiles


        % LWP is calculated by integrating from cloud bottom to
        % cloud top. We are going to calculated the total LWP along a
        % horizontal profile to determine whether or not precipitation was
        % present

        % Compute the total LWP
        horz_profs.lwp{nn} = trapz(horz_profs.horz_dist{nn}, horz_profs.lwc{nn});            % g/m^2

        % ------ Compute the CDP LWP ---------
        horz_profs.lwp_CDP{nn} = trapz(horz_profs.horz_dist{nn}, horz_profs.lwc_CDP{nn});


        % ------ Compute the 2DC LWP ---------
        horz_profs.lwp_2DC{nn} = trapz(horz_profs.horz_dist{nn}, horz_profs.lwc_2DC{nn});



    end




else

    % if vocalsRex is not a cell, then there is only one profile



    % Compute the total LWP
    horz_profs.lwp = trapz(horz_profs.horz_dist, horz_profs.lwc);            % g/m^2

    % ------ Compute the CDP LWP ---------
    % compute the CDP LWP by integration over the cloud depth
    horz_profs.lwp_CDP = trapz(horz_profs.horz_dist, horz_profs.lwc_CDP);


    % ------ Compute the 2DC LWP ---------
    % compute the 2DC LWP by integration over the cloud depth
    horz_profs.lwp_2DC = trapz(horz_profs.horz_dist, horz_profs.lwc_2DC);




end








end