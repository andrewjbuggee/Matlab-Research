%% Finding VOCALS-REx data where the plane flies cleanly through a cloud layer
% Save these vertical profiles


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


%       (3) stop_at_max_lwc - this is a logical flag, whose values can either
%       be true or false. If true, each vertical profile found will be
%       truncated at the peak value of the liquid water content. Most
%       vertical droplet profiles will have a LWC that grows with altitude.
%       If this is true, the profile usually peters out shortly after the
%       peak LWC value. IF this is false, the data will be kept until the
%       the LWC is below the LWC threshold defined above.


%       (4) Nc_threshold - (droplets/cm^3) this is a threshold that helps
%       the function find profiles with some tangible physical meaning. At
%       times there are confounding measurements where the LWC is greater
%       than the defined threshold but the total number concentration is
%       less than 1. Typically this coincides with erroneous droplet size
%       estimates. This value will be used to ensure only contiguous data
%       above this threshold will satisfy the profile search.



function [vert_profs] = find_verticalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, stop_at_max_lwc, Nc_threshold)

% Define the length of consecutive values found above the liquid water
% content threshold that is required to deem it a vertical profile.
consecutive_length_threshold = 5;

% Store whether or not the 2DC was conforming
vert_profs.flag_2DC_data_is_conforming = vocalsRex.flag_2DC_data_is_conforming;

% Let's also store the time vector in UTC format

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

% IF the CDP data was sampled at 10Hz, make a time vector to use for
% indexing
if length(vocalsRex.re_CDP)>length(vocalsRex.time)
    time_sps10 = linspace(vocalsRex.time(1), vocalsRex.time(end), length(vocalsRex.re_CDP));
end


% ---- Vertical Profile Requirements ----
% dz/dt must be non-zero.
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


% The plane's vertical velocity exceeds 2 m/s on average when it ascends or
% descends. Use this the window the data. Find

index_ascend_or_descend = find(abs(dz_dt_mean)>2);


% Find consecutive logical true values, which represent stand alone
% profiles where the plane is climbing or descending.
index_consec = find(diff(index_ascend_or_descend)~=1);

% include a 1 to start
index_consec = [0, index_consec];


% store the total number of droplets in each bin to calculate the LWC for
% each instrument
Nc_per_bin = cell(1, length(index_consec)-1);

% -----------------------------------------------------------------------
% For each break in the data, create a profile
for ii = 1:length(index_consec)-1

    % read in the total number of droplets per unit volume
    vert_profs.Nc{ii} = vocalsRex.total_Nc(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.Nc_CDP{ii} = vocalsRex.total_Nc_CDP(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.Nc_2DC{ii} = vocalsRex.total_Nc_2DC(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';

    % read in the total number of droplets for each size bin, but there is
    % no need to store this variable
    Nc_per_bin{ii} = vocalsRex.Nc(:,index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))));

    % read in the number of droplets for each size bin, and divide by the
    % width of the size bin to estimate the droplet distribution
    vert_profs.nr{ii} = vocalsRex.Nc(:,index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))./...
        repmat(diff(double(vocalsRex.drop_radius_bin_edges))', 1, length(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1)))));

    vert_profs.time{ii} = vocalsRex.time(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.time_utc{ii} = UTC_starttime + double(vert_profs.time{ii})./3600;
    vert_profs.latitude{ii} = vocalsRex.latitude(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.longitude{ii} = vocalsRex.longitude(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.altitude{ii} = vocalsRex.altitude(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';

    if vocalsRex.flag_2DC_data_is_conforming==true
        vert_profs.re{ii} = vocalsRex.re(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
        vert_profs.re_2DC{ii} = vocalsRex.re_2DC(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    else
        vert_profs.mean_r_2DC{ii} = vocalsRex.mean_r_2DC(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    end

    % check to see if the CDP data was sampled at 1Hz or 10 Hz
    if length(vocalsRex.re_CDP)>length(vocalsRex.time)
        % CDP data sampled at 10Hz
        indices_sps10 = time_sps10>=vert_profs.time{ii}(1) & time_sps10<=vert_profs.time{ii}(end);
        vert_profs.re_CDP{ii} = vocalsRex.re_CDP(indices_sps10)';
        vert_profs.lwc_CDP{ii} = vocalsRex.lwc_CDP(indices_sps10)';
    else
        % CDP data sampled at 1Hz
        vert_profs.re_CDP{ii} = vocalsRex.re_CDP(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
        vert_profs.lwc_CDP{ii} = vocalsRex.lwc_CDP(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    end

    vert_profs.lwc{ii} = vocalsRex.lwc(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.lwc_2DC{ii} = vocalsRex.lwc_2DC(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';

    vert_profs.horz_wind_speed{ii} = vocalsRex.horz_wind_speed(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))));
    vert_profs.horz_wind_direction{ii} = vocalsRex.horz_wind_direction(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))));


    % ----------- These variables are not being used -----------------
    %     vert_profs.SWT{ii} = vocalsRex.SWT(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    %     vert_profs.SWB{ii} = vocalsRex.SWB(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    %     vert_profs.LWT{ii} = vocalsRex.LWT(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    %     vert_profs.LWB{ii} = vocalsRex.LWB(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';

end

% Grab the non-time-stamped data as well
vert_profs.drop_radius_bin_edges = double(vocalsRex.drop_radius_bin_edges);
vert_profs.drop_radius_bin_center = double(vocalsRex.drop_radius_bin_center);
vert_profs.startTime = vocalsRex.startTime;                                                                    % We have to assume that this is in UTC time as well




% -----------------------------------------------------------------------
% Now find which of these profiles start and end with Nc<1 AND have some values above 1

index2delete = [];

for ii = 1:length(vert_profs.Nc)

    %[vert_profs.Nc{ii}(1), vert_profs.Nc{ii}(end), any(vert_profs.Nc{ii}>10^7), all(vert_profs.Nc{ii}<1)]

    % Check to see if any of these statements are met. If so, delete the
    % profile
    if any(vert_profs.Nc{ii}>10^7) || all(vert_profs.Nc{ii}<Nc_threshold)

        % if this is true, mark the index for deletion
        index2delete = [index2delete, ii];

    end




end




% Let's delete all cells that met the above conditions
vert_profs.Nc(index2delete) = [];
Nc_per_bin(index2delete) = [];
vert_profs.nr(index2delete) = [];
vert_profs.time(index2delete) = [];
vert_profs.time_utc(index2delete) = [];
vert_profs.lwc(index2delete) = [];
vert_profs.latitude(index2delete) = [];
vert_profs.longitude(index2delete) = [];
vert_profs.altitude(index2delete) = [];

if vocalsRex.flag_2DC_data_is_conforming==true
    vert_profs.re(index2delete) = [];
    vert_profs.re_2DC(index2delete) = [];
else
    vert_profs.mean_r_2DC(index2delete) = [];
end

vert_profs.re_CDP(index2delete) = [];
vert_profs.lwc_CDP(index2delete) = [];
vert_profs.lwc_2DC(index2delete) = [];
vert_profs.Nc_CDP(index2delete) = [];
vert_profs.Nc_2DC(index2delete) = [];
vert_profs.horz_wind_speed(index2delete) = [];
vert_profs.horz_wind_direction(index2delete) = [];

% ----------- These variables are not being used -----------------
% vert_profs.SWT(index2delete) = [];
% vert_profs.SWB(index2delete) = [];
% vert_profs.LWT(index2delete) = [];
% vert_profs.LWB(index2delete) = [];




% -----------------------------------------------------------------------
% Usually at the end of each profile there are several zeros. Let's delete
% all 0 values at the end of the vector ONLY when all remaining values are
% 0

% for ii = 1:length(vert_profs.Nc)
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
max_lwc = zeros(1, length(vert_profs.lwc));


% create a zero vector that defines which profiles don't meet requirements
% and need to be delted
%profile_idx_2delete = zeros(1, length(vert_profs.lwc));


for ii = 1:length(vert_profs.lwc)

    % Find the first data point and the last data point where the LWC
    % threshold is exceeded. That is, before the first index, several
    % values should be below the threshold. And several data points
    % after the last index should also be below the threshold

    indexes_above_LWC_Nc_threshold = vert_profs.lwc{ii}>=LWC_threshold & vert_profs.Nc{ii}>=Nc_threshold;


    % find the first 1, the first data point that exceeds the LWC
    % threshold. The first value in the index below is where the
    % profile will start. The last index below will be the end of the
    % profile
    indexes2keep = find(indexes_above_LWC_Nc_threshold);

    % We also need to check every profile to ensure the liquid water
    % content stays above our threshold for the entire profile. It should
    % only drop below our threshold at before and after the profile found
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
                % content threshold
                if all(vert_profs.lwc{ii}(firstIndex:lastIndex)>=LWC_threshold & vert_profs.Nc{ii}(firstIndex:lastIndex)>=Nc_threshold)==true

                    % compute the length of the vector
                    consecutive_length(firstIndex, lastIndex) = length(vert_profs.lwc{ii}(firstIndex:lastIndex));

                end

            end
        end

        % Now we find the first and last index that corresponds to the longest
        % consecutive string of values above the LWC threshold
        [~, max_idx] = max(consecutive_length, [], 'all');

        % the row and column correspond the the first and last index,
        % respectively
        [firstIndex, lastIndex] = ind2sub(size(consecutive_length), max_idx);



        % find the max LWC and the index between the first and last index.
        % The max value wil be used to remove profiles that don't have a max
        % LWC above the minimum threshold. The max_index will be used if the
        % profile needs to stop at the max LWC value as well
        [max_lwc(ii), index_max_lwc] = max(vert_profs.lwc{ii}(firstIndex:lastIndex));

        [~, index_absolute_max_lwc] = max(vert_profs.lwc{ii});


    end



    if stop_at_max_lwc == true


        if isempty(indexes2keep)==true

            % do nothing if this is true. The values will simply be set to
            % zero.

        else

            % if the plane is ascending, and the LWC is growing with
            % altitude within the cloud, the max LWC value will the final
            % index. If the plane is descending, and the LWC is growing
            % with altitude within the cloud, then the LWC will be the
            % first index.

            % The data will start at the first index above the minimum
            % threshold, found above, and end at the index of maximum LWC,
            % found in this if statement

            dz_dt = diff(vert_profs.altitude{ii}(firstIndex:lastIndex))./diff(double(vert_profs.time{ii}(firstIndex:lastIndex)))';

            % first let's smooth out the LWC profile
            lwc_smooth = movmean(double(vert_profs.lwc{ii}(firstIndex:lastIndex)), 12);
            dLWC_dz = diff(lwc_smooth)./diff(vert_profs.altitude{ii}(firstIndex:lastIndex));

            % if the liquid water content is increasing with altitude, then
            % the max value is the last index. If the liquid water content
            % decreases with altitude, then the first index is the max
            % value

            % first determine if the plane is ascending or descending

            if mean(dz_dt)>0
                % if true then the plane is ascending
                if mean(dLWC_dz)>0
                    % and the lwc is increasing, make the last index the
                    % max LWC value

                    lastIndex = firstIndex + index_max_lwc -1;

                elseif mean(dLWC_dz)<0
                    % then the lwc is decreasing with altitude and we
                    % should set the first index with the max value

                    firstIndex = firstIndex + index_max_lwc -1;

                end

            elseif mean(dz_dt)<0
                % then the plane is descending
                if mean(dLWC_dz)>0
                    % then the lwc increases with decreasing height, which
                    % means the lwc decreases with increasing height, and
                    % so the first index should be set as the max index.

                    firstIndex = firstIndex + index_max_lwc -1;

                    % and we keep the last index

                elseif mean(dLWC_dz)<0
                    % then the lwc is descreasing as we descend, which
                    % means the LWC is increasing with increasing altitude,
                    % and the last index should be set as the max value

                    lastIndex = firstIndex + index_max_lwc -1;

                end


            end




        end


    end


    % if none of the values within the vertical profile are above the LWC
    % threshold, delete this profile. But for now, just set it to be a
    % vector of 0's
    if isempty(indexes2keep)==true

        vert_profs.Nc{ii} = 0;
        Nc_per_bin{ii} = 0;
        vert_profs.nr{ii} = 0;
        vert_profs.time{ii} = 0;
        vert_profs.time_utc{ii} = 0;
        vert_profs.lwc{ii} = 0;
        vert_profs.latitude{ii} = 0;
        vert_profs.longitude{ii} = 0;
        vert_profs.altitude{ii} = 0;

        if vocalsRex.flag_2DC_data_is_conforming==true
            vert_profs.re{ii} = 0;
            vert_profs.re_2DC{ii} = 0;
        else
            vert_profs.mean_r_2DC{ii} = 0;
        end


        vert_profs.re_CDP{ii} = 0;
        vert_profs.lwc_CDP{ii} = 0;
        vert_profs.lwc_2DC{ii} = 0;
        vert_profs.Nc_CDP{ii} = 0;
        vert_profs.Nc_2DC{ii} = 0;
        vert_profs.horz_wind_speed{ii} = 0;
        vert_profs.horz_wind_direction{ii} = 0;

        % ----------- These variables are not being used -----------------
        %         vert_profs.SWT{ii} = 0;
        %         vert_profs.SWB{ii} = 0;
        %         vert_profs.LWT{ii} = 0;
        %         vert_profs.LWB{ii} = 0;


        % if the longest consecutive vector of measurements above the
        % defined liquid water content threshold is less than the
        % consecutive length threshold, we won't keep this vertical
        % profile.

    elseif length(firstIndex:lastIndex)<consecutive_length_threshold

        vert_profs.Nc{ii} = 0;
        Nc_per_bin{ii} = 0;
        vert_profs.nr{ii} = 0;
        vert_profs.time{ii} = 0;
        vert_profs.time_utc{ii} = 0;
        vert_profs.lwc{ii} = 0;
        vert_profs.latitude{ii} = 0;
        vert_profs.longitude{ii} = 0;
        vert_profs.altitude{ii} = 0;

        if vocalsRex.flag_2DC_data_is_conforming==true
            vert_profs.re{ii} = 0;
            vert_profs.re_2DC{ii} = 0;
        else
            vert_profs.mean_r_2DC{ii} = 0;
        end

        vert_profs.re_CDP{ii} = 0;
        vert_profs.lwc_CDP{ii} = 0;
        vert_profs.lwc_2DC{ii} = 0;
        vert_profs.Nc_CDP{ii} = 0;
        vert_profs.Nc_2DC{ii} = 0;
        vert_profs.horz_wind_speed{ii} = 0;
        vert_profs.horz_wind_direction{ii} = 0;

        % ----------- These variables are not being used -----------------
        %         vert_profs.SWT{ii} = 0;
        %         vert_profs.SWB{ii} = 0;
        %         vert_profs.LWT{ii} = 0;
        %         vert_profs.LWB{ii} = 0;




    else

        % delete all time-stamped data that meet the above logical statement

        vert_profs.Nc{ii} = vert_profs.Nc{ii}(firstIndex:lastIndex);
        Nc_per_bin{ii} = Nc_per_bin{ii}(:,firstIndex:lastIndex);
        vert_profs.nr{ii} = vert_profs.nr{ii}(:, firstIndex:lastIndex);
        vert_profs.time{ii} = vert_profs.time{ii}(firstIndex:lastIndex);
        vert_profs.time_utc{ii} = vert_profs.time_utc{ii}(firstIndex:lastIndex);
        vert_profs.lwc{ii} = vert_profs.lwc{ii}(firstIndex:lastIndex);
        vert_profs.latitude{ii} = vert_profs.latitude{ii}(firstIndex:lastIndex);
        vert_profs.longitude{ii} = vert_profs.longitude{ii}(firstIndex:lastIndex);
        vert_profs.altitude{ii} = vert_profs.altitude{ii}(firstIndex:lastIndex);

        if vocalsRex.flag_2DC_data_is_conforming==true
            vert_profs.re{ii} = vert_profs.re{ii}(firstIndex:lastIndex);
            vert_profs.re_2DC{ii} = vert_profs.re_2DC{ii}(firstIndex:lastIndex);
        else
            vert_profs.mean_r_2DC{ii} = vert_profs.mean_r_2DC{ii}(firstIndex:lastIndex);
        end

        % check to see if were dealing with SPS1 of SPS10 data
        if length(vocalsRex.re_CDP)>length(vocalsRex.time)
            error([newline, 'I dont know how to handle SPS10 data.', newline])
        else
            vert_profs.re_CDP{ii} = vert_profs.re_CDP{ii}(firstIndex:lastIndex);
            vert_profs.lwc_CDP{ii} = vert_profs.lwc_CDP{ii}(firstIndex:lastIndex);
        end

        vert_profs.lwc_2DC{ii} = vert_profs.lwc_2DC{ii}(firstIndex:lastIndex);
        vert_profs.Nc_CDP{ii} = vert_profs.Nc_CDP{ii}(firstIndex:lastIndex);
        vert_profs.Nc_2DC{ii} = vert_profs.Nc_2DC{ii}(firstIndex:lastIndex);
        vert_profs.horz_wind_speed{ii} = vert_profs.horz_wind_speed{ii}(firstIndex:lastIndex);
        vert_profs.horz_wind_direction{ii} = vert_profs.horz_wind_direction{ii}(firstIndex:lastIndex);

        % ----------- These variables are not being used -----------------
        %         vert_profs.SWT{ii} = vert_profs.SWT{ii}(firstIndex:lastIndex);
        %         vert_profs.SWB{ii} = vert_profs.SWB{ii}(firstIndex:lastIndex);
        %         vert_profs.LWT{ii} = vert_profs.LWT{ii}(firstIndex:lastIndex);
        %         vert_profs.LWB{ii} = vert_profs.LWB{ii}(firstIndex:lastIndex);


    end


    % Flag profiles where the LWC has been set to a single value comprised of a 0
    % These need to be deleted.

    if length(vert_profs.lwc{ii})==1 && vert_profs.lwc{ii}==0

        profile_idx_2delete(ii) = true;


    % Also check to see if the total profile found is less than 50 meters
    % deep. If it is, flag it for deletion
    elseif abs(vert_profs.altitude{ii}(1) - vert_profs.altitude{ii}(end))<50

        profile_idx_2delete(ii) = true;

    else

        profile_idx_2delete(ii) = false;

    end



end






% the is_profile_zero logical array defines which profiles to delete


% if this is true, delete all entries for the index
if sum(profile_idx_2delete)>0

    % delete all time-stamped data that meet the above logical statement

    vert_profs.Nc(profile_idx_2delete) = [];
    Nc_per_bin(profile_idx_2delete) = [];
    vert_profs.nr(profile_idx_2delete) = [];
    vert_profs.time(profile_idx_2delete) = [];
    vert_profs.time_utc(profile_idx_2delete) = [];
    vert_profs.lwc(profile_idx_2delete) = [];
    vert_profs.latitude(profile_idx_2delete) = [];
    vert_profs.longitude(profile_idx_2delete) = [];
    vert_profs.altitude(profile_idx_2delete) = [];

    if vocalsRex.flag_2DC_data_is_conforming==true
        vert_profs.re(profile_idx_2delete) = [];
        vert_profs.re_2DC(profile_idx_2delete) = [];
    else
        vert_profs.mean_r_2DC(profile_idx_2delete) = [];
    end


    vert_profs.re_CDP(profile_idx_2delete) = [];
    vert_profs.lwc_CDP(profile_idx_2delete) = [];
    vert_profs.lwc_2DC(profile_idx_2delete) = [];
    vert_profs.Nc_CDP(profile_idx_2delete) = [];
    vert_profs.Nc_2DC(profile_idx_2delete) = [];
    vert_profs.horz_wind_speed(profile_idx_2delete) = [];
    vert_profs.horz_wind_direction(profile_idx_2delete) = [];

    % ----------- These variables are not being used -----------------
    %     vert_profs.SWT(profile_idx_2delete) = [];
    %     vert_profs.SWB(profile_idx_2delete) = [];
    %     vert_profs.LWT(profile_idx_2delete) = [];
    %     vert_profs.LWB(profile_idx_2delete) = [];

end


% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Store the LWC threshold used
vert_profs.LWC_threshold = LWC_threshold;            % g/m^3

% Store the Nc threshold used
vert_profs.Nc_threshold = Nc_threshold;            % g/m^3






% ----------------------------------------------------------------------
% ------------------ Compute Liquid Water Path -------------------------
% ----------------------------------------------------------------------



% we want to compute the total LWP, and the LWP for the two instruments
% used to measure droplets.
% The CDP probe measures radii between 0.875 and 25.055 microns.
% The 2DC probe measures radii between 31.25 and 793 microns

if iscell(vert_profs.altitude)==true

    num_profiles = length(vert_profs.altitude);


    % step through each profile
    for nn = 1:num_profiles


        % LWP is calculated by integrating from cloud bottom to
        % cloud top. If the plane decreases in altitude, we need to
        % integrate from the end of the profile (cloud bottom) to the
        % begining (cloud top). If the plane is ascending, we do
        % the opposite.

        dz_dt = diff(reshape(vert_profs.altitude{nn}, [],1))./diff(reshape(vert_profs.time{nn}, [], 1));

        if mean(dz_dt)>0
            % then the plane is ascending!

            % Compute the total LWP
            vert_profs.lwp{nn} = trapz(vert_profs.altitude{nn}, vert_profs.lwc{nn});            % g/m^2

            % ------ Compute the CDP LWP ---------
            vert_profs.lwp_CDP{nn} = trapz(vert_profs.altitude{nn}, vert_profs.lwc_CDP{nn});


            % ------ Compute the 2DC LWP ---------
            vert_profs.lwp_2DC{nn} = trapz(vert_profs.altitude{nn}, vert_profs.lwc_2DC{nn});


        elseif mean(dz_dt)<0
            % then the plane is descending!

            % Compute the total LWP
            vert_profs.lwp{nn} = trapz(flipud(vert_profs.altitude{nn}), flipud(vert_profs.lwc{nn}));            % g/m^2

            % ------ Compute the CDP LWP ---------
            vert_profs.lwp_CDP{nn} = trapz(flipud(vert_profs.altitude{nn}), flipud(vert_profs.lwc_CDP{nn}));


            % ------ Compute the 2DC LWP ---------
            % compute the 2DC LWP by integration over the cloud depth
            vert_profs.lwp_2DC{nn} = trapz(flipud(vert_profs.altitude{nn}), flipud(vert_profs.lwc_2DC{nn}));





        end



    end




else

    % if vocalsRex is not a cell, then there is only one profile


    % LWP is calculated by integrating from cloud bottom to
    % cloud top. If the plane decreases in altitude, we need to
    % integrate from the end of the profile (cloud bottom) to the
    % begining (cloud top). If the plane is ascending, we do
    % the opposite.

    dz_dt = diff(reshape(vert_profs.altitude, [],1))./diff(reshape(vert_profs.time, [], 1));

    if mean(dz_dt)>0
        % then the plane is ascending!

        % Compute the total LWP
        vert_profs.lwp = trapz(vert_profs.altitude, vert_profs.lwc);            % g/m^2

        % ------ Compute the CDP LWP ---------
        % compute the CDP LWP by integration over the cloud depth
        vert_profs.lwp_CDP = trapz(vert_profs.altitude, vert_profs.lwc_CDP);


        % ------ Compute the 2DC LWP ---------
        % compute the 2DC LWP by integration over the cloud depth
        vert_profs.lwp_2DC = trapz(vert_profs.altitude, vert_profs.lwc_2DC);


    elseif mean(dz_dt)<0
        % then the plane is descending!

        % Compute the total LWP
        vert_profs.lwp = trapz(flipud(vert_profs.altitude), flipud(vert_profs.lwc));            % g/m^2

        % ------ Compute the CDP LWP ---------
        % compute the CDP LWP by integration over the cloud depth
        vert_profs.lwp_CDP = trapz(flipud(vert_profs.altitude), flipud(vert_profs.lwc_CDP));


        % ------ Compute the 2DC LWP ---------
        % compute the 2DC LWP by integration over the cloud depth
        vert_profs.lwp_2DC = trapz(flipud(vert_profs.altitude), flipud(vert_profs.lwc_2DC));





    end




end



% ----------------------------------------------------------------------
% ----------- COMPUTE THE HORIZONTAL DISTANCE TRAVELLED ----------------
% ----------------------------------------------------------------------

% Also compute the cloud depth and then the slant path travelled in the
% cloud by the plane us pythagreous' theorem. The compute the zenith angle
% of the slant path with repsect to the cloud base.


% why do we need this? Optical depth is a path dependent quantity. The
% properties of the medium that make it attenuating must be integrated
% along the entire line of sight. Cloud optical depth is defind as the
% vertically integrated optical depth from cloud top to cloud bottom.
% Assuming a plane-parallel cloud, the equation is: 
% int( sec(theta) * Qe * pi * r^2(z) * Nc(z)) dz
% The term sec(theta) accounts for the zenith angle off from a true
% vertical path (one that travels against the vector of the gravitational
% force). To estimate what sec(theta) is, we need to compute the horizontal
% distance travelled during the vertical profile. Cos(theta) is the ratio
% of the the vertical distance travelled to the hypotenuse. And sec(theta)
% is 1/cos(theta).


% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with a length unit of meters.
wgs84 = wgs84Ellipsoid("m");        % units of meters


for nn = 1:length(vert_profs.Nc)

    % Step through each point to get the linear distance travelled as a
    % vector
    horz_distance_travelled = zeros(1, length(vert_profs.latitude{nn}));

    for xx = 2:length(vert_profs.latitude{nn})

        % Find the linear distance between the start and end of the horizontal profile.
        % When you specify a reference ellipsoid as input to the distance function,
        % the function returns linear distances in the units of the semimajor axis of the ellipsoid.
        horz_distance_travelled(xx) = distance(vert_profs.latitude{nn}(1), vert_profs.longitude{nn}(1),...
            vert_profs.latitude{nn}(xx), vert_profs.longitude{nn}(xx),wgs84);

    end

    vert_profs.horz_dist{nn} = horz_distance_travelled;

    % Using the horizontal distance travelled, and the vertical depth of
    % the cloud, use pythareous's theorem to estimate the slant path
    % travelled within the cloud
    vert_profs.cloud_depth{nn} = max(vert_profs.altitude{nn}) - min(vert_profs.altitude{nn});
    vert_profs.slath_path{nn} = sqrt(vert_profs.horz_dist{nn}(end)^2 + vert_profs.cloud_depth{nn}^2);

    % Compute the zenith angle of the slant path with respect to the cloud
    % base, assuming a plane-parallel cloud and a straight line for the
    % planes trajectory

    vert_profs.VR_zenith_angle{nn} = atand(vert_profs.horz_dist{nn}(end)/vert_profs.cloud_depth{nn});


end




% ------------------------------------------------------------------
% ------------------ Compute optical depth -------------------------
% ------------------------------------------------------------------

% compute the optical depth for each vertical profile and store it


% optical depth is defined to be 0 at cloud top and increasing towards
% cloud bottom

% We do NOT include the zenith angle of the plane as it travels through the
% cloud because we are concerned with the cloud optical depth and not the
% optical depth along the slant path. Cloud optical depth is defined soley
% as the optical depth over the vertical dimension of a cloud.

% if there is more than one profile, the data will be stored in a cell
% array
if iscell(vert_profs.altitude)==true

    num_profiles = length(vert_profs.altitude);


    % step through each profile
    for nn = 1:num_profiles


        vector_length = length(vert_profs.altitude{nn});
        vert_profs.tau{nn} = zeros(1,vector_length-1);

        % compute sec(theta) by first computing the hypotenuse
        


        % Optical thickness is defined by integrating from cloud top to
        % cloud bottom. If the plane increases in altitude, we need to
        % integrating from the end of the profile (cloud top) to the
        % begining (cloud bottom). If the plane is descending, we do
        % the opposite.

        dz_dt = diff(vert_profs.altitude{nn})./diff(vert_profs.time{nn})';


        if mean(dz_dt)>0
            % then the plane is ascending!


            % step through the altitude array
            for ii = 1:vector_length-1


                % we have to convert Nc and re to have the same units as the alitude,
                % which is in meters

                if vocalsRex.flag_2DC_data_is_conforming==true
                    re_meters = vert_profs.re{nn}(vector_length-ii:vector_length)./1e6;                      % meters
                else
                    % What choice do we have? I guess we will just use the
                    % effevtive radius from the 2DC data, but this will
                    % underestimate the optical depth
                    % if the re CDP data is taken at 10 samples per second,
                    % only take every 10th value. 
                    if length(vocalsRex.re_CDP)>length(vocalsRex.time)
                        error([newline, 'I dont know how to handle SPS10 data!', newline])
                    
                    else
                        re_meters = vert_profs.re_CDP{nn}(vector_length-ii:vector_length)./1e6;                      % meters
                    end

                end

                total_Nc_meters = vert_profs.Nc{nn}(vector_length-ii:vector_length).*1e6;                           % #/m^3
                altitude = vert_profs.altitude{nn}(end) -  vert_profs.altitude{nn}(vector_length-ii:vector_length); % meters


                % We assume the droplet size is appreciably larger than the
                % incident wavelength (something in the visible, like 550 nm)
                % so we can assume the extinction efficiency is 2. This leads
                % to the following equation for optical depth:
                % tau = integral( Q_e * pi * r^2(z) * N_c(z) ) dz
                % where Q_e is set to 2
                %vert_profs.tau{nn}(ii) = 2*pi* trapz(flipud(altitude), flipud(re_meters.^2 .* total_Nc_meters));

                % Or we can use a pre-computed mie table to use a more
                % accurate value for the extinction efficiency.
                %Q_e = interp_mie_computed_tables([linspace(550, 550, length(re_meters))', re_meters.*1e6], 'mono', true);

                % Or we could compute the average extinction efficiency
                % over a droplet size distrubution
                [~, Qe_avg, ~] = average_mie_over_size_distribution(re_meters.*1e6, linspace(10,10,length(re_meters)),...
                    550, 'water', 'gamma', ii);

                vert_profs.tau{nn}(ii) = pi* trapz(flipud(altitude), flipud(Qe_avg' .* re_meters.^2 .* total_Nc_meters));

            end

        elseif mean(dz_dt)<0
            % then the plane is descending!

            % step through the altitude array
            for ii = 1:vector_length-1


                % we have to convert Nc and re to have the same units as the alitude,
                % which is in meters
                if vocalsRex.flag_2DC_data_is_conforming==true
                    re_meters = vert_profs.re{nn}(1:ii+1)./1e6;                      % meters
                else
                    % What choice do we have? I guess we will just use the
                    % effevtive radius from the 2DC data, but this will
                    % underestimate the optical depth
                    re_meters = vert_profs.re_CDP{nn}(1:ii+1)./1e6;                      % meters
                end

                total_Nc_meters = vert_profs.Nc{nn}(1:ii+1).*1e6;                           % #/m^3
                altitude = vert_profs.altitude{nn}(1) -  vert_profs.altitude{nn}(1:ii+1);   % meters


                % We assume the droplet size is appreciably larger than the
                % incident wavelength (something in the visible, like 550 nm)
                % so we can assume the extinction efficiency is 2. This leads
                % to the following equation for optical depth:
                % tau = integral( Q_e * pi * r^2(z) * N_c(z) ) dz
                % where Q_e is set to 2
                %vert_profs.tau{nn}(ii) = 2*pi* trapz(altitude, re_meters.^2 .* total_Nc_meters);


                % Or we can use a pre-computed mie table to use a more
                % accurate value for the extinction efficiency.
                %Q_e = interp_mie_computed_tables([linspace(550, 550, length(re_meters))', re_meters'.*1e6], 'mono', true);


                % Or we could compute the average extinction efficiency
                % over a droplet size distrubution
                [~, Qe_avg, ~] = average_mie_over_size_distribution(re_meters.*1e6, linspace(10,10,length(re_meters)),...
                    550, 'water', 'gamma', ii);


                vert_profs.tau{nn}(ii) = pi* trapz(altitude, Qe_avg(:,end) .* re_meters.^2 .* total_Nc_meters);



            end

        end



        % add a zero at the begining!
        vert_profs.tau{nn} = [0,vert_profs.tau{nn}];

    end




else

    % if vocalsRex is not a cell, then there is only one profile

    vector_length = length(vert_profs.altitude);
    vert_profs.tau = zeros(1,vector_length-1);

    % Optical thickness is defined by integrating from cloud top to
    % cloud bottom. If the plane increases in altitude, we need to
    % integrating from the end of the profile (cloud top) to the
    % begining (cloud bottom). If the plane is descending, we do
    % the opposite.

    dz_dt = diff(vert_profs.altitude{nn})./diff(vert_profs.time{nn})';

    if mean(dz_dt)>0


        % step through the altitude array
        for ii = 1:vector_length-1

            % we have to convert Nc and re to have the same units as the alitude,
            % which is in meters
            if vocalsRex.flag_2DC_data_is_conforming==true
                re_meters = vert_profs.re(vector_length-ii:vector_length)./1e6;                      % meters
            else
                % What choice do we have? I guess we will just use the
                % effevtive radius from the 2DC data, but this will
                % underestimate the optical depth
                re_meters = vert_profs.re_CDP(vector_length-ii:vector_length)./1e6;                      % meters
            end

            total_Nc_meters = vert_profs.Nc(vector_length-ii:vector_length).*1e6;                           % #/m^3
            altitude = vert_profs.altitude(end) -  vert_profs.altitude(vector_length-ii:vector_length);


            % We assume the droplet size is appreciably larger than the
            % incident wavelength (something in the visible, like 550 nm)
            % so we can assume the extinction efficiency is 2. This leads
            % to the following equation for optical depth:
            % tau = integral( Q_e * pi * r^2(z) * N_c(z) ) dz
            % where Q_e is set to 2
            %vert_profs.tau(ii) = 2*pi* trapz(flipud(altitude), flipud(re_meters.^2 .* total_Nc_meters));

            % Or we can use a pre-computed mie table to use a more
            % accurate value for the extinction efficiency.
            %Q_e = interp_mie_computed_tables([linspace(550, 550, length(re_meters))', re_meters.*1e6], 'mono', true);

            % Or we could compute the average extinction efficiency
            % over a droplet size distrubution
            [~, Qe_avg, ~] = average_mie_over_size_distribution(re_meters.*1e6, linspace(10,10,length(re_meters)),...
                550, 'water', 'gamma', ii);

            vert_profs.tau(ii) = pi* trapz(flipud(altitude), flipud(Qe_avg' .* re_meters.^2 .* total_Nc_meters));





        end

    elseif mean(dz_dt)<0

        % The plane is descending! The cloud top values come first


        % we have to convert Nc and re to have the same units as the alitude,
        % which is in meters
        if vocalsRex.flag_2DC_data_is_conforming==true
            re_meters = vert_profs.re(1:ii+1)./1e6;                      % meters
        else
            % What choice do we have? I guess we will just use the
            % effevtive radius from the 2DC data, but this will
            % underestimate the optical depth
            re_meters = vert_profs.re_CDP(1:ii+1)./1e6;                     % meters
        end

        total_Nc_meters = vert_profs.Nc(1:ii+1).*1e6;                           % #/m^3
        altitude = vert_profs.altitude(1) -  vert_profs.altitude(1:ii+1);


        % We assume the droplet size is appreciably larger than the
        % incident wavelength (something in the visible, like 550 nm)
        % so we can assume the extinction efficiency is 2. This leads
        % to the following equation for optical depth:
        % tau = integral( Q_e * pi * r^2(z) * N_c(z) ) dz
        % where Q_e is set to 2
        %vert_profs.tau(ii) = 2*pi* trapz(altitude, re_meters.^2 .* total_Nc_meters);


        % Or we can use a pre-computed mie table to use a more
        % accurate value for the extinction efficiency.
        %Q_e = interp_mie_computed_tables([linspace(550, 550, length(re_meters))', re_meters'.*1e6], 'mono', true);


        % Or we could compute the average extinction efficiency
        % over a droplet size distrubution
        [~, Qe_avg, ~] = average_mie_over_size_distribution(re_meters.*1e6, linspace(10,10,length(re_meters)),...
            550, 'water', 'gamma', ii);


        vert_profs.tau(ii) = pi* trapz(altitude, Qe_avg(:,end) .* re_meters.^2 .* total_Nc_meters);



    else

        error([newline, 'Something is wrong with the calculation of vertical velocity'], newline)

    end


    % add a zero at the begining!
    vert_profs.tau = [0,vert_profs.tau];


end










end