%% Analyze composite soundings from VOCALS-REx


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    folderpath_radioSonde = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/allSoundings_composite_5mb/'];

    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

elseif strcmp(which_computer,'andrewbuggee')==true


    folderpath_radioSonde = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/radiosonde/allSoundings_composite_5mb/'];


    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


end



% ---------------------
% define the ensemble filename
% profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat';
profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_04-Dec-2025.mat';



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

sounding_files = dir(folderpath_radioSonde);

% extract the names
% sounding_file_names = {sounding_files.name};

%% Step through each profile and find the radiosonde closest in space and time
% The radiosonde should be over the ocean, rather than over land


% the composite radiosonde data is saved with the day it was recorded
% on. So first find the radiosonde file for the day the current
% ensemble profile was recorded


total_column_pw = zeros(length(ensemble_profiles), 1);
above_cloud_pw = zeros(length(ensemble_profiles), 1);

cloud_topHeight = zeros(length(ensemble_profiles), 1);
cloud_baseHeight = zeros(length(ensemble_profiles), 1);

cloud_topPressure = zeros(length(ensemble_profiles), 1);
cloud_basePressure = zeros(length(ensemble_profiles), 1);

temp_prof = cell(length(ensemble_profiles), 1);
watVap_prof = cell(length(ensemble_profiles), 1);
pressure_prof = cell(length(ensemble_profiles), 1);
altitude = cell(length(ensemble_profiles), 1);


for nn = 1:length(ensemble_profiles)

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
    radiosonde = read_cls_file([folderpath_radioSonde, sounding_files(mm).name]);

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
    [radiosonde_min, idx_min] = min(abs(time_diff_C130_radiosonde));            % duration object - minimum time

    % Check that the radiosonde profile selected is over ocean
    % But some days don't have any radiosondes over the ocean. In that
    % case, just use the radiosonde closest in time
    % Also, the soundings taken at Paposo seem lousy. Don't use them
    if strcmp(radiosonde(idx_min).release_site(1:3), 'R/V') ~= true ||...
            strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true

        release_sites = {radiosonde.release_site};
        onShip = zeros(1, length(release_sites));

        for rr = 1:length(release_sites)

            if strcmp(release_sites{rr}(1:3), 'R/V')==true

                onShip(rr) = true;

            end

        end

        % if there are radiosondes released over the ocean, then I will find
        % the radiosonde launched closest in time to the measurement time of
        % the vertical profile

        if sum(onShip)>0

            while strcmp(radiosonde(idx_min).release_site(1:3), 'R/V') ~= true ||...
                    strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true

                % this radiosonde is over land, find the radiosonde closest in time
                % to the vertical profile over ocean
                time_diff_C130_radiosonde(idx_min) = duration(48,0,0);                     % set to 24 hours
                % Find the minimum
                [radiosonde_min, idx_min] = min(abs(time_diff_C130_radiosonde));            % m - minimum distance

            end

        elseif sum(onShip)==0 && strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true

            % Find the closest time the isn't at the Paposo site
            while strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true

                % this radiosonde is over land, find the radiosonde closest in time
                % to the vertical profile over ocean
                time_diff_C130_radiosonde(idx_min) = duration(48,0,0);                     % set to 24 hours
                % Find the minimum
                [radiosonde_min, idx_min] = min(abs(time_diff_C130_radiosonde));            % m - minimum distance

            end


        end


    end



    %% Find the radiosonde closest in space to the C130 sampled profile



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
    % [radiosonde_min, idx_min] = min(dist_btwn_C130_radiosonde, [], 'all');            % m - minimum distance
    %
    %
    % % Check that the radiosonde profile selected is over ocean
    % % But some days don't have any radiosondes over the ocean. In that
    % % case, just use the radiosonde closest in time
    % % Also, the soundings taken at Paposo seem lousy. Don't use them
    % if strcmp(radiosonde(idx_min).release_site(1:3), 'R/V') ~= true ||...
    %         strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true
    %
    %     release_sites = {radiosonde.release_site};
    %     onShip = zeros(1, length(release_sites));
    %
    %     for rr = 1:length(release_sites)
    %
    %         if strcmp(release_sites{rr}(1:3), 'R/V')==true
    %
    %             onShip(rr) = true;
    %
    %         end
    %
    %     end
    %
    %     % if there are radiosondes released over the ocean, then I will find
    %     % the radiosonde launched closest in time to the measurement time of
    %     % the vertical profile
    %
    %     if sum(onShip)>0
    %
    %         while strcmp(radiosonde(idx_min).release_site(1:3), 'R/V') ~= true ||...
    %                 strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true
    %
    %             % this radiosonde is over land, find the radiosonde closest
    %             % in space
    %             % to the vertical profile over ocean
    %             dist_btwn_C130_radiosonde(idx_min) = 1e10;                     % set to 10 million meters
    %             % Find the minimum
    %             [radiosonde_min, idx_min] = min(dist_btwn_C130_radiosonde, [], 'all');            % m - minimum distance
    %
    %         end
    %
    %     elseif sum(onShip)==0 && strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true
    %
    %         % Find the closest time the isn't at the Paposo site
    %         while strcmp(radiosonde(idx_min).release_site(1:4), 'SCQN') == true
    %
    %             % this radiosonde is over land, find the radiosonde closest
    %             % in space
    %             % to the vertical profile over ocean
    %             dist_btwn_C130_radiosonde(idx_min) = 1e10;                     % set to 10 million meters
    %             % Find the minimum
    %             [radiosonde_min, idx_min] = min(dist_btwn_C130_radiosonde, [], 'all');            % m - minimum distance
    %
    %         end
    %
    %
    %     end
    %
    %
    % end




    %%
    % check that the radiosonde selected has a cloud
    % Using Wang et al. 1995, (actually Wang and Rossow 1995) to detect
    % cloud base and top
    % relative humidity is the 5th column, altitude of the balloon is the
    % 15th column.

    % -------------------------------------------------------------------
    % --------------------- DETECTING CLOUD BASE ------------------------
    % -------------------------------------------------------------------
    % (a) find values of relative humidity above 87%
    idx_87 = radiosonde(idx_min).sounding_data(:,5) >= 87;


    cloud_base = [];

    if sum(idx_87) > 0

        lvls_87 = find(idx_87);

        % Find the cloud base level


        for ll = 1:length(lvls_87)
            % the first lvl will always count as a cloud layer
            if ll==1

                % Check to see that this isn't the surface level and that there is
                % atleast 2 levels with a RH of at least 87%
                if lvls_87(ll) ~= 1 || lvls_87(ll) ~= 2 &&...
                        radiosonde(idx_min).sounding_data((lvls_87(ll)+1),5) >= 87

                    cloud_base = [cloud_base; lvls_87(ll)];

                    % mark cloud-base as found
                    cloud_base_determined_radSonde = true;

                end





            else

                % check to see if this level is adjacent to the one before
                if lvls_87(ll)-lvls_87(ll-1) <= 3

                    % store as moist layer


                    % continue to the next level

                else

                    % Then this is a different moist layer to check

                    % Check to see that this isn't the surface level and that there is
                    % atleast 2 levels with a RH of at least 87%
                    if lvls_87(ll) ~= 1 || lvls_87(ll) ~= 2 &&...
                            radiosonde(idx_min).sounding_data((lvls_87(ll)+1),5) >= 87

                        cloud_base = [cloud_base; lvls_87(ll)];

                        % mark cloud-base as found
                        cloud_base_determined_radSonde = true;

                    end



                end


            end


        end


    else

        % (b) if the level is not the surface (above the first few data points)
        % RH is at least 84% but less than 87%, and RH increases by at least 3%
        % from the previous level
        idx_84_to_87 = radiosonde(idx_min).sounding_data(:,5) >= 84 &...
            radiosonde(idx_min).sounding_data(:,5) < 87;
        % check the levels just below the levels satisfying this criteria
        % change by at least 3%
        lvls_84 = find(idx_84_to_87);



        if numel(lvls_84) > 0

            for ll = 1:length(lvls_84)
                % the first lvl will always count as a cloud layer
                if ll==1

                    % compute the difference between the RH at the current level
                    % and the mean RH of the previous 3 levels. If this is a moist
                    % layer, the difference is positive and at least 3
                    if lvls_84(ll)>3

                        change_in_RH = radiosonde(idx_min).sounding_data(lvls_84(ll), 5) - ...
                            mean(radiosonde(idx_min).sounding_data((lvls_84(ll)-3):(lvls_84(ll)-1), 5));

                        if change_in_RH > 3
                            % Mark this level as a cloud base
                            cloud_base = [cloud_base; lvls_84(ll)];

                            % mark cloud-base as found
                            cloud_base_determined_radSonde = true;

                        end


                    else

                        % (c) the RH is at least 84 if this is a surface level
                        % So if this is at level 1 2 or 3, is the RH at least 84?

                        % Mark this level as a cloud base
                        cloud_base = [cloud_base; lvls_84(ll)];

                        % mark cloud-base as found
                        cloud_base_determined_radSonde = true;


                    end





                else

                    % check to see if this level is adjacent to the one before
                    if lvls_84(ll)-lvls_84(ll-1) <= 2

                        % store as moist layer


                        % continue to the next level

                    else

                        % Then this is a different moist layer to check
                        if lvls_84(ll)>3

                            change_in_RH = radiosonde(idx_min).sounding_data(lvls_84(ll), 5) - ...
                                mean(radiosonde(idx_min).sounding_data((lvls_84(ll)-3):(lvls_84(ll)-1), 5));

                            if change_in_RH > 3
                                % Mark this level as a cloud base
                                cloud_base = [cloud_base; lvls_84(ll)];

                                % mark cloud-base as found
                                cloud_base_determined_radSonde = true;

                            end


                        else

                            % (c) the RH is at least 84 if this is a surface level
                            % So if this is at level 1 2 or 3, is the RH at least 84?

                            % Mark this level as a cloud base
                            cloud_base = [cloud_base; lvls_84(ll)];

                            % mark cloud-base as found
                            cloud_base_determined_radSonde = true;


                        end







                    end


                end


            end


        else

            % how else can we define cloud top?
            % for now, let's use the vocals-rex CDP measurement to define
            % cloud base in the the radiosonde data
            [~, cloud_base] = min( abs( min(ensemble_profiles{nn}.altitude) -...
                radiosonde(idx_min).sounding_data(:,15)));


            % mark cloud-base as NOT found
            cloud_base_determined_radSonde = false;


        end
        % -------------------------------------------------------------------
        % -------------------------------------------------------------------






    end










    % ------------------------------------------------------------------
    % --------------------- DETECTING CLOUD TOP ------------------------
    % ------------------------------------------------------------------
    % The top of the moist layer is detected as the level where any one of
    % the following three criteria are met

    % Find the first layer below 84%
    % stepping through each cloud base...
    idx_top2Check = zeros(length(cloud_base), 1);

    cloud_top = zeros(length(cloud_base), 1);

    if cloud_base_determined_radSonde == true

        for bb = 1:length(cloud_base)

            idx_2Check = (cloud_base(bb)+1):size(radiosonde(idx_min).sounding_data,1);

            idx_above_base_lessThan_84 = radiosonde(idx_min).sounding_data(idx_2Check,5) < 84;

            % numbered indicies
            % num_idx = find(idx_above_base_lessThan_84);

            % find the first index above cloud bottom that is below a relative
            % humidity of 84 for at least 4 levels
            cc = 1;
            while (idx_above_base_lessThan_84(cc)==0 &&...
                    all(idx_above_base_lessThan_84((cc+1):(cc+3))==1)) ~= true


                cc = cc + 1;

            end




            idx_top2Check(bb) = cloud_base(bb) + cc;




            % (a)  where the RH is at least 87%
            % is the topmost layer at least 87?
            if radiosonde(idx_min).sounding_data(idx_top2Check(bb), 5) >= 87

                % Mark this level as the cloud top
                cloud_top(bb) = idx_top2Check(bb);








                % (b) RH is at least 84 but less than 87%, and there is at least a 3%
                % between the current layer and the higher layer (3% decrease)
            elseif (idx_top2Check(bb)+3) < size(radiosonde(idx_min).sounding_data,1)

                idx_aboveCloud = (idx_top2Check(bb)+1):(idx_top2Check(bb)+3);

                change_in_RH = radiosonde(idx_min).sounding_data(idx_top2Check(bb), 5) - ...
                    mean(radiosonde(idx_min).sounding_data(idx_aboveCloud, 5));

                if change_in_RH > 3
                    % Mark this level as a cloud base
                    cloud_top(bb) = idx_top2Check(bb);


                else

                    % set cloud top with NaN
                    cloud_top(bb) = NaN;


                end


            else

                % (c) the RH is at least 84 if this is a surface level
                % So if this is at level 1 2 or 3, is the RH at least 84?

                % Mark this level as a cloud base
                cloud_top(bb) = idx_top2Check(bb);


            end



        end

    else


        % how else can we define cloud top?
        % for now, let's use the vocals-rex CDP measurement to define
        % cloud base in the the radiosonde data
        [~, cloud_top] = min( abs( max(ensemble_profiles{nn}.altitude) -...
            radiosonde(idx_min).sounding_data(:,15)));

    end
    % ------------------------------------------------------------------
    % ------------------------------------------------------------------

    %% Determine the above cloud precipitable water amount

    % If there are multiple layers, find the cloud top height closest in
    % altitude to the cloud top height measured by the CDP
    if length(cloud_top)>1

        [~, idx_cldTop_layer] = min(radiosonde(idx_min).sounding_data(cloud_top, 15) - ...
            max(ensemble_profiles{nn}.altitude));

        % Do the same for cloud base, if needed
        if length(cloud_base)>1

            [~, idx_cldBase_layer] = min(radiosonde(idx_min).sounding_data(cloud_base, 15) - ...
                min(ensemble_profiles{nn}.altitude));

        end

        % Compute the mass of water vapor per unit volume at all altitudes
        % above the cloud top
        [rho, waterVapor_concentration_cm3] = compute_waterVapor_concentration_from_radiosonde(radiosonde(idx_min).sounding_data(:,4),...
            radiosonde(idx_min).sounding_data(:,3));



        % check total column amount
        total_column_pw(nn) = trapz(radiosonde(idx_min).sounding_data(:,15), rho);            % kg / m^2
        % total_column_pw = trapz(radiosonde(idx_minTime).sounding_data(:,15), waterVapor_concentration_cm3.*1e6 .* (con.Mol_mass_h2o_vap/con.N_A))            % kg / m^2


        % compute the above cloud precipitable water amount
        above_cloud_pw(nn) = trapz(radiosonde(idx_min).sounding_data(cloud_top(idx_cldTop_layer):end,15), rho(cloud_top(idx_cldTop_layer):end));            % kg / m^2



        % ---------------------------------------------------------------
        % ------- Save the cloud base and top height and pressure -------
        % ---------------------------------------------------------------
        cloud_topHeight(nn) = radiosonde(idx_min).sounding_data(cloud_top(idx_cldTop_layer), 15);   % meters
        cloud_topPressure(nn) = radiosonde(idx_min).sounding_data(cloud_top(idx_cldTop_layer), 15);   % meters

        cloud_baseHeight(nn) = radiosonde(idx_min).sounding_data(cloud_base(idx_cldBase_layer), 15);   % meters
        cloud_basePressure(nn) = radiosonde(idx_min).sounding_data(cloud_base(idx_cldBase_layer), 15);   % meters

        % --------------------------------------------------


    else


        % there is only one layer, so just compute the total and above
        % cloud precipitable water amount

        % Compute the mass of water vapor per unit volume at all altitudes
        % above the cloud top
        [rho, waterVapor_concentration_cm3] = compute_waterVapor_concentration_from_radiosonde(radiosonde(idx_min).sounding_data(:,4),...
            radiosonde(idx_min).sounding_data(:,3));


        % check total column amount
        total_column_pw(nn) = trapz(radiosonde(idx_min).sounding_data(:,15), rho);            % kg / m^2
        % total_column_pw = trapz(radiosonde(idx_minTime).sounding_data(:,15), waterVapor_concentration_cm3.*1e6 .* (con.Mol_mass_h2o_vap/con.N_A))            % kg / m^2


        % compute the above cloud precipitable water amount
        above_cloud_pw(nn) = trapz(radiosonde(idx_min).sounding_data(cloud_top:end,15), rho(cloud_top:end));            % kg / m^2



        % ---------------------------------------------------------------
        % ------- Save the cloud base and top height and pressure -------
        % ---------------------------------------------------------------
        cloud_topHeight(nn) = radiosonde(idx_min).sounding_data(cloud_top, 15);   % meters
        cloud_topPressure(nn) = radiosonde(idx_min).sounding_data(cloud_top, 15);   % meters

        cloud_baseHeight(nn) = radiosonde(idx_min).sounding_data(cloud_base, 15);   % meters
        cloud_basePressure(nn) = radiosonde(idx_min).sounding_data(cloud_base, 15);   % meters

        % --------------------------------------------------



    end


    % --------------------------------------------------
    % --- Store the vertical profiles of T, P and RH ---
    % --------------------------------------------------
    pressure_prof{nn} = radiosonde(idx_min).sounding_data(:, 2);   % mb
    temp_prof{nn} = radiosonde(idx_min).sounding_data(:, 3);   % C
    watVap_prof{nn} = radiosonde(idx_min).sounding_data(:, 5);   % % - Relative humidity
    altitude{nn} = radiosonde(idx_min).sounding_data(:, 15);   % meters

    % --------------------------------------------------


    % --------------------------------------------------------------
    % --- Store the name, date and start time of this radiosonde ---
    % --------------------------------------------------------------
    anicllary_info.lat_release(nn) = radiosonde(idx_min).latitude; % degrees latitude
    anicllary_info.long_release(nn) = radiosonde(idx_min).longitude; % degrees latitude

    anicllary_info.dateTime_release{nn} = radiosonde(idx_min).release_time; % degrees latitude

    % --------------------------------------------------------------




    % plot results
    % plot_radiosonde_wvConcentration_with_US_STD_ATM(waterVapor_concentration_cm3, radiosonde(idx_min).sounding_data(:,15))



    % clear a few variables
    clear cloud_base cloud_top waterVapor_concentration_cm3 radiosonde idx_min






end


%% Save the total and above cloud precipitable water vapor amounts


if strcmp(which_computer,'anbu8374')==true


    folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/paper2_prior_stats/'];



elseif strcmp(which_computer,'andrewbuggee')==true


    folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/radiosonde/paper2_prior_stats/'];


end



% ----- Save Paper 2 prior stats -------
% if exist("dist_btwn_C130_radiosonde", "var")==1
%
%     save([folderpath_2save,'precipitable_water_stats_for_paper2_closest_radiosonde_in_space_',...
%         char(datetime("today")),'.mat'],...
%         'total_column_pw', 'above_cloud_pw')
%
%
% elseif exist("time_diff_C130_radiosonde", "var")==1
%
%     save([folderpath_2save,'precipitable_water_stats_for_paper2_closest_radiosonde_in_time_',...
%         char(datetime("today")),'.mat'],...
%         'total_column_pw', 'above_cloud_pw')
%
%
% end




% ----- Save Radiosonde profiles and prior stats -------
if exist("dist_btwn_C130_radiosonde", "var")==1

    save([folderpath_2save,'radiosonde_profiles_for_paper2_measurements_closest_radiosonde_in_space_',...
        char(datetime("today")),'.mat'],...
        'total_column_pw', 'above_cloud_pw', 'watVap_prof', "temp_prof", "pressure_prof",...
        "cloud_topHeight", "cloud_topPressure", "cloud_baseHeight", "cloud_basePressure",...
        "altitude", 'anicllary_info')


elseif exist("time_diff_C130_radiosonde", "var")==1

    save([folderpath_2save,'radiosonde_profiles_for_paper2_measurements_closest_radiosonde_in_time_',...
        char(datetime("today")),'.mat'],...
        'total_column_pw', 'above_cloud_pw', 'watVap_prof', "temp_prof", "pressure_prof",...
        "cloud_topHeight", "cloud_topPressure", "cloud_baseHeight", "cloud_basePressure",...
        "altitude", 'anicllary_info')


end






%% What distribution does the ACPW follow?


% --------------------------------------------------------
% ------- ABOVE CLOUD PRECIPITABLE WATER FITTING ---------
% --------------------------------------------------------
% Let's also fit these three distributions to the optical depth data

% fit the number concentration data to a normal distribution
acpw_fit_normal = fitdist(combined_aboveCloud_pw_timeAndSpace, 'normal');
[acpw_reject_normal, acpw_p_normal] = chi2gof(tau_c, 'CDF', acpw_fit_normal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the number concentration content data to a log-normal distribution
acpw_fit_lognormal = fitdist(tau_c, 'lognormal');
[acpw_reject_lognormal, acpw_p_lognormal] = chi2gof(tau_c, 'CDF', acpw_fit_lognormal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the total number concentration data to a gamma distribution - use my custom
% libRadtran gamma distribution
acpw_fit_gamma = prob.GammaDistribution_libRadtran.fit(tau_c);
[acpw_reject_gamma, acpw_p_gamma] = chi2gof(tau_c, 'CDF', acpw_fit_gamma,...
    'alpha', significance_lvl, 'NParams', 2);

% Plot results
lgnd_fnt = 20;

figure; histogram(tau_c,'NumBins', 30, 'Normalization','pdf')
hold on
xVals = linspace(min(tau_c), max(tau_c), 1000);
plot(xVals, pdf(acpw_fit_normal, xVals))
plot(xVals, pdf(acpw_fit_lognormal, xVals))
plot(xVals, pdf(acpw_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
title('$acpw$ statistics and fits', ...
    'FontSize', 20, 'Interpreter', 'latex')
