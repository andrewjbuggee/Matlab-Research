%% Calculate ensemble statistics of LWC, re, and Nc along horizontal profiles from VOCALS-REx


% Andrew John Buggee

clear variables

%% Loop through all VOCALS-REx files and load all horizontal profiles


% Read all the file names

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    % Mac folder name
    foldername = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/';

elseif strcmp(whatComputer, 'andrewbuggee')==true

    % Macbook folder name
    foldername = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
        'VOCALS_REx/vocals_rex_data/SPS_1/'];


end




% read all file names in the above listed folder
folder_contents = dir(foldername);

% grab just the vocals filenames

filename = cell(1, length(folder_contents));
backspace = 0;

for nn = 1:length(folder_contents)

    if length(folder_contents(nn).name)>1

        if strcmp(folder_contents(nn).name(1:2), 'RF')==true

            filename{nn - backspace} = folder_contents(nn).name;

        else

            % delete this cell array
            filename(nn - backspace) = [];
            % add one to backspace
            backspace = backspace+1;

        end

    else

        % delete this cell array
        filename(nn - backspace) = [];

        % add one to backspace
        backspace = backspace+1;

    end

end



%% Read each data set and find horizontal profiles

% --------------------------------------------------------
% ---------- Define the profile thresholds --------------
% --------------------------------------------------------

% define the LWC threshold
ensemble_profiles.inputs.LWC_threshold = 0.03;          % g/m^3

% what is the maximum vertical displacement allowed for a single horizontal
% profile?
ensemble_profiles.inputs.max_vert_displacement = 20;     % meters

% define the total number concentration threshold
ensemble_profiles.inputs.Nc_threshold = 1;       %  #-droplets/cm^3

% if sorting for precipitation, provide a drizzle/precip threshold.
ensemble_profiles.sort_for_precip_driz = true;

% the logic flag below tells the code to use either profiles with
% precipitation or those without
ensemble_profiles.keep_precip_drizzle_profiles = false;             % if false, keep non-precip profiles only

% The threshold is defined as the total 2DC LWP
ensemble_profiles.precip_driz_threshold = 35;         % g/m^2

% Load data

for nn = 1:length(filename)

    disp(['nn = ',num2str(nn)])

    vocalsRex = readVocalsRex([foldername,filename{nn}]);

    % find the horizontal profiles
    horz_profs = find_horizontalProfiles_VOCALS_REx(vocalsRex, ensemble_profiles.inputs.LWC_threshold,...
        ensemble_profiles.inputs.Nc_threshold, ensemble_profiles.inputs.max_vert_displacement);

    if ensemble_profiles.sort_for_precip_driz == true
        % sort profiles into those with drizzle/precipitaiton and those without
        % The index below indicates with profiles meet the threshold
        % requirements for precipitation
        index_precip_drizzle = sort_horz_profs_for_precipitation(horz_profs, ensemble_profiles.precip_driz_threshold);
        % Do you wish to keep the profiles with precipitation or those
        % without?
        if ensemble_profiles.keep_precip_drizzle_profiles==false
            % We want the index values that only pertain to
            % non-precipitating profiles
            index_precip_drizzle = setxor((1:length(horz_profs.lwc)), index_precip_drizzle);
        end


        % grab just the variables of interest for from each horizontal profile
        if nn==1

            % we only have droplet effective radius from both instruments if the 2DC
            % data is non-zero
            if horz_profs.flag_2DC_data_is_conforming==true

                % then we have an effective radius that uses data from both instruments
                ensemble_profiles.re = horz_profs.re(index_precip_drizzle);                          % both instruments
                % and we have an effective radius from just the 2DC data
                ensemble_profiles.re_2DC = horz_profs.re_2DC(index_precip_drizzle);                   % from the 2DC instrument only

            else

                % we don't have an effective radius for the 2DC data
                % What we we have is the first moment
                ensemble_profiles.mean_r_2DC = horz_profs.mean_r_2DC(index_precip_drizzle);            % from the 2DC instrument only

            end

            ensemble_profiles.altitude = horz_profs.altitude(index_precip_drizzle);
            ensemble_profiles.horz_dist = horz_profs.horz_dist(index_precip_drizzle);
            ensemble_profiles.re_CDP = horz_profs.re_CDP(index_precip_drizzle);
            ensemble_profiles.lwc = horz_profs.lwc(index_precip_drizzle);
            ensemble_profiles.lwc_CDP = horz_profs.lwc_CDP(index_precip_drizzle);
            ensemble_profiles.lwc_2DC = horz_profs.lwc_2DC(index_precip_drizzle);
            ensemble_profiles.Nc = horz_profs.Nc(index_precip_drizzle);
            ensemble_profiles.time = horz_profs.time(index_precip_drizzle);

        else


            % we only have droplet effective radius from both instruments if the 2DC
            % data is non-zero
            if horz_profs.flag_2DC_data_is_conforming==true

                % then we have an effective radius that uses data from both instruments
                ensemble_profiles.re = [ensemble_profiles.re, horz_profs.re(index_precip_drizzle)];                          % both instruments
                % and we have an effective radius from just the 2DC data
                ensemble_profiles.re_2DC = [ensemble_profiles.re_2DC, horz_profs.re_2DC(index_precip_drizzle)];        % from the 2DC instrument only

                % Store all the other variables

                ensemble_profiles.altitude = [ensemble_profiles.altitude, horz_profs.altitude(index_precip_drizzle)];
                ensemble_profiles.horz_dist = [ensemble_profiles.horz_dist, horz_profs.horz_dist(index_precip_drizzle)];
                ensemble_profiles.re_CDP = [ensemble_profiles.re_CDP, horz_profs.re_CDP(index_precip_drizzle)];
                ensemble_profiles.lwc = [ensemble_profiles.lwc, horz_profs.lwc(index_precip_drizzle)];
                ensemble_profiles.lwc_CDP = [ensemble_profiles.lwc_CDP, horz_profs.lwc_CDP(index_precip_drizzle)];
                ensemble_profiles.lwc_2DC = [ensemble_profiles.lwc_2DC, horz_profs.lwc_2DC(index_precip_drizzle)];
                ensemble_profiles.Nc = [ensemble_profiles.Nc, horz_profs.Nc(index_precip_drizzle)];
                ensemble_profiles.time = [ensemble_profiles.time, horz_profs.time(index_precip_drizzle)];

            else

                % we don't have an effective radius for the 2DC data
                % What we we have is the first moment
                % FOR NOW - STORE THE CDP RE DATA AS THE RE DATA
                % THIS IS JUSTIFIED ONLY IF THE 2DC THRESHOLD IS SET TO A
                % VERY LOW VALUE

                % check to see if this field exists
                if isfield(ensemble_profiles, 'mean_r_2DC')==true

                    ensemble_profiles.mean_r_2DC = [ensemble_profiles.mean_r_2DC, horz_profs.mean_r_2DC(index_precip_drizzle)];            % from the 2DC instrument only

                else

                    ensemble_profiles.mean_r_2DC = horz_profs.mean_r_2DC(index_precip_drizzle);

                end

                % set the overall effective radius to be only the CDP
                % calculated effective radius
                ensemble_profiles.re = [ensemble_profiles.re, horz_profs.re_CDP(index_precip_drizzle)];                          % both instruments

                ensemble_profiles.altitude = [ensemble_profiles.altitude, horz_profs.altitude(index_precip_drizzle)];
                ensemble_profiles.horz_dist = [ensemble_profiles.horz_dist, horz_profs.horz_dist(index_precip_drizzle)];
                ensemble_profiles.re_CDP = [ensemble_profiles.re_CDP, horz_profs.re_CDP(index_precip_drizzle)];
                ensemble_profiles.lwc = [ensemble_profiles.lwc, horz_profs.lwc(index_precip_drizzle)];
                ensemble_profiles.lwc_CDP = [ensemble_profiles.lwc_CDP, horz_profs.lwc_CDP(index_precip_drizzle)];
                ensemble_profiles.lwc_2DC = [ensemble_profiles.lwc_2DC, horz_profs.lwc_2DC(index_precip_drizzle)];
                ensemble_profiles.Nc = [ensemble_profiles.Nc, horz_profs.Nc(index_precip_drizzle)];
                ensemble_profiles.time = [ensemble_profiles.time, horz_profs.time(index_precip_drizzle)];

            end



        end

    else

        % grab just the variables of interest for from each horizontal profile
        if nn==1

            % we only have droplet effective radius from both instruments if the 2DC
            % data is non-zero
            if horz_profs.flag_2DC_data_is_conforming==true

                % then we have an effective radius that uses data from both instruments
                ensemble_profiles.re = horz_profs.re;                          % both instruments
                % and we have an effective radius from just the 2DC data
                ensemble_profiles.re_2DC = horz_profs.re_2DC;                   % from the 2DC instrument only

            else

                % we don't have an effective radius for the 2DC data
                % What we we have is the first moment
                ensemble_profiles.mean_r_2DC = horz_profs.mean_r_2DC;            % from the 2DC instrument only

            end

            ensemble_profiles.altitude = horz_profs.altitude;
            ensemble_profiles.horz_dist = horz_profs.horz_dist;
            ensemble_profiles.re_CDP = horz_profs.re_CDP;
            ensemble_profiles.lwc = horz_profs.lwc;
            ensemble_profiles.lwc_CDP = horz_profs.lwc_CDP;
            ensemble_profiles.lwc_2DC = horz_profs.lwc_2DC;
            ensemble_profiles.Nc = horz_profs.Nc;
            ensemble_profiles.time = horz_profs.time;

        else



            % we only have droplet effective radius from both instruments if the 2DC
            % data is non-zero
            if horz_profs.flag_2DC_data_is_conforming==true

                % then we have an effective radius that uses data from both instruments
                ensemble_profiles.re = [ensemble_profiles.re, horz_profs.re];                          % both instruments
                % and we have an effective radius from just the 2DC data
                ensemble_profiles.re_2DC = [ensemble_profiles.re_2DC, horz_profs.re_2DC];        % from the 2DC instrument only

                % Store all the other variables

                ensemble_profiles.altitude = [ensemble_profiles.altitude, horz_profs.altitude];
                ensemble_profiles.horz_dist = [ensemble_profiles.horz_dist, horz_profs.horz_dist];
                ensemble_profiles.re_CDP = [ensemble_profiles.re_CDP, horz_profs.re_CDP];
                ensemble_profiles.lwc = [ensemble_profiles.lwc, horz_profs.lwc];
                ensemble_profiles.lwc_CDP = [ensemble_profiles.lwc_CDP, horz_profs.lwc_CDP];
                ensemble_profiles.lwc_2DC = [ensemble_profiles.lwc_2DC, horz_profs.lwc_2DC];
                ensemble_profiles.Nc = [ensemble_profiles.Nc, horz_profs.Nc];
                ensemble_profiles.time = [ensemble_profiles.time, horz_profs.time];

            else

                % we don't have an effective radius for the 2DC data
                % What we we have is the first moment
                % FOR NOW - STORE THE CDP RE DATA AS THE RE DATA
                % THIS IS JUSTIFIED ONLY IF THE 2DC THRESHOLD IS SET TO A
                % VERY LOW VALUE

                % check to see if this field exists
                if isfield(ensemble_profiles, 'mean_r_2DC')==true

                    ensemble_profiles.mean_r_2DC = [ensemble_profiles.mean_r_2DC, horz_profs.mean_r_2DC];            % from the 2DC instrument only

                else

                    ensemble_profiles.mean_r_2DC = horz_profs.mean_r_2DC;

                end

                % set the overall effective radius to be only the CDP
                % calculated effective radius
                ensemble_profiles.re = [ensemble_profiles.re, horz_profs.re_CDP];                          % both instruments

                ensemble_profiles.altitude = [ensemble_profiles.altitude, horz_profs.altitude];
                ensemble_profiles.horz_dist = [ensemble_profiles.horz_dist, horz_profs.horz_dist];
                ensemble_profiles.re_CDP = [ensemble_profiles.re_CDP, horz_profs.re_CDP];
                ensemble_profiles.lwc = [ensemble_profiles.lwc, horz_profs.lwc];
                ensemble_profiles.lwc_CDP = [ensemble_profiles.lwc_CDP, horz_profs.lwc_CDP];
                ensemble_profiles.lwc_2DC = [ensemble_profiles.lwc_2DC, horz_profs.lwc_2DC];
                ensemble_profiles.Nc = [ensemble_profiles.Nc, horz_profs.Nc];
                ensemble_profiles.time = [ensemble_profiles.time, horz_profs.time];

            end





        end

    end


end



% save the ensemble profiles
if ensemble_profiles.sort_for_precip_driz==true

    if ensemble_profiles.keep_precip_drizzle_profiles==true



        save([foldername,'ensemble_horizontal_profiles_with_precip_from_',num2str(length(filename)),...
            '_files_LWC-threshold_',num2str(ensemble_profiles.inputs.LWC_threshold), '_Nc-threshold_',...
            num2str(ensemble_profiles.inputs.Nc_threshold), '_',char(datetime("today")),'.mat'],...
            'ensemble_profiles', 'filename')

    else

        save([foldername,'ensemble_horizontal_profiles_without_precip_from_',num2str(length(filename)),...
            '_files_LWC-threshold_',num2str(ensemble_profiles.inputs.LWC_threshold), '_Nc-threshold_',...
            num2str(ensemble_profiles.inputs.Nc_threshold), '_',char(datetime("today")),'.mat'],...
            'ensemble_profiles', 'filename')

    end

else


    save([foldername,'ensemble_horizontal_profiles_from_',num2str(length(filename)), '_files_LWC-threshold_',...
        num2str(ensemble_profiles.inputs.LWC_threshold), '_Nc-threshold_',...
        num2str(ensemble_profiles.inputs.Nc_threshold), '_maxVertDisplacement_',...
        num2str(ensemble_profiles.inputs.max_vert_displacement), 'm_', char(datetime("today")),'.mat'],...
        'ensemble_profiles', 'filename')

end






%% Plot 