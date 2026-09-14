%% Calculate ensemble statistics vertical profiles measured by VOCALS-REx


% Andrew John Buggee

clear variables

%% Loop through all VOCALS-REx files and load all vertical profiles


% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    foldername = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

elseif strcmp(which_computer,'andrewbuggee')==true

    % Macbook folder name
    foldername = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

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



%% Read each data set and find vertical profiles

% --------------------------------------------------------
% ---------- Define the profile thresholds --------------
% --------------------------------------------------------

% define the LWC threshold
inputs.LWC_threshold = 0.05;       % g/m^3

% do you want to cut off the profile at the maximum LWC value?
inputs.stop_at_max_LWC = false;

% define the total number concentration threshold
inputs.Nc_threshold = 10;       %  #-droplets/cm^3

% if sorting for precipitation, provide a drizzle/precip threshold.
inputs.sort_for_precip_driz = true;

% the logic flag below tells the code to use either profiles with
% precipitation or those without
inputs.keep_precip_drizzle_profiles = true;             % if false, keep non-precip profiles only

% The threshold is defined as the total 2DC LWP
inputs.precip_driz_threshold = 5;         % g/m^2


% start an empty array
ensemble_profiles = cell([]);

% Load data

for nn = 1:length(filename)

    disp(['nn = ',num2str(nn)])

    vocalsRex = readVocalsRex([foldername,filename{nn}]);

    % find the vertical profiles
    vert_profs = find_verticalProfiles_VOCALS_REx_ver2(vocalsRex, inputs.LWC_threshold,...
        inputs.stop_at_max_LWC, inputs.Nc_threshold, which_computer);

    if inputs.sort_for_precip_driz == true
        % sort profiles into those with drizzle/precipitaiton and those without
        % The index below indicates which profiles met the threshold
        % requirements for precipitation
        indexes_with_precip = sort_vert_profs_for_precipitation(vert_profs, inputs.precip_driz_threshold);
        % Do you wish to keep the profiles with precipitation or those
        % without?
        if inputs.keep_precip_drizzle_profiles==false
            % We want the index values that only pertain to
            % non-precipitating profiles
            indexes_2keep = setxor((1:length(vert_profs)), indexes_with_precip);

        else
            indexes_2keep = 1:length(vert_profs);
        end

        % what is the current length of the ensemble profile structure?
        l = length(ensemble_profiles);

        % grab just the variables of interest from each vertical profile
        for mm = 1:length(indexes_2keep)


            % Grab all variables!
            if l==0 & mm==1
                ensemble_profiles{1} = vert_profs(indexes_2keep(mm));
            else

                ensemble_profiles{l+mm} = vert_profs(indexes_2keep(mm));
            end


            % we only have droplet effective radius from both instruments if the 2DC
            % data is non-zero
            % if vert_profs(indexes_2keep(mm)).flag_2DC_data_is_conforming
            % 
            %     % then we have an effective radius that uses data from both instruments
            %     ensemble_profiles(l+mm).re = vert_profs(indexes_2keep(mm)).re;                          % both instruments
            %     % and we have an effective radius from just the 2DC data
            %     ensemble_profiles(l+mm).re_2DC = vert_profs(indexes_2keep(mm)).re_2DC;                   % from the 2DC instrument only
            % 
            % else
            % 
            %     % we don't have an effective radius for the 2DC data
            %     % What we we have is the first moment
            %     % FOR NOW - STORE THE CDP RE DATA AS THE RE DATA
            %     % THIS IS JUSTIFIED ONLY IF THE 2DC THRESHOLD IS SET TO A
            %     % VERY LOW VALUE
            % 
            %     % check to see if this field exists
            %     if isfield(ensemble_profiles, 'mean_r_2DC')==true
            % 
            %         % we don't have an effective radius for the 2DC data
            %         % What we we have is the first moment
            %         ensemble_profiles(l+mm).mean_r_2DC = vert_profs(indexes_2keep(mm)).mean_r_2DC;            % from the 2DC instrument only
            % 
            % 
            %     end
            % 
            % 
            % end


            % ensemble_profiles(l+mm).altitude = vert_profs(indexes_2keep(mm)).altitude;
            % ensemble_profiles(l+mm).tau = vert_profs(indexes_2keep(mm)).tau;
            % ensemble_profiles(l+mm).re_CDP = vert_profs(indexes_2keep(mm)).re_CDP;
            % ensemble_profiles(l+mm).lwc = vert_profs(indexes_2keep(mm)).lwc;
            % ensemble_profiles(l+mm).lwc_CDP = vert_profs(indexes_2keep(mm)).lwc_CDP;
            % ensemble_profiles(l+mm).lwc_2DC = vert_profs(indexes_2keep(mm)).lwc_2DC;
            % ensemble_profiles(l+mm).total_Nc = vert_profs(indexes_2keep(mm)).total_Nc;
            % ensemble_profiles(l+mm).time = vert_profs(indexes_2keep(mm)).time;
            % ensemble_profiles(l+mm).Nc_vapor = vert_profs(indexes_2keep(mm)).Nc_vapor;
            % ensemble_profiles(l+mm,:).Nc_vapor_partPres = vert_profs(indexes_2keep(mm)).Nc_vapor_partPres;



        end






    else


        % what is the current length of the ensemble profile structure?
        l = length(ensemble_profiles);

        % grab just the variables of interest for from each vertical profile
        for mm = 1:length(indexes_2keep)


            % Grab all variables!
            % Grab all variables!
            if l==0 & mm==1
                ensemble_profiles{1} = vert_profs(indexes_2keep(mm));
            else

                ensemble_profiles{l+mm} = vert_profs(indexes_2keep(mm));
            end



            % % we only have droplet effective radius from both instruments if the 2DC
            % % data is non-zero
            % if vert_profs(indexes_2keep(mm)).flag_2DC_data_is_conforming
            % 
            %     % then we have an effective radius that uses data from both instruments
            %     ensemble_profiles(l+mm).re = vert_profs(indexes_2keep(mm)).re;                          % both instruments
            %     % and we have an effective radius from just the 2DC data
            %     ensemble_profiles(l+mm).re_2DC = vert_profs(indexes_2keep(mm)).re_2DC;                   % from the 2DC instrument only
            % 
            % else
            % 
            %     % we don't have an effective radius for the 2DC data
            %     % What we we have is the first moment
            %     % FOR NOW - STORE THE CDP RE DATA AS THE RE DATA
            %     % THIS IS JUSTIFIED ONLY IF THE 2DC THRESHOLD IS SET TO A
            %     % VERY LOW VALUE
            % 
            %     % check to see if this field exists
            %     if isfield(ensemble_profiles, 'mean_r_2DC')==true
            % 
            %         % we don't have an effective radius for the 2DC data
            %         % What we we have is the first moment
            %         ensemble_profiles(l+mm).mean_r_2DC = vert_profs(indexes_2keep(mm)).mean_r_2DC;            % from the 2DC instrument only
            % 
            % 
            %     end
            % 
            % 
            % end
            % 
            % 
            % ensemble_profiles(l+mm).altitude = vert_profs(indexes_2keep(mm)).altitude;
            % ensemble_profiles(l+mm).tau = vert_profs(indexes_2keep(mm)).tau;
            % ensemble_profiles(l+mm).re_CDP = vert_profs(indexes_2keep(mm)).re_CDP;
            % ensemble_profiles(l+mm).lwc = vert_profs(indexes_2keep(mm)).lwc;
            % ensemble_profiles(l+mm).lwc_CDP = vert_profs(indexes_2keep(mm)).lwc_CDP;
            % ensemble_profiles(l+mm).lwc_2DC = vert_profs(indexes_2keep(mm)).lwc_2DC;
            % ensemble_profiles(l+mm).total_Nc = vert_profs(indexes_2keep(mm)).total_Nc;
            % ensemble_profiles(l+mm).time = vert_profs(indexes_2keep(mm)).time;
            % ensemble_profiles(l+mm).Nc_vapor = vert_profs(indexes_2keep(mm)).Nc_vapor;
            % ensemble_profiles(l+mm).Nc_vapor_partPres = vert_profs(indexes_2keep(mm)).Nc_vapor_partPres;


        end





    end


end



% save the ensemble profiles
if inputs.sort_for_precip_driz==true

    if inputs.keep_precip_drizzle_profiles==true



        save([foldername,'ensemble_profiles_with_precip_from_',num2str(length(filename)), '_files_LWC-threshold_',...
            num2str(inputs.LWC_threshold), '_Nc-threshold_',...
            num2str(inputs.Nc_threshold), '_',char(datetime("today")),'.mat'],...
            'ensemble_profiles', 'filename')

    else

        save([foldername,'ensemble_profiles_without_precip_from_',num2str(length(filename)), '_files_LWC-threshold_',...
            num2str(inputs.LWC_threshold), '_Nc-threshold_',...
            num2str(inputs.Nc_threshold), '_drizzleLWP-threshold_',...
            num2str(inputs.precip_driz_threshold),'_',char(datetime("today")),'.mat'],...
            'ensemble_profiles', 'filename')

    end

else


    save([foldername,'ensemble_profiles_from_',num2str(length(filename)), '_files_LWC-threshold_',...
        num2str(inputs.LWC_threshold), '_Nc-threshold_',...
        num2str(inputs.Nc_threshold), '_',char(datetime("today")),'.mat'],...
        'ensemble_profiles', 'filename')

end




%% Switch to the struct-of-cell-arrays layout for the analysis below
% The collection loop above builds -- and the save above writes -- a CELL ARRAY
% of per-profile structs (ensemble_profiles{nn}.re). Everything from here down
% was written against the older layout, a struct of cell arrays
% (ensemble_profiles.re{nn}). Convert once here rather than rewriting the
% analysis. The file already written to disk keeps the cell-array layout.
%
% inputs is attached as ensemble_profiles.inputs, which the figure titles read
% for the LWC and Nc thresholds.

ensemble_profiles = ensembleCell2StructOfCells(ensemble_profiles, inputs);


%% Ensemble profiles contains ALL vertical profiles whether they increase or decrease
% Let's split them up into increasing droplet size or decreasing droplet
% size

% Let's step through each profile found
for nn = 1:length(ensemble_profiles.re)

    % Not every profile has a combined CDP + 2DC effective radius: when the
    % 2DC data did not conform, only re_CDP exists. Fall back to it, the same
    % way the optical-depth segmenting further down does.
    if isempty(ensemble_profiles.re{nn})
        re_profile = ensemble_profiles.re_CDP{nn};
    else
        re_profile = ensemble_profiles.re{nn};
    end

    % compute a moving average every 5 data points which will smooth out
    % the small scale variability
    mean_re_profile = movmean(re_profile, 5);

    % Now compute the average derivative of this smooth droplet profile with
    % respect to altitude
    dre_dz = mean(diff(mean_re_profile)./diff(ensemble_profiles.altitude{nn}));

    % sort out profiles between those where droplet size increase with
    % altitude, and those that decrease
    if dre_dz>0

        % grab just the variables of interest for from each vertical profile
        if exist("increasing_profiles", "var")==false

            increasing_profiles.altitude = ensemble_profiles.altitude(nn);
            increasing_profiles.tau = ensemble_profiles.tau(nn);
            increasing_profiles.re = {re_profile};
            increasing_profiles.lwc = ensemble_profiles.lwc(nn);
            increasing_profiles.total_Nc = ensemble_profiles.total_Nc(nn);
            increasing_profiles.time = ensemble_profiles.time(nn);
            % NOTE: the two lines below were copy-paste leftovers from the
            % collection loop above -- they referenced l, mm, vert_profs and
            % indexes_2keep, none of which belong to this loop, and indexing
            % increasing_profiles(l+mm) would have turned it into a struct
            % array. The identical lines are already commented out where they
            % came from. Nothing consumes Nc_vapor downstream.
            % increasing_profiles.Nc_vapor = ensemble_profiles.Nc_vapor(nn);
            % increasing_profiles.Nc_vapor_partPres = ensemble_profiles.Nc_vapor_partPres(nn);

        else

            increasing_profiles.altitude = [increasing_profiles.altitude, ensemble_profiles.altitude(nn)];
            increasing_profiles.tau = [increasing_profiles.tau, ensemble_profiles.tau(nn)];
            increasing_profiles.re = [increasing_profiles.re, {re_profile}];
            increasing_profiles.lwc = [increasing_profiles.lwc, ensemble_profiles.lwc(nn)];
            increasing_profiles.total_Nc = [increasing_profiles.total_Nc, ensemble_profiles.total_Nc(nn)];
            increasing_profiles.time = [increasing_profiles.time, ensemble_profiles.time(nn)];

        end

    elseif dre_dz<0
        % then the profile decreases with increasing altitude, on average

        % grab just the variables of interest for from each vertical profile
        if exist("decreasing_profiles", "var")==false

            decreasing_profiles.altitude = ensemble_profiles.altitude(nn);
            decreasing_profiles.tau = ensemble_profiles.tau(nn);
            decreasing_profiles.re = {re_profile};
            decreasing_profiles.lwc = ensemble_profiles.lwc(nn);
            decreasing_profiles.total_Nc = ensemble_profiles.total_Nc(nn);
            decreasing_profiles.time = ensemble_profiles.time(nn);

        else

            decreasing_profiles.altitude = [decreasing_profiles.altitude, ensemble_profiles.altitude(nn)];
            decreasing_profiles.tau = [decreasing_profiles.tau, ensemble_profiles.tau(nn)];
            decreasing_profiles.re = [decreasing_profiles.re, {re_profile}];
            decreasing_profiles.lwc = [decreasing_profiles.lwc, ensemble_profiles.lwc(nn)];
            decreasing_profiles.total_Nc = [decreasing_profiles.total_Nc, ensemble_profiles.total_Nc(nn)];
            decreasing_profiles.time = [decreasing_profiles.time, ensemble_profiles.time(nn)];

        end

    end



end


% store the LWC threshold

increasing_profiles.lwc_threshold = inputs.LWC_threshold;
decreasing_profiles.lwc_threshold = inputs.LWC_threshold;




%% Plot a histogram of the optical depths and the geometric depths

cloud_optical_depth = zeros(1, length(ensemble_profiles));
cloud_geometric_thickness = zeros(1, length(ensemble_profiles));


for nn = 1:length(ensemble_profiles)

    cloud_optical_depth(nn) = ensemble_profiles{nn}.tau(end);
    cloud_geometric_thickness(nn) = abs(ensemble_profiles{nn}.altitude(1) - ensemble_profiles{nn}.altitude(end));

end


figure;

subplot(1,2,1)
histogram(cloud_geometric_thickness, 'NumBins',20)
xlabel('Cloud Geometric Thickness (m)')
ylabel('Counts')
title([num2str(length(cloud_optical_depth)), ' vertical profiles'])

subplot(1,2,2)
histogram(cloud_optical_depth, 'NumBins',20)
xlabel('Cloud Optical Thickness')
ylabel('Counts')

set(gcf, 'Position', [0 0 1000 550])


%% Plot cloud adiabaticity, Wood (2005)
% Each profile carries an .adiabaticity structure from
% compute_adiabaticity_Wood2005 (attached in
% find_verticalProfiles_VOCALS_REx_ver2). Two estimates are shown:
%   profile method  f_ad(z) = LWC(z) / [Gamma_ad (z - z_cb)]   (Wood Fig. 4c)
%   path method     f_ad = LWP_obs / [(1/2) Gamma_ad h^2]      (Wood Table 3)
% with Gamma_ad evaluated at the mean in-cloud temperature and pressure
% (ATX, PSXC). The plotting function accepts either ensemble layout.
% The returned summary vectors can be used to segment the ensemble by
% adiabaticity, the same way LWP and optical depth are used below.

adiab_summary = plot_ensemble_adiabaticity(ensemble_profiles, 'VOCALS-REx');


%%  Segment re, LWC, Nc into N bins along optical depth

% In order to compute a median vertical profile, we have to first normalize
% the vertical extent so that all profiles lie between values [0,1]. Then
% we break up the vertical component in n discrete bins. Within each bin we
% can compute the mean, median and standard deviation

n_bins = 30; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
vertically_segmented_attributes = cell(n_bins, 3);


normalized_tau = cell(1, length(ensemble_profiles));


for nn = 1:length(ensemble_profiles)

    % first we need to normalize the vertical component of all profiles
    normalized_tau{nn} = ensemble_profiles{nn}.tau./ensemble_profiles{nn}.tau(end);

    % the data is stored in altitude space. If we wish to have the data
    % oriented in optical depth space, we need to reverse the order,
    % because optical depth is meausred from cloud top to cloud bottom
    % so we need to check if the plane is ascending or descending so we
    % know whether the data starts at cloud top or cloud bottom
    % We want all data oriented in tau space. So look for ascending
    % profiles, and flip the vector. If the plane is descending, we don't
    % need to do anything
    if mean(diff(ensemble_profiles{nn}.altitude)./diff(ensemble_profiles{nn}.time))>0
        % if true, then the plane is ascending, grab and flip the variables
        % of interest
        if isfield(ensemble_profiles{nn}, 're')==false
            re = fliplr(ensemble_profiles{nn}.re_CDP);

        else
            re = fliplr(ensemble_profiles{nn}.re);
        end

        lwc = fliplr(ensemble_profiles{nn}.lwc);
        Nc = fliplr(ensemble_profiles{nn}.total_Nc);

    else
        % if false, then the plane is descending, just grab the variables
        % of interest
        if isfield(ensemble_profiles{nn}, 're')==false
            re = ensemble_profiles{nn}.re_CDP;
        else
            re = ensemble_profiles{nn}.re;
        end

        lwc = ensemble_profiles{nn}.lwc;
        Nc = ensemble_profiles{nn}.total_Nc;

    end





    % for each profile, we need to segment the variables of interest into n
    % bins.

    for bb = 1:length(bin_edges)-1

        % grab all re, LWC, and Nc values within each bin. Segment them
        % accordingly
        if bb==1
            index_segment = normalized_tau{nn}>=bin_edges(bb) & normalized_tau{nn}<=bin_edges(bb+1);

        else
            index_segment = normalized_tau{nn}>bin_edges(bb) & normalized_tau{nn}<=bin_edges(bb+1);
        end

        % store the effective radius values
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}, re(index_segment)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}, lwc(index_segment)];

        % store the total droplet number concentration values
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}, Nc(index_segment)];



    end



end



%% Make errorbar plot of the mean and standard deviation of effective radius at each vertical bin


re_mean = zeros(n_bins, 1);
re_std = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the mean value for the current bin
    re_mean(bb) = mean(vertically_segmented_attributes{bb,1});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    re_std(bb) = std(vertically_segmented_attributes{bb,1});         % microns - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(re_mean, bin_center, re_std, re_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Mean r_e')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])


%% Make errorbar plot of the median and standard deviation of effective radius at each vertical bin


re_median = zeros(n_bins, 1);
re_std = zeros(n_bins, 1);
re_lowP = zeros(n_bins, 1);
re_highP = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);


lowP = 10;
highP = 90;

for bb = 1:n_bins

    % compute the mean value for the current bin
    re_median(bb) = median(vertically_segmented_attributes{bb,1});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    re_std(bb) = std(vertically_segmented_attributes{bb,1});         % microns - standard deviation

    % compute the 10th percentile effective radius from the distribution at
    % this level
    re_lowP(bb) = prctile(vertically_segmented_attributes{bb, 1}, lowP);         % microns - standard deviation

    % compute the 90th percentile effective radius from the distribution at
    % this level
    re_highP(bb) = prctile(vertically_segmented_attributes{bb, 1}, highP);         % microns - standard deviation



    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(re_median, bin_center, re_std, re_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Median r_e')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])


%% Plot the mean and median effective radius profiles

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_mean-re_std; flipud(re_mean + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(re_mean, bin_center, 'Color', mySavedColors(1, 'fixed'));
hold on

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_median-re_std; flipud(re_median + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(re_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<r_e(z)>$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])



%% Plot the mean effective radius profile with a shaded area representing the standard deviation

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_mean-re_std; flipud(re_mean + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(re_mean, bin_center, 'Color', mySavedColors(1, 'fixed'));

% plot a droplet profile fit over to see if it is adiabatic or not
profile_type = 'subadiabatic_aloft';
re_fit = create_droplet_profile2([re_mean(1), re_mean(end)], bin_center, 'optical_depth', profile_type);

hold on
plot(re_fit, bin_center, 'Color', 'k', 'LineWidth',1);

legend('','', 'subadiabatic', 'location', 'best')


set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('<r_e>')
ylabel('Normalized Optical Depth')
title('Mean r_e')

set(gcf, 'Position', [0 0 650 600])


%% Make errorbar plot of the mean and standard deviation of liquid water content at each vertical bin


lwc_mean = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);
lwc_lowP = zeros(n_bins, 1);
lwc_highP = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the mean value for the current bin
    lwc_mean(bb) = mean(vertically_segmented_attributes{bb,2});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % microns - standard deviation

    % compute the 10th percentile from the distribution at
    % this level
    lwc_lowP(bb) = prctile(vertically_segmented_attributes{bb, 2}, lowP);         % microns - standard deviation

    % compute the 90th percentile from the distribution at
    % this level
    lwc_highP(bb) = prctile(vertically_segmented_attributes{bb, 2}, highP);         % microns - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(lwc_mean, bin_center, lwc_std, lwc_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Mean LWC')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])



%% Plot the mean liquid water content with a shaded area representing the standard deviation

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_mean-lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(3,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(lwc_mean, bin_center, 'Color', mySavedColors(3, 'fixed'));



set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('<LWC(z)>')
ylabel('Normalized Optical Depth')
title('Mean LWC')

set(gcf, 'Position', [0 0 650 600])


%% Make errorbar plot of the median and standard deviation of liquid water content at each vertical bin


lwc_median = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the mean value for the current bin
    lwc_median(bb) = median(vertically_segmented_attributes{bb,2});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % microns - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(lwc_median, bin_center, lwc_std, lwc_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Median LWC')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])



%% Plot the median liquid water content with a shaded area representing the standard deviation

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(4,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median profile
plot(lwc_median, bin_center, 'Color', mySavedColors(4, 'fixed'));

% Plot a line between the first and last value of LWC
hold on
plot([lwc_median(2), lwc_median(end)], [bin_center(2), bin_center(end)], 'k', 'Linewidth', 1)

legend('','','Linear Fit', 'Location', 'best')



set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('<LWC(z)>')
ylabel('Normalized Optical Depth')
title('Median LWC')

set(gcf, 'Position', [0 0 650 600])




%% Plot the mean and median Liquid Water content profiles

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_mean-lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(6,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mead profile
plot(lwc_mean, bin_center, 'Color', mySavedColors(6, 'fixed'));
hold on

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(4,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(lwc_median, bin_center, 'Color', mySavedColors(4, 'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<LWC(z)>$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])




%% Plot the mean and median Number Concentration profiles



Nc_median = zeros(n_bins, 1);
Nc_mean = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);
Nc_lowP = zeros(n_bins, 1);
Nc_highP = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the median value for the current bin
    Nc_median(bb) = median(vertically_segmented_attributes{bb,3});       % microns - mean droplet radius value

    % compute the mean value for the current bin
    Nc_mean(bb) = mean(vertically_segmented_attributes{bb,3});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    Nc_std(bb) = std(vertically_segmented_attributes{bb,3});         % microns - standard deviation

    % compute the 10th percentile from the distribution at
    % this level
    Nc_lowP(bb) = prctile(vertically_segmented_attributes{bb, 3}, lowP);         % microns - standard deviation

    % compute the 90th percentile from the distribution at
    % this level
    Nc_highP(bb) = prctile(vertically_segmented_attributes{bb, 3}, highP);         % microns - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end


figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [Nc_mean-Nc_std; flipud(Nc_mean + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(Nc_mean, bin_center, 'Color', mySavedColors(1, 'fixed'));
hold on

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [Nc_median-Nc_std; flipud(Nc_median + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(Nc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])


%% Make a subplot showing the median profile of droplet size, liquid water content and number concentration

figure;

% plot the median effective radius
subplot(1,3,1)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
% x = [re_median-re_std; flipud(re_median + re_std)];
x = [re_lowP; flipud(re_highP)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(re_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')


% plot the median liquid water content profile
subplot(1,3,2)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
% x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
x = [lwc_lowP; flipud(lwc_highP)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(lwc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

% set the figure title
title(['Median Profiles:  $LWC \geq$', num2str(ensemble_profiles{1}.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles{1}.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex')



% plot the median droplet number concentration
subplot(1,3,3)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
% x = [Nc_median-Nc_std; flipud(Nc_median + Nc_std)];
x = [Nc_lowP; flipud(Nc_highP)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(Nc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

% set the size of the figure
set(gcf, 'Position', [0 0 1200 625])


%% What is the average difference between cloud top and cloud bottom radii?

% I've broken up all of the data into vertical segments.
% What is the average radius value for the bottom most segment?
re_top_median = mean(vertically_segmented_attributes{1,1});
re_bottom_median = median(vertically_segmented_attributes{end,1});

% NOTE: nothing downstream consumes diff_top_bottom_1/2 at present. The lines
% below were missing their profile index entirely, and indexed the profile in
% SAMPLING order, which puts index 1 at cloud base for ascending profiles. They
% are oriented cloud top first here, to match the name and the rest of this
% script.

n_profiles = length(ensemble_profiles.re);
diff_top_bottom_1 = zeros(1, n_profiles);
diff_top_bottom_2 = zeros(1, n_profiles);

for nn = 1:n_profiles

    % combined CDP + 2DC effective radius, or CDP alone when 2DC did not conform
    if isempty(ensemble_profiles.re{nn})
        re_profile = ensemble_profiles.re_CDP{nn};
    else
        re_profile = ensemble_profiles.re{nn};
    end

    % orient cloud top first
    if mean(diff(ensemble_profiles.altitude{nn}) ./ diff(ensemble_profiles.time{nn})) > 0
        re_profile = fliplr(reshape(re_profile, 1, []));
    else
        re_profile = reshape(re_profile, 1, []);
    end

    diff_top_bottom_1(nn) = re_profile(1) - re_profile(end);
    diff_top_bottom_2(nn) = re_profile(2) - re_profile(end);

end



%% Look at vertical profiles with different LWP regimes


% First plot all of LWP to determine what regimes to plot

LWP = zeros(1, length(ensemble_profiles.lwc));

for nn = 1:length(ensemble_profiles.lwc)

    LWP(nn) = abs(trapz(ensemble_profiles.altitude{nn}, ensemble_profiles.lwc{nn}));

end

figure; plot(LWP, '.')
grid on; grid minor
xlabel('Profile Index')
ylabel('LWP (g/m^2)')


% Let's find the profiles with a total LWP < 20 g/m^2

index_20 = LWP<=20;

plot_median_ensemble_LWC_re_NC_vs_norm_altitude(ensemble_profiles, index_20)

% Create textbox
annotation('textbox',[0.866 0.9376 0.1115 0.0527],'String',{'LWP<20 g/m^2'},...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');


% Let's find the profiles with a total LWP >20 and <50 g/m^2

index_50 = LWP>20 & LWP<=50;

plot_median_ensemble_LWC_re_NC_vs_norm_altitude(ensemble_profiles, index_50)

annotation('textbox',[0.866 0.9376 0.1115 0.0527],'String',{'20<LWP<50 g/m^2'},...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');



% Let's find the profiles with a total LWP >50 and <100 g/m^2

index_100 = LWP>50 & LWP<=100;

plot_median_ensemble_LWC_re_NC_vs_norm_altitude(ensemble_profiles, index_100)

annotation('textbox',[0.866 0.9376 0.1115 0.0527],'String',{'50<LWP<100 g/m^2'},...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');


%% Plot the median droplet size profile for several different LWP regimes


% Compute the LWP for each profile
LWP = zeros(1, length(ensemble_profiles.lwc));

for nn = 1:length(ensemble_profiles.lwc)

    LWP(nn) = abs(trapz(ensemble_profiles.altitude{nn}, ensemble_profiles.lwc{nn}));

end


% Define the LWP regimes
% LWP_bin_edges = [22, 50, 90, 140];      % g/m^2
LWP_bin_edges = [0, 20, 40, 62, 100, 140];      % g/m^2

index_regime = false(length(LWP_bin_edges), length(LWP));
legend_str = cell(1, length(LWP_bin_edges));

for nn = 1:length(LWP_bin_edges)


    if nn>=1 && nn<length(LWP_bin_edges)

        index_regime(nn,:) = LWP>=LWP_bin_edges(nn) & LWP<LWP_bin_edges(nn+1);

        % Define the legend string
        legend_str{nn} = [num2str(LWP_bin_edges(nn)),'$\leq LWP<$', num2str(LWP_bin_edges(nn+1))];

    else
        index_regime(nn,:) = LWP>LWP_bin_edges(nn);

        % define the legend string
        legend_str{nn} = ['$LWP>$', num2str(LWP_bin_edges(nn)), ' $g/m^{2}$'];

    end

end


% define the changing variable for the plot
changing_variable = 'LWP';


% Combine all these indices and plot the indexed median profiles
% NOTE: plot_multiple_median_ensemble_LWC_re_NC_vs_altitude does not exist -- it
% has never been written, here or anywhere else in the repo. Only the optical
% depth version, plot_multiple_median_ensemble_LWC_re_NC_vs_norm_optical_depth,
% exists, and it is already called with these same regimes. This call is
% commented out so the script runs; write the altitude version and uncomment it
% if you want that figure.
% plot_multiple_median_ensemble_LWC_re_NC_vs_altitude(ensemble_profiles, index_regime,...
%     changing_variable, legend_str)




%% Plot the median droplet size profile for several different optical depth regimes


% find the cloud optical depth for each profile
tau = zeros(1, length(ensemble_profiles.lwc));

for nn = 1:length(ensemble_profiles.tau)

    tau(nn) = ensemble_profiles.tau{nn}(end);

end


% Define the optical depth regimes
tau_bin_edges = [0, 5, 9, 14, 21, 30];      % g/m^2

index_regime = false(length(tau_bin_edges), length(tau));
legend_str = cell(1, length(tau_bin_edges));

for nn = 1:length(tau_bin_edges)


    if nn>=1 && nn<length(tau_bin_edges)

        index_regime(nn,:) = tau>=tau_bin_edges(nn) & tau<tau_bin_edges(nn+1);

        % Define the legend string
        legend_str{nn} = [num2str(tau_bin_edges(nn)),'$\leq \tau <$', num2str(tau_bin_edges(nn+1))];

    else
        index_regime(nn,:) = tau>tau_bin_edges(nn);

        % define the legend string
        legend_str{nn} = ['$\tau >$', num2str(tau_bin_edges(nn))];

    end

end


% define the changing variable for the plot
changing_variable = '$\tau$';


% Combine all these indices and plot the indexed median profiles
% NOTE: plot_multiple_median_ensemble_LWC_re_NC_vs_altitude does not exist -- it
% has never been written, here or anywhere else in the repo. Only the optical
% depth version, plot_multiple_median_ensemble_LWC_re_NC_vs_norm_optical_depth,
% exists, and it is already called with these same regimes. This call is
% commented out so the script runs; write the altitude version and uncomment it
% if you want that figure.
% plot_multiple_median_ensemble_LWC_re_NC_vs_altitude(ensemble_profiles, index_regime,...
%     changing_variable, legend_str)
