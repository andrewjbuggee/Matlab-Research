%% Plot the medain profile for LWC, re, and Nc versus normalized optical depth

% only use the profiles defined by the logical index provided


function plot_median_ensemble_LWC_re_NC_vs_norm_optical_depth(ensemble_profiles, index)

% delete the profiles that are not desired for this plot

if length(index)<length(ensemble_profiles.re)

    indices2delete = setxor(1:length(ensemble_profiles.lwc), index);

    ensemble_profiles.altitude(indices2delete) = [];
    ensemble_profiles.tau(indices2delete) = [];
    ensemble_profiles.re(indices2delete) = [];
    ensemble_profiles.lwc(indices2delete) = [];
    ensemble_profiles.Nc(indices2delete) = [];
    ensemble_profiles.time(indices2delete) = [];

end


    


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


normalized_tau = cell(1, length(index));


for nn = 1:length(index)

    % first we need to normalize the vertical component of all profiles
    normalized_tau{nn} = ensemble_profiles.tau{nn}./ensemble_profiles.tau{nn}(end);

    % the data is stored in altitude space. If we wish to have the data
    % oriented in optical depth space, we need to reverse the order,
    % because optical depth is meausred from cloud top to cloud bottom
    % so we need to check if the plane is ascending or descending so we
    % know whether the data starts at cloud top or cloud bottom
    % We want all data oriented in tau space. So look for ascending
    % profiles, and flip the vector. If the plane is descending, we don't
    % need to do anything
    if mean(diff(ensemble_profiles.altitude{nn})./diff(ensemble_profiles.time{nn}))>0
        % if true, then the plane is ascending, grab and flip the variables
        % of interest
        re = flipud(ensemble_profiles.re{nn});
        lwc = flipud(ensemble_profiles.lwc{nn});
        Nc = flipud(ensemble_profiles.Nc{nn});

    else
        % if false, then the plane is descending, just grab the variables
        % of interest
        re = ensemble_profiles.re{nn};
        lwc = ensemble_profiles.lwc{nn};
        Nc = ensemble_profiles.Nc{nn};

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
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; re(index_segment)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; lwc(index_segment)];

        % store the total droplet number concentration values
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; Nc(index_segment)];



    end



end


%%  Compute the Median LWC, re, and Nc of each layer

re_median = zeros(n_bins, 1);
re_std = zeros(n_bins, 1);
re_avg_deviation_from_median = zeros(n_bins, 1);
re_min = zeros(n_bins, 1);

lwc_median = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);

Nc_median = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);

bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the median value for the current bin
    re_median(bb) = median(vertically_segmented_attributes{bb,1});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    re_std(bb) = std(vertically_segmented_attributes{bb,1});         % microns - standard deviation

    % The std is the average deviation from the mean value. Compute the
    % average deviation from the median value
    re_avg_deviation_from_median(bb) = sqrt(mean((vertically_segmented_attributes{bb,1} - re_median(bb)).^2));

    % find the smallest droplet size value at each bin
    re_min(bb) = min(vertically_segmented_attributes{bb,1});

    % compute the median value for the current bin
    lwc_median(bb) = median(vertically_segmented_attributes{bb,2});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % microns - standard deviation

    % compute the median value for the current bin
    Nc_median(bb) = median(vertically_segmented_attributes{bb,3});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    Nc_std(bb) = std(vertically_segmented_attributes{bb,3});         % microns - standard deviation


    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end


%% Make a subplot of all 3 median profiles


figure; 

% plot the median effective radius
subplot(1,3,1)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
% x = [re_min; flipud(re_median + re_avg_deviation_from_median)];
x = [re_median - re_avg_deviation_from_median; flipud(re_median + re_avg_deviation_from_median)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(re_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Optical Depth', 'Interpreter', 'latex')


% plot the median liquid water content profile
subplot(1,3,2)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(lwc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth', 'Interpreter','latex')

% set the figure title
title(['Median Profiles:  $LWC \geq$', num2str(ensemble_profiles.inputs.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles.inputs.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex')



% plot the median droplet number concentration
subplot(1,3,3)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
x = [Nc_median-Nc_std; flipud(Nc_median + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(Nc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth', 'Interpreter', 'latex')

% set the size of the figure
set(gcf, 'Position', [0 0 1200 625])








end

