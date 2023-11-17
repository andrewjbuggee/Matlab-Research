%% Plot the medain profile for LWC, re, and Nc versus normalized altitude

% only use the profiles defined by the logical index provided


function plot_median_ensemble_LWC_re_NC_vs_norm_altitude(ensemble_profiles, index)

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


normalized_altitude = cell(1, length(index));


for nn = 1:length(index)

    % first we need to normalize the vertical component of all profiles
    normalized_altitude{nn} = (ensemble_profiles.altitude{nn} - min(ensemble_profiles.altitude{nn}))./...
                               (max(ensemble_profiles.altitude{nn}) - min(ensemble_profiles.altitude{nn}));

    % the data is stored in altitude space. 

    re = ensemble_profiles.re{nn};
    lwc = ensemble_profiles.lwc{nn};
    Nc = ensemble_profiles.Nc{nn};



    % for each profile, we need to segment the variables of interest into n
    % bins.

    for bb = 1:length(bin_edges)-1

        % grab all re, LWC, and Nc values within each bin. Segment them
        % accordingly
        if bb==1
            index_segment = normalized_altitude{nn}>=bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);

        else
            index_segment = normalized_altitude{nn}>bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);
        end

        % store the effective radius values
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; re(index_segment)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; lwc(index_segment)];

        % store the total droplet number concentration values
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; Nc(index_segment)];



    end



end


%% Create a PDF object at each level in the cloud and fit a distribution to this PDF

% store the refection of each null hypothesis and the p-value for each
% chi-squared test

re_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
re_p_normal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
re_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
re_p_gamma = zeros(1, size(vertically_segmented_attributes,1));


lwc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));



Nc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));

for bb = 1:size(vertically_segmented_attributes, 1)


    % -----------------------------------------------
    % ------- EFFECTIVE DROPLET RADIUS FITTING ------
    % -----------------------------------------------


    % fit the effective radius data to a normal distribution
    re_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'normal');
    [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb));

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'lognormal');
    [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb));

    % fit the effective radius data to a gamma distribution
    re_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'gamma');
    [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(re_fit_normal(bb)); title('Normal Fit'); xlabel('r_e (\mum)')
%     subplot(1,3,2); plot(re_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('r_e (\mum)')
%     subplot(1,3,3); plot(re_fit_gamma(bb)); title('Gamma Fit'); xlabel('r_e (\mum)')
%     set(gcf, 'Position', [0 0 1200 500])




    % -------------------------------------------
    % ------- LIQUID WATER CONTENT FITTING ------
    % -------------------------------------------


    % fit the liquid water content data to a normal distribution
    lwc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'normal');
    [lwc_reject_normal(bb), lwc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_normal(bb));

    % fit the liquid water content data to a log-normal distribution
    lwc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'lognormal');
    [lwc_reject_lognormal(bb), lwc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_lognormal(bb));

    % fit the liquid water content data to a gamma distribution
    lwc_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'gamma');
    [lwc_reject_gamma(bb), lwc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF', lwc_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(lwc_fit_normal(bb)); title('Normal Fit'); xlabel('LWC (g/m^{2})')
%     subplot(1,3,2); plot(lwc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (g/m^{2})')
%     subplot(1,3,3); plot(lwc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (g/m^{2})')
%     set(gcf, 'Position', [0 0 1200 500])




    % -------------------------------------------
    % ------- NUMBER CONCENTRATION FITTING ------
    % -------------------------------------------


    % fit the number concentration data to a normal distribution
    Nc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'normal');
    [Nc_reject_normal(bb), Nc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_normal(bb));

    % fit the number concentration content data to a log-normal distribution
    Nc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'lognormal');
    [Nc_reject_lognormal(bb), Nc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_lognormal(bb));

    % fit the number concentration content data to a gamma distribution
    Nc_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'gamma');
    [Nc_reject_gamma(bb), Nc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF', Nc_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(Nc_fit_normal(bb)); title('Normal Fit'); xlabel('N_c (cm^{-3})')
%     subplot(1,3,2); plot(Nc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (cm^{-3})')
%     subplot(1,3,3); plot(Nc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (cm^{-3})')
%     set(gcf, 'Position', [0 0 1200 500])
    

end


% Now let's find the where the hypothesis was not rejected (reject_ = 0)
% which means the chi-squared test is confident in the choice of
% distribution to within 5% uncertainty

bin_names = {'Normal', 'Log-Normal', 'Gamma'};
% -----------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS FITTING ------
% -----------------------------------------------
 [max__re_p, idx_re_p] = max([re_p_normal; re_p_lognormal; re_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_re_p==1), sum(idx_re_p==2), sum(idx_re_p==3)]); 
% title('r_e best distribution fit'); ylabel('Counts')



% -------------------------------------------
% ------- LIQUID WATER CONTENT FITTING ------
% -------------------------------------------
[max__lwc_p, idx_lwc_p] = max([lwc_p_normal; lwc_p_lognormal; lwc_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_lwc_p==1), sum(idx_lwc_p==2), sum(idx_lwc_p==3)]); 
% title('LWC best distribution fit'); ylabel('Counts')


% -------------------------------------------
% ------- NUMBER CONCENTRATION FITTING ------
% -------------------------------------------

[max__Nc_p, idx_Nc_p] = max([Nc_p_normal; Nc_p_lognormal; Nc_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_Nc_p==1), sum(idx_Nc_p==2), sum(idx_Nc_p==3)]); 
% title('N_c best distribution fit'); ylabel('Counts')


%%  Compute the Median LWC, re, and Nc of each layer

% ---- most common best fit distribution for r_e was is the log-normal dist ---
re_median = zeros(n_bins, 1);
re_std = zeros(n_bins, 1);
re_avg_deviation_from_median = zeros(n_bins, 1);
re_logNormal_std = zeros(n_bins, 1);
re_logNormal_mean = zeros(n_bins, 1);


% ---- most common best fit distribution for LWC was is the normal dist ---
lwc_median = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);


% ---- most common best fit distribution for N_c was is the gamma dist ---
Nc_median = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);

bin_center = zeros(n_bins, 1);



for bb = 1:n_bins
    

    % ----- COMPUTE STATISTICS FOR DROPLET SIZE -----

    % compute the median value for the current bin
    re_median(bb) = median(vertically_segmented_attributes{bb,1});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    re_std(bb) = std(vertically_segmented_attributes{bb,1});         % microns - standard deviation

    % The std is the average deviation from the mean value. Compute the
    % average deviation from the median value
    re_avg_deviation_from_median(bb) = sqrt(mean((vertically_segmented_attributes{bb,1} - re_median(bb)).^2));

    % find the mean of the log normal distribution
    re_logNormal_mean(bb) = exp(re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2 /2);

    % find squareroot of the variance of the lognormal distribution
    re_logNormal_std(bb) = sqrt(exp(2*re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2)*(exp(re_fit_lognormal(bb).sigma^2) - 1));


    % ----- COMPUTE STATISTICS FOR LIQUID WATER CONTENT -----

    % compute the median value for the current bin
    lwc_median(bb) = median(vertically_segmented_attributes{bb,2});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % microns - standard deviation




    % ----- COMPUTE STATISTICS FOR DROPLET NUMBER CONCENTRATION -----

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
% x = [re_logNormal_mean - re_logNormal_std; flipud(re_logNormal_mean + re_logNormal_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(re_median, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')


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


grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter','latex')

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


grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set the size of the figure
set(gcf, 'Position', [0 0 1250 625])








end

