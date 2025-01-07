%% This script creates all figures used in my first paper


% By Andrew John Buggee



%% FIGURE 1 - Plot the ensemble MEAN of droplet size, liquid water content and
% number concentration for non-precipitating clouds. Add an adiabatic
% profile for the liquid water content and effective radius to show the
% mean profiles are close to adiabatic, supporting my assumption.
% ----- For Vertical Profiles ----


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    % --- non-precip profiles only, LWC>0.03, Nc>1  ----
    load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
        '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_19-Sep-2023'])


    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
    %     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
    %         '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.005, Nc>1  ----
    %     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023.mat'])

elseif strcmp(whatComputer, 'andrewbuggee')==true

    % --- non-precip profiles only, LWC>0.005, Nc>1 ----
    % load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %  '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_19-Sep-2023'])

    % --- all profiles, LWC>0.005, Nc>1 ----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
        '/ensemble_profiles_from_14_files_LWC-threshold_0.005_Nc-threshold_1_14-Sep-2023'])


end




% using the mean ensemble function to plot the mean vertical profile of
% the ensemble




% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Segment re, LWC, Nc into N bins along optical depth

% In order to compute a mean vertical profile, we have to first normalize
% the vertical extent so that all profiles lie between values [0,1]. Then
% we break up the vertical component in n discrete bins. Within each bin we
% can compute the mean, median and standard deviation

n_bins = 30; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
vertically_segmented_attributes = cell(n_bins, 3);


normalized_altitude = cell(1, length(ensemble_profiles.lwc));


for nn = 1:length(ensemble_profiles.lwc)

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



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Create a PDF object at each level in the cloud and fit a distribution to this PDF

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



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Compute the mean LWC, re, and Nc of each layer

% ---- most common best fit distribution for r_e was is the log-normal dist ---
re_logNormal_std = zeros(n_bins, 1);
re_logNormal_mean = zeros(n_bins, 1);


% ---- most common best fit distribution for LWC was is the normal dist ---
lwc_mean = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);


% ---- most common best fit distribution for N_c was is the normal dist ---
Nc_mean = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);

bin_center = zeros(n_bins, 1);




for bb = 1:n_bins


    % ----- COMPUTE STATISTICS FOR DROPLET SIZE -----

    % find the mean of the log normal distribution
    re_logNormal_mean(bb) = exp(re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2 /2);

    % find squareroot of the variance of the lognormal distribution
    re_logNormal_std(bb) = sqrt(exp(2*re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2)*(exp(re_fit_lognormal(bb).sigma^2) - 1));



    % ----- COMPUTE STATISTICS FOR LIQUID WATER CONTENT -----

    % compute the mean value for the current bin
    % the mean of the distribution (the standard way of computing the expected value)
    % is also the mean of the normal distribution. They are identical.
    lwc_mean(bb) = mean(vertically_segmented_attributes{bb,2});       % g/cm^3 - mean liqiud water content
    %lwc_mean(bb) = lwc_fit_gamma(bb).mean;       % g/cm^3 - mean liqiud water content

    % compute the standard deviation of the current bin
    % the std of the distribution (the standard way of computing the squareroot of the variance)
    % is also the std of the normal distribution. They are identical.
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % g/cm^3 - standard deviation
    %lwc_std(bb) = lwc_fit_gamma(bb).std;         % g/cm^3 - standard deviation



    % ----- COMPUTE STATISTICS FOR DROPLET NUMBER CONCENTRATION -----

    % compute the mean value for the current bin
    %Nc_mean(bb) = Nc_fit_gamma(bb).mean;       % cm^(-3) - mean number concentration
    Nc_mean(bb) = Nc_fit_normal(bb).mean;       % cm^(-3) - mean number concentration

    % compute the standard deviation of the current bin
    %Nc_std(bb) = Nc_fit_gamma(bb).std;         % cm^(-3) - standard deviation
    Nc_std(bb) = Nc_fit_normal(bb).std;         % cm^(-3) - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end



% ----------------------------------------------------------------------
% -------------- Make a subplot of all 3 mean profiles -----------------
% ----------------------------------------------------------------------




figure;

% plot the mean effective radius
subplot(1,3,1)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_logNormal_mean - re_logNormal_std; flipud(re_logNormal_mean + re_logNormal_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(re_logNormal_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([4, 10])                   % microns



% plot the mean liquid water content profile
subplot(1,3,2)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_mean-lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(lwc_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter','latex')

% set x axis boundaries
xlim([0, 0.6])                   % g/cm^3

% set the figure title
title(['Mean Profiles:  $LWC \geq$', num2str(ensemble_profiles.inputs.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles.inputs.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex')





% plot the mean droplet number concentration
subplot(1,3,3)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [Nc_mean-Nc_std; flipud(Nc_mean + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(Nc_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([0, 320])                   % cm^(-3)

% set the size of the figure
set(gcf, 'Position', [0 0 1255 625])




% ----------------------------------------------------------------------
% ----------------------- Adiabatic Curve Fits -------------------------
% ----------------------------------------------------------------------

% ------------------------- LIQUID WATER CONTENT ---------------------------
% ----- Fit an adiabatic curve to the mean liquid water content profile -----
nudge_from_top =2;
lwc_slope = (lwc_mean(end-nudge_from_top) - lwc_mean(1))/(bin_center(end-nudge_from_top) - bin_center(1));
lwc_intercept = lwc_mean(1) - lwc_slope*bin_center(1);
lwc_adiabatic_fit = lwc_slope*bin_center + lwc_intercept;

% add to subplot(1,3,1)
subplot(1,3,2); hold on
plot(lwc_adiabatic_fit(1:end-nudge_from_top), bin_center(1:end-nudge_from_top), 'k', 'LineWidth',2)

% -- Include a legend in the lower right-hand corner of the 2nd subplot --
legend({'Standard Deviation', 'Mean Profile', 'Adiabatic Fit'}, 'Interpreter','latex',...
    'Location','southeast', 'FontSize', 17)



% -------------------- EFFECTIVE DROPLET RADIUS -----------------------
% Plot an adiabatic curve fit to the mean droplet radius profile
nudge_from_bottom = 3;

% ----- Fit an adiabatic curve to the mean droplet profile -----
% use droplet profile function to create adiabatic fit
re_adiabatic_fit = create_droplet_profile2([re_logNormal_mean(end), re_logNormal_mean(1 + nudge_from_bottom)],...
    bin_center(1+nudge_from_bottom:end),'altitude', 'adiabatic');

% derive droplet profile from adiabatic fit of LWC


% add to subplot(1,3,1)
subplot(1,3,1); hold on
plot(re_adiabatic_fit, bin_center(1+nudge_from_bottom:end), 'k', 'LineWidth',2)





% ---------------------- DROPLET NUMBER CONCENTRATION ---------------------------
% ----- Fit an adiabatic curve to Number Concentration based on Adiabatic fits above -----
% Nc_adiabatic_fit = 3*(lwc_adiabatic_fit./1e6)./...
%                 (4*pi*(re_adiabatic_fit*1e-4).^3);




% --- Create a Textbox stating these are only non-precipitating clouds ---
% Create textbox
annotation('textbox',[0.134 0.802 0.142 0.114],...
    'String',{'Non-Precipitating clouds only ($LWP_{2DC}<1 \,g/m^{2}$)'},...
    'LineWidth',2,...
    'Interpreter','latex',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');





%% FIGURE 2 - Plot the ensemble MEDIAN of droplet size, liquid water content and
% number concentration for non-precipitating clouds. Add an adiabatic
% profile for the liquid water content and effective radius to show the
% mean profiles are close to adiabatic, supporting my assumption


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    % --- non-precip profiles only, LWC>0.03, Nc>1  ----
    %     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
    %         '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_19-Sep-2023'])


    % ***** Using new VOCALS READ function with LWC adjustment *****
    % --- non-precip profiles only, LWC>0.03, Nc>1  ----
    load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
        '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_05-Nov-2023.mat'])


    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
    %     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
    %         '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.005, Nc>1  ----
    %     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023.mat'])

elseif strcmp(whatComputer, 'andrewbuggee')==true

    % --- non-precip profiles only, LWC>0.005, Nc>1 ----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    load(['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
        'vocals_rex_data/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_05-Nov-2023'])

    % --- all profiles, LWC>0.005, Nc>1 ----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/ensemble_profiles_from_14_files_LWC-threshold_0.005_Nc-threshold_1_14-Sep-2023'])


end




% using the median ensemble function to plot the median vertical profile of
% the ensemble




% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Segment re, LWC, Nc into N bins along optical depth

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


normalized_altitude = cell(1, length(ensemble_profiles.lwc));


for nn = 1:length(ensemble_profiles.lwc)

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



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Create a PDF object at each level in the cloud and fit a distribution to this PDF

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

figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_re_p==1), sum(idx_re_p==2), sum(idx_re_p==3)]);
title('r_e best distribution fit'); ylabel('Counts')



% -------------------------------------------
% ------- LIQUID WATER CONTENT FITTING ------
% -------------------------------------------
[max__lwc_p, idx_lwc_p] = max([lwc_p_normal; lwc_p_lognormal; lwc_p_gamma],[], 1);

figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_lwc_p==1), sum(idx_lwc_p==2), sum(idx_lwc_p==3)]);
title('LWC best distribution fit'); ylabel('Counts')


% -------------------------------------------
% ------- NUMBER CONCENTRATION FITTING ------
% -------------------------------------------

[max__Nc_p, idx_Nc_p] = max([Nc_p_normal; Nc_p_lognormal; Nc_p_gamma],[], 1);

figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_Nc_p==1), sum(idx_Nc_p==2), sum(idx_Nc_p==3)]);
title('N_c best distribution fit'); ylabel('Counts')



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Compute the Median LWC, re, and Nc of each layer

% ---- most common best fit distribution for r_e was is the log-normal dist ---
re_logNormal_std = zeros(n_bins, 1);
re_logNormal_median = zeros(n_bins, 1);
re_custom_logNormal_std_larger = zeros(n_bins, 1);
re_custom_logNormal_std_smaller = zeros(n_bins, 1);

% ---- most common best fit distribution for LWC was is the normal dist ---
% lwc_median = zeros(n_bins, 1);
% lwc_std = zeros(n_bins, 1);
lwc_logNormal_std = zeros(n_bins, 1);
lwc_logNormal_median = zeros(n_bins, 1);
lwc_custom_logNormal_std_larger = zeros(n_bins, 1);
lwc_custom_logNormal_std_smaller = zeros(n_bins, 1);

% ---- most common best fit distribution for N_c was is the normal dist ---
Nc_median = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);

bin_center = zeros(n_bins, 1);




for bb = 1:n_bins


    % ----- COMPUTE STATISTICS FOR DROPLET SIZE -----

    % find the mean of the log normal distribution
    re_logNormal_median(bb) = re_fit_lognormal(bb).median;

    % find squareroot of the variance of the lognormal distribution
    re_logNormal_std(bb) = sqrt(exp(2*re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2)*(exp(re_fit_lognormal(bb).sigma^2) - 1));

    % Let's also compute the average deviation from the median value when
    % radii are larger and smaller than the median.
    % For radii larger than the median...
    idx_larger = vertically_segmented_attributes{bb,1}>re_logNormal_median(bb);
    re_larger = vertically_segmented_attributes{bb,1}(idx_larger);
    re_custom_logNormal_std_larger(bb) = sqrt(mean((re_larger - re_logNormal_median(bb)).^2));
    %re_custom_logNormal_std_larger(bb) = mean(re_larger - re_logNormal_median(bb));

    % For radii smaller than the median...
    idx_smaller = vertically_segmented_attributes{bb,1}<re_logNormal_median(bb);
    re_smaller = vertically_segmented_attributes{bb,1}(idx_smaller);
    re_custom_logNormal_std_smaller(bb) = sqrt(mean((re_smaller - re_logNormal_median(bb)).^2));
    %re_custom_logNormal_std_smaller(bb) = mean(re_logNormal_median(bb) - re_smaller);



    % ----- COMPUTE STATISTICS FOR LIQUID WATER CONTENT -----

    % find the mean of the log normal distribution
    lwc_logNormal_median(bb) = lwc_fit_lognormal(bb).median;

    % find squareroot of the variance of the lognormal distribution
    lwc_logNormal_std(bb) = sqrt(exp(2*lwc_fit_lognormal(bb).mu + lwc_fit_lognormal(bb).sigma^2)*(exp(lwc_fit_lognormal(bb).sigma^2) - 1));

    % Let's also compute the average deviation from the median value when
    % lwc is larger and smaller than the median.
    % For lwc larger than the median...
    idx_larger = vertically_segmented_attributes{bb,2}>lwc_logNormal_median(bb);
    lwc_larger = vertically_segmented_attributes{bb,2}(idx_larger);
    lwc_custom_logNormal_std_larger(bb) = sqrt(mean((lwc_larger - lwc_logNormal_median(bb)).^2));

    % For lwc smaller than the median...
    idx_smaller = vertically_segmented_attributes{bb,2}<lwc_logNormal_median(bb);
    lwc_smaller = vertically_segmented_attributes{bb,2}(idx_smaller);
    lwc_custom_logNormal_std_smaller(bb) = sqrt(mean((lwc_smaller - lwc_logNormal_median(bb)).^2));
    %re_custom_logNormal_std_smaller(bb) = mean(re_logNormal_median(bb) - re_smaller);

    % --- for normal distribution statistics ---
    % compute the mean value for the current bin
    % the mean of the distribution (the standard way of computing the expected value)
    % is also the mean of the normal distribution. They are identical.
    %lwc_median(bb) = median(vertically_segmented_attributes{bb,2});       % g/cm^3 - mean liqiud water content

    % compute the standard deviation of the current bin
    % the std of the distribution (the standard way of computing the squareroot of the variance)
    % is also the std of the normal distribution. They are identical.
    %lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % g/cm^3 - standard deviation




    % ----- COMPUTE STATISTICS FOR DROPLET NUMBER CONCENTRATION -----

    % --- Use Normal distribution statistics ---
    % compute the mean value for the current bin
    Nc_median(bb) = Nc_fit_normal(bb).median;       % cm^(-3) - mean number concentration

    % compute the standard deviation of the current bin
    %Nc_std(bb) = Nc_fit_gamma(bb).std;         % cm^(-3) - standard deviation
    Nc_std(bb) = Nc_fit_normal(bb).std;         % cm^(-3) - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end



% ----------------------------------------------------------------------
% -------------- Make a subplot of all 3 mean profiles -----------------
% ----------------------------------------------------------------------




figure;

% --------------------------------
% plot the median effective radius
% --------------------------------

subplot(1,3,1)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
%x = [re_logNormal_median - re_logNormal_std; flipud(re_logNormal_median + re_logNormal_std)];
x = [re_logNormal_median - re_custom_logNormal_std_smaller; flipud(re_logNormal_median + re_custom_logNormal_std_larger)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(re_logNormal_median, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([4, 12])                   % microns


% --------------------------------------------
% plot the median liquid water content profile
% --------------------------------------------

subplot(1,3,2)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
% x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
x = [lwc_logNormal_median - lwc_custom_logNormal_std_smaller; flipud(lwc_logNormal_median + lwc_custom_logNormal_std_larger)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(lwc_logNormal_median, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter','latex')

% set x axis boundaries
xlim([0, 0.6])                   % g/cm^3

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

% set x axis boundaries
xlim([0, 320])                   % cm^(-3)

% set the size of the figure
set(gcf, 'Position', [0 0 1255 625])




% ----------------------------------------------------------------------
% ----------------------- Adiabatic Curve Fits -------------------------
% ----------------------------------------------------------------------

% ------------------------- LIQUID WATER CONTENT ---------------------------
% ----- Fit an adiabatic curve to the mean liquid water content profile -----
nudge_from_top = 2;
% lwc_slope = (lwc_median(end-nudge_from_top) - lwc_median(1))/(bin_center(end-nudge_from_top) - bin_center(1));
% lwc_intercept = lwc_median(1) - lwc_slope*bin_center(1);
lwc_slope = (lwc_logNormal_median(end-nudge_from_top) - lwc_logNormal_median(1))/(bin_center(end-nudge_from_top) - bin_center(1));
lwc_intercept = lwc_logNormal_median(1) - lwc_slope*bin_center(1);
lwc_adiabatic_fit = lwc_slope*bin_center + lwc_intercept;

% add to subplot(1,3,1)
subplot(1,3,2); hold on
plot(lwc_adiabatic_fit(1:end-nudge_from_top), bin_center(1:end-nudge_from_top), 'k', 'LineWidth',2)

% -- Include a legend in the lower right-hand corner of the 2nd subplot --
legend({'Standard Deviation', 'Median Profile', 'Adiabatic Fit'}, 'Interpreter','latex',...
    'Location','southeast', 'FontSize', 17)



% -------------------- EFFECTIVE DROPLET RADIUS -----------------------
% Plot an adiabatic curve fit to the mean droplet radius profile
nudge_from_bottom = 0;

% ----- Fit an adiabatic curve to the mean droplet profile -----
% use droplet profile function to create adiabatic fit
re_adiabatic_fit = create_droplet_profile2([re_logNormal_median(end), re_logNormal_median(1 + nudge_from_bottom)],...
    bin_center(1+nudge_from_bottom:end),'altitude', 'adiabatic');

% derive droplet profile from adiabatic fit of LWC


% add to subplot(1,3,1)
subplot(1,3,1); hold on
plot(re_adiabatic_fit, bin_center(1+nudge_from_bottom:end), 'k', 'LineWidth',2)





% ---------------------- DROPLET NUMBER CONCENTRATION ---------------------------
% ----- Fit an adiabatic curve to Number Concentration based on Adiabatic fits above -----
% density = 1;            % g/cm^3
% Nc_adiabatic_fit = 3*(lwc_adiabatic_fit./1e6)./...
%                 (4*pi*density*(re_adiabatic_fit*1e-4).^3);          % #/cm^3



subplot(1,3,3)
hold on;
%plot(Nc_adiabatic_fit, bin_center, 'k', 'LineWidth', 2)
xline(mean(Nc_median), 'Color', 'black', 'LineWidth', 2, 'Alpha',1,...
    'Label', [num2str(round(mean(Nc_median))), ' $cm^{-3}$'], 'Interpreter','latex',...
    'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','middle',...
    'Fontsize', 20)



% --- Create a Textbox stating these are only non-precipitating clouds ---
% Create textbox
% annotation('textbox',[0.134 0.802 0.07 0.06],...
%     'String',{'All Profiles'},...
%     'LineWidth',2,...
%     'Interpreter','latex',...
%     'FontSize',17,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off');


annotation('textbox',[0.134 0.802 0.142 0.114],...
    'String',{'Non-Precipitating clouds only ($LWP_{2DC}<1 \,g/m^{2}$)'},...
    'LineWidth',2,...
    'Interpreter','latex',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');







%% FIGURE 3

% Plot weighting functions of the first 7 MODIS spectral channels
% These weighting functions were created using the VOCALS-REx data set from
% Nov-9-2023. The droplet profile and optical depth were modeled after the
% vertical profile sampled at 1.734 hrs after the plane took off. The SZA
% was set as the value measured by MODIS for the median pixel, the pixel
% found closest the C130 aircraft in the middle of it's ascent through the
% cloud.

filenames = {'2D_MC_05-Sep-2023_Wavelength_469_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_555_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_645_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_858.5_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_1240_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_1640_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_2130_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat'};

% plot as a pdf
probability_str = 'pdf';

% define the wavelengths as the changing variables
wavelength = modisBands(1:7);
changing_variable = wavelength(:,1);        % center wavelenghts in nm

plot_probability_scat_top_maxDepth_with_changing_variable(filenames, probability_str ,changing_variable)





%% TAKE 2 - Figure 3 - Weighting Functions


clear variables


%wavelength = [555, 1240, 1640, 2000, 2130];

% filenames = {'2D_MC_04-Dec-2022_Wavelength_578_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
%     '2D_MC_04-Dec-2022_Wavelength_1650_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
%     '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
%     '2D_MC_04-Dec-2022_Wavelength_3700_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat'};

% filenames = {'2D_MC_07-Dec-2023_Wavelength_555_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_0.mat',...
%             '2D_MC_07-Dec-2023_Wavelength_1240_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_0.mat',...
%             '2D_MC_07-Dec-2023_Wavelength_1640_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_0.mat',...
%             '2D_MC_07-Dec-2023_Wavelength_2000_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_0.mat',...
%             '2D_MC_07-Dec-2023_Wavelength_2130_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_0.mat'};


% ---------------------------------------------------------------------------------------
filenames = {'2D_MC_05-Sep-2023_Wavelength_469_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_555_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_645_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_858.5_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_1240_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_1640_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
    '2D_MC_05-Sep-2023_Wavelength_2130_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat'};

% define the wavelengths as the changing variables
wavelength = modisBands(1:7);
% ---------------------------------------------------------------------------------------

% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';



% Do you want to smooth the raw PDF's?
smooth_curves = true;


% Define a set of colors based on the number of files
C = mySavedColors(1:length(filenames), 'fixed');


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));


% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results'];



elseif strcmp(whatComputer,'andrewbuggee')

    error([newline, 'Where is the new folder?', newline])

else
    error('I dont recognize this computer user name')
end


if strcmp(pwd,saved_simulations)==false
    cd(saved_simulations)
end



% Start figure
figure;

if smooth_curves==false

    % Plot the raw PDF's

    for nn = 1:length(filenames)


        % Load a simulation
        load(filenames{nn})



        % First select those photons that were scattered out the top

        index_scatter_out_top = final_state.scatter_out_top_INDEX;

        [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
            histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);



        % Plot the conditional probability
        plot(scatter_out_top_maxDepth_PDF,...
            scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2, 'Color',C(nn,:))
        hold on



        % Create legend string
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end



else

    % If this is true, we smooth each PDF to make a nice pretty plot, but
    % at the expense of loosing the PDF (the smoothed functions likely
    % won't integrate to 0)


    for nn = 1:length(filenames)


        % Load a simulation
        load(filenames{nn})



        % First select those photons that were scattered out the top

        index_scatter_out_top = final_state.scatter_out_top_INDEX;

        [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
            histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);



        % -------------------------------------------------------------
        % Integrate the drolet profile with the weighting function to
        % get an average effective radius measured, and thus an average
        % optical depth.
        % -------------------------------------------------------------
        if nn~=0
            % create an re vector that is the same length as our weighting
            % function
            new_tau = linspace(inputs.dropletProfile.tau_layer_mid_points(1), inputs.dropletProfile.tau_layer_mid_points(end), length(scatter_out_top_maxDepth_PDF));
            re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, new_tau);
            re_avg = trapz(new_tau, re .* scatter_out_top_maxDepth_PDF);
            tau_avg(nn) = interp1(re, new_tau,re_avg);


        end
        % -------------------------------------------------------------
        % -------------------------------------------------------------

        % Create smooth spline function
        f=fit((scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2)',scatter_out_top_maxDepth_PDF', 'smoothingspline','SmoothingParam',0.95);

        % Plot the conditional probability
        plot(f(scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2),...
            scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2, 'Color',C(nn,:))
        hold on



        % Create legend string
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end


    % horizontal line width
    horizontal_linewidth = 4;
    line_font_size = 23;

    for nn = 1:length(filenames)
        % Plot line of constant tau associated with retrieval depth
        yline(tau_avg(nn),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color',C(nn,:),'Label',...
            ['Depth of retrieved $r_e$ for ',num2str(wavelength(nn)/1e3),' $\mu m$'], 'Interpreter','latex',...
            'FontSize',line_font_size,'LabelVerticalAlignment','middle')

    end

    %     % Plot line of constant average tau for 0.66 microns
    %
    %     yline(tau_avg(1),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
    %         ['Depth of retrieved $r_e$ for ',num2str(wavelength(1)/1e3),' $\mu m$'], 'Interpreter','latex',...
    %         'FontSize',line_font_size,'LabelVerticalAlignment','bottom')
    %
    %     % Plot line of constant average tau for 1.6 microns
    %
    %     yline(tau_avg(2),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
    %         ['Depth of retrieved $r_e$ for ',num2str(wavelength(2)/1e3),' $\mu m$'], 'Interpreter','latex',...
    %         'FontSize',line_font_size,'LabelVerticalAlignment','bottom')
    %
    %     % Plot line of constant average tau for 2.2 microns
    %
    %     yline(tau_avg(3),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
    %         ['Depth of retrieved $r_e$ for ',num2str(wavelength(3)/1e3),' $\mu m$'], 'Interpreter','latex',...
    %         'FontSize',line_font_size,'LabelVerticalAlignment','top')
    %
    %     % Plot line of constant average tau for 3.7 microns
    %
    %     yline(tau_avg(4),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
    %         ['Depth of retrieved $r_e$ for ',num2str(wavelength(4)/1e3),' $\mu m$'], 'Interpreter','latex',...
    %         'FontSize',line_font_size,'LabelVerticalAlignment','top')







end



% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$P(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% Create title
title({'Conditional probability of photons that scatter out cloud top',...
    'reaching a max depth of $\tau$'},'Interpreter','latex')


% Create textbox with simulation properties

% Textbox
dim = [0.685 0.5 0 0];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';


% Create Legend
legend(legend_str,'Interpreter','latex','Location','northwest','FontSize',22)


set(gcf, 'Position',[0 0 1000 630])


%clear variables




%% Compute the effective radius retrieval using a single NIR wavelength and the median droplet profile found above


% ----- Load a weighting function -----

filenames = {'2D_MC_05-Sep-2023_Wavelength_2130_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat'};

% define the wavelengths as the changing variables
wavelength = modisBands(7);

% Load a simulation
load(filenames{1})


% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';


% First select those photons that were scattered out the top

index_scatter_out_top = final_state.scatter_out_top_INDEX;

[scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
    histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);

normalized_tau = scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2;
normalized_tau = normalized_tau./normalized_tau(end);

%normalized_tau = scatter_out_top_maxDepth_PDF_tau_edges./max(scatter_out_top_maxDepth_PDF_tau_edges);

% Create smooth spline function
f=fit(normalized_tau',scatter_out_top_maxDepth_PDF', 'smoothingspline','SmoothingParam',0.95);

tau_vector = scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2;
f=fit(tau_vector',scatter_out_top_maxDepth_PDF', 'smoothingspline','SmoothingParam',0.95);


% -------------------------------------------------------------
% Integrate the drolet profile with the weighting function to
% get an average effective radius measured, and thus an average
% optical depth.
% -------------------------------------------------------------


% re_logNormal_median is normalized, so simply multiply the bin_center
% values by the optical depth of scatter_out_top and integrate over this
% region

tau_2_integrate = bin_center * max(scatter_out_top_maxDepth_PDF_tau_edges);

re_retrieved = trapz(tau_2_integrate, f(tau_2_integrate) .* re_logNormal_meidan);

if nn~=0
    % create an re vector that is the same length as our weighting
    % function
    new_tau = linspace(inputs.dropletProfile.tau_layer_mid_points(1), inputs.dropletProfile.tau_layer_mid_points(end), length(scatter_out_top_maxDepth_PDF));
    re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, new_tau);
    re_avg = trapz(new_tau, re .* scatter_out_top_maxDepth_PDF);
    tau_avg(nn) = interp1(re, new_tau,re_avg);


end
% -------------------------------------------------------------
% -------------------------------------------------------------









%% FIGURE 4

% Plot the single wavelength retrieved droplet radius for each of the 7 wavelengths

figure;
plot(inputs.dropletProfile.re, inputs.dropletProfile.tau_layer_mid_points)
% flip the y axis so 0 is at the depth
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$r_e(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

hold on
for nn = 1:length(filenames)

    % Load a simulation
    load(filenames{nn})

    % compute the PDF of photons scattering out the top after reaching a
    % max depth of tau
    % First select those photons that were scattered out the top

    index_scatter_out_top = final_state.scatter_out_top_INDEX;

    [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_edges] = ...
        histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);

    % independent TAU variable for the PDF
    tau_pdf = scatter_out_top_maxDepth_PDF_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_edges)/2;

    % interpolate to get a droplet profile the same lenght as the PDF
    interp_re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, tau_pdf, 'spline');

    retrieved_re = trapz(tau_pdf, interp_re .* scatter_out_top_maxDepth_PDF);

    % add a vertical line to the plot
    xline(retrieved_re, 'LineStyle','--', 'LineWidth',2, 'Color', 'k',...
        'Label',['$\lambda = $',  num2str(round(inputs.mie.wavelength(1))), ' $\mu m$'],...
        'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','bottom', ...
        'Interpreter','latex', 'FontSize', 15)



end

set(gcf, 'Position',[0 0 1000 630])





%% FIGURE 5

% ---- Standard Deviation of Horizontal Profiles
% Compute the standard deviation of droplet size over some identified
% length scale. Slide this length scale across the horizontal profile and
% compute the standard deviation for each window. Do this for each profile
% in a data set


clear variables

% Plot the liquid wtaer content, effective radius and number concentration
% as a function of horizontal distance travelled for a single day of data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    folder_path = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/';

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
        'VOCALS_REx/vocals_rex_data/SPS_1/'];


end

% Oct-15-2008 Data
%filename = 'RF01.20081015.164800_201200.PNI.nc';

% Oct-18-2008 Data
%filename = 'RF02.20081018.130300_213000.PNI.nc';

% Oct-21-2008 Data
%filename = 'RF03.20081021.060000_142400.PNI.nc';

% Oct-25-2008 Data
% filename = 'RF05.20081025.062900_152500.PNI.nc';

% ----- November 9 data -----
filename = 'RF11.20081109.125700_213600.PNI.nc';

% ------ November 11 data -----
%filename = 'RF12.20081111.125000_214500.PNI.nc';

% ----- November 13 data -----
%filename = 'RF13.20081113.125700_215700.PNI.nc';

% ----- November 15 data -----
%filename = 'RF14.20081115.125800_220400.PNI.nc';


vocalsRex = readVocalsRex([folder_path, filename]);


% ---- set thresholds for the LWC and Nc ---
LWC_threshold = 0.03;       % g/m^3
Nc_threshold = 1;           % cm^{-3}
max_vertical_displacement = 10;     % meters

% ---- Find all Horizontal Profiles ---
horz_profs = find_horizontalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, Nc_threshold, max_vertical_displacement);

% --- Define the horizontal length scale to compute statistics over ---
length_scale = 1000;        % meters

% Loop through each profile

for nn = 1:length(horz_profs.lwc)

    % step trough each individual droplet size profile
    mean_lengthScale{nn} = [];
    std_lengthScale{nn} = [];

    for rr = 2:length(horz_profs.re{nn})

        % First, check to see if the remaining distance between rr and the
        % end of our horizontal profile is greater than our length scale.
        % If not, break the for loop
        if (horz_profs.horz_dist{nn}(end) - horz_profs.horz_dist{nn}(rr))>length_scale

            idx_displacement = 0;
            % find data points that make up 1 km
            dist_1km = horz_profs.horz_dist{nn}(rr+idx_displacement) - horz_profs.horz_dist{nn}(rr-1);

            while dist_1km<length_scale

                % step to the next data point
                idx_displacement = idx_displacement + 1;
                dist_1km = horz_profs.horz_dist{nn}(rr+idx_displacement) - horz_profs.horz_dist{nn}(rr-1);

            end

            % when the distance between two data points reaches 1000 meters,
            % stop and calculate the mean and standard deviation
            mean_lengthScale{nn} = [mean_lengthScale{nn}, mean(horz_profs.re{nn}(rr-1 : rr+idx_displacement))];         % microns
            std_lengthScale{nn} = [std_lengthScale{nn}, std(horz_profs.re{nn}(rr-1 : rr+idx_displacement))];           % microns

        else

            break

        end

    end

end



% ----- PLOT 1 -----
% Plot the mean versus the standard deviation for each profile
figure;
for nn = 1:length(mean_lengthScale)
    plot(std_lengthScale{nn}, mean_lengthScale{nn}, '.-', 'Linewidth', 1.75,...
        'MarkerSize', 15)
    hold on

end


% Include an x axis label on the middle plot
xlabel('$\sigma(r_e)$ ($\mu m$)', 'Interpreter','latex');

grid on; grid minor;
ylabel('$\left<r_e\right>$ ($\mu m$)', 'Interpreter','latex')

title(['Mean and STD over length scale of ', num2str(length_scale), ' meters'], 'Interpreter','latex')

% set plot size
set(gcf, 'Position', [0 0 1200 625])



% ----- PLOT 2 -----
% Try plotting the standard deviation for every length scale as a histogram
legend_str = cell(1, length(std_lengthScale));

figure;
for nn = 1:length(mean_lengthScale)

    histogram(std_lengthScale{nn}, 10, 'FaceAlpha', 0.5)
    hold on

    legend_str{nn} = ['index = ', num2str(nn), ' $\left<\sigma(r_e)\right> = \,$', num2str(mean(std_lengthScale{nn})), ' $\mu m$'];



end

% Include an x axis label on the middle plot
xlabel('$\sigma(r_e)$ ($\mu m$)', 'Interpreter','latex');

grid on; grid minor;
ylabel('Counts', 'Interpreter','latex')

title(['STD over length scale of ', num2str(length_scale), ' meters'], 'Interpreter','latex')

legend(legend_str, 'Interpreter','latex', 'Location','best', 'FontSize', 19)

% set plot size
set(gcf, 'Position', [0 0 1200 625])







%% FIGURE 6

% ----- Plot Horizontal Profiles -----
% Plot the liquid wtaer content, effective radius and number concentration
% as a function of horizontal distance travelled for a single day of data


clear variables


% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    folder_path = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/';

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
        'VOCALS_REx/vocals_rex_data/SPS_1/'];


end

% Oct-15-2008 Data
%filename = 'RF01.20081015.164800_201200.PNI.nc';

% Oct-18-2008 Data
%filename = 'RF02.20081018.130300_213000.PNI.nc';

% Oct-21-2008 Data
%filename = 'RF03.20081021.060000_142400.PNI.nc';

% Oct-25-2008 Data
% filename = 'RF05.20081025.062900_152500.PNI.nc';

% ----- November 9 data -----
filename = 'RF11.20081109.125700_213600.PNI.nc';

% ------ November 11 data -----
%filename = 'RF12.20081111.125000_214500.PNI.nc';

% ----- November 13 data -----
%filename = 'RF13.20081113.125700_215700.PNI.nc';

% ----- November 15 data -----
%filename = 'RF14.20081115.125800_220400.PNI.nc';


vocalsRex = readVocalsRex([folder_path, filename]);


% ---- set thresholds for the LWC and Nc ---
LWC_threshold = 0.03;       % g/m^3
Nc_threshold = 1;           % cm^{-3}
max_vertical_displacement = 15;     % meters

% ---- Find all Horizontal Profiles ---
horz_profs = find_horizontalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, Nc_threshold, max_vertical_displacement);


% ---- Plot the properties of all horizontal profiles collected ---
normalize_distance = false;
plot_horiztonal_profiles_LWC_and_re_and_Nc(horz_profs, 1:length(horz_profs.re), normalize_distance)


% ----- Plot JUST the effective radius as a function of distance ----

% Define the indices you'd like to plot
% indices_2_plot = [3, 8, 10];
indices_2_plot = [3, 4, 10, 12];

% define the colors of each curve
C = mySavedColors([6,7,9], 'fixed');

N_cuvres = length(indices_2_plot);

legend_str = cell(1, 2*N_cuvres);

std_val = zeros(1, N_cuvres);

mean_val = zeros(1, N_cuvres);

figure;
for nn = 1:N_cuvres

    % if normalize distance is true, all distance vectors will be
    % normalized between 0 and 1

    if normalize_distance==true

        norm_dist = (horz_profs.horz_dist{indices_2_plot(nn)} - min(horz_profs.horz_dist{indices_2_plot(nn)}))./...
            (max(horz_profs.horz_dist{indices_2_plot(nn)}) - min(horz_profs.horz_dist{indices_2_plot(nn)}));

        % plot the effective radius
        % if the 2DC data is compliant, plot the effective radius computed
        % using both instruments
        if horz_profs.flag_2DC_data_is_conforming==true

            plot(norm_dist, horz_profs.re{indices_2_plot(nn)}, 'Color',C(nn,:));

            hold on;

            % plot the average value as a dashed line
            mean_val(nn) = mean(horz_profs.re{indices_2_plot(nn)});
            constant_y_val = linspace(mean_val(nn), mean_val(nn), length(horz_profs.horz_dist{indices_2_plot(nn)}));
            plot(norm_dist, constant_y_val,'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

            % compute the standard deviation
            std_val(nn) = std(horz_profs.re{indices_2_plot(nn)});


        else

            % if the 2DC data is non-conforming, use only the CDP data and
            % make a note of it
            plot(norm_dist, horz_profs.re_CDP{indices_2_plot(nn)}, 'Color',C(nn,:));

            % plot the average value as a dashed line
            mean_val(nn) = mean(horz_profs.re_CDP{indices_2_plot(nn)});
            constant_y_val = linspace(mean_val(nn), mean_val(nn), length(horz_profs.horz_dist{indices_2_plot(nn)}));
            plot(norm_dist, constant_y_val,'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

            % compute the standard deviation
            std_val(nn) = std(horz_profs.re_CDP{indices_2_plot(nn)});

        end
        hold on



    else

        % --- DATA IS IN METERS - PLOT IN KILOMETERS ---

        % plot the effective radius
        % if the 2DC data is compliant, plot the effective radius computed
        % using both instruments

        if horz_profs.flag_2DC_data_is_conforming==true

            plot(horz_profs.horz_dist{indices_2_plot(nn)}./1e3, horz_profs.re{indices_2_plot(nn)},...
                'Color',C(nn,:));

            hold on;

            % plot the average value as a dashed line
            mean_val(nn) = mean(horz_profs.re{indices_2_plot(nn)});
            constant_y_val = linspace(mean_val(nn), mean_val(nn), length(horz_profs.horz_dist{indices_2_plot(nn)}));
            plot(horz_profs.horz_dist{indices_2_plot(nn)}./1e3, constant_y_val,...
                'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

            % compute the standard deviation
            std_val(nn) = std(horz_profs.re{indices_2_plot(nn)});

        else

            % if the 2DC data is non-conforming, use only the CDP data and
            % make a note of it
            plot(horz_profs.horz_dist{indices_2_plot(nn)}./1e3, horz_profs.re_CDP{indices_2_plot(nn)},...
                'Color',C(nn,:));

            % plot the average value as a dashed line
            mean_val(nn) = mean(horz_profs.re_CDP{indices_2_plot(nn)});
            constant_y_val = linspace(mean_val(nn), mean_val(nn), length(horz_profs.horz_dist{indices_2_plot(nn)}));
            plot(horz_profs.horz_dist{indices_2_plot(nn)}./1e3, constant_y_val,...
                'LineStyle','--','LineWidth',1.5,'Color', C(nn,:));

            % compute the standard deviation
            std_val(nn) = std(horz_profs.re_CDP{indices_2_plot(nn)});

        end
        hold on


    end

    legend_str{2*nn - 1} = ['$\sigma$ = ', num2str(round(std_val(nn), 2)), ' $\mu m$'];
    % skip one because of the mean value
    legend_str{2*nn} = ['$\left<r_e \right>$ = ', num2str(round(mean_val(nn), 2)), ' $\mu m$'];


end


grid on; grid minor;
% if the 2DC data is compliant, plot the effective radius computed
% using both instruments
if horz_profs.flag_2DC_data_is_conforming==true
    ylabel('$r_e$ ($\mu m$)', 'Interpreter','latex')
else
    % if the 2DC data is non-conforming, use only the CDP data and
    % make a note of it
    ylabel('$r_e$ ($\mu m$) - (CDP only)', 'Interpreter','latex')
end

% include a title in the middle plot
if isfield(horz_profs, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(horz_profs.LWC_threshold),' $g/m^{3}$',...
        '     $N_c \geq$ ', num2str(horz_profs.Nc_threshold), ' $cm^{-3}$',...
        '     Max vert displacement: ', num2str(horz_profs.max_vert_displacement), ' $m$'], 'interpreter', 'latex')

elseif isfield(horz_profs.inputs, 'LWC_threshold')==true
    title(['$LWC \geq$ ', num2str(horz_profs.inputs.LWC_threshold),' $g/m^{3}$',...
        '     $N_c \geq$ ', num2str(horz_profs.inputs.Nc_threshold), ' $cm^{-3}$',...
        '     Max vert displacement: ', num2str(horz_profs.inputs.max_vert_displacement), ' $m$'], 'interpreter', 'latex')

end

% Include an x axis label on the middle plot
if normalize_distance==true

    xlabel('Normalized Horizontal Distance Travelled', 'Interpreter','latex');
else

    xlabel('Horizontal Distance Travelled ($km$)', 'Interpreter','latex');
end




% in the third subplot, define the indices_2_plot being plotted
legend(legend_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20)

% set plot size
set(gcf, 'Position', [0 0 1200 625])






%% FIGURE 7 - Plot the ensemble MEAN of droplet size, liquid water content and
% number concentration for non-precipitating clouds. Add an adiabatic
% profile for the liquid water content and effective radius to show the
% mean profiles are close to adiabatic, supporting my assumption.
% ----- For Horizontal Profiles ----


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023.mat'])



elseif strcmp(whatComputer, 'andrewbuggee')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 10 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 50 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_50_30-Oct-2023'])



end




% using the mean ensemble function to plot the mean horizontal profile of
% the ensemble




% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Segment re into N bins along optical depth

% In order to compute a mean horizontal profile, we have to first normalize
% the horizontal extent so that all profiles lie between values [0,1]. Then
% we break up the horizontal component in n discrete bins. Within each bin we
% can compute the mean, median and standard deviation

n_bins = 30; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
horizontally_segmented_attributes = cell(n_bins, 1);


normalized_dist = cell(1, length(ensemble_profiles.lwc));


for nn = 1:length(ensemble_profiles.lwc)

    % first we need to normalize the vertical component of all profiles
    normalized_dist{nn} = ensemble_profiles.horz_dist{nn}./max(ensemble_profiles.horz_dist{nn});

    % the data is stored in horizontal distance space.

    re = ensemble_profiles.re{nn};




    % for each profile, we need to segment the variables of interest into n
    % bins.

    for bb = 1:length(bin_edges)-1

        % grab all re values within each bin. Segment them
        % accordingly
        if bb==1
            index_segment = normalized_dist{nn}>=bin_edges(bb) & normalized_dist{nn}<=bin_edges(bb+1);

        else
            index_segment = normalized_dist{nn}>bin_edges(bb) & normalized_dist{nn}<=bin_edges(bb+1);
        end

        % store the effective radius values
        horizontally_segmented_attributes{bb, 1} = [horizontally_segmented_attributes{bb, 1}; re(index_segment)];




    end



end



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Create a PDF object at each level in the cloud and fit a distribution to this PDF

% store the rejection of each null hypothesis and the p-value for each
% chi-squared test

re_reject_normal = zeros(1, size(horizontally_segmented_attributes,1));
re_p_normal = zeros(1, size(horizontally_segmented_attributes,1));

re_reject_lognormal = zeros(1, size(horizontally_segmented_attributes,1));
re_p_lognormal = zeros(1, size(horizontally_segmented_attributes,1));

re_reject_gamma = zeros(1, size(horizontally_segmented_attributes,1));
re_p_gamma = zeros(1, size(horizontally_segmented_attributes,1));



for bb = 1:size(horizontally_segmented_attributes, 1)


    % -----------------------------------------------
    % ------- EFFECTIVE DROPLET RADIUS FITTING ------
    % -----------------------------------------------


    % fit the effective radius data to a normal distribution
    re_fit_normal(bb) = fitdist(horizontally_segmented_attributes{bb,1}, 'normal');
    [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(horizontally_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb));

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal(bb) = fitdist(horizontally_segmented_attributes{bb,1}, 'lognormal');
    [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(horizontally_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb));

    % fit the effective radius data to a gamma distribution
    re_fit_gamma(bb) = fitdist(horizontally_segmented_attributes{bb,1}, 'gamma');
    [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(horizontally_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
    %         figure; subplot(1,3,1); plot(re_fit_normal(bb)); title('Normal Fit'); xlabel('r_e (\mum)')
    %         subplot(1,3,2); plot(re_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('r_e (\mum)')
    %         subplot(1,3,3); plot(re_fit_gamma(bb)); title('Gamma Fit'); xlabel('r_e (\mum)')
    %         set(gcf, 'Position', [0 0 1200 500])



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



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Compute the mean re of each layer

% ---- most common best fit distribution for r_e was is the log-normal dist ---
re_logNormal_std = zeros(n_bins, 1);
re_logNormal_mean = zeros(n_bins, 1);



bin_center = zeros(n_bins, 1);




for bb = 1:n_bins


    % ----- COMPUTE STATISTICS FOR DROPLET SIZE -----

    % find the mean of the log normal distribution
    re_logNormal_mean(bb) = exp(re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2 /2);

    % find squareroot of the variance of the lognormal distribution
    re_logNormal_std(bb) = sqrt(exp(2*re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2)*(exp(re_fit_lognormal(bb).sigma^2) - 1));

    % ---------- COMPUTE BIN CENTER VALUES ------
    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end



% ----------------------------------------------------------------------
% ------------------- Make plot of mean re profile ---------------------
% ----------------------------------------------------------------------




figure;

% plot the mean effective radius

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
y = [re_logNormal_mean - re_logNormal_std; flipud(re_logNormal_mean + re_logNormal_std)];
x = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(bin_center, re_logNormal_mean, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
ylabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
xlabel('Normalized Distance', 'Interpreter', 'latex')




% set the figure title
title(['Mean Profiles:  $LWC \geq$', num2str(ensemble_profiles.inputs.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles.inputs.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex')


% set the size of the figure
set(gcf, 'Position', [0 0 1255 625])



% --- Create a Textbox stating these are only non-precipitating clouds ---
% Create textbox
annotation('textbox',[0.134 0.802 0.142 0.114],...
    'String',{'Non-Precipitating clouds only ($LWP_{2DC}<25 \,g/m^{2}$)'},...
    'LineWidth',2,...
    'Interpreter','latex',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');







%% FIGURE 8 - Standard Deviation of all ensemble profiles for different length scales

% ---- Standard Deviation of Horizontal Profiles ---------
% Compute the standard deviation of droplet size over some identified
% length scale. Slide this length scale across the horizontal profile and
% compute the standard deviation for each window. Do this for each profile
% in a data set


clear variables


% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
%     load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
%         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023.mat'])


    % ****** Using LWC Adjustment in readVocals ******
    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 35 g/m^2         -----
    load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_35_06-Nov-2023'])



elseif strcmp(whatComputer, 'andrewbuggee')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 10 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023'])


    % ****** Using LWC Adjustment in readVocals ******
    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 35 g/m^2         -----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_35_06-Nov-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 50 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_50_30-Oct-2023'])



end



% --- Define the horizontal length scale to compute statistics over ---
length_scale = [1000, 5000];        % meters

% Loop through each profile in the ensemble

for ll = 1:length(length_scale)

    % Set up an empty array for each length scale
    % step trough each individual droplet size profile
    mean_lengthScale{ll} = [];
    std_lengthScale{ll} = [];

    for nn = 1:length(ensemble_profiles.altitude)



        for rr = 2:length(ensemble_profiles.re{nn})

            % First, check to see if the remaining distance between rr and the
            % end of our horizontal profile is greater than our length scale.
            % If not, break the for loop
            if (ensemble_profiles.horz_dist{nn}(end) - ensemble_profiles.horz_dist{nn}(rr))>length_scale(ll)

                idx_displacement = 0;
                % find data points that make up the length scale
                dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                while dist_km<length_scale(ll)

                    % step to the next data point
                    idx_displacement = idx_displacement + 1;
                    dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                end

                % when the distance between two data points reaches 1000 meters,
                % stop and calculate the mean and standard deviation
                mean_lengthScale{ll} = [mean_lengthScale{ll}, mean(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];         % microns
                std_lengthScale{ll} = [std_lengthScale{ll}, std(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];           % microns

            else

                break

            end

        end

    end

end





% ----- PLOT 1 -----
% Try plotting the standard deviation for every length scale as a histogram

% compute
figure;

for ll = 1:length(length_scale)



    % ------- Plot Histogram and Mean Value ------

    histogram(std_lengthScale{ll}, 100, 'FaceAlpha', 0.5, 'FaceColor', mySavedColors(ll+2, 'fixed'))
    hold on


    % Plot a vertical line showing the average value of the distribution
    xline(mean(std_lengthScale{ll}), 'LineWidth', 3, 'LineStyle', '--',...
        'Color', mySavedColors(ll+2, 'fixed'), 'FontSize', 23, 'FontWeight','bold')

    legend_str{2*ll - 1} = ['Length Scale = ', num2str(length_scale(ll)/1e3), ' $km$'];
    legend_str{2*ll} = ['$\left< \sigma_{r_e} \right> $ = ', num2str(round(mean(std_lengthScale{ll}), 2)), '$\mu m$'];

    ylabel('Counts', 'Interpreter','latex')





    % ------- Plot CDF and 0.5 yline ------

    %     histogram(std_lengthScale{ll}, 100, 'FaceAlpha', 0.5, 'Normalization', 'cdf')
    %     hold on
    %
    %
    %     % Plot a vertical line showing the average value of the distribution
    %     yline(0.5, 'LineWidth', 2, 'Label','Half Total Counts',...
    %         'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
    %         'Interpreter', 'latex', 'LineStyle', '--',...
    %         'Color', 'k', 'FontSize', 23, 'FontWeight','bold')
    %
    %     ylabel('CDF', 'Interpreter','latex')


end


% ------ Pretty Plot Stuff -------

% Include an x axis label on the middle plot
xlabel('$\sigma(r_e)$ ($\mu m$)', 'Interpreter','latex');

grid on; grid minor;
ylabel('Counts', 'Interpreter','latex')

title(['STD over different length scales'], 'Interpreter','latex')

legend(legend_str, 'Location', 'best', 'interpreter', 'latex', 'Fontsize', 22)


% set plot size
set(gcf, 'Position', [0 0 1200 625])








%% FIGURE 9 - Standard Deviation of all ensemble profiles for different length scales and for different mean droplet sizes
% Sort by mean droplet size over the entire horizontal profile

clear variables


% Let's sort the horizontal profiles into 3 regimes:
% The mean value of the horizontal profile will be broken up into values
% less than 6 microns, between 6 and 7 microns, and greater than 7 microns.
re_regimes = [0, 6.25, 7.25, 8.25, inf];



% ---- Standard Deviation of Horizontal Profiles ---------
% Compute the standard deviation of droplet size over some identified
% length scale. Slide this length scale across the horizontal profile and
% compute the standard deviation for each window. Do this for each profile
% in a data set




% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true






elseif strcmp(whatComputer, 'andrewbuggee')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 10 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 50 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_50_30-Oct-2023'])



end



% --- Define the horizontal length scale to compute statistics over ---
length_scale = [1000];        % meters


% Let's compute the mean value of each horizontal profile
mean_overProfile = zeros(length(ensemble_profiles.re), 1);


% Loop through each profile in the ensemble

for ll = 1:length(length_scale)

    % Set up an empty array for each length scale
    % step trough each individual droplet size profile
    mean_lengthScale{ll, length(re_regimes)-1} = [];
    std_lengthScale{ll, length(re_regimes)-1} = [];




    for nn = 1:length(ensemble_profiles.altitude)

        % first, compute the average value over the entire profile
        mean_overProfile(nn) = mean(ensemble_profiles.re{nn});

        % Sort based on the defined boundaries
        for mm = 1:length(re_regimes)-1

            idx_regime(mm) = mean_overProfile(nn)>=re_regimes(mm) & mean_overProfile(nn)<re_regimes(mm+1);

        end


        for rr = 2:length(ensemble_profiles.re{nn})

            % First, check to see if the remaining distance between rr and the
            % end of our horizontal profile is greater than our length scale.
            % If not, break the for loop
            if (ensemble_profiles.horz_dist{nn}(end) - ensemble_profiles.horz_dist{nn}(rr))>length_scale(ll)

                idx_displacement = 0;
                % find data points that make up the length scale
                dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                while dist_km<length_scale(ll)

                    % step to the next data point
                    idx_displacement = idx_displacement + 1;
                    dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                end


                % when the distance between two data points reaches 1000 meters,
                % stop and calculate the mean and standard deviation


                mean_lengthScale{ll, idx_regime} = [mean_lengthScale{ll, idx_regime}, mean(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];         % microns
                std_lengthScale{ll, idx_regime} = [std_lengthScale{ll, idx_regime}, std(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];           % microns

            else

                break

            end

        end


    end

end





% ----- PLOT 1 -----
% Try plotting the standard deviation for every length scale as a histogram

% compute

for ll = 1:length(length_scale)

    figure;

    for mm = 1:length(re_regimes)-1

        % ------- Plot Histogram and Mean Value ------

        h(mm) = histogram(std_lengthScale{ll,mm}, 100, 'FaceAlpha', 0.5, ...
            'FaceColor', mySavedColors(mm, 'fixed'));
        hold on


        % Plot a vertical line showing the average value of the distribution
        legend_str{mm} = [num2str(re_regimes(mm)), '$ < \; \left<r_e\right> \; \leq$ ', num2str(re_regimes(mm+1)),...
            ' $\mu m$; Avg. = ',num2str(round(mean(std_lengthScale{ll,mm}), 2)), ' $\mu m$'];


    end


    ylabel('Counts', 'Interpreter','latex')





    % ------- Plot CDF and 0.5 yline ------

    %     histogram(std_lengthScale{ll}, 100, 'FaceAlpha', 0.5, 'Normalization', 'cdf')
    %     hold on
    %
    %
    %     % Plot a vertical line showing the average value of the distribution
    %     yline(0.5, 'LineWidth', 2, 'Label','Half Total Counts',...
    %         'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
    %         'Interpreter', 'latex', 'LineStyle', '--',...
    %         'Color', 'k', 'FontSize', 23, 'FontWeight','bold')
    %
    %     ylabel('CDF', 'Interpreter','latex')





    % ------ Pretty Plot Stuff -------

    % Include an x axis label on the middle plot
    xlabel('$\sigma(r_e)$ ($\mu m$)', 'Interpreter','latex');

    grid on; grid minor;
    ylabel('Counts', 'Interpreter','latex')

    title(['STD over length scale of ', num2str(length_scale(ll)), ' meters'], 'Interpreter','latex')

    legend(legend_str, 'Interpreter', 'latex', 'FontSize', 23, 'FontWeight','bold',...
        'Location', 'best')


    % set plot size
    set(gcf, 'Position', [0 0 1200 625])

    % Set the order of the different histograms
    uistack(h(1), 'top')
    uistack(h(3), 'top')


end






%% FIGURE 10 - Standard Deviation of all ensemble profiles for different length scales and for different mean droplet sizes
% Sort by mean droplet size of the length segment

clear variables


% Let's sort every length segment into 3 regimes:
% The mean value of the horizontal profile will be broken up into values
% less than 6 microns, between 6 and 7 microns, and greater than 7 microns.
re_regimes = [0, 6.25, 7.25, inf];



% ---- Standard Deviation of Horizontal Profiles ---------
% Compute the standard deviation of droplet size over some identified
% length scale. Slide this length scale across the horizontal profile and
% compute the standard deviation for each window. Do this for each profile
% in a data set




% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023.mat'])




elseif strcmp(whatComputer, 'andrewbuggee')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 10 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 50 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_50_30-Oct-2023'])



end



% --- Define the horizontal length scale to compute statistics over ---
length_scale = [1000];        % meters



% Loop through each profile in the ensemble

for ll = 1:length(length_scale)

    % Set up an empty array for each length scale
    % step trough each individual droplet size profile
    mean_lengthScale{ll} = [];
    std_lengthScale{ll, length(re_regimes)-1} = [];




    for nn = 1:length(ensemble_profiles.altitude)



        for rr = 2:length(ensemble_profiles.re{nn})

            % First, check to see if the remaining distance between rr and the
            % end of our horizontal profile is greater than our length scale.
            % If not, break the for loop
            if (ensemble_profiles.horz_dist{nn}(end) - ensemble_profiles.horz_dist{nn}(rr))>length_scale(ll)

                idx_displacement = 0;
                % find data points that make up the length scale
                dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                while dist_km<length_scale(ll)

                    % step to the next data point
                    idx_displacement = idx_displacement + 1;
                    dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                end


                % when the distance between two data points reaches 1000 meters,
                % stop and calculate the mean and standard deviation


                mean_lengthScale{ll} = [mean_lengthScale{ll}, mean(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];         % microns

                % Sort based on the defined boundaries
                for mm = 1:length(re_regimes)-1

                    idx_regime(mm) = mean_lengthScale{ll}(end)>=re_regimes(mm) & mean_lengthScale{ll}(end)<re_regimes(mm+1);

                end

                std_lengthScale{ll, idx_regime} = [std_lengthScale{ll, idx_regime}, std(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];           % microns

            else

                break

            end

        end


    end

end





% ----- PLOT 1 -----
% Try plotting the standard deviation for every length scale as a histogram


for ll = 1:length(length_scale)

    figure;

    for mm = 1:length(re_regimes)-1

        % ------- Plot Histogram and Mean Value ------

        h(mm) = histogram(std_lengthScale{ll,mm}, 100, 'FaceAlpha', 0.5, ...
            'FaceColor', mySavedColors(mm, 'fixed'));
        hold on


        % Plot a vertical line showing the average value of the distribution
        legend_str{mm} = [num2str(re_regimes(mm)), '$ < \; \left<r_e\right> \; \leq$ ', num2str(re_regimes(mm+1)),...
            ' $\mu m$; Avg. = ',num2str(round(mean(std_lengthScale{ll,mm}), 2)), ' $\mu m$'];


    end


    ylabel('Counts', 'Interpreter','latex')





    % ------- Plot CDF and 0.5 yline ------

    %     histogram(std_lengthScale{ll}, 100, 'FaceAlpha', 0.5, 'Normalization', 'cdf')
    %     hold on
    %
    %
    %     % Plot a vertical line showing the average value of the distribution
    %     yline(0.5, 'LineWidth', 2, 'Label','Half Total Counts',...
    %         'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
    %         'Interpreter', 'latex', 'LineStyle', '--',...
    %         'Color', 'k', 'FontSize', 23, 'FontWeight','bold')
    %
    %     ylabel('CDF', 'Interpreter','latex')





    % ------ Pretty Plot Stuff -------

    % Include an x axis label on the middle plot
    xlabel('$\sigma(r_e)$ ($\mu m$)', 'Interpreter','latex');

    grid on; grid minor;
    ylabel('Counts', 'Interpreter','latex')

    title(['STD over length scale of ', num2str(length_scale(ll)), ' meters'], 'Interpreter','latex')

    legend(legend_str, 'Interpreter', 'latex', 'FontSize', 23, 'FontWeight','bold',...
        'Location', 'best')


    % set plot size
    set(gcf, 'Position', [0 0 1200 625])

    % Set the order of the different histograms
    %uistack(h(end), 'top')
    uistack(h(2), 'bottom')


end







%% FIGURE 11 - Standard Deviation of all ensemble profiles as a function of length

% ---- Standard Deviation of Horizontal Profiles ---------
% Compute the standard deviation of droplet size over some identified
% length scale. Slide this length scale across the horizontal profile and
% compute the standard deviation for each window. Do this for each profile
% in a data set


clear variables


% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023.mat'])




elseif strcmp(whatComputer, 'andrewbuggee')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 10 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 50 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_50_30-Oct-2023'])



end



% --- Define the horizontal length scale to compute statistics over ---
length_scale = 1000:100:5000;        % meters

% Loop through each profile in the ensemble

for ll = 1:length(length_scale)

    % Set up an empty array for each length scale
    % step trough each individual droplet size profile
    mean_lengthScale{ll} = [];
    std_lengthScale{ll} = [];

    for nn = 1:length(ensemble_profiles.altitude)



        for rr = 2:length(ensemble_profiles.re{nn})

            % First, check to see if the remaining distance between rr and the
            % end of our horizontal profile is greater than our length scale.
            % If not, break the for loop
            if (ensemble_profiles.horz_dist{nn}(end) - ensemble_profiles.horz_dist{nn}(rr))>length_scale(ll)

                idx_displacement = 0;
                % find data points that make up the length scale
                dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                while dist_km<length_scale(ll)

                    % step to the next data point
                    idx_displacement = idx_displacement + 1;
                    dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                end

                % when the distance between two data points reaches 1000 meters,
                % stop and calculate the mean and standard deviation
                mean_lengthScale{ll} = [mean_lengthScale{ll}, mean(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];         % microns
                std_lengthScale{ll} = [std_lengthScale{ll}, std(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];           % microns

            else

                break

            end

        end

    end

end





% ----- PLOT 1 -----
% Try plotting the standard deviation for every length scale as a histogram

figure;

for ll = 1:length(length_scale)



    plot(length_scale(ll)./1e3, mean(std_lengthScale{ll}), '.', 'MarkerSize', 25, ...
        'Color', mySavedColors(1, 'fixed'))
    hold on


end

% Include an x axis label on the middle plot
xlabel('Pixel Length Scale (km)', 'Interpreter','latex');

grid on; grid minor;
ylabel('$\left< \sigma_{r_e} \right>$  ($\mu m$)', 'Interpreter','latex')

title('Mean Standard Deviation for different Length Scales', 'Interpreter','latex')


% set plot size
set(gcf, 'Position', [0 0 1200 625])






%% FIGURE 12 - Compute the spread of each STD distribution for all ensemble profiles as a function of length

% ---- Standard Deviation of Horizontal Profiles ---------
% Compute the standard deviation of droplet size over some identified
% length scale. Slide this length scale across each horizontal profile and
% compute the standard deviation for each window. Do this for each profile
% in a data set and do this for a range of length scales.

% Then, let's fit a distribution and find the width parameter of that
% distribution


clear variables


% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023.mat'])



elseif strcmp(whatComputer, 'andrewbuggee')==true


    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 10 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %         '/SPS_1/ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 25 g/m^2         -----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
        'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_25_30-Oct-2023'])



    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
    % ------           2DC LWP < 50 g/m^2         -----
    %     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/SPS_1/',...
    %         'ensemble_horizontal_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_drizzleLWP-threshold_50_30-Oct-2023'])



end



% --- Define the horizontal length scale to compute statistics over ---
length_scale = 1000:100:5000;        % meters


% store the refection of each null hypothesis and the p-value for each
% chi-squared test

dist_reject_normal = zeros(1, length(length_scale));
dist_p_normal = zeros(1, length(length_scale));

dist_reject_lognormal = zeros(1, length(length_scale));
dist_p_lognormal = zeros(1, length(length_scale));

dist_reject_gamma = zeros(1, length(length_scale));
dist_p_gamma = zeros(1, length(length_scale));




% Loop through each profile in the ensemble

for ll = 1:length(length_scale)

    % Set up an empty array for each length scale
    % step trough each individual droplet size profile
    std_lengthScale{ll} = [];

    for nn = 1:length(ensemble_profiles.altitude)



        for rr = 2:length(ensemble_profiles.re{nn})

            % First, check to see if the remaining distance between rr and the
            % end of our horizontal profile is greater than our length scale.
            % If not, break the for loop
            if (ensemble_profiles.horz_dist{nn}(end) - ensemble_profiles.horz_dist{nn}(rr))>length_scale(ll)

                idx_displacement = 0;
                % find data points that make up the length scale
                dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                while dist_km<length_scale(ll)

                    % step to the next data point
                    idx_displacement = idx_displacement + 1;
                    dist_km = ensemble_profiles.horz_dist{nn}(rr+idx_displacement) - ensemble_profiles.horz_dist{nn}(rr-1);

                end

                % when the distance between two data points reaches 1000 meters,
                % stop and calculate the mean and standard deviation
                std_lengthScale{ll} = [std_lengthScale{ll}, std(ensemble_profiles.re{nn}(rr-1 : rr+idx_displacement))];           % microns

            else

                break

            end

        end

    end

    % After we have the distribution for a given length scale, compute the
    % distribution fit parameters

    % fit the effective radius data to a normal distribution
    dist_fit_normal(ll) = fitdist(std_lengthScale{ll}', 'normal');
    [dist_reject_normal(ll), dist_p_normal(ll)] = chi2gof(std_lengthScale{ll}', 'CDF',dist_fit_normal(ll));

    % fit the effective radius data to a log-normal distribution
    dist_fit_lognormal(ll) = fitdist(std_lengthScale{ll}', 'lognormal');
    [dist_reject_lognormal(ll), dist_p_lognormal(ll)] = chi2gof(std_lengthScale{ll}', 'CDF',dist_fit_lognormal(ll));

    % fit the effective radius data to a gamma distribution
    dist_fit_gamma(ll) = fitdist(std_lengthScale{ll}', 'gamma');
    [dist_reject_gamma(ll), dist_p_gamma(ll)] = chi2gof(std_lengthScale{ll}', 'CDF',dist_fit_gamma(ll));



end


bin_names = {'Normal', 'Log-Normal', 'Gamma'};
% -----------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS FITTING ------
% -----------------------------------------------
[max_dist_p, idx_dist_p] = max([dist_p_normal; dist_p_lognormal; dist_p_gamma],[], 1);

figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_dist_p==1), sum(idx_dist_p==2), sum(idx_dist_p==3)]);
title('best distribution fit'); ylabel('Counts')


% ----- PLOT -----
% The gamma distribution is the best fit for the majority of length scales
% Compute the length of two standard deviations (+/- 1 sigma) for each
% distribution using the gamma fit




for ll = 1:length(length_scale)

    pm_1sigma(ll) = 2*dist_fit_gamma(ll).std;       % microns

    % grab the median of the distribution
    % does this value change with pixel length scale?
    dist_median(ll) = dist_fit_gamma(ll).median;    %   microns


    % grab the mean of the distribution
    % does this value change with pixel length scale?
    dist_mean(ll) = dist_fit_gamma(ll).mean;    %   microns

end


% ----- PLOT +/- 1SIGMA (2*SIGMA) -----
figure;
plot(length_scale./1e3, pm_1sigma, '.', 'MarkerSize', 25, ...
    'Color', mySavedColors(1, 'fixed'))
hold on

% Include an x axis label on the middle plot
xlabel('Pixel Length Scale (km)', 'Interpreter','latex');

grid on; grid minor;
ylabel('$\pm \, 1 \sigma$  ($\mu m$)', 'Interpreter','latex')

title('$\pm \, 1 \sigma$ of the distribution of STDs for each length scale', 'Interpreter','latex')


% set plot size
set(gcf, 'Position', [0 0 1200 625])



% ----- PLOT MEDIAN -----
figure;
plot(length_scale./1e3, dist_median, '.', 'MarkerSize', 25, ...
    'Color', mySavedColors(1, 'fixed'))
hold on

% Include an x axis label on the middle plot
xlabel('Pixel Length Scale (km)', 'Interpreter','latex');

grid on; grid minor;
ylabel('Median  ($\mu m$)', 'Interpreter','latex')

title('Median of the distribution of STDs for each length scale', 'Interpreter','latex')


% set plot size
set(gcf, 'Position', [0 0 1200 625])




% ----- PLOT MEAN -----
figure;
plot(length_scale./1e3, dist_mean, '.', 'MarkerSize', 25, ...
    'Color', mySavedColors(1, 'fixed'))
hold on

% Include an x axis label on the middle plot
xlabel('Pixel Length Scale (km)', 'Interpreter','latex');

grid on; grid minor;
ylabel('Mean  ($\mu m$)', 'Interpreter','latex')

title('Mean of the distribution of STDs for each length scale', 'Interpreter','latex')


% set plot size
set(gcf, 'Position', [0 0 1200 625])


%% Subplot comparing l2 norm residual for synthetic data bewteen different measurement uncertainty scenarios

clear variables

% define axes label font size
axes_label_font_size = 40;

% define axes tick label font size
axes_tick_label_font_size = 25;

% define colorbar font size
cb_font_size = 40;

% define contour label size
contour_label_size = 25;


% LOAD DATA SET

load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-11-Dec-2024_rev1.mat')

% Create mesh grid
[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);


% Define synthetic model data

r_top_truth = 10.17;
r_bot_truth = 4.74;
tau_c_truth = 5.96;

[~, idx_r_top] = min(abs(r_top_fine - r_top_truth));
[~, idx_r_bot] = min(abs(r_bot_fine - r_bot_truth));
[~, idx_tau_c] = min(abs(tau_c_fine - tau_c_truth));

synthetic_measurement = reshape(Refl_model_fine(idx_r_top, idx_r_bot, idx_tau_c, :), [],1);


% --- Create synthetic measurements with 1% uncertinaty ---
% -------------------------------------
% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
measurement_uncert = 0.05;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve

% Compute the new synethtic measurement with gaussian noise
% *** Gaussian noise can be either positive or negative. Meaning, an
% uncertainty of 5% implies the true value can lie anywhere between the
% measured value +/- 5% of the measured value
% define a
synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert/3) .*...
    randn(length(inputs.bands2run), 1);

% define the synthetic relfectance uncertainty
synthetic_measurement_uncert = measurement_uncert .* synthetic_measurement_with_noise;






% Using an exact modeled estimate without noise
use_l2_norm = true;

rms_residual = [];
rms_uncert = [];

if use_l2_norm==false

    % Compute the rms difference between the measurements and the modeled
    % reflectances
    rms_residual = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));


    % Compute the RMS of the synthetic measurement uncertainty
    rms_uncert = sqrt( mean( synthetic_measurement_uncert.^2));

else

    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rms_residual = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));

    % Compute the l2 norm of the synthetic measurement uncertainty
    rms_uncert = sqrt( sum( synthetic_measurement_uncert.^2));

end



% Find the states with the lowest rms residul

% find the smallest rms residual value, omitting nans
[~, idx_min] = min(rms_residual, [], 'all', 'omitnan');

r_top_min = R_top_fine(idx_min);




% Create Contour plot of rms residual between true EMIT measurements and the libRadTran modeled measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----

% define the optical depth slice you'd like to plot
idx_rTop_5pct = r_top_fine == r_top_min(1);

% Create figure
figure;

subplot(1,2,1)


% rms residual values to plot
lvls = [0, 1:24];


% Create filled contour
[c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop_5pct,:, :)./rms_uncert, length(r_bot_fine),...
    length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', 'k');
%clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold');

% % Create contour
% [c1,h1] = contour(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop_5pct,:, :)./rms_uncert, length(r_bot_fine),...
%     length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(20, 'fixed'));
% clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(20, 'fixed'));



% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create title
if use_l2_norm==false

    title(['RMS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_5pct)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

else

    title(['RSS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_5pct)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

end

idx_uncert = rms_residual./rms_uncert <= 1;
percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(rms_residual);
disp([newline,'Percent of state space within the convergence region using 5% measurement uncertainty: ',...
    num2str(100*percent_states_less_than_rms_uncert),'%', newline])


ylim([min(r_top_min(1)-r_bot_fine), max(r_top_min(1)-r_bot_fine)])

grid on; grid minor

% --- Create synthetic measurements with 1% uncertinaty ---
% -------------------------------------
% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
measurement_uncert = 0.01;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve

% Compute the new synethtic measurement with gaussian noise
% *** Gaussian noise can be either positive or negative. Meaning, an
% uncertainty of 5% implies the true value can lie anywhere between the
% measured value +/- 5% of the measured value
% define a
synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert/3) .*...
    randn(length(inputs.bands2run), 1);

% define the synthetic relfectance uncertainty
synthetic_measurement_uncert = measurement_uncert .* synthetic_measurement_with_noise;

% Using an exact modeled estimate without noise
use_l2_norm = true;

rms_residual = [];
rms_uncert = [];

if use_l2_norm==false

    % Compute the rms difference between the measurements and the modeled
    % reflectances
    rms_residual = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));


    % Compute the RMS of the synthetic measurement uncertainty
    rms_uncert = sqrt( mean( synthetic_measurement_uncert.^2));

else

    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rms_residual = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));

    % Compute the l2 norm of the synthetic measurement uncertainty
    rms_uncert = sqrt( sum( synthetic_measurement_uncert.^2));

end



% Find the states with the lowest rms residul

% find the smallest rms residual value, omitting nans
[~, idx_min] = min(rms_residual, [], 'all', 'omitnan');

r_top_min = R_top_fine(idx_min);


% define the optical depth slice you'd like to plot
idx_rTop_1pct = r_top_fine == r_top_min(1);



s2 = subplot(1,2,2);

% set subplot position
set(s2, 'Position', [0.51 0.11 0.334659090909091 0.815])




% Create filled contour
[c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop_1pct,:, :)./rms_uncert, length(r_bot_fine),...
    length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', 'k');
%clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold');

% Create colorbar
cb = colorbar();
% create colorbar label
ylabel(cb, '$\sqrt{ \Sigma{ \left(R(\vec{x}) - \vec{m} \right)^{2} }} / \sqrt{ \Sigma{ \left(\delta \vec{m} \right)^{2}}}$',...
    'FontSize', cb_font_size, 'Interpreter', 'latex')
clim([lvls(1), lvls(end)])



% % Create contour
% [c1,h1] = contour(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop_1pct,:, :)./rms_uncert, length(r_bot_fine),...
%     length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(20, 'fixed'));
% clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(20, 'fixed'));




% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create title
if use_l2_norm==false

    title(['RMS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_1pct)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

else

    title(['RSS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop_1pct)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 16);

end


ylim([min(r_top_min(1)-r_bot_fine), max(r_top_min(1)-r_bot_fine)])

grid on; grid minor



% set the figure size to be proportional to the length of the r_top and
% r_bot vectors

set(gcf, 'Position', [0 0 2400 1200])

idx_uncert = rms_residual./rms_uncert <= 1;
percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(rms_residual);
disp([newline,'Percent of state space within the convergence region using 1% measurement uncertainty: ',...
    num2str(100*percent_states_less_than_rms_uncert),'%', newline])



clear variables





%% Subplot comparing l2 norm residual for synthetic data bewteen 7 and 35 wavelengths employed in the retrieval

clear variables


% define axes label font size
axes_label_font_size = 40;

% define axes tick label font size
axes_tick_label_font_size = 25;

% define colorbar font size
cb_font_size = 40;

% define contour label size
contour_label_size = 25;

% LOAD DATA SET

load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-11-Dec-2024_rev1.mat')

% Create mesh grid
[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);


% Define synthetic model data

r_top_truth = 10.17;
r_bot_truth = 4.74;
tau_c_truth = 5.96;

[~, idx_r_top] = min(abs(r_top_fine - r_top_truth));
[~, idx_r_bot] = min(abs(r_bot_fine - r_bot_truth));
[~, idx_tau_c] = min(abs(tau_c_fine - tau_c_truth));

synthetic_measurement = reshape(Refl_model_fine(idx_r_top, idx_r_bot, idx_tau_c, :), [],1);


% --- Create synthetic measurements with 1% uncertinaty ---
% -------------------------------------
% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
measurement_uncert = 0.01;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve

% Compute the new synethtic measurement with gaussian noise
% *** Gaussian noise can be either positive or negative. Meaning, an
% uncertainty of 5% implies the true value can lie anywhere between the
% measured value +/- 5% of the measured value
% define a
synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert/3) .*...
    randn(length(inputs.bands2run), 1);

% define the synthetic relfectance uncertainty
synthetic_measurement_uncert = measurement_uncert .* synthetic_measurement_with_noise;


% Define 7 MODIS channels
% Let's now seperate out the interpolated relfectance at the seven MODIS
% wavelengths
wl_MODIS7_idx = [1, 4, 6, 7, 19, 23, 29];

% define the synthetic measurement
synthetic_measurement_with_noise_MODIS7 = synthetic_measurement_with_noise(wl_MODIS7_idx);
synthetic_measurement_uncert_MODIS7 = synthetic_measurement_uncert(wl_MODIS7_idx);

% Grab the modeled data at just the 7 MODIS bands
Refl_model_fine_MODIS7 = Refl_model_fine(:,:,:, wl_MODIS7_idx);



% Using an exact modeled estimate without noise
use_l2_norm = true;


if use_l2_norm==false

    % Compute the rms difference between the measurements and the modeled
    % reflectances
    rms_residual = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));

    rms_residual_MODIS7 = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise_MODIS7, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine_MODIS7).^2, 4));

    % Compute the RMS of the synthetic measurement uncertainty
    rms_uncert = sqrt( mean( synthetic_measurement_uncert.^2));

    % Compute the RMS of the synthetic measurement uncertainty using 7
    % MODIS channels
    rms_uncert_MODIS7 = sqrt( mean( synthetic_measurement_uncert_MODIS7.^2));


else

    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rms_residual = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));

    rms_residual_MODIS7 = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise_MODIS7, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine_MODIS7).^2, 4));

    % Compute the l2 norm of the synthetic measurement uncertainty
    rms_uncert = sqrt( sum( synthetic_measurement_uncert.^2));

    % Compute the l2 norm of the synthetic measurement uncertainty using 7
    % MODIS channels
    rms_uncert_MODIS7 = sqrt( sum( synthetic_measurement_uncert_MODIS7.^2));


end


% Find the states with the lowest rms residul
% find the smallest rms residual value, omitting nans
[~, idx_min] = min(rms_residual, [], 'all', 'omitnan');

r_top_min = R_top_fine(idx_min);



% Find the states with the lowest rms residul for the 7 MODIS wavelengths
% find the smallest rms residual value, omitting nans
[~, idx_min_MODIS7] = min(rms_residual_MODIS7, [], 'all', 'omitnan');

r_top_min_MODIS7 = R_top_fine(idx_min_MODIS7);




% ------------- PLOT l2 norm using 7 MODIS wavelengths ------------
% Create Contour plot of rms residual between 7 MODIS wavelengths and
% Synthetic measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----

% define the slice at some cloud top radii you'd like to plot
idx_rTop_MODIS7 = r_top_fine == r_top_min_MODIS7(1);

% Create figure
f = figure;

s1 = subplot(1,2,1);


% rms residual values to plot
lvls = [0, 1:24];


% % Create contour plot
% [c1,h1] = contour(tau_c_fine, r_top_min_MODIS7(1)-r_bot_fine, reshape(rms_residual_MODIS7(idx_rTop_MODIS7,:, :)./rms_uncert_MODIS7, length(r_bot_fine),...
%     length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(20, 'fixed'));
% clabel(c1,h1,'FontSize',contour_label_size,'FontWeight','bold', 'Color', mySavedColors(20, 'fixed'));


% Create filled contour plot
[c1,h1] = contourf(tau_c_fine, r_top_min_MODIS7(1)-r_bot_fine, reshape(rms_residual_MODIS7(idx_rTop_MODIS7,:, :)./rms_uncert_MODIS7, length(r_bot_fine),...
    length(tau_c_fine)),  lvls, 'LineWidth', 3, 'EdgeColor', 'k');
%clabel(c1,h1,'FontSize',contour_label_size,'FontWeight','bold', 'Color', mySavedColors(9, 'fixed'));


% Set tick label font size
ax = gca(f);
ax.FontSize = axes_tick_label_font_size;

% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create title
if use_l2_norm==false

    title(['RMS Residual at $r_{top} = $', num2str(r_top_fine(idx_rTop_MODIS7)),...
        ' - 7 Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty'],'Interpreter','latex', 'FontSize', 16);

else

    title(['RSS Residual at $r_{top} = $', num2str(r_top_fine(idx_rTop_MODIS7)),...
        ' - 7 Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty'],'Interpreter','latex', 'FontSize', 16);

end

grid on; grid minor



ylim([min(r_top_min(1)-r_bot_fine), max(r_top_min(1)-r_bot_fine)])


% ------------- PLOT l2 norm using 35 MODIS wavelengths ------------
% Create Contour plot of rms residual between 35 MODIS wavelengths and
% Synthetic measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----

% define the slice at some cloud top radii you'd like to plot
idx_rTop = r_top_fine == r_top_min(1);

s2 = subplot(1,2,2);

% set subplot position
set(s2, 'Position', [0.51 0.11 0.334659090909091 0.815])





% % Create contour
% [c1,h1] = contour(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop,:, :)./rms_uncert, length(r_bot_fine),...
%     length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', mySavedColors(20, 'fixed'));
% clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(20, 'fixed'));


% Create filled contour
[c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop,:, :)./rms_uncert, length(r_bot_fine),...
    length(tau_c_fine)),  lvls, 'LineWidth', 3, 'EdgeColor', 'k');
%clabel(c1,h1,'FontSize', contour_label_size,'FontWeight','bold', 'Color', mySavedColors(9, 'fixed'));

% Create colorbar
cb = colorbar();
% create colorbar label
ylabel(cb, '$\sqrt{ \Sigma{ \left(R(\vec{x}) - \vec{m} \right)^{2} }} / \sqrt{ \Sigma{ \left(\delta \vec{m} \right)^{2}}}$',...
    'FontSize', cb_font_size, 'Interpreter', 'latex')
clim([lvls(1), lvls(end)])




% Set tick label font size
ax = gca(f);
ax.FontSize = axes_tick_label_font_size;

% Create ylabel
ylabel('$r_{top}^{*} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', axes_label_font_size);

% Create title
if use_l2_norm==false

    title(['RMS Residual at $r_{top} = $', num2str(r_top_fine(idx_rTop)),...
        ' - 35 Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty'],'Interpreter','latex', 'FontSize', 16);

else

    title(['RSS Residual at $r_{top} = $', num2str(r_top_fine(idx_rTop)),...
        ' - 35 Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty'],'Interpreter','latex', 'FontSize', 16);

end

ylim([min(r_top_min(1)-r_bot_fine), max(r_top_min(1)-r_bot_fine)])

grid on; grid minor




%set the figure size
set(gcf, 'Position', [0 0 2400 1200])


idx_uncert_MODIS7 = rms_residual_MODIS7./rms_uncert_MODIS7 <= 1;
percent_states_less_than_rms_uncert_MODIS7 = sum(idx_uncert_MODIS7, 'all')/numel(rms_residual_MODIS7);
disp([newline,'Percent of state space within the convergence region using 7 MODIS channels: ',...
    num2str(100*percent_states_less_than_rms_uncert_MODIS7),'%', newline])



idx_uncert = rms_residual./rms_uncert <= 1;
percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(rms_residual);
disp([newline,'Percent of state space within the convergence region using 35 EMIT channels: ',...
    num2str(100*percent_states_less_than_rms_uncert),'%', newline])



clear variables


