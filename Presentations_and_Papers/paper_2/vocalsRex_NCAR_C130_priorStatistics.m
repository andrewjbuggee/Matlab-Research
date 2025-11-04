% Analysis of all vertical profiles from non-precipitating clouds measured
% during the VOCALS-REx campaign on-board the NCAR C130 aircract

% By Andrew John Buggee


%% Plot the ensemble distribution of effective radius and optical depth
% for non-precipitating clouds.
% ----- For Vertical Profiles -----


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true


    % Location of ensemble data
    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    % 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_30-Oct-2025_rev1.mat'])

    load([folderpath,...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_02-Nov-2025.mat'])



elseif strcmp(whatComputer, 'andrewbuggee')==true

    % Location of ensemble data
    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
        'vocals_rex_data/NCAR_C130/SPS_1/'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    % 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_30-Oct-2025_rev1.mat'])

    load([folderpath,...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_02-Nov-2025.mat'])


end





%%
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Segment re, LWC, Nc into N bins along optical depth to normalize each
%  profile

% In order to compute a mean vertical profile, we have to first normalize
% the vertical extent so that all profiles lie between values [0,1]. Then
% we break up the vertical component in n discrete bins.

n_bins = 30; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
vertically_segmented_attributes = cell(n_bins, 3);

% Store all cloud optical depths
tau_c = zeros(length(ensemble_profiles),1);


normalized_altitude = cell(1, length(ensemble_profiles));




for nn = 1:length(ensemble_profiles)

    % first we need to normalize the vertical component of all profiles
    normalized_altitude{nn} = (ensemble_profiles{nn}.altitude - min(ensemble_profiles{nn}.altitude))./...
        (max(ensemble_profiles{nn}.altitude) - min(ensemble_profiles{nn}.altitude));

    % the data is stored in altitude space.

    if ensemble_profiles{nn}.flag_2DC_data_is_conforming==true

        re = ensemble_profiles{nn}.re;
    else

        % Since we are analyzing non-precipitating clouds, as long as the
        % 2DC LWP is small, we can simply take the CDP data as the
        % effective radius data
        if ensemble_profiles{nn}.lwp_2DC<5

            re = ensemble_profiles{nn}.re_CDP;

        end

    end

    lwc = ensemble_profiles{nn}.lwc;
    Nc_total = sum(ensemble_profiles{nn}.Nc, 1);

    % check to see if any values are 0
    if sum(re==0)>=1

        % This happens sometimes in multi layer clouds. Set it to be a
        % small non-zero value
        re(re==0) = 0.01; % microns;

    end

    if sum(lwc==0)>=1

        % This happens sometimes in multi layer clouds. Set it to be a
        % small non-zero value
        lwc(lwc==0) = 0.001; % microns;

    end

    if sum(Nc_total==0)>=1

        % This happens sometimes in multi layer clouds. Set it to be a
        % small non-zero value
        Nc_total(Nc_total==0) = 0.1; % microns;

    end



    % Also store the cloud optical thickness for each profile
    tau_c(nn,1) = max(ensemble_profiles{nn}.tau);

    % Check if there is ever a 0
    if tau_c(nn,1)==0

        error('Why?')

    end




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
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; reshape(re(index_segment),[],1)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; reshape(lwc(index_segment), [],1)];

        % store the total number concentration
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; reshape(Nc_total(index_segment), [],1)];



    end



end



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% We know want to determine the statistics of the effective radius, lwc and
% total number concentration as a function of height within the cloud. At
% each normalized level, let's fit a normal, lognormal and gamma function
% to the distribution and see which is the best fit by check the
% chi-squared test

significance_lvl = 0.05;

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

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [re_reject_normal(bb), re_p_normal(bb)] = kstest(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
        'alpha', significance_lvl);

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'lognormal');

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [re_reject_lognormal(bb), re_p_lognormal(bb)] = kstest(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
        'alpha', significance_lvl);

    % fit the effective radius data to a gamma distribution - use my custom
    % libRadtran gamma distribution
    re_fit_gamma(bb) = prob.GammaDistribution_libRadtran.fit(vertically_segmented_attributes{bb,1});

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [re_reject_gamma(bb), re_p_gamma(bb)] = kstest(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
        'alpha', significance_lvl);

    % make plot of normal log-normal and gamma fits
    % figure; subplot(1,3,1); plot(re_fit_normal(bb)); title('Normal Fit'); xlabel('r_e (\mum)')
    % subplot(1,3,2); plot(re_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('r_e (\mum)')
    % subplot(1,3,3); plot(re_fit_gamma(bb)); title('Gamma Fit'); xlabel('r_e (\mum)')
    % set(gcf, 'Position', [0 0 1200 500])











end



% -------------------------------------------
% ----------- OPTICAL DEPTH FITTING ---------
% -------------------------------------------
% Let's also fit these three distributions to the optical depth data

% fit the number concentration data to a normal distribution
tauC_fit_normal = fitdist(tau_c, 'normal');
[tauC_reject_normal, tauC_p_normal] = chi2gof(tau_c, 'CDF', tauC_fit_normal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the number concentration content data to a log-normal distribution
tauC_fit_lognormal = fitdist(tau_c, 'lognormal');
[tauC_reject_lognormal, tauC_p_lognormal] = chi2gof(tau_c, 'CDF', tauC_fit_lognormal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the total number concentration data to a gamma distribution - use my custom
% libRadtran gamma distribution
tauC_fit_gamma = prob.GammaDistribution_libRadtran.fit(tau_c);
[tauC_reject_gamma, tauC_p_gamma] = chi2gof(tau_c, 'CDF', tauC_fit_gamma,...
    'alpha', significance_lvl, 'NParams', 2);

% Plot results
figure; histogram(tau_c,'NumBins', 30, 'Normalization','pdf')
hold on
xVals = linspace(min(tau_c), max(tau_c), 1000);
plot(xVals, pdf(tauC_fit_normal, xVals))
plot(xVals, pdf(tauC_fit_lognormal, xVals))
plot(xVals, pdf(tauC_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')
title('$\tau_{C}$ statistics and fits', ...
        'FontSize', 20, 'Interpreter', 'latex')












% Now let's find the where the hypothesis was not rejected (reject_ = 0)
% which means the chi-squared test is confident in the choice of
% distribution to within 5% uncertainty

bin_names = {'Normal', 'Log-Normal', 'Gamma'};
cloudTop_idx = 21:30;
cloudBot_idx = 1:10;
% ------------------------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS AT CLOUD TOP FITTING ------
% ------------------------------------------------------------
num_rejects_reTop = sum([re_reject_normal(cloudTop_idx); re_reject_lognormal(cloudTop_idx);...
    re_reject_gamma(cloudTop_idx)], 2, "omitnan");

figure; histogram('Categories', bin_names, 'BinCounts', num_rejects_reTop);
title('Number of distribution rejections of $r_{top}$ data', 'interpreter', 'latex'); ylabel('Counts')



% ---------------------------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS AT CLOUD BOTTOM FITTING ------
% ---------------------------------------------------------------
num_rejects_reBot = sum([re_reject_normal(cloudBot_idx); re_reject_lognormal(cloudBot_idx);...
    re_reject_gamma(cloudBot_idx)], 2, "omitnan");

figure; histogram('Categories', bin_names, 'BinCounts', num_rejects_reBot);
title('Number of distribution rejections of $r_{bot}$ data', 'interpreter', 'latex'); ylabel('Counts')



% ------------------------------------------
% ------- CLOUD OPTICAL DEPTH FITTING ------
% ------------------------------------------

figure; histogram('Categories', bin_names, 'BinCounts', [tauC_reject_normal, tauC_reject_lognormal, tauC_reject_gamma]);
title('Number of distribution rejections of $r_{bot}$ data', 'interpreter', 'latex'); ylabel('Counts')


% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------



%% What if we fit statistics for cloud the ensemble of all effective radius values at the
% upper regions of the cloud for cloud top, and the lower regions of the
% cloud for cloud bottom?

lgnd_fnt = 20;

% *** COMBINING DIFFERENT VERTICAL BINS LEADS TO A REJECTION OF ALL FIT
% TYPES ***

cloudTop_idx = 26:30;
cloudBot_idx = 6:10;
significance_lvl = 0.1;



re_top_ensemble = [];
re_bot_ensemble = [];

for bb = 1:length(cloudTop_idx)

    re_top_ensemble = [re_top_ensemble; vertically_segmented_attributes{cloudTop_idx(bb), 1}];

    re_bot_ensemble = [re_bot_ensemble; vertically_segmented_attributes{cloudBot_idx(bb), 1}];

end

% ------------------------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS AT CLOUD TOP FITTING ------
% ------------------------------------------------------------

% fit the effective radius data to a normal distribution
re_top_fit_normal = fitdist(re_top_ensemble, 'normal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_top_reject_normal, re_top_p_normal] = kstest(re_top_ensemble, 'CDF', re_top_fit_normal,...
    'alpha', significance_lvl);


% fit the effective radius data to a log-normal distribution
re_top_fit_lognormal = fitdist(re_top_ensemble, 'lognormal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_top_reject_lognormal, re_top_p_lognormal] = kstest(re_top_ensemble, 'CDF',re_top_fit_lognormal,...
    'alpha', significance_lvl);


% fit the effective radius data to a gamma distribution - use my custom
% libRadtran gamma distribution
re_top_fit_gamma = prob.GammaDistribution_libRadtran.fit(re_top_ensemble);

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_top_reject_gamma, re_top_p_gamma] = kstest(re_top_ensemble, 'CDF', re_top_fit_gamma,...
    'alpha', significance_lvl);


% Plot results
figure; histogram(re_top_ensemble,'NumBins', 100, 'Normalization','pdf')
hold on
xVals = linspace(min(re_top_ensemble), max(re_top_ensemble), 1000);
plot(xVals, pdf(re_top_fit_normal, xVals))
plot(xVals, pdf(re_top_fit_lognormal, xVals))
plot(xVals, pdf(re_top_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')
title('$r_{top}$ statistics and fits', ...
        'FontSize', 20, 'Interpreter', 'latex')




% ------------------------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS AT CLOUD BOT FITTING ------
% ------------------------------------------------------------

% fit the effective radius data to a normal distribution
re_bot_fit_normal = fitdist(re_bot_ensemble, 'normal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_bot_reject_normal, re_bot_p_normal] = kstest(re_bot_ensemble, 'CDF', re_bot_fit_normal,...
    'alpha', significance_lvl);


% fit the effective radius data to a log-normal distribution
re_bot_fit_lognormal = fitdist(re_bot_ensemble, 'lognormal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_bot_reject_lognormal, re_bot_p_lognormal] = kstest(re_bot_ensemble, 'CDF',re_bot_fit_lognormal,...
    'alpha', significance_lvl);


% fit the effective radius data to a gamma distribution - use my custom
% libRadtran gamma distribution
re_bot_fit_gamma = prob.GammaDistribution_libRadtran.fit(re_bot_ensemble);

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_bot_reject_gamma, re_bot_p_gamma] = kstest(re_bot_ensemble, 'CDF', re_bot_fit_gamma,...
    'alpha', significance_lvl);


% Plot results
figure; histogram(re_bot_ensemble,'NumBins', 100, 'Normalization','pdf')
hold on
xVals = linspace(min(re_bot_ensemble), max(re_bot_ensemble), 1000);
plot(xVals, pdf(re_bot_fit_normal, xVals))
plot(xVals, pdf(re_bot_fit_lognormal, xVals))
plot(xVals, pdf(re_bot_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')
title('$r_{bot}$ statistics and fits', ...
        'FontSize', 20, 'Interpreter', 'latex')








%% Try the same analysis as above by trim away the outliers
% outliers are defined as counts of 1
% [N, edges] = histcounts(re_top_ensemble, 'BinEdges', [0:1:round(max(re_top_ensemble))+1]);
% re_top_median = median(re_top_ensemble, 'all');

% for re at cloud top
% remove all values beyond 30 microns. There are only a few and they all
% have a single count
re_top_ensemble_trimmed = re_top_ensemble;
re_top_ensemble_trimmed(re_top_ensemble_trimmed>30) = [];
% remove all values below 5 microns. There are only a few and they all
% have a single count
re_top_ensemble_trimmed(re_top_ensemble_trimmed<5) = [];



% remove all values beyond 15 microns. There are only a few and they all
% have a single count
re_bot_ensemble_trimmed = re_bot_ensemble;
re_bot_ensemble_trimmed(re_bot_ensemble_trimmed>15) = [];



significance_lvl = 0.04;



% ------------------------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS AT CLOUD TOP FITTING ------
% ------------------------------------------------------------

% fit the effective radius data to a normal distribution
re_top_fit_normal = fitdist(re_top_ensemble_trimmed, 'normal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_top_reject_normal, re_top_p_normal] = kstest(re_top_ensemble_trimmed, 'CDF', re_top_fit_normal,...
    'alpha', significance_lvl);


% fit the effective radius data to a log-normal distribution
re_top_fit_lognormal = fitdist(re_top_ensemble_trimmed, 'lognormal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_top_reject_lognormal, re_top_p_lognormal] = kstest(re_top_ensemble_trimmed, 'CDF',re_top_fit_lognormal,...
    'alpha', significance_lvl);


% fit the effective radius data to a gamma distribution - use my custom
% libRadtran gamma distribution
re_top_fit_gamma = prob.GammaDistribution_libRadtran.fit(re_top_ensemble_trimmed);

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_top_reject_gamma, re_top_p_gamma] = kstest(re_top_ensemble_trimmed, 'CDF', re_top_fit_gamma,...
    'alpha', significance_lvl);


% Plot results
figure; histogram(re_top_ensemble_trimmed,'NumBins', 100, 'Normalization','pdf')
hold on
xVals = linspace(min(re_top_ensemble_trimmed), max(re_top_ensemble_trimmed), 1000);
plot(xVals, pdf(re_top_fit_normal, xVals))
plot(xVals, pdf(re_top_fit_lognormal, xVals))
plot(xVals, pdf(re_top_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')
title('$r_{top}$ trimmed statistics and fits', ...
        'FontSize', 20, 'Interpreter', 'latex')




% ------------------------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS AT CLOUD BOT FITTING ------
% ------------------------------------------------------------

% fit the effective radius data to a normal distribution
re_bot_fit_normal = fitdist(re_bot_ensemble_trimmed, 'normal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_bot_reject_normal, re_bot_p_normal] = kstest(re_bot_ensemble_trimmed, 'CDF', re_bot_fit_normal,...
    'alpha', significance_lvl);


% fit the effective radius data to a log-normal distribution
re_bot_fit_lognormal = fitdist(re_bot_ensemble_trimmed, 'lognormal');

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_bot_reject_lognormal, re_bot_p_lognormal] = kstest(re_bot_ensemble_trimmed, 'CDF',re_bot_fit_lognormal,...
    'alpha', significance_lvl);


% fit the effective radius data to a gamma distribution - use my custom
% libRadtran gamma distribution
re_bot_fit_gamma = prob.GammaDistribution_libRadtran.fit(re_bot_ensemble_trimmed);

% *** Chi-squared doesn't work well ***
% chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
% which is more robust to outliers
% [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
%     'alpha', significance_lvl, 'NParams', 2);
[re_bot_reject_gamma, re_bot_p_gamma] = kstest(re_bot_ensemble_trimmed, 'CDF', re_bot_fit_gamma,...
    'alpha', significance_lvl);


% Plot results
figure; histogram(re_bot_ensemble_trimmed,'NumBins', 100, 'Normalization','pdf')
hold on
xVals = linspace(min(re_bot_ensemble_trimmed), max(re_bot_ensemble_trimmed), 1000);
plot(xVals, pdf(re_bot_fit_normal, xVals))
plot(xVals, pdf(re_bot_fit_lognormal, xVals))
plot(xVals, pdf(re_bot_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
             'Color', 'white', 'TextColor', 'k')
title('$r_{bot}$ trimmed statistics and fits', ...
        'FontSize', 20, 'Interpreter', 'latex')







%%

% Compile statistics using the best fit distributions


% ---- most common best fit distribution for r_top ---
if num_rejects_reTop(1)<num_rejects_reTop(2) && num_rejects_reTop(1)<num_rejects_reTop(3)

    % normal fit was best


elseif num_rejects_reTop(2)<num_rejects_reTop(1) && num_rejects_reTop(2)<num_rejects_reTop(3)

    % lognormal fit was best


elseif num_rejects_reTop(3)<num_rejects_reTop(1) && num_rejects_reTop(3)<num_rejects_reTop(2)

    % gamma fit was best


elseif num_rejects_reTop(2)<num_rejects_reTop(1) && num_rejects_reTop(2)==num_rejects_reTop(3)

    % lognormal and gamma fits had same number of rejections. Use
    % lognormal






end










%%
% -------------------------------------------------------------
% ----- Make a subplot of r_top, r_bot and tau_c profiles -----
% -------------------------------------------------------------


axes_fnt_sz = 20;
title_fnt_sz = 25;


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
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex', 'Fontsize',axes_fnt_sz)
ylabel('Normalized Altitude', 'Interpreter', 'latex', 'FontSize',axes_fnt_sz)

% set x axis boundaries
xlim([2, 20])                   % microns






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
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex', 'FontSize', axes_fnt_sz)
ylabel('Normalized Altitude', 'Interpreter','latex', 'FontSize', axes_fnt_sz)

% set x axis boundaries
xlim([0, 0.6])                   % g/cm^3

% set the figure title
title(['Mean Profiles:  $LWC \geq$', num2str(ensemble_profiles{nn}.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles{nn}.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex', 'FontSize', title_fnt_sz)





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
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex', "FontSize", axes_fnt_sz)
ylabel('Normalized Altitude', 'Interpreter', 'latex', 'FontSize', axes_fnt_sz)

% set x axis boundaries
xlim([0, 320])                   % cm^(-3)

% set the size of the figure
set(gcf, 'Position', [0 0 1255 700])


%%

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



