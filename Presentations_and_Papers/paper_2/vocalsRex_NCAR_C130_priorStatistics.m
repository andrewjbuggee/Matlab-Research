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

    % --- all profiles, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    %     'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_02-Nov-2025.mat']);



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


% -----------------------------------------------------------------
% ------------------- For the covariance matrix -------------------
% -----------------------------------------------------------------
% We need 1 observation of cloud effective radius at top and bottom
% for each cloud optical depth and each above cloud precipitable
% water
re_bot_sample = zeros(length(ensemble_profiles), 1);
re_top_sample = zeros(length(ensemble_profiles), 1);

% define which levels you wish to sample 'cloud top', and which you
% wish to sample for 'cloud bottom'
cloud_top_lvl_sample = round(0.9 * n_bins);
cloud_bot_lvl_sample = round(0.1 * n_bins);


% remove all values beyond 30 microns. There are only a few and they all
% have a single count
re_top_upper_trim_limit = 30;

% remove all values below 4 microns. There are only a few and they all
% have a single count
re_top_lower_trim_limit = 4;


% remove all values beyond 15 microns. There are only a few and they all
% have a single count
re_bot_upper_trim_limit = 16;

% -----------------------------------------------------------------
% -----------------------------------------------------------------

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


        % ------------------- For the covariance matrix -------------------
        % We need 1 observation of cloud effective radius at top and bottom
        % for each cloud optical depth and each above cloud precipitable
        % water

        if bb == cloud_top_lvl_sample

            re_top_sample(nn) = vertically_segmented_attributes{bb, 1}(end);

            % check to see if this as the same value as the previous entry
            if nn>1 && re_bot_sample(nn) == re_top_sample(nn-1)
                % we want to sample another value
                % Sample another value for the cloud bottom effective radius
                if re_bot_sample(nn-1) ~= vertically_segmented_attributes{bb, 1}(end-1)

                    error([newline, 'There was no sample from this cloud region for this profile!', newline])

                    re_bot_sample(nn) = vertically_segmented_attributes{bb, 1}(end-1);

                end

            end


            % We don't want to include outliers. Check to make sure these
            % are within the predetermined bounds
            if re_top_sample(nn) < re_top_lower_trim_limit || re_top_sample(nn) > re_top_upper_trim_limit

                idx_used = find(index_segment);
                % take the last value and try to move towards more data
                % points, not less
                if (ensemble_profiles{nn}.altitude(end) - ensemble_profiles{nn}.altitude(1))>0
                    idx = idx_used(end)+1;
                else
                    idx = idx_used(end)-1;
                end
                re_top_sample(nn) = re(idx);

                while re_top_sample(nn) < re_top_lower_trim_limit || re_top_sample(nn) > re_top_upper_trim_limit

                    %  move towards more data
                    % points, not less
                    if (ensemble_profiles{nn}.altitude(end) - ensemble_profiles{nn}.altitude(1))>0
                        idx = idx+1;
                    else
                        idx = idx-1;
                    end

                    re_top_sample(nn) = re(idx);
                end


            end



        elseif bb == cloud_bot_lvl_sample

            re_bot_sample(nn) = vertically_segmented_attributes{bb, 1}(end);

            % check to see if this as the same value as the previous entry
            if nn>1 && re_bot_sample(nn) == re_top_sample(nn-1)
                % we want to sample another value
                % Sample another value for the cloud bottom effective radius
                if re_bot_sample(nn-1) ~= vertically_segmented_attributes{bb, 1}(end-1)

                    error([newline, 'There was no sample from this cloud region for this profile!', newline])

                    re_bot_sample(nn) = vertically_segmented_attributes{bb, 1}(end-1);

                end

            end

            % We don't want to include outliers. Check to make sure these
            % are within the predetermined bounds
            if re_bot_sample(nn) > re_bot_upper_trim_limit

                idx_used = find(index_segment);
                % take the last value and try to move towards more data
                % points, not less
                if (ensemble_profiles{nn}.altitude(end) - ensemble_profiles{nn}.altitude(1))>0
                    idx = idx_used(end)+1;
                else
                    idx = idx_used(end)-1;
                end
                re_bot_sample(nn) = re(idx);

                while re_bot_sample(nn) > re_bot_upper_trim_limit

                    %  move towards more data
                    % points, not less
                    if (ensemble_profiles{nn}.altitude(end) - ensemble_profiles{nn}.altitude(1))>0
                        idx = idx+1;
                    else
                        idx = idx-1;
                    end

                    re_bot_sample(nn) = re(idx);
                end
            end


        end




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
lgnd_fnt = 20;

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

% # of bins = 30
% cloudTop_idx = 26:30;
% cloudBot_idx = 6:10;

cloudTop_idx = 26:30;
cloudBot_idx = 1:5;

% cloudTop_idx = 27;
% cloudBot_idx = 3;


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
re_top_upper_trim_limit = 30;
re_top_ensemble_trimmed = re_top_ensemble;
re_top_ensemble_trimmed(re_top_ensemble_trimmed > re_top_upper_trim_limit) = [];

% remove all values below 4 microns. There are only a few and they all
% have a single count
re_top_lower_trim_limit = 4;
re_top_ensemble_trimmed(re_top_ensemble_trimmed < re_top_lower_trim_limit) = [];



% remove all values beyond 15 microns. There are only a few and they all
% have a single count
re_bot_upper_trim_limit = 15;
re_bot_ensemble_trimmed = re_bot_ensemble;
re_bot_ensemble_trimmed(re_bot_ensemble_trimmed > re_bot_upper_trim_limit) = [];



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




















%% Load Above cloud precipitable water data drom VOCALS-REx radisonde data


if strcmp(whatComputer,'anbu8374')==true


    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/paper2_prior_stats/'];


elseif strcmp(whatComputer,'andrewbuggee')==true


    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/radiosonde/paper2_prior_stats/'];

end

radiosonde = load([folderpath,...
    'precipitable_water_stats_for_paper2_combined_12-Nov-2025.mat']);

% --------------------------------------------------------
% ------- ABOVE CLOUD PRECIPITABLE WATER FITTING ---------
% --------------------------------------------------------
% Let's also fit these three distributions to the optical depth data
significance_lvl = 0.05;    % 5% risk of false rejection


% Kolmogorov test is better for this data set because of the outliers

% if p < significance_lvl:
%     h = 1 (true)  → REJECT the null hypothesis → Distribution is a BAD fit
% else:
%     h = 0 (false) → FAIL TO REJECT → Distribution is a GOOD fit (or at least acceptable)


% fit the number concentration data to a normal distribution
acpw_fit_normal = fitdist(radiosonde.combined_aboveCloud_pw_timeAndSpace, 'normal');
% [acpw_reject_normal, acpw_p_normal] = chi2gof(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
%     'CDF', acpw_fit_normal,'alpha', significance_lvl, 'NParams', 2);
[acpw_reject_normal, acpw_p_normal] = kstest(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
    'CDF', acpw_fit_normal,'alpha', significance_lvl);


% fit the number concentration content data to a log-normal distribution
acpw_fit_lognormal = fitdist(radiosonde.combined_aboveCloud_pw_timeAndSpace, 'lognormal');
% [acpw_reject_lognormal, acpw_p_lognormal] = chi2gof(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
%     'CDF', acpw_fit_lognormal,'alpha', significance_lvl, 'NParams', 2);
[acpw_reject_lognormal, acpw_p_lognormal] = kstest(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
    'CDF', acpw_fit_lognormal,'alpha', significance_lvl);


% fit the total number concentration data to a gamma distribution - use my custom
% libRadtran gamma distribution
acpw_fit_gamma = prob.GammaDistribution_libRadtran.fit(radiosonde.combined_aboveCloud_pw_timeAndSpace);
% [acpw_reject_gamma, acpw_p_gamma] = chi2gof(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
%     'CDF', acpw_fit_gamma,'alpha', significance_lvl, 'NParams', 2);
[acpw_reject_gamma, acpw_p_gamma] = kstest(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
    'CDF', acpw_fit_gamma,'alpha', significance_lvl);

% Plot results
lgnd_fnt = 20;

figure; histogram(radiosonde.combined_aboveCloud_pw_timeAndSpace,'NumBins', 30, 'Normalization','pdf')
hold on
xVals = linspace(min(radiosonde.combined_aboveCloud_pw_timeAndSpace),...
    max(radiosonde.combined_aboveCloud_pw_timeAndSpace), 1000);
plot(xVals, pdf(acpw_fit_normal, xVals))
plot(xVals, pdf(acpw_fit_lognormal, xVals))
plot(xVals, pdf(acpw_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
title('$acpw$ statistics and fits', ...
    'FontSize', 20, 'Interpreter', 'latex')







%%
% -------------------------------------------------------------
% ----- Make quantile-quantile plots of the data -----
% -------------------------------------------------------------

fnt_sz = 15;
mrkr_sz = 10;
line_width = 1.5;
line_width_2 = 2.5;
lgnd_fnt = 15;

figure;

% -- Linear state vector --

% --- Plot for the droplet size at cloud top ---
subplot(4,2,1)
qp1 = qqplot(re_top_ensemble_trimmed);
set(qp1(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp1(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('Effective radius at cloud top', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp1))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

% --- Plot for the droplet size at cloud bottom ---
subplot(4,2,3)
qp2 = qqplot(re_bot_ensemble_trimmed);
set(qp2(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp2(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('Effective radius at cloud bottom', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp2))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

% --- Plot for cloud optical depth ---
subplot(4,2,5)
qp3 = qqplot(tau_c);
set(qp3(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp3(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('Cloud optical depth', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp3))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

% --- Plot for the above cloud precipitable water ---
subplot(4,2,7)
qp4 = qqplot(radiosonde.combined_aboveCloud_pw_timeAndSpace);
set(qp4(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp4(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('above cloud precipitable water', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp4))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')






% -- log state vector --

% --- Plot for the droplet size at cloud top ---
subplot(4,2,2)
qp5 = qqplot(log(re_top_ensemble_trimmed));
set(qp5(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp5(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('$\ln($Effective radius at cloud top$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp5))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')


% --- Plot for the droplet size at cloud bottom ---
subplot(4,2,4)
qp6 = qqplot(log(re_bot_ensemble_trimmed));
set(qp6(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp6(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('$\ln($Effective radius at cloud bottom$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp6))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

% --- Plot for cloud optical depth ---
subplot(4,2,6)
qp7 = qqplot(log(tau_c));
set(qp7(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp7(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('$\ln($Cloud optical depth$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp7))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')


% --- Plot for the above cloud precipitable water ---
subplot(4,2,8)
qp8 = qqplot(log(radiosonde.combined_aboveCloud_pw_timeAndSpace));
set(qp8(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp8(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('$\ln($above cloud precipitable water$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp8))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

set(gcf, 'Position', [0,0, 1700, 950])




%% Compute the covariance matrix

% Each column represents the samples of a random variable, and each row are
% the observations

% For the covariance matrix, the number of samples need to be the same for
% each varaible. And if they have a relationship, it would be best to
% sample the different variables at the same time (or location, or whatever
% the independent variable is).

% I have 73 vertical profiles. How should I arrange the data to take 1
% value of r_top r_bot and optical depth?

prior_cov = cov([re_top_sample, re_bot_sample, tau_c, radiosonde.combined_aboveCloud_pw_timeAndSpace]);

% The prior covariance must be symmetric positive definite
try chol(prior_cov)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end


%% Save the prior covariance matrix and supporting variables to the paper_2 folder



if strcmp(whatComputer,'anbu8374')==true


    folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/'];



elseif strcmp(whatComputer,'andrewbuggee')==true


    folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/'];


end


combined_aboveCloud_pw_timeAndSpace = radiosonde.combined_aboveCloud_pw_timeAndSpace;

save([folderpath_2save,'prior_covarance_matrix_', char(datetime("today")),'.mat'],...
    'prior_cov', 're_top_sample', 're_bot_sample', 'tau_c', 'combined_aboveCloud_pw_timeAndSpace')


%% Clear variables

% there are too many
clear variables