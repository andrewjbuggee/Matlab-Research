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

% n_bins = 30; % number of segments the noramlized vertical component is broken up into
% ** USE 20 BINS FOR SAME NUMBER OF VERTICAL LAYERS AS RETRIEVAL FORWARD MODEL **
n_bins = 20; % number of segments the noramlized vertical component is broken up into


bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for:
%   (1) droplet size
%   (2) liquid water content
%   (3) total number concentration
%   (4) drop distribution alpha paramter (related to effective variance)

% ** Vertically_segmented_attributes are defined base to top **
% The first cell is at and near cloud base. 
% The final cell is at or near cloud top
vertically_segmented_attributes = cell(n_bins, 4);

% Store all cloud optical depths
tau_c = zeros(length(ensemble_profiles),1);

% Store all cloud top heights and geometric depths
cloudTopHeight = zeros(length(ensemble_profiles),1);
cloudDepth = zeros(length(ensemble_profiles),1);

% store the effective variance
eff_variance = cell(n_bins, 1);



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
    alpha_param = ensemble_profiles{nn}.gammaFit.alpha;

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

    if sum(alpha_param==0)>=1 

        % Sometimes there isn't a good fit. Set it to the default
        % libRadtran value of 7
        alpha_param(alpha_param==0) = 7;

    end

    if sum(isnan(alpha_param))>=1

        % Sometimes there isn't a good fit. Use the previous value
        idx_nan = find(isnan(alpha_param));
        idx_nonNan = find(~isnan(alpha_param));

        for aa = 1:numel(idx_nan)

            % find closest non-nan index
            [~, min_idx] = min( abs( idx_nan(aa) - idx_nonNan ));
                
             alpha_param(idx_nan(aa)) = alpha_param(idx_nonNan(min_idx));
             
        end

    end








    % store the cloud optical thickness for each profile
    tau_c(nn,1) = max(ensemble_profiles{nn}.tau);


    % store the cloud optical thickness for each profile
    cloudTopHeight(nn,1) = max(ensemble_profiles{nn}.altitude);  % meters


    % store the cloud optical thickness for each profile
    cloudDepth(nn,1) = ensemble_profiles{nn}.cloud_depth;        % meters

    % Check if there is ever a 0
    if tau_c(nn,1)==0

        error('Why?')

    elseif cloudTopHeight(nn,1)==0

        error('Why?')

    elseif cloudDepth(nn,1)==0

        error('Why?')

    end

    % -------------------------------------------------------------------
    % *********************** IMPORTANT!!!! *****************************
    % -------------------------------------------------------------------
    % The vertically segmented profiles organize data from cloud base
    % to cloud top. It starts with bin edge at 0 and goes to 1
    % -------------------------------------------------------------------
    % -------------------------------------------------------------------
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

        % store the gamma droplet distirbution alpha parameter (related to
        % the effective variance)
        vertically_segmented_attributes{bb, 4} = [vertically_segmented_attributes{bb, 4}; reshape(alpha_param(index_segment), [],1)];



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




alpha_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
alpha_p_normal = zeros(1, size(vertically_segmented_attributes,1));

alpha_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
alpha_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

alpha_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
alpha_p_gamma = zeros(1, size(vertically_segmented_attributes,1));




effVariance_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
effVariance_p_normal = zeros(1, size(vertically_segmented_attributes,1));

effVariance_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
effVariance_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

effVariance_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
effVariance_p_gamma = zeros(1, size(vertically_segmented_attributes,1));



% Step through vertically_segmented_attributes from cloud base to cloud top
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











    % -----------------------------------------------------------
    % ------- ALPHA PARAMETER (EFFECTIVE VARIANCE) FITTING ------
    % -----------------------------------------------------------


    % fit the effective radius data to a normal distribution
    alpha_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,4}, 'normal');

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [alpha_reject_normal(bb), alpha_p_normal(bb)] = kstest(vertically_segmented_attributes{bb,4}, 'CDF', alpha_fit_normal(bb),...
        'alpha', significance_lvl);


    % fit the effective radius data to a log-normal distribution
    alpha_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,4}, 'lognormal');

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [alpha_reject_lognormal(bb), alpha_p_lognormal(bb)] = kstest(vertically_segmented_attributes{bb,4}, 'CDF', alpha_fit_lognormal(bb),...
        'alpha', significance_lvl);


    % fit the effective radius data to a gamma distribution - use my custom
    % libRadtran gamma distribution
    alpha_fit_gamma(bb) = prob.GammaDistribution_libRadtran.fit(vertically_segmented_attributes{bb,4});

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [alpha_reject_gamma(bb), alpha_p_gamma(bb)] = kstest(vertically_segmented_attributes{bb,4}, 'CDF', alpha_fit_gamma(bb),...
        'alpha', significance_lvl);








    % -----------------------------------------------------------
    % --------------- TRUE EFFECTIVE VARIANCE FITTING -----------
    % -----------------------------------------------------------
    

    % ** store the alpha parameter as the effective variance! **
    eff_variance{bb} = 1./(vertically_segmented_attributes{bb,4} + 3);


    % fit the effective radius data to a normal distribution
    effVariance_fit_normal(bb) = fitdist(eff_variance{bb}, 'normal');

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    [effVariance_reject_normal(bb), effVariance_p_normal(bb)] = kstest(eff_variance{bb}, 'CDF', effVariance_fit_normal(bb),...
        'alpha', significance_lvl);


    % fit the effective radius data to a log-normal distribution
    effVariance_fit_lognormal(bb) = fitdist(eff_variance{bb}, 'lognormal');

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    [effVariance_reject_lognormal(bb), effVariance_p_lognormal(bb)] = kstest(eff_variance{bb}, 'CDF', effVariance_fit_lognormal(bb),...
        'alpha', significance_lvl);


    % fit the effective radius data to a gamma distribution - use my custom
    % libRadtran gamma distribution
    effVariance_fit_gamma(bb) = prob.GammaDistribution_libRadtran.fit(eff_variance{bb});

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    [effVariance_reject_gamma(bb), effVariance_p_gamma(bb)] = kstest(eff_variance{bb}, 'CDF', effVariance_fit_gamma(bb),...
        'alpha', significance_lvl);


        






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





% -------------------------------------------
% -------------- CLOUD TOP HEIGHT -----------
% -------------------------------------------
% Let's also fit these three distributions to the optical depth data

% fit the number concentration data to a normal distribution
cloudTopHeight_fit_normal = fitdist(cloudTopHeight, 'normal');
[cloudTopHeight_reject_normal, cloudTopHeight_p_normal] = chi2gof(cloudTopHeight, 'CDF', cloudTopHeight_fit_normal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the number concentration content data to a log-normal distribution
cloudTopHeight_fit_lognormal = fitdist(cloudTopHeight, 'lognormal');
[cloudTopHeight_reject_lognormal, cloudTopHeight_p_lognormal] = chi2gof(cloudTopHeight, 'CDF', cloudTopHeight_fit_lognormal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the total number concentration data to a gamma distribution - use my custom
% libRadtran gamma distribution
cloudTopHeight_fit_gamma = prob.GammaDistribution_libRadtran.fit(cloudTopHeight);
[cloudTopHeight_reject_gamma, cloudTopHeight_p_gamma] = chi2gof(cloudTopHeight, 'CDF', cloudTopHeight_fit_gamma,...
    'alpha', significance_lvl, 'NParams', 2);

% Plot results
lgnd_fnt = 20;

figure; histogram(cloudTopHeight,'NumBins', 30, 'Normalization','pdf')
hold on
xVals = linspace(min(cloudTopHeight), max(cloudTopHeight), 1000);
plot(xVals, pdf(cloudTopHeight_fit_normal, xVals))
plot(xVals, pdf(cloudTopHeight_fit_lognormal, xVals))
plot(xVals, pdf(cloudTopHeight_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
title('Cloud Top Height statistics and fits', ...
    'FontSize', 20, 'Interpreter', 'latex')





% -------------------------------------------
% ----------- CLOUD DEPTH FITTING ---------
% -------------------------------------------
% Let's also fit these three distributions to the cloud depth data

% fit the number concentration data to a normal distribution
cloudDepth_fit_normal = fitdist(cloudDepth, 'normal');
[cloudDepth_reject_normal, cloudDepth_p_normal] = chi2gof(cloudDepth, 'CDF', cloudDepth_fit_normal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the number concentration content data to a log-normal distribution
cloudDepth_fit_lognormal = fitdist(cloudDepth, 'lognormal');
[cloudDepth_reject_lognormal, cloudDepth_p_lognormal] = chi2gof(cloudDepth, 'CDF', cloudDepth_fit_lognormal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the total number concentration data to a gamma distribution - use my custom
% libRadtran gamma distribution
cloudDepth_fit_gamma = prob.GammaDistribution_libRadtran.fit(cloudDepth);
[cloudDepth_reject_gamma, cloudDepth_p_gamma] = chi2gof(cloudDepth, 'CDF', cloudDepth_fit_gamma,...
    'alpha', significance_lvl, 'NParams', 2);

% Plot results
lgnd_fnt = 20;

figure; histogram(cloudDepth,'NumBins', 30, 'Normalization','pdf')
hold on
xVals = linspace(min(cloudDepth), max(cloudDepth), 1000);
plot(xVals, pdf(cloudDepth_fit_normal, xVals))
plot(xVals, pdf(cloudDepth_fit_lognormal, xVals))
plot(xVals, pdf(cloudDepth_fit_gamma, xVals))
grid on; grid minor
legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
title('Cloud Depth statistics and fits', ...
    'FontSize', 20, 'Interpreter', 'latex')












% Now let's find the where the hypothesis was not rejected (reject_ = 0)
% which means the chi-squared test is confident in the choice of
% distribution to within 5% uncertainty

% ** Cloud base is the first idx **
% ** Cloud top is the final idx **

bin_names = {'Normal', 'Log-Normal', 'Gamma'};
if n_bins==30

    cloudTop_idx = 21:30;
    cloudBot_idx = 1:10;

elseif n_bins == 20
    cloudTop_idx = 16:20;
    cloudBot_idx = 1:5;

end

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
% 
% figure; histogram('Categories', bin_names, 'BinCounts', [tauC_reject_normal, tauC_reject_lognormal, tauC_reject_gamma]);
% title('Number of distribution rejections of $\tau_{C}$ data', 'interpreter', 'latex'); ylabel('Counts')


% ------------------------------------------
% ---------- CLOUD TOP HEIGHT FITTING ------
% ------------------------------------------
% 
% figure; histogram('Categories', bin_names, 'BinCounts', [cloudTopHeight_reject_normal, cloudTopHeight_reject_lognormal,...
%     cloudTopHeight_reject_gamma]);
% title('Number of distribution rejections of cloud top height data', 'interpreter', 'latex'); ylabel('Counts')


% ------------------------------------------
% ----------- CLOUD DEPTH FITTING ----------
% ------------------------------------------

% figure; histogram('Categories', bin_names, 'BinCounts', [cloudDepth_reject_normal, cloudDepth_reject_lognormal, cloudDepth_reject_gamma]);
% title('Number of distribution rejections of cloud depth data', 'interpreter', 'latex'); ylabel('Counts')
% 


% ---------------------------------------------------------------
% ------- ALPHA PARAMETER (EFFECTIVE VARIANCE) FITTING ------
% ---------------------------------------------------------------
num_rejects_alpha = sum([alpha_reject_normal; alpha_reject_lognormal;...
    alpha_reject_gamma], 2, "omitnan");

figure; histogram('Categories', bin_names, 'BinCounts', num_rejects_alpha);
title('Number of distribution rejections of $\alpha$ data', 'interpreter', 'latex'); ylabel('Counts')





% ---------------------------------------------------------------
% ------------------- EFFECTIVE VARIANCE FITTING ----------------
% ---------------------------------------------------------------
num_rejects_effVar = sum([effVariance_reject_normal; effVariance_reject_lognormal;...
    effVariance_reject_gamma], 2, "omitnan");

figure; histogram('Categories', bin_names, 'BinCounts', num_rejects_effVar);
title('Number of distribution rejections of $v_e$ data', 'interpreter', 'latex'); ylabel('Counts')


% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------




%% fit a PDF to all the alpha parameter data and determine what range, from 1 to X, results in 75% of the measurements

all_alpha = [];
for nn = 1:size(vertically_segmented_attributes, 1)
    all_alpha = [all_alpha; vertically_segmented_attributes{nn, 4}];
end

all_alpha_fit = fitdist(all_alpha, 'lognormal');
[all_alpha_reject, all_alpha_p] = kstest(all_alpha, 'CDF', all_alpha_fit,...
        'alpha', significance_lvl);

figure; histogram(all_alpha, 100);
title('All Alpha Parameter measurements', 'interpreter', 'latex'); ylabel('Counts')
figure; plot(all_alpha_fit)
title('PDF from all alpha parameter measurements', 'interpreter', 'latex'); ylabel('PDF')

% Plot the CDF
% create a range of x values from the minimum to the maximum of your data
x = linspace(min(all_alpha), max(all_alpha), 1000);
y_cdf = cdf(all_alpha_fit, x); % Calculate the CDF values

plot(x, y_cdf);
xlabel('Alpha Parameter');
ylabel('Cumulative Probability (F(x))');
title('Cumulative Distribution Function');
grid on; grid minor

% At alpha = 50, we've covered ~85% of all observations.
% At alpha = 30, we've covered ~68% of all observations.





%% What if we fit statistics for cloud the ensemble of all effective radius values at the
% upper regions of the cloud for cloud top, and the lower regions of the
% cloud for cloud bottom?

lgnd_fnt = 20;

% *** COMBINING DIFFERENT VERTICAL BINS LEADS TO A REJECTION OF ALL FIT
% TYPES ***

% # of bins = 30
% cloudTop_idx = 26:30;
% cloudBot_idx = 6:10;


if n_bins==30

    cloudTop_idx = 26:30;
    cloudBot_idx = 1:5;

elseif n_bins == 20

    cloudTop_idx = 18:20;
    cloudBot_idx = 1:3;

end

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
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *

%% Load Above cloud precipitable water data from VOCALS-REx radisonde data
% ------------------------------------------------------------------------

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

xlabel('$pw_{ac}$ ($mm$)', 'Interpreter','latex', 'FontSize', lgnd_fnt+3)
ylabel('Counts', 'Interpreter','latex', 'FontSize', lgnd_fnt+3)

legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
title('$acpw$ statistics and fits', ...
    'FontSize', 20, 'Interpreter', 'latex')








%%
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *
% *


%% Load Above cloud precipitable water data from ERA5 Reanalysis Data
% ------------------------------------------------------------------------

if strcmp(whatComputer,'anbu8374')==true


    folderpath_era5 = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/'];


elseif strcmp(whatComputer,'andrewbuggee')==true

    folderpath_era5 = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/'];

end

% load ERA5 data
era5 = load([folderpath_era5,...
    'ERA5_profiles_closest_to_73-VR_profiles_01-Feb-2026.mat']);



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
acpw_fit_normal = fitdist(era5.above_cloud_pw_usingVR, 'normal');
% [acpw_reject_normal, acpw_p_normal] = chi2gof(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
%     'CDF', acpw_fit_normal,'alpha', significance_lvl, 'NParams', 2);
[acpw_reject_normal, acpw_p_normal] = kstest(era5.above_cloud_pw_usingVR,...
    'CDF', acpw_fit_normal,'alpha', significance_lvl);


% fit the number concentration content data to a log-normal distribution
acpw_fit_lognormal = fitdist(era5.above_cloud_pw_usingVR, 'lognormal');
% [acpw_reject_lognormal, acpw_p_lognormal] = chi2gof(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
%     'CDF', acpw_fit_lognormal,'alpha', significance_lvl, 'NParams', 2);
[acpw_reject_lognormal, acpw_p_lognormal] = kstest(era5.above_cloud_pw_usingVR,...
    'CDF', acpw_fit_lognormal,'alpha', significance_lvl);


% fit the total number concentration data to a gamma distribution - use my custom
% libRadtran gamma distribution
acpw_fit_gamma = prob.GammaDistribution_libRadtran.fit(era5.above_cloud_pw_usingVR);
% [acpw_reject_gamma, acpw_p_gamma] = chi2gof(radiosonde.combined_aboveCloud_pw_timeAndSpace,...
%     'CDF', acpw_fit_gamma,'alpha', significance_lvl, 'NParams', 2);
[acpw_reject_gamma, acpw_p_gamma] = kstest(era5.above_cloud_pw_usingVR,...
    'CDF', acpw_fit_gamma,'alpha', significance_lvl);

% Plot results
lgnd_fnt = 20;

figure; histogram(era5.above_cloud_pw_usingVR,'NumBins', 30, 'Normalization','pdf')
hold on
xVals = linspace(min(era5.above_cloud_pw_usingVR),...
    max(era5.above_cloud_pw_usingVR), 1000);
plot(xVals, pdf(acpw_fit_normal, xVals))
plot(xVals, pdf(acpw_fit_lognormal, xVals))
plot(xVals, pdf(acpw_fit_gamma, xVals))
grid on; grid minor

xlabel('$pw_{ac}$ ($mm$)', 'Interpreter','latex', 'FontSize', lgnd_fnt+3)
ylabel('Counts', 'Interpreter','latex', 'FontSize', lgnd_fnt+3)

legend('data', 'normal fit', 'lognormal fit', 'gamma fit', 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
title('$acpw$ from ERA5 w/ VR Adjustment', ...
    'FontSize', 20, 'Interpreter', 'latex')












%%
% -------------------------------------------------------------
% ----- Make quantile-quantile plots of the data -----
% -------------------------------------------------------------

fnt_sz = 25;
mrkr_sz = 10;
line_width = 1.5;
line_width_2 = 2.5;
lgnd_fnt = 25;

fig1 = figure;

% -- Linear state vector --

% --- Plot for the droplet size at cloud top ---
subplot(4,2,1)
qp1 = qqplot(re_top_ensemble_trimmed);
set(qp1(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp1(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
xlabel('')
% ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('')
title('Effective radius at cloud top', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp1))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------

% --- Plot for the droplet size at cloud bottom ---
subplot(4,2,3)
qp2 = qqplot(re_bot_ensemble_trimmed);
set(qp2(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp2(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
xlabel('')
% ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('')
title('Effective radius at cloud bottom', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp2))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------

% --- Plot for cloud optical depth ---
subplot(4,2,5)
qp3 = qqplot(tau_c);
set(qp3(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp3(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
xlabel('')
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('Cloud optical depth', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp3))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------

% --- Plot for the above cloud precipitable water ---
% % ** using radiosonde data **
% subplot(4,2,7)
% qp4 = qqplot(radiosonde.combined_aboveCloud_pw_timeAndSpace);
% set(qp4(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
% set(qp4(3), 'LineStyle', '--', 'LineWidth', line_width_2);
% grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
% % ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('')
% title('Above-Cloud Precipitable Water', 'Interpreter','latex', 'FontSize', fnt_sz)
% % compute the R^2 value from the figure handle and print this in the legend
% legend(['$R^2 = $', num2str(compute_qqplot_R2(qp4))], 'location',...
%     'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
%     'Color', 'white', 'TextColor', 'k')

% ** using ERA5 data with VR adjustment **
subplot(4,2,7)
qp4 = qqplot(era5.above_cloud_pw_usingVR);
set(qp4(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp4(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('')
title('Above Cloud Precipitable Water', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp4))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------





% -- log state vector --

% --- Plot for the droplet size at cloud top ---
subplot(4,2,2)
qp5 = qqplot(log(re_top_ensemble_trimmed));
set(qp5(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp5(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
xlabel('')
% ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('')
title('$\ln($Effective radius at cloud top$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp5))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------


% --- Plot for the droplet size at cloud bottom ---
subplot(4,2,4)
qp6 = qqplot(log(re_bot_ensemble_trimmed));
set(qp6(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp6(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
xlabel('')
% ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('')
title('$\ln($Effective radius at cloud bottom$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp6))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------

% --- Plot for cloud optical depth ---
subplot(4,2,6)
qp7 = qqplot(log(tau_c));
set(qp7(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp7(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
xlabel('')
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('$\ln($Cloud optical depth$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp7))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------


% --- Plot for the above cloud precipitable water ---
% % ** using radiosonde data **
% subplot(4,2,8)
% qp4 = qqplot(log(radiosonde.combined_aboveCloud_pw_timeAndSpace));
% set(qp4(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
% set(qp4(3), 'LineStyle', '--', 'LineWidth', line_width_2);
% grid on; grid minor;
% xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
% % ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('')
% title('$\ln($Above Cloud Precipitable Water$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% % compute the R^2 value from the figure handle and print this in the legend
% legend(['$R^2 = $', num2str(compute_qqplot_R2(qp4))], 'location',...
%     'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
%     'Color', 'white', 'TextColor', 'k')

% ** using ERA5 data with VR adjustment **
subplot(4,2,8)
qp4 = qqplot(log(era5.above_cloud_pw_usingVR));
set(qp4(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp4(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
% ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('')
title('$\ln($Above Cloud Precipitable Water$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp4))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')
% ------------------------------------------------------------------------


set(gcf, 'Position', [0,0, 1700, 950])


% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
        error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig1,[folderpath_figs,'Quantile-Quantile plot for all 4 variables.fig']);


% save .png with 400 DPI resolution
% remove title
exportgraphics(fig1,[folderpath_figs,'Quantile-Quantile plot for all 4 variables.jpg'],'Resolution', 500);
% -------------------------------------
% -------------------------------------





% create a seperate plot for cloud top height and cloud depth
figure; 

% -- linear state vector --

% --- Plot for cloud top height ---
subplot(2,2,1)
qp3 = qqplot(cloudTopHeight);
set(qp3(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp3(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('Cloud top height', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp3))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

% --- Plot for cloud Depth ---
subplot(2,2,2)
qp3 = qqplot(cloudDepth);
set(qp3(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp3(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('Cloud Depth', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp3))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

% -- log state vector --

% --- Plot for cloud top height ---
subplot(2,2,3)
qp3 = qqplot(log(cloudTopHeight));
set(qp3(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp3(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title(['$\ln($Cloud Top Height$)$'], 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp3))], 'location',...
    'best','Interpreter','latex', 'Location','best', 'FontSize', lgnd_fnt,...
    'Color', 'white', 'TextColor', 'k')

% --- Plot for cloud Depth ---
subplot(2,2,4)
qp3 = qqplot(log(cloudDepth));
set(qp3(1), 'MarkerSize', mrkr_sz, 'LineWidth', line_width); % Update for top ensemble
set(qp3(3), 'LineStyle', '--', 'LineWidth', line_width_2);
grid on; grid minor;
xlabel('Standard Normal Quantiles', 'Interpreter','latex', 'FontSize', fnt_sz)
ylabel('Quantiles of Input Sample', 'Interpreter','latex', 'FontSize', fnt_sz)
title('$\ln($Cloud Depth$)$', 'Interpreter','latex', 'FontSize', fnt_sz)
% compute the R^2 value from the figure handle and print this in the legend
legend(['$R^2 = $', num2str(compute_qqplot_R2(qp3))], 'location',...
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

% ****!!!! Using Radiosonde data !!!!***
% prior_cov_lin = cov([re_top_sample, re_bot_sample, tau_c, radiosonde.combined_aboveCloud_pw_timeAndSpace]);
% prior_cov_log = cov(log([re_top_sample, re_bot_sample, tau_c, radiosonde.combined_aboveCloud_pw_timeAndSpace]));


% ****!!!! Using ERA5 data w/ VR adjustment !!!!***
prior_cov_lin = cov([re_top_sample, re_bot_sample, tau_c, era5.above_cloud_pw_usingVR]);
prior_cov_log = cov(log([re_top_sample, re_bot_sample, tau_c, era5.above_cloud_pw_usingVR]));


% prior_cov_lin_noACPW = cov([re_top_sample, re_bot_sample, tau_c]);
% 
% prior_cov_log_noACPW = cov(log([re_top_sample, re_bot_sample, tau_c]));


fnt_sz = 25;

% The prior covariance must be symmetric positive definite
try chol(prior_cov_lin)
    disp('Matrix is symmetric positive definite.')

    % plot heat maps of the covariance matrix
    % Plot heat map for linear covariance matrix
    fig2 = figure;
    imagesc(prior_cov_lin);
    colormap("hot")
    cb = colorbar;
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = fnt_sz;
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fnt_sz;


    title('Heat Map of Linear Covariance Matrix', 'Interpreter', 'latex', 'FontSize', fnt_sz);
    % define the variable names along the x and y axis
    % Define variable names for the heat map axes
    % variableNames = {'Effective Radius Top', 'Effective Radius Bottom', 'Optical Depth', 'Above Cloud PW'};
    variableNames = {'$r_{top}$', '$r_{bot}$', '$\tau_c$', '$pw_{ac}$'};
    set(gca, 'XTick', 1:length(variableNames), 'XTickLabel', variableNames, 'YTick', 1:length(variableNames),...
        'YTickLabel', variableNames, 'TickLabelInterpreter', 'latex');
    
    % display to covariance matrix values on the heat map
    % Display covariance matrix values on the heat map
    textStrings = num2str(prior_cov_lin(:), '%.2f'); % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings)); % Remove any space padding
    [xPos, yPos] = meshgrid(1:size(prior_cov_lin, 2), 1:size(prior_cov_lin, 1)); % Create x and y coordinates
    hStrings = text(xPos(:), yPos(:), textStrings(:), 'HorizontalAlignment', 'center'); % Create text objects
    midValue = mean(get(gca, 'CLim')); % Get the middle value of the color range
    textColors = repmat(prior_cov_lin(:) > midValue, 1, 3); % Choose white or black for the text color
    set(hStrings, {'Color'}, num2cell(textColors, 2)); % Change the color of the text
    % update the font size of the covaraince values displayed on the heat
    % map
    % Update the font size of the covariance values displayed on the heat map
    set(hStrings, 'FontSize', fnt_sz);
    
    % the color of the covaraince matrix values in the last two columns of
    % the last two rows needs to be white to show up better against the
    % background color
    % The color of the covariance matrix values in the last two columns of the last two rows needs to be white to show up better against the background color
    
    % set(hStrings(end-1:end, :), 'Color', 'white');
    % set(hStrings(end-4, :), 'Color', 'white');
    set(hStrings(:, :), 'Color', 'white');
   
    % The color of the covaraince matrix value for the 3rd column and 3rd
    % row needs to be black to show up better against the background
    set(hStrings(end-5, :), 'Color', 'black');

    % ** Paper Worthy **
    % -------------------------------------
    % ---------- Save figure --------------
    % save .fig file
    if strcmp(whatComputer,'anbu8374')==true
        error(['Where do I save the figure?'])
    elseif strcmp(whatComputer,'andrewbuggee')==true
        folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
    end
    saveas(fig2,[folderpath_figs,'Linear a prioiri covariance matrix using acpw from ERA5 with VR adjustment.fig']);


    % save .png with 400 DPI resolution
    % remove title
    title('');
    exportgraphics(fig2,[folderpath_figs,'Linear a prioiri covariance matrix using acpw from ERA5 with VR adjustment.jpg'],...
        'Resolution', 500);
    % -------------------------------------
    % -------------------------------------




    

    % Plot heat map for logarithmic covariance matrix
    fig3 = figure;
    imagesc(prior_cov_log);
    colormap("hot")
    cb = colorbar;
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = fnt_sz;
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fnt_sz;


    title('Heat Map of Logarithmic Covariance Matrix', 'Interpreter', 'latex', 'FontSize', fnt_sz);
    % define the variable names along the x and y axis
    % Define variable names for the heat map axes
    % variableNames = {'log(Effective Radius Top)', 'log(Effective Radius Bottom)',...
    %     'log(Optical Depth)', 'log(Above Cloud PW)'};
    variableNames = {'$\ln{(r_{top})}$', '$\ln{(r_{bot})}$', '$\ln{(\tau_c)}$', '$\ln{(pw_{ac})}$'};
    set(gca, 'XTick', 1:length(variableNames), 'XTickLabel', variableNames, 'YTick', 1:length(variableNames),...
        'YTickLabel', variableNames, 'TickLabelInterpreter', 'latex');
    
    % display to covariance matrix values on the heat map
    % Display covariance matrix values on the heat map
    textStrings = num2str(prior_cov_log(:), '%.3f'); % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings)); % Remove any space padding
    [xPos, yPos] = meshgrid(1:size(prior_cov_log, 2), 1:size(prior_cov_log, 1)); % Create x and y coordinates
    hStrings = text(xPos(:), yPos(:), textStrings(:), 'HorizontalAlignment', 'center'); % Create text objects
    midValue = mean(get(gca, 'CLim')); % Get the middle value of the color range
    textColors = repmat(prior_cov_log(:) > midValue, 1, 3); % Choose white or black for the text color
    set(hStrings, {'Color'}, num2cell(textColors, 2)); % Change the color of the text
    % update the font size of the covaraince values displayed on the heat
    % map
    % Update the font size of the covariance values displayed on the heat map
    set(hStrings, 'FontSize', fnt_sz);
    % the color of the covaraince matrix values in the last two columns of
    % the last two rows needs to be white to show up better against the
    % background color
    % The color of the covariance matrix values in the last two columns of the last two rows needs to be white to show up better against the background color
    
    % set(hStrings(end-1, :), 'Color', 'white');
    % set(hStrings(end-4, :), 'Color', 'white');
    % set(hStrings(end, :), 'Color', 'black');

    set(hStrings(:), 'Color', 'white');
    set(hStrings(11), 'Color', 'black');
    set(hStrings(end), 'Color', 'black');

    % The color of the covaraince matrix value for the 3rd column and 3rd
    % row needs to be black to show up better against the background
    set(hStrings(end-5, :), 'Color', 'black');

    % ** Paper Worthy **
    % -------------------------------------
    % ---------- Save figure --------------
    % save .fig file
    if strcmp(whatComputer,'anbu8374')==true
        error(['Where do I save the figure?'])
    elseif strcmp(whatComputer,'andrewbuggee')==true
        folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
    end
    saveas(fig3,[folderpath_figs,'Logarithmic a prioiri covariance matrix using acpw from ERA5 with VR adjustment.fig']);


    % save .png with 400 DPI resolution
    % remove title
    title('');
    exportgraphics(fig3,[folderpath_figs,'Logarithmic a prioiri covariance matrix using acpw from ERA5 with VR adjustment.jpg'],'Resolution', 500);
    % % -------------------------------------
    % -------------------------------------





catch ME
    disp('Matrix is not symmetric positive definite')
end


%% Save the prior covariance matrix and supporting variables to the paper_2 folder



if strcmp(whatComputer,'anbu8374')==true


    folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/'];



elseif strcmp(whatComputer,'andrewbuggee')==true


    folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/'];


end


% acpw = radiosonde.combined_aboveCloud_pw_timeAndSpace;
acpw = era5.above_cloud_pw_usingVR;

% save([folderpath_2save,'prior_covarance_matrix_', char(datetime("today")),'.mat'],...
%     'prior_cov_lin', 'prior_cov_log', 'prior_cov_lin_noACPW', "prior_cov_log_noACPW",'re_top_sample', 're_bot_sample',...
%     'tau_c', 'acpw')

save([folderpath_2save,'prior_covarance_matrix_', char(datetime("today")),'.mat'],...
    'prior_cov_lin', 'prior_cov_log','re_top_sample', 're_bot_sample',...
    'tau_c', 'acpw')

%% Save variables for the forward model covariance matrix

% ** save the cloud top height **
% save([folderpath_2save,'VR_cloud_top_height_obs_', char(datetime("today")),'.mat'],...
%     'cloudTopHeight', 'cloudDepth')

% ** save the effective variance (alpha parameter) values for each prof **
% alpha parameter (effective variance) best fits a log-normal distribution
% (0 rejects) compared to the gamma dist (9 rejects)

% Effective variance best fits a log-normal distribution or a gamma
% distribution (1 reject and 0 rejects respectively)
% save([folderpath_2save,'VR_effective_variance_at_normalized_altitudes_20-levels_', char(datetime("today")),'.mat'],...
%     'alpha_fit_lognormal', 'alpha_fit_normal', 'effVariance_fit_lognormal', 'effVariance_fit_normal')


%% Clear variables

% there are too many
clear variables