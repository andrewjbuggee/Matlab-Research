% Find the best distribution fit for the histogram of droplet distribution
% data

% Try three different fits: normal, gamma and lognormal

% input the droplet distribution data for a single vertical profile. Code
% will step through each measurement and fit the three distributions above
% to the data and compute the R^2 values

% INPUTS:
% -------
%       (1) droplet_dist -

%       (2) bin edges -

%       (3) binCenters -

%       (3) significance_lvl - a value of 0.05 means you're willing to accept a 5%
%           chance of incorrectly rejecting a distribution that's actually a good fit (Type I error)


% OUTPUTS:
% --------
%       (1) h_horm - % h = 0 means fail to reject (good fit)
%                   h = 1 means reject (poor fit)


% By Andrew John Buggee

%%


function [normFit, logNormFit, gammaFit] = find_bestFitDist_dropDist(droplet_dist,...
    binEdges, binCenters, significance_lvl)


%% Make double

binCenters = double(binCenters)';

binEdges = double(binEdges)';

%%

% length of the profile - number of distribution measurements
data_length = size(droplet_dist,2);

% set up output vectors
h_norm = zeros(data_length, 1); 
p_norm = zeros(data_length, 1);
norm_mean = zeros(data_length, 1); 
norm_std = zeros(data_length, 1);



h_logNorm = zeros(data_length, 1); 
p_logNorm = zeros(data_length, 1); 
logNorm_mean = zeros(data_length, 1); 
logNorm_std = zeros(data_length, 1); 



h_gamma = zeros(data_length, 1); 
p_gamma = zeros(data_length, 1); 
% gamma_shape = zeros(data_length, 1); 
% gamma_scale = zeros(data_length, 1); 
gamma_rEff = zeros(data_length, 1); 
gamma_alpha = zeros(data_length, 1); 


% Values of 0 cause the lognoraml fit to fail. Set these to some small
% value greater than 0
droplet_dist(droplet_dist==0) = 1e-10;

for zz = 1:data_length


    % First, recreate the raw samples from the histogram data

    % Reconstruct individual samples by replicating each bin center
    samples = [];
    for ii = 1:length(binCenters)
        % Round counts to integers and replicate bin centers
        numSamples = round(droplet_dist(ii, zz));
        samples = [samples; repmat(binCenters(ii), numSamples, 1)];
    end

    % sum_samples = sum(samples);
    num_total_droplets = sum(round(droplet_dist(:,zz)));


    % make sure the number of observations is greater than 10

    if length(samples)>10

        % The chi-squared method checks the differences between the observed
        % counts in each bin with the expected counts in each bin. The p-value
        % is the probability of getting the chi-squared value we got by chance.
        % Then, we check to see if the p-value is less than or equal to the
        % significance value. The significance value is the probability that we
        % reject the fit, but it's actually a good one:

        % If p-value < α: Reject the null hypothesis (the distribution is NOT a good fit)
        % If p-value ≥ α: Fail to reject (the distribution is a reasonable fit)


        %% Fit a normal distribution and compute the chi^2 goodness of fit


        % ----- Using the reconstructed samples -----
        gauss_fit = [];
        gauss_fit = fitdist(samples, 'normal');

        % expected_counts_gauss = num_total_droplets * pdf(gauss_fit, binCenters);

        % There are two free parameters that define the normal distribution
        % [h_norm(zz), p_norm(zz)] = chi2gof(binCenters,'Ctrs',binCenters,'Frequency', droplet_dist(:,zz), ...
        %     'Expected',expected_counts_gauss, 'Alpha', significance_lvl,...
        %     'NParams',2);
        [h_norm(zz), p_norm(zz)] = chi2gof(samples, 'CDF', gauss_fit, 'Alpha', significance_lvl,...
            'NParams',2);

        % store all normal fit output
        norm_mean(zz) = gauss_fit.mu;
        norm_std(zz) = gauss_fit.std; 




        %% Fit a lognormal distribution and compute the R^2 value

        % ----- Using the reconstructed samples -----
        logNormal_fit = [];
        logNormal_fit = fitdist(samples, 'Lognormal');

        % expected_counts_logNorm = num_total_droplets * pdf(logNormal_fit, binCenters);

        % There are two free parameters that define the lognormal distribution
        % [h_logNorm(zz), p_logNorm(zz)] = chi2gof(binCenters,'Ctrs',binCenters,'Frequency', droplet_dist(:,zz), ...
        %     'Expected',expected_counts_logNorm, 'Alpha', significance_lvl,...
        %     'NParams',2);
        [h_logNorm(zz), p_logNorm(zz)] = chi2gof(samples, 'CDF', logNormal_fit, 'Alpha', significance_lvl,...
            'NParams',2);

        % Store all lognormal fit output
        logNorm_mean(zz) = logNormal_fit.mu;
        logNorm_std(zz) = logNormal_fit.sigma;




        %% Fit a gamma distribution and compute the R^2 value

        % ----- Using the reconstructed samples -----
        gamma_fit = [];
        % gamma_fit = fitdist(samples, 'Gamma');
        gamma_fit = prob.GammaDistribution_libRadtran.fit(samples);

        % expected_counts_gamma = num_total_droplets * pdf(gamma_fit, binCenters);

        % There are two free parameters that define the Gamma distribution
        % [h_gamma(zz), p_gamma(zz)] = chi2gof(binCenters,'Ctrs',binCenters,'Frequency', droplet_dist(:,zz), ...
        %     'Expected',expected_counts_gamma, 'Alpha', significance_lvl,...
        %     'NParams',2);

        [h_gamma(zz), p_gamma(zz)] = chi2gof(samples, 'CDF', gamma_fit, 'Alpha', significance_lvl,...
            'NParams',2);


        % Store all gamma fit output
        % gamma_shape(zz) = gamma_fit.a;
        % gamma_scale(zz) = gamma_fit.b;
        gamma_rEff(zz) = gamma_fit.r_eff;
        gamma_alpha(zz) = gamma_fit.alpha;

        % alpha values higher than 150 cause issues because of the
        % enormously high value of the normalizing coefficient
        % Gamma(150) = 3.8e260
        % Limit alpha values to 150
        if gamma_fit.alpha>150

            gamma_alpha(zz) = 150;

        end


        %% Plot the three fits against the histogram data

        % plot the data as a histogram and overlay the distribution fit
        % figure; histogram('BinEdges',binEdges', 'BinCounts',droplet_dist(:,zz), 'Normalization','pdf')
        % hold on;
        % 
        % x_values = linspace(min(samples), max(samples), 500);
        % 
        % y_values_gauss = pdf(gauss_fit, x_values);
        % plot(x_values, y_values_gauss, 'r-', 'LineWidth', 2);
        % 
        % y_values_logNorm = pdf(logNormal_fit, x_values);
        % plot(x_values, y_values_logNorm, 'b-', 'LineWidth', 2);
        % 
        % y_values_gamma = pdf(gamma_fit, x_values);
        % plot(x_values, y_values_gamma, 'g-', 'LineWidth', 2);
        % 
        % title('Data and Fitted Normal Distribution');
        % xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex');
        % ylabel('Probability Density');
        % 
        % legend('observations', 'normal', 'log-normal', 'gamma',...
        %     'Interpreter','latex', 'Location','northwest', 'FontSize', 20,...
        %     'Color', 'white', 'TextColor', 'k')
        % 
        % xlim([0, 30])


    else

        % There are not enough observations to fit a distribution
        % set output values to nan
        h_norm(zz) = NaN;
        p_norm(zz) = NaN;

        h_logNorm(zz) = NaN;
        p_logNorm(zz) = NaN;

        h_gamma(zz) = NaN;
        p_gamma(zz) = NaN;

        norm_mean(zz) = NaN;
        norm_std(zz) = NaN;

        logNorm_mean(zz) = NaN;
        logNorm_std(zz) = NaN;

        gamma_alpha(zz) = NaN;
        gamma_rEff(zz) = NaN;


    end


end


%% Collect all output into 3 structures



% Store all normal fit output
normFit.mean = norm_mean;
normFit.std = norm_std;

% store the h-value and p-value for the normal fit
normFit.h_test = h_norm;
normFit.p_value = p_norm;


% Store all lognormal fit output
logNormFit.mean = logNorm_mean;
logNormFit.std = logNorm_std;

% store the h-value and p-value for the lognormal fit
logNormFit.h_test = h_logNorm;
logNormFit.p_value = p_logNorm;



% Store all gamma fit output
% gammaFit.shape = gamma_shape;
% gammaFit.scale = gamma_scale;
gammaFit.rEff = gamma_rEff;
gammaFit.alpha = gamma_alpha;

% store the h-value and p-value for the gamma fit
gammaFit.h_test = h_gamma;
gammaFit.p_value = p_gamma;





end

