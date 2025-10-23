% Find the best distribution fit for the histogram of droplet distribution
% data

% Try three different fits: normal, gamma and lognormal

% input the droplet distribution data for a single vertical profile. Code
% will step through each measurement and fit the three distributions above
% to the data and compute the R^2 values

% INPUTS:
% -------
%       (1) vert_profs - vertical profiles cell array found from the
%       find_verticalProfiles_VOCALS_REx function

%       (2) indexes2plot - the indexes associated with the vertical
%       profiles you wish to plot

%       (3) index_altitude - this is the index within the vertical profile
%       defining the altitude to plot. ** This must be the same length as
%       indexes2plot. There should be a unique index_altitude for each
%       index2plot

%       (3) radius_limits - the range of radii you wish to display on this
%       histrogram. This will trim the x-axis. Input as [r_min, r_max] in
%       units of microns


% By Andrew John Buggee

%%


function [h_norm, p_norm, h_logNorm, p_logNorm, h_gamma, p_gamma] = find_bestFitDist_dropDist(droplet_dist,...
    binEdges, binCenters, significance_lvl)


%% Make double

binCenters = double(binCenters);

binEdges = double(binEdges);

%%

% length of the profile - number of distribution measurements
data_length = size(droplet_dist,2);

% set up hypothesis testing vectors
h_norm = zeros(1, data_length);
h_logNorm = zeros(1, data_length);
h_gamma = zeros(1, data_length);

p_norm = zeros(1, data_length);
p_logNorm = zeros(1, data_length);
p_gamma = zeros(1, data_length);

% Values of 0 cause the lognoraml fit to fail. Set these to some small
% value greater than 0
droplet_dist(droplet_dist==0) = 1e-10;

for zz = 1:data_length


    % normalize distribution
    
    % The chi-squared method checks the differences between the observed
    % counts in each bin with the expected counts in each bin. The p-value
    % is the probability of getting the chi-squared value we got by chance.
    % Then, we check to see if the p-value is less than or equal to the
    % significance value. The significance value is the probability that we
    % reject the fit, but it's actually a good one:

        % If p-value < α: Reject the null hypothesis (the distribution is NOT a good fit)
        % If p-value ≥ α: Fail to reject (the distribution is a reasonable fit)


    %% Fit a normal distribution and compute the chi^2 goodness of fit



    % gauss_fit = fitdist(droplet_dist(:,zz), 'normal');
    %
    % [h, p, st] = chi2gof(droplet_dist(:,zz),'CDF',gauss_fit);


    gauss_fit = fitdist(droplet_dist(:,zz), 'normal');

    expected_counts_gauss = sum(droplet_dist(:,zz)) * pdf(gauss_fit, binCenters);

    % There are two free parameters that define the normal distribution
    [h_norm(zz), p_norm(zz)] = chi2gof(binCenters,'Ctrs',binCenters,'Frequency', droplet_dist(:,zz), ...
        'Expected',expected_counts_gauss, 'Alpha', significance_lvl,...
        'NParams',2);

    %% Fit a lognormal distribution and compute the R^2 value


    logNormal_fit = fitdist(droplet_dist(:,zz), 'Lognormal');

    expected_counts_logNorm = sum(droplet_dist(:,zz)) * pdf(logNormal_fit, binCenters);

    % There are two free parameters that define the lognormal distribution
    [h_logNorm(zz), p_logNorm(zz)] = chi2gof(binCenters,'Ctrs',binCenters,'Frequency', droplet_dist(:,zz), ...
        'Expected',expected_counts_logNorm, 'Alpha', significance_lvl,...
        'NParams',2);



    %% Fit a gamma distribution and compute the R^2 value

    gamma_fit = fitdist(droplet_dist(:,zz), 'Gamma');

    expected_counts_gamma = sum(droplet_dist(:,zz)) * pdf(gamma_fit, binCenters);

    % There are two free parameters that define the Gamma distribution
    [h_gamma(zz), p_gamma(zz)] = chi2gof(binCenters,'Ctrs',binCenters,'Frequency', droplet_dist(:,zz), ...
        'Expected',expected_counts_gamma, 'Alpha', significance_lvl,...
        'NParams',2);



end






end

