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


function [goodness_of_fits] = find_bestFitDist_dropDist_2(droplet_dist, binEdges, binCenters)


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

    % Use the R^2 value to determine the best fit


    %% Fit a normal distribution and compute the R^2 goodness of fit


    [gauss_fit, gof_gauss] = fit(binCenters', droplet_dist(:,zz), 'gauss1');


    %% Fit a lognormal distribution and compute the R^2 value


    [logNormal_fit, gof_logNormal] = fit(binCenters', droplet_dist(:,zz), 'Lognormal');




    %% Fit a gamma distribution and compute the R^2 value

    [gamma_fit, gof_gauss] = fit(binCenters', droplet_dist(:,zz), 'Gamma');




end






end

