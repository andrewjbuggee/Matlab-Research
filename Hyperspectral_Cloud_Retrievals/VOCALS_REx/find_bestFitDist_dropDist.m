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


function [fits, goodness] = find_bestFitDist_dropDist(droplet_dist)


% length of the profile - number of distribution measurements
data_length = size(droplet_dist,2);

for zz = 1:data_length



    %% Fit a normal distribution and compute the chi^2 goodness of fit 

    gauss_fit = fitdist(droplet_dist(:,zz), 'normal');

    
    %% Fit a lognormal distribution and compute the R^2 value




    %% Fit a gamma distribution and compute the R^2 value




end




end

