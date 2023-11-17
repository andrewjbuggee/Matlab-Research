% Fit a lognormal particle size distribution to a droplet distribution data set

% INPUTS:
%   (1) nr_data - this is a vector of the droplet distribution data. The
%   units must by #/cm^(-3)/micron. It must be the same length as the r
%   vector

%   (2) r - this is the independent vector associated with the distribution
%   data nr_data. It must be in units of microns.

%   

% OUTPUTS:
%   (1) nr_fit - this is a vector of the best fit found. It represents the
%   droplet distribution data using the best fit parameters

%   (2) fit_parameters - the modal radius, the effective varaiance, and the RMS 
%   difference between the data and the best fit


% By Andrew John Buggee

%%

function [nr_fit, fit_parameters] = fit_lognormal_size_distribution_kokhanovsky(nr_data, r)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=2
    error([newline,'Not enough inputs. Need 3: droplet effective radius, total number concentration, and the alpha parameter value.', newline])
end

% Check to make sure r is the same length as nr_data

if length(r)~=length(nr_data)

    error([newline,'The inputs r and nr_data must be the same length', newline])           % says who?
end



%% using the independent vector provide, compute the gamma distribution for many alpha and r0 values
% This definition follows the definition of A. Kokhonovsky in 'Cloud
% Optics' on page 5

% make sure r and nr_data are column vectors
r = reshape(r, [],1);
nr_data = reshape(nr_data, [], 1);

% Compute the total number of droplets
N0 = trapz(r, nr_data);



% % lets find the modal radius
% [~,idx] = max(nr_data);
% r_modal = r(idx);
% 
% 
% std_dev_fit = logspace(-2,0,100);            % values of mu to test
% 
% N = 1./(sqrt(2*pi) * std_dev_fit);  % normalization constant
% 
% 
% % Compute the fit for the first 100 values of mu
% nr_fits = zeros(numel(r), numel(std_dev_fit));
% 
% rms_diff = zeros(1, numel(std_dev_fit));
% 
% for nn = 1:numel(std_dev_fit)
%     
%     nr_fits(:, nn) = N0 * N(nn) * 1./r .* exp(-log(r./r_modal).^2 ./(2*std_dev_fit(nn)^2));       % lognormal droplet distribution
% 
%     % Find the least squares difference between the fit true data
% 
%     rms_diff(nn) = sqrt(sum((nr_data - nr_fits(:,nn)).^2, 1));
% 
% end
% 
% 
% 
% % find the mu value associated witht he minimum rms difference
% [min_rms, min_idx] = min(rms_diff);
% 
% % Grab the mu value and the nr_fit associated with the minimum rms 
% std_dev = std_dev_fit(min_idx);
% nr_fit = nr_fits(:, min_idx);













% ----- For varrying r_modal and mu_fit --------

% define a course vector of sigma values, we will interpolate on a finer grid later
std_dev_fit = 1:100;            % values of sigma to test

% compute the log_normal distribution for different values of sigma and r_modal
r_modal = 1:60;        % microns

% use this to run a loop with changing r_modal values
[Sigma, R_modal] = meshgrid(std_dev_fit, r_modal);


rms_diff = zeros(numel(r_modal), numel(std_dev_fit));

for nn = 1:numel(r_modal)
    
    N = 1./(sqrt(2*pi) * std_dev_fit);  % normalization constant

    nr_fits = N0 * N(nn) * 1./r .* exp(-log(r./r_modal(nn)).^2 ./(2*std_dev_fit.^2));       % lognormal droplet distribution

    % Find the least squares difference between the fit true data

    rms_diff(nn,:) = sqrt(sum(( repmat(nr_data, 1, numel(std_dev_fit)) - nr_fits).^2, 1));

end



% find the std_dev and the modal radius associated with the minimum rms difference
[~, first_min_idx] = min(rms_diff, [], 'all');
first_min_std_dev = Sigma(first_min_idx);
first_min_r_modal = R_modal(first_min_idx);


% ---------------------------------------------------------------
% ------- REFINE THE SOLUTION BY RUNNING ANOTHER SEARCH ---------
% ---------------------------------------------------------------

% Let's run a more refined search with the bounds surrounding the minimum
% value found above
% define a range of values around the minimum point found above
std_dev_fit = linspace(0.8*first_min_std_dev, 1.2*first_min_std_dev, 1000);            % values of mu to test

% compute the gamma distribution for different values of mu and r_modal
r_modal = linspace(0.5 * first_min_r_modal, 1.5 * first_min_r_modal, 1000);        % microns

% use this to run a loop with changing r_modal values
[Sigma, R_modal] = meshgrid(std_dev_fit, r_modal);


rms_diff = zeros(numel(r_modal), numel(std_dev_fit));

for nn = 1:numel(r_modal)
    
    N = 1./(sqrt(2*pi) * std_dev_fit);  % normalization constant

    nr_fits = N0 * N(nn) * 1./r .* exp(-log(r./r_modal(nn)).^2 ./(2*std_dev_fit.^2));       % lognormal droplet distribution

    % Find the least squares difference between the fit true data

    rms_diff(nn,:) = sqrt(sum(( repmat(nr_data, 1, numel(std_dev_fit)) - nr_fits).^2, 1));

end



% find the mu value associated with the minimum rms difference
[fit_parameters.min_rms, min_idx] = min(rms_diff, [], 'all');


% Grab the mu value and the nr_fit associated with the minimum rms 
fit_parameters.std_dev = Sigma(min_idx);                % my value (effective varaince?)
fit_parameters.r_modal = R_modal(min_idx);

% Now make the fit with the parameters found above
N = 1./(sqrt(2*pi) * fit_parameters.std_dev);  % normalization constant

nr_fit = N0 * N * 1./r .* exp(-log(r./fit_parameters.r_modal).^2 ./(2*fit_parameters.std_dev.^2));              % lognoraml droplet distribution


% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------






end