% Fit a gamma particle size distribution to a droplet distribution data set

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

%   (2) fit_parameters - the modal radius, the effective varaiance, the RMS 
%   difference between the data and the best fit, and the variance
%   coefficient.


% By Andrew John Buggee

%%

function [nr_fit, fit_parameters] = fit_gamma_size_distribution_kokhanovsky(nr_data, r)

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



%% using the independent vector provided, compute the gamma distribution for many alpha and r0 values
% This definition follows the definition of A. Kokhonovsky in 'Cloud
% Optics' on page 2

% make sure r and nr_data are column vectors
r = reshape(r, [],1);
nr_data = reshape(nr_data, [], 1);

% Compute the total number of droplets
N0 = trapz(r, nr_data);




% ----------------------------------------------------------
% -------------- For varrying just mu_fit ------------------
% ----------------------------------------------------------

% % define the range of mu values to check
% mu_fit = 1:0.01:50;            % values of mu to test
% 
% % Compute the gamma distribution for varrying values of mu
% 
% % lets find the modal radius
% [~,idx] = max(nr_data);
%  r_modal = r(idx);
% 
% 
% N = mu_fit.^(mu_fit+1)./(gamma(mu_fit+1) .* r_modal.^(mu_fit+1));  % normalization constant
% 
% 
% nr_fits = N0 * N .* r.^mu_fit .* exp(-mu_fit .*r ./r_modal);                            % gamma droplet distribution
% 
% rms_diff = sqrt(sum(( repmat(nr_data, 1, numel(mu_fit)) - nr_fits).^2, 1));
% 
% % find the mu value associated with the minimum rms difference
% [fit_parameters.min_rms, min_idx] = min(rms_diff, [], 'all');
% 
% 
% % Grab the mu value and the nr_fit associated with the minimum rms 
% 
% fit_parameters.mu = mu_fit(min_idx);
% fit_parameters.var_coef = 1/sqrt(1+ fit_parameters.mu);
% nr_fit = nr_fits(:, min_idx);
% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------





% ----- For varrying r_modal and mu_fit --------

% define a short range a mu values, we will interpolate later
mu_fit = 1:100;            % values of mu to test

% compute the gamma distribution for different values of mu and r_modal
r_modal = 1:60;        % microns

% use this to run a loop with changing r_modal values
[Mu, R_modal] = meshgrid(mu_fit, r_modal);


%nr_fits = zeros(numel(r), numel(mu_fit), numel(r_modal));

rms_diff = zeros(numel(r_modal), numel(mu_fit));

for nn = 1:numel(r_modal)
    
    N = mu_fit.^(mu_fit+1)./(gamma(mu_fit+1) .* r_modal(nn).^(mu_fit+1));  % normalization constant

    nr_fits = N0 * N .* r.^mu_fit .* exp(-mu_fit .*r./r_modal(nn));                            % gamma droplet distribution

    % Find the least squares difference between the fit true data

    rms_diff(nn,:) = sqrt(sum(( repmat(nr_data, 1, numel(mu_fit)) - nr_fits).^2, 1));

end



% find the mu value associated with the minimum rms difference
[~, first_min_idx] = min(rms_diff, [], 'all');
first_min_mu = Mu(first_min_idx);
first_min_r_modal = R_modal(first_min_idx);


% ---------------------------------------------------------------
% ------- REFINE THE SOLUTION BY RUNNING ANOTHER SEARCH ---------
% ---------------------------------------------------------------

% Let's run a more refined search with the bounds surrounding the minimum
% value found above
% define a range of values around the minimum point found above
mu_fit = linspace(0.9*first_min_mu, 1.1*first_min_mu, 1000);            % values of mu to test

% compute the gamma distribution for different values of mu and r_modal
r_modal = linspace(0.5 * first_min_r_modal, 1.5 * first_min_r_modal, 1000);        % microns

% use this to run a loop with changing r_modal values
[Mu, R_modal] = meshgrid(mu_fit, r_modal);


rms_diff = zeros(numel(r_modal), numel(mu_fit));

for nn = 1:numel(r_modal)
    
    N = mu_fit.^(mu_fit+1)./(gamma(mu_fit+1) .* r_modal(nn).^(mu_fit+1));  % normalization constant

    nr_fits = N0 * N .* r.^mu_fit .* exp(-mu_fit .*r./r_modal(nn));                            % gamma droplet distribution

    % Find the least squares difference between the fit true data

    rms_diff(nn,:) = sqrt(sum(( repmat(nr_data, 1, numel(mu_fit)) - nr_fits).^2, 1));

end



% find the mu value associated with the minimum rms difference
[fit_parameters.min_rms, min_idx] = min(rms_diff, [], 'all');


% Grab the mu value and the nr_fit associated with the minimum rms 
fit_parameters.mu = Mu(min_idx);                % my value (effective varaince?)
fit_parameters.r_modal = R_modal(min_idx);
fit_parameters.var_coef = 1/sqrt(1+ fit_parameters.mu);

% Now make the fit with the parameters found above
N = fit_parameters.mu.^(fit_parameters.mu+1)./(gamma(fit_parameters.mu+1) .*...
    fit_parameters.r_modal.^(fit_parameters.mu+1));  % normalization constant

nr_fit = N0 * N .* r.^fit_parameters.mu .* exp(-fit_parameters.mu .*r./fit_parameters.r_modal);                            % gamma droplet distribution


% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------




end