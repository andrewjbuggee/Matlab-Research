%% Create forward model covariance matrix

% *** CURRENT FORWARD MODEL UNCERTAINTIES CONSIDERED ***
% (1) Adiabatic droplet profile assumption
% (2) Cloud top height assumption
% (3) Droplet distribution effective variance assumption


% There are many components that could be taken into account when
% considering the forward model uncertainty. Here is an incomplete list of
% the more important set of forward model parameters that may cause
% uncertainty
%   (1) assumed effective variance of the droplet size distribution
%   (2) assumed shape of the droplet size distribution?
%   (3) assumed vertical profile of water vapor concentration
%   (4) assumed cloud top height
%   (5) deviations between the true droplet profile and the theoretical
%   profile retrieved
%   (6) vertical profile of carbon dioxide?
%   (7) vertical profile of aerosols?
%   (8) vertical profile of temperature?






% By Andrew John Buggee

%%

function [GN_inputs] = create_HySICS_forMod_cov_ver4_log_reP_CTH_Ve(GN_inputs)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE FORWARD MODEL PARAMETERS ARE INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------






%----------------------------------------------------------------
% --- Define the mean and std of each forward model parameter ---
%----------------------------------------------------------------




%----------------------------------------------------------------
% ------------ Uncertainty due to adiabatic profile -------------
%----------------------------------------------------------------

% Define the uncertainty due to assuming a theoretical adiabatic
% droplet profile
% create a droplet profile using the a priori values and the number of
% layers assumed in our forward model
% re is a column vector that starts at cloud base and grows towards cloud top
GN_inputs.model.forward_model.re.mean{1} = create_droplet_profile2(exp([GN_inputs.model.apriori(1),...
    GN_inputs.model.apriori(2)]), sort(GN_inputs.RT.z),...
    GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector

% These values are the mean, define the standard deviation. This
% reflects the uncertainty between the theoretical value and the true
% profile
% Values at the top and bottom have more uncerainty due to low LWC and
% entrainment at cloud top
GN_inputs.model.forward_model.re.std = zeros(size(GN_inputs.model.forward_model.re.mean{1}));
GN_inputs.model.forward_model.re.std(1:2) = 2; % microns
GN_inputs.model.forward_model.re.std(end-1:end) = 2; % microns
GN_inputs.model.forward_model.re.std(3:end-2) = 1; % microns


% *** For computing the covariance of the logarithm of the parameters ***
% Assume each assumed droplet effective radius along the profile follows a
% lognormal distribution, where the mean is the assumed value and the std
% is the value defined above

% Take the log of these and compute the variance
% The main diagonal is var(log(y))
n_samples = 10000;
samples_effective_radius = zeros(length(GN_inputs.model.forward_model.re.mean{1}), n_samples);

for nn = 1:length(GN_inputs.model.forward_model.re.mean{1})

    % convert mean and std from a normal distribution to a lognormal
    % distribution
    mu_log = log(GN_inputs.model.forward_model.re.mean{1}(nn)) -...
        (1/2 * log( (GN_inputs.model.forward_model.re.std(nn)^2 / GN_inputs.model.forward_model.re.mean{1}(nn)^2) +1));

    std_log = sqrt( log( (GN_inputs.model.forward_model.re.std(nn)^2 / GN_inputs.model.forward_model.re.mean{1}(nn)^2) +1));


    logNorm_dist = makedist("Lognormal", mu_log, std_log);

    samples_effective_radius(nn, :) = logNorm_dist.random(1, n_samples);


end








%----------------------------------------------------------------
% --------------- Uncertainty of Cloud Top Height ---------------
%----------------------------------------------------------------

% ** VOCALS-REx in-situ measurements result in a mean cloud top height
% of 1203 meters and a mean cloud depth of about 230 meters
% ** testing the retrieval when I lack knowledge of cloud top precisely **

% load the set of VOCALS-REx in-situ observations
if strcmp(GN_inputs.which_computer, 'anbu8374')==true

    cloud_top_obs = load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/VR_cloud_top_height_obs_19-Jan-2026.mat']);

elseif strcmp(GN_inputs.which_computer, 'andrewbuggee')==true

    cloud_top_obs = load(['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/VR_cloud_top_height_obs_19-Jan-2026.mat']);

elseif strcmp(GN_inputs.which_computer, 'curc')==true

    cloud_top_obs = load(['/projects/anbu8374/Matlab-Research/Presentations_and_Papers/',...
        'paper_2/VR_cloud_top_height_obs_19-Jan-2026.mat']);

end


% Define the uncertainty due incorrect knowledge of cloud top height
% create a droplet profile using the a priori values and the number of
% layers assumed in our forward model
% re is a column vector that starts at cloud base and grows towards cloud top
% GN_inputs.model.forward_model.cloudTopHeight.mean = GN_inputs.RT.z_topBottom(1);    % km - cloud top height
GN_inputs.model.forward_model.cloudTopHeight.mean = mean(cloud_top_obs.cloudTopHeight)/1e3;                  % km


% Define the standard deviation. This
% reflects the uncertainty between the assumed value and the true value
% Set the standard deviation to 218 meters, which is the value determined from the
% VOCALS-REx in-situ data
% GN_inputs.model.forward_model.cloudTopHeight.std = 0.218;  % km
GN_inputs.model.forward_model.cloudTopHeight.std = std( cloud_top_obs.cloudTopHeight)/1e3;  % km

% check to make sure the cloud top height is larger than the std. If not,
% reduce the cloud top height std to a value less than the cloud top height
if GN_inputs.model.forward_model.cloudTopHeight.mean < GN_inputs.model.forward_model.cloudTopHeight.std 

    GN_inputs.model.forward_model.cloudTopHeight.std = 0.5 * GN_inputs.model.forward_model.cloudTopHeight.mean;

end


% *** For computing the covariance of the logarithm of the parameters ***
% Assume the cloud top height follows a lognormal distribution,
% where the mean is the assumed value and the std
% is the value defined above

% Take the log of these and compute the variance
% The main diagonal is var(log(y))
n_samples = 10000;

% convert mean and std from a normal distribution to a lognormal
% distribution
mu_log = log(GN_inputs.model.forward_model.cloudTopHeight.mean) -...
    (1/2 * log( (GN_inputs.model.forward_model.cloudTopHeight.std^2 / GN_inputs.model.forward_model.cloudTopHeight.mean^2) +1));

std_log = sqrt( log( (GN_inputs.model.forward_model.cloudTopHeight.std^2 / GN_inputs.model.forward_model.cloudTopHeight.mean^2) +1));


logNorm_dist = makedist("Lognormal", mu_log, std_log);

samples_cloudTopHeight = logNorm_dist.random(1, n_samples);








%----------------------------------------------------------------
% ------------ Uncertainty due to effective variance ------------
%----------------------------------------------------------------

% Define the uncertainty due to assuming a particular value of the
% effective variance of the droplet distribution


% ** VOCALS-REx in-situ measurements **
% load the set of VOCALS-REx in-situ observations
if strcmp(GN_inputs.which_computer, 'anbu8374')==true

    alpha_effVar = load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/VR_effective_variance_at_normalized_altitudes_20-levels_19-Jan-2026.mat']);

elseif strcmp(GN_inputs.which_computer, 'andrewbuggee')==true

    alpha_effVar = load(['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/VR_effective_variance_at_normalized_altitudes_20-levels_19-Jan-2026.mat']);

elseif strcmp(GN_inputs.which_computer, 'curc')==true

    alpha_effVar = load(['/projects/anbu8374/Matlab-Research/Presentations_and_Papers/',...
        'paper_2/VR_effective_variance_at_normalized_altitudes_20-levels_19-Jan-2026.mat']);

end




% First, define the vector of effective variances, which is defined as the
% alpha parameter for libRadtran
GN_inputs.model.forward_model.alpha.mean{1} = GN_inputs.RT.distribution_var;   % alpha vector - the alpha parameter at each cloud layer for a gamma droplet distribution
GN_inputs.model.forward_model.alpha.std = GN_inputs.RT.distribution_var_std;   % the standard deviation of alpha values fournd at each level from the VOCALS-REx in-situ data



% *** For computing the covariance of the logarithm of the parameters ***
% Assume each assumed droplet effective radius along the profile follows a
% lognormal distribution, where the mean is the assumed value and the std
% is the value defined above

% Take the log of these and compute the variance
% The main diagonal is var(log(y))
n_samples = 10000;
samples_effective_variance = zeros(length(GN_inputs.model.forward_model.alpha.mean{1}), n_samples);

for nn = 1:length(GN_inputs.model.forward_model.re.mean{1})

    % convert mean and std from a normal distribution to a lognormal
    % distribution
    mu_log = log(GN_inputs.model.forward_model.alpha.mean{1}(nn)) -...
        (1/2 * log( (GN_inputs.model.forward_model.alpha.std(nn)^2 /...
        GN_inputs.model.forward_model.alpha.mean{1}(nn)^2) +1));

    std_log = sqrt( log( (GN_inputs.model.forward_model.alpha.std(nn)^2 /...
        GN_inputs.model.forward_model.alpha.mean{1}(nn)^2) +1));


    logNorm_dist = makedist("Lognormal", mu_log, std_log);

    samples_effective_variance(nn, :) = logNorm_dist.random(1, n_samples);


end











%----------------------------------------------------------
% ----------- Define the Covariance Matrix ----------------
%----------------------------------------------------------

% For now lets claim the desired variables are independent
% GN_inputs.model.forward_model.covariance = diag((GN_inputs.model.forward_model.re.std).^2);


% ** ensure the magnitude of cloud top height is similar to effective radius **

% you need to take the log of the data!
% ** Combine both data sets! **
GN_inputs.model.forward_model.covariance = diag(var( log( [samples_effective_radius; samples_cloudTopHeight;...
    samples_effective_variance] ), [], 2));     %  log space covaraince

GN_inputs.model.forward_model.covariance_lin = diag( var( [samples_effective_radius; samples_cloudTopHeight;...
    samples_effective_variance], [], 2));     %  log space covaraince

%----------------------------------------------------
% ------ Define the Variance of each Variable  ------
%----------------------------------------------------
GN_inputs.model.forward_model.variance = diag(GN_inputs.model.forward_model.covariance);

GN_inputs.model.forward_model.variance_lin = diag(GN_inputs.model.forward_model.covariance_lin);















end










