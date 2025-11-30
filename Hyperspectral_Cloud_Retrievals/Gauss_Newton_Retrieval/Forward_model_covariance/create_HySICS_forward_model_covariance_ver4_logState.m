%% Create forward model covariance matrix


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

function [GN_inputs] = create_HySICS_forward_model_covariance_ver4_logState(GN_inputs)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE FORWARD MODEL PARAMETERS ARE INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------






%----------------------------------------------------------------
% --- Define the mean and std of each forward model parameter ---
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





%----------------------------------------------------------
% ----------- Define the Covariance Matrix ----------------
%----------------------------------------------------------

% For now lets claim the desired variables are independent
GN_inputs.model.forward_model.covariance = diag((GN_inputs.model.forward_model.re.std).^2);


%----------------------------------------------------
% ------ Define the Variance of each Variable  ------
%----------------------------------------------------
GN_inputs.model.forward_model.variance = diag(GN_inputs.model.forward_model.covariance);















end










