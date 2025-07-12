%% create measurement prior


% By Andrew J. Buggee
%%

function [GN_inputs] = create_HySICS_measurement_covariance(GN_inputs, simulated_measurements)


covariance_type = GN_inputs.measurement.covariance_type;



% ---**--- Important Quantity ---**---
% According To "VALIDATION OF MODIS-DERIVED TOP-OF-ATMOSPHERE SPECTRAL RADIANCES BY MEANS OF VICARIOUS CALIBRATION"
%GN_inputs.measurement.uncertainty = 0.02; % percentage of measurement uncertainty for reflectance




% --------------------------------------------------------
% Create the covariance matrix by taking the cross
% correlation between spectral channels with the modis data
% ---------------------------------------------------------

if strcmp(covariance_type,'computed') == true

    error([newline, 'Not capable of doing this calcualtion right now.', newline])


elseif strcmp(covariance_type,'independent') == true
    
    % -------------------------------------------------
    % *** Compute the Measurement Covariance Matrix ***
    % -------------------------------------------------
    % The measurement covariance matrix is the sum of the measurement
    % uncertainty and the forward model uncertainty (Poulsen et al. 2012,
    % pg 1893)

    % create the covaraince matrix of the model parameters
    % if the covariance matrix is diagonal, then we are assuming each
    % measurement (spectral channel) is independent of one another



    % if each uncertainty represents the standard deviation, the
    % variance is the square of each value.
    % the refelctance uncertanties are listed in percentages. So we
    % multiply these percentages with the modis reflectance values to
    % get the uncertainty in reflectance.

    if isfield(simulated_measurements.inputs, 'measurement')==true &&...
            isfield(simulated_measurements.inputs.measurement, 'uncert')==true

        % use the definine measurement uncertainty. Convert to a percent
        measurement_uncertainty_percent = simulated_measurements.inputs.measurement.uncert * 100;
    
    else

        % If you're not modeling measurement uncertainty, assume a near zero value for all
        % HySICS spectral bands
        measurement_uncertainty_percent = 0.01;      % percent!

    end

    % Define the forward model uncertainty
    forward_model_uncertainty = 0.01;               % percent

    GN_inputs.measurement.uncertainty = linspace((measurement_uncertainty_percent + forward_model_uncertainty)/100,...
        (measurement_uncertainty_percent + forward_model_uncertainty)/100,...
        length(GN_inputs.bands2run))';        % fraction

    GN_inputs.measurement.rss_uncert = sqrt(sum( (simulated_measurements.Refl_model.* GN_inputs.measurement.uncertainty).^2))';   % 1/sr - reflectance

    % Lets assume the percentage given is the standard deviation
    % According to King and Vaughn (2012): 'the values along the main
    % diagonal correspond to the square of the uncertainty estimate for
    % each wavelength channel'

    GN_inputs.measurement.variance = (simulated_measurements.Refl_model.* GN_inputs.measurement.uncertainty).^2;     % 1/sr^2 - reflectance squared




    % Create a diagonal matrix where each entry is the variance of that
    % spectral channel for reflectance measurements
    GN_inputs.measurement.covariance = diag(GN_inputs.measurement.variance);






end



% Define the convergence limit. Convergence is defined using the residual,
% which is the difference between the true and estimated measurements.
% We take the RMS of the residual using all spectral channels. This is how
% we define the convergence limit. If the residual is the difference
% between the true measurement and the estimated measurement, and the true
% measurement has an uncertainty of 10%, then our estimate measurement
% should be within this uncertainty. Using MODIS, we can compute the
% RMS uncertainty vector and set this as the convergence limit.

% the convergence limit should be the RSS of the measurement uncertainty,
% in units of reflectance!
GN_inputs.convergence_limit = sqrt(sum((simulated_measurements.Refl_model .* GN_inputs.measurement.uncertainty).^2));           % 1/sr - Root-sum-square of the reflectance uncertainty
%GN_inputs.convergence_limit = linspace(0.01, 0.01, length(pixels2use.res1km.linearIndex));  % generic convergence limit






end



