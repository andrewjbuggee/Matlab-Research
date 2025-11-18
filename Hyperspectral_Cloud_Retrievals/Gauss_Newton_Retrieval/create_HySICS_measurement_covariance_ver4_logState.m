%% create measurement prior


% By Andrew J. Buggee
%%

function [GN_inputs] = create_HySICS_measurement_covariance_ver4_logState(GN_inputs, simulated_measurements)


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

    if isfield(simulated_measurements, 'Refl_model_with_noise')==true 

        % use the definine measurement uncertainty. Convert to a percent
        measurement_uncertainty_percent = simulated_measurements.inputs.measurement.uncert * 100;  % percent

    else

        % If you're not modeling measurement uncertainty, assume a near zero value for all
        % HySICS spectral bands
        measurement_uncertainty_percent = 0.0001;      % percent!

    end

    % Define the forward model uncertainty
    forward_model_uncertainty_percent = 0.01;               % percent

    GN_inputs.measurement.uncertainty = linspace((measurement_uncertainty_percent + forward_model_uncertainty_percent)/100,...
        (measurement_uncertainty_percent + forward_model_uncertainty_percent)/100,...
        length(GN_inputs.bands2run))';        % fraction


    % *** For computing the covariance of the logarithm of the measurements ***
    % Assume each spectral channel is follows a gaussian distribution with
    % the mean value as the measurement and the standard deviation as the
    % total uncertainty.
    % create a set of synthetic measurements for each channel
    % Take the log of these and compute the variance
    % The main diagonal is var(log(y))
    n_samples = 10000;
    measurement_samples_synthetic = zeros(n_samples, length(simulated_measurements.Refl_model));
    for nn = 1:length(simulated_measurements.Refl_model)

        measurement_samples_synthetic(:,nn) = (GN_inputs.measurement.uncertainty(nn) .* randn(n_samples, 1)) +...
            simulated_measurements.Refl_model(nn);

    end




    % Compute the uncertainty of each spectral measurement in linear space
    % first
    meas_uncertainty_absolute = simulated_measurements.Refl_model .* GN_inputs.measurement.uncertainty;

    % Does it make sense to compute the root-sum-square uncertainty using
    % the measurements in log space? Values less than 1 approach 1 infinity
    % when you take the log. So, smaller absolute unceratinties leads to a
    % larger RSS uncert. Let's keep this in linear space
    GN_inputs.measurement.rss_uncert_linear = sqrt(sum( (meas_uncertainty_absolute ).^2));   % 1/sr - reflectance



    % Lets assume the percentage given is the standard deviation
    % According to King and Vaughn (2012): 'the values along the main
    % diagonal correspond to the square of the uncertainty estimate for
    % each wavelength channel'

    GN_inputs.measurement.variance_noLog = var(measurement_samples_synthetic, [], 1);     % 1/sr^2 - reflectance squared




    % Create a diagonal matrix where each entry is the variance of that
    % spectral channel for reflectance measurements
    GN_inputs.measurement.covariance_noLog = diag(GN_inputs.measurement.variance_noLog);

    % % To avoid values of -infinity, set the zeros to some small non-zero
    % % value
    % GN_inputs.measurement.covariance_noLog(GN_inputs.measurement.covariance_noLog==0) = 1e-10;
    % 
    % % Define the covariance of the log of the prior
    % GN_inputs.measurement.covariance = log(GN_inputs.measurement.covariance_noLog);

    GN_inputs.measurement.covariance = diag(var(log(measurement_samples_synthetic),[], 1));





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
GN_inputs.convergence_limit = GN_inputs.measurement.rss_uncert_linear;           % 1/sr - Root-sum-square of the reflectance uncertainty
%GN_inputs.convergence_limit = linspace(0.01, 0.01, length(pixels2use.res1km.linearIndex));  % generic convergence limit






end



