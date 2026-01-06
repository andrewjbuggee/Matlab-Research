%% create measurement prior with EMIT measurements


% By Andrew J. Buggee
%%

function [GN_inputs] = create_EMIT_measurement_cov_ver4_log_no_FM_uncert_perPixel(GN_inputs, emit, pixel_num)


% define the covariance type
covariance_type = GN_inputs.measurement.covariance_type;




% --------------------------------------------------------
% Create the covariance matrix by taking the cross
% correlation between spectral channels with the modis data
% ---------------------------------------------------------

if strcmp(covariance_type,'computed') == true

    data = cat(3,emit.EV.m250.reflectance(:, pixel_num),emit.EV.m500.reflectance(:, pixel_num));
    data = data(:,:,GN_inputs.spectral_bins);
    for bb = 1:length(GN_inputs.spectral_bins)
        for ii = 1:length(modisInputs.pixels2use.res500m.row)
            data2run(ii,bb) = data(modisInputs.pixels2use.res500m.row(ii),modisInputs.pixels2use.res500m.col(ii),bb);
        end
    end

    GN_inputs.measurement.covariance = cov(data2run);



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


    % ---------------------------------------------------------------
    % Create the measurement covariance using the radiance uncertainty
    % ----------------------------------------------------------------


    % Let's assume gaussian measurement uncertainty. An instrument with 1%
    % measurement uncertainty means that 95% of identical measurements fall
    % within 1% of the 'true' value. For a normal distribution, 95% of all
    % measurements fall with [mean - 2*std, mean + 2*std].


    % if each uncertainty represents the standard deviation, the
    % variance is the square of each value.

    % the main diagonal of the measurement covariance matrix is the
    % variance of the measurement uncertainty. Since we assume Gaussian
    % statistics, the uncertainty is assumed to be the standard deviaiton.
    % The values should be in units of reflectance


    % use the definine measurement uncertainty. Convert to a percent
    GN_inputs.measurement.uncert_frac = emit.reflectance.uncertainty(:, pixel_num) ./...
        emit.reflectance.value(:, pixel_num);
    reflec_uncert = emit.reflectance.uncertainty(:, pixel_num);  % 1/sr
    reflec_mean = emit.reflectance.value(:, pixel_num);           % 1/sr


    % *** For computing the covariance of the logarithm of the measurements ***
    % Assume each spectral channel follows a gaussian distribution with
    % the mean value as the measurement and the standard deviation as the
    % uncertainty.
    % create a set of synthetic measurements for each channel
    % Take the log of these and compute the variance
    % The main diagonal is var(log(y))
    n_samples = 10000;
    measurement_samples_synthetic = zeros(n_samples, length(reflec_uncert));
    
    for nn = 1:length(reflec_uncert)

        measurement_samples_synthetic(:,nn) = (reflec_uncert(nn) .* randn(n_samples, 1)) +...
            reflec_mean(nn);

    end
   




    % Does it make sense to compute the root-sum-square uncertainty using
    % the measurements in log space? Values less than 1 approach 1 infinity
    % when you take the log. So, smaller absolute unceratinties leads to a
    % larger RSS uncert. Let's keep this in linear space
    GN_inputs.measurement.rss_uncert_linear = sqrt(sum( (reflec_uncert).^2));   % 1/sr - reflectance



    % Lets assume the percentage given is the standard deviation
    % According to King and Vaughn (2012): 'the values along the main
    % diagonal correspond to the square of the uncertainty estimate for
    % each wavelength channel'

    GN_inputs.measurement.variance_noLog = var(measurement_samples_synthetic, [], 1);     % 1/sr^2 - reflectance squared




    % Create a diagonal matrix where each entry is the variance of that
    % spectral channel for reflectance measurements
    GN_inputs.measurement.covariance_lin = diag(GN_inputs.measurement.variance_noLog);

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











