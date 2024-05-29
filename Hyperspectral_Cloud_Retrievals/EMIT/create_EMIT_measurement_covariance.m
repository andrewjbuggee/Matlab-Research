%% create measurement prior with EMIT measurements


% By Andrew J. Buggee
%%

function [inputs] = create_EMIT_measurement_covariance(inputs,modis,modisInputs,pixels2use)


covariance_type = inputs.measurement.covariance_type;



% ---**--- Important Quantity ---**---
% According To "VALIDATION OF MODIS-DERIVED TOP-OF-ATMOSPHERE SPECTRAL RADIANCES BY MEANS OF VICARIOUS CALIBRATION"
%GN_inputs.measurement.uncertainty = 0.02; % percentage of measurement uncertainty for reflectance

% Define the number of spectral channels
num_bands_2run = length(inputs.bands2run);


% --------------------------------------------------------
% Create the covariance matrix by taking the cross
% correlation between spectral channels with the modis data
% ---------------------------------------------------------

if strcmp(covariance_type,'computed') == true

    data = cat(3,modis.EV.m250.reflectance,modis.EV.m500.reflectance);
    data = data(:,:,inputs.spectral_bins);
    for bb = 1:length(inputs.spectral_bins)
        for ii = 1:length(modisInputs.pixels2use.res500m.row)
            data2run(ii,bb) = data(modisInputs.pixels2use.res500m.row(ii),modisInputs.pixels2use.res500m.col(ii),bb);
        end
    end

    inputs.measurement.covariance = cov(data2run);



elseif strcmp(covariance_type,'independent') == true

    % create the covaraince matrix of the model parameters
    % if the covariance matrix is diagonal, then we are assuming each
    % measurement (spectral channel) is independent of one another
    

    % ---------------------------------------------------------------
    % Create the measurement covariance using the radiance uncertainty
    % product
    % ----------------------------------------------------------------

    inputs.measurement.variance = zeros(num_bands_2run,length(pixels2use.res1km.linearIndex));
    inputs.measurement.covariance = zeros(num_bands_2run,num_bands_2run,length(pixels2use.res1km.linearIndex));

    % Step through each pixel being used
    for pp = 1:length(pixels2use.res1km.linearIndex)

        % Grab the row and column
        r = pixels2use.res1km.row(pp);
        c = pixels2use.res1km.col(pp);

        % Step through each band
        for bb = 1:num_bands_2run

            band_num = inputs.bands2use(bb);
            % if each uncertainty represents the standard deviation, the
            % variance is the square of each value.
            % the refelctance uncertanties are listed in percentages. So we
            % multiply these percentages with the modis reflectance values to
            % get the uncertainty in reflectance.


            % Lets start by converting the percentage to a decimal
            inputs.measurement.uncertainty(bb,pp) = 0.01* double(modis.EV1km.reflectanceUncert(r,c,band_num));        % uncertainty as a decimal


            % Lets assume the percentage given is the standard deviation
            % According to King and Vaughn (2012): 'the values along the main
            % diagonal correspond to the square of the uncertainty estimate for
            % each wavelength channel'

            inputs.measurement.variance(bb,pp) = (modis.EV1km.reflectance(r,c,band_num).* inputs.measurement.uncertainty(bb,pp)).^2;


        end

        % Create a diagonal matrix where each entry is the variance of that
        % spectral channel for reflectance measurements
        inputs.measurement.covariance(:,:,pp) = diag(inputs.measurement.variance(:,pp));

    end





end



% Define the convergence limit. Convergence is defined using the residual,
% which is the difference between the true and estimated measurements.
% We take the RMS of the residual using all spectral channels. This is how
% we define the convergence limit. If the residual is the difference
% between the true measurement and the estimated measurement, and the true
% measurement has an uncertainty of 10%, then our estimate measurement
% should be within this uncertainty. Using MODIS, we can compute the
% RMS uncertainty vector and set this as the convergence limit.

%GN_inputs.convergence_limit = sqrt(sum(GN_inputs.measurement.uncertainty.^2, 1));        % uncertainty as a decimal
inputs.convergence_limit = linspace(0.01, 0.01, length(pixels2use.res1km.linearIndex));  % generic convergence limit






end



