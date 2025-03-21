%% create measurement prior with EMIT measurements


% By Andrew J. Buggee
%%

function [inputs] = create_EMIT_measurement_covariance(inputs, emit, pixels2use)


% define the covariance type
covariance_type = inputs.measurement.covariance_type;


% Define the number of spectral channels
num_bands_2run = length(inputs.bands2run);

% define the number of pixels to run
num_pixels = length(pixels2use);


% --------------------------------------------------------
% Create the covariance matrix by taking the cross
% correlation between spectral channels with the modis data
% ---------------------------------------------------------

if strcmp(covariance_type,'computed') == true

    data = cat(3,emit.EV.m250.reflectance,emit.EV.m500.reflectance);
    data = data(:,:,inputs.spectral_bins);
    for bb = 1:length(inputs.spectral_bins)
        for ii = 1:length(modisInputs.pixels2use.res500m.row)
            data2run(ii,bb) = data(modisInputs.pixels2use.res500m.row(ii),modisInputs.pixels2use.res500m.col(ii),bb);
        end
    end

    inputs.measurement.covariance = cov(data2run);



elseif strcmp(covariance_type,'independent') == true

    % create the measurement covariance matrix using EMIT uncertainty
    % estimates

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
    inputs.measurement.variance = (emit.reflectance.uncertainty(inputs.bands2run,:)).^2;      % variance in reflectance


    inputs.measurement.covariance = zeros(num_bands_2run, num_bands_2run, num_pixels);

    % Step through each pixel being used
    for pp = 1:num_pixels


        % Lets assume the percentage given is the standard deviation
        % According to King and Vaughn (2012): 'the values along the main
        % diagonal correspond to the square of the uncertainty estimate for
        % each wavelength channel'



        % Create a diagonal matrix where each entry is the variance of that
        % spectral channel for reflectance measurements
        inputs.measurement.covariance(:,:,pp) = diag(inputs.measurement.variance(:,pp));

    end





end







end



