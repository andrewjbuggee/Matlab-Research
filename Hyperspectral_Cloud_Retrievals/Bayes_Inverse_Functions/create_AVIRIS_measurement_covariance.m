%% create measurement prior


% By Andrew J. Buggee
%%

function [bayes_inputs] = create_AVIRIS_measurement_covariance(bayes_inputs,data_struct,data_inputs)


covariance_type = bayes_inputs.measurement.covariance_type;

    
    bayes_inputs.spectral_bins = length(data_struct.wavelengths); % number of spectral channels, for this application
    
    
    
    % ---**--- Important Quantity ---**---
    bayes_inputs.measurement.uncertainty = 0.03; % percentage of measurement uncertainty according To Coddington et al. (2012) paragraph 24
    
    % measurement (spectral channel) is independent of one another
    bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
    
    % lets create our measurement pdf at each wavelength. We do this by
    % sampling a gaussian pdf with a  mean and variance as defined above
    
    
    
   




end



