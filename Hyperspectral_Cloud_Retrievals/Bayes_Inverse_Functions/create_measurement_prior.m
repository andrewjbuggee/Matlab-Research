% create measurement prior


% By Andrew J. Buggee
%%


function [bayes_inputs] = create_measurement_prior(bayes_inputs,data_struct,data_inputs)

bayes_inputs.data_type = inputname(2); % data set to be read. tells function what sort of measurement device we are using
measurement_prior = bayes_inputs.measurement.prior;

% ---------------------------------------
% --- Stuff for the Measurement Prior ---
% ---------------------------------------


% the measurement priors depend on what measurement device (MODIS or
% AVIRIS) we are using

% the measurement prior could be quite sophisticated. We could use
% knowledge of the spectrometer in question to estimate things like
% measurement uncertainty, and spectral resolution. We can use the average
% atmospheric composition to compute a spectral transmittance. Then, using
% planks law, we can estimate where our measurements might look like using
% this spectral transmittance. We can vary this prior by changing the
% transmittance for different cloud types and scence types. We also need to
% decide how correlated each pixel is from one another.

if strcmp(bayes_inputs.data_type,'aviris')==true
    
    bayes_inputs.spectral_bins = length(data_struct.wavelengths); % number of spectral channels, for this application
    
    
    if strcmp(measurement_prior,'custom')
        
        % lets create the variance and mean for each model parameter
        bayes_inputs.measurement.variance = [4,10]; % variance for the effective radius (microns squared) and optical thickness respectively
        bayes_inputs.measurement.mean = [10,15]; % expected values for the effective radius (microns) and the optical depth
        bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
        
        error('mvnpdf is not properly set up yet!')
        bayes_inputs.measurement.prior = mvnpdf(bayes_inputs.measurement.mean,bayes_inputs.measurement.covariance);
        
    elseif strcmp(measurement_prior,'standard_normal')
        % if the standard normal option is chosen, the code assumes the model
        % parameters take on the standard normal gaussian, where the the
        % expected value for each parameter is 0, and the variance is 1
        
        bayes_inputs.measurement.varaince = ones(1,bayes_inputs.spectral_bins);
        bayes_inputs.measurement.mean = zeros(1,bayes_inputs.spectral_bins);
        
        % create the covaraince matrix of the model parameters
        bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
        
        error('standard normal option is not properly set up yet!')
        
    elseif strcmp(measurement_prior,'gaussian')
        % if the 'gaussian' option is chosen, a multivariate gaussian pdf with
        % custom variance and mean will be implemented
        
        % we use the aviris measurement uncertainty according. Radiance is
        % measured in units of microwatts per square centimeter per
        % nanometer per steradian. AVIRIS radiometric calibration factors are
        % calculated by measuring the response of AVIRIS to an integrating
        % sphere (a known target illuminated by a known light source). This
        % calibration is accurate to within 7%, absolute, over time.
        
        % ---**--- Important Quantity ---**---
        bayes_inputs.measurement.uncertainty = 0.03; % percentage of measurement uncertainty according To Coddington et al. (2012) paragraph 24
        
        % we will calculate the plank function to estimate the measurement
        % expected values across the AVIRIS spectrum. The AVIRIS data is
        % reported in micro-watts/cm^2/nm/sr. So we need to convert the plank
        % function calculation, which is in watts/m^2/nm/sr
        
        T_sun = 5800; % Kelvin - temperature of the photosphere of the sun
        wavelengths = data_struct.wavelengths;
        radiance = planks_function(wavelengths,T_sun,'nanometers'); % watts/m^2/nm/sr
        
        % now lets compute the spectral radiance at the top of Earths
        % atmosphere
        r_sun = 696340000; % meters - radius of the sun
        dist_sun2earth = 1.5e11; % meters - average distance from sun to earth
        
        radiance = radiance * (r_sun/dist_sun2earth)^2; % watts/m^2/nm/sr - top of atmophsere spectral radiance
        
        % guess at the transmittance
        transmittance = 0.5; % constant over all wavelengths
        % convert radiance to aviris units. The radiance calculation below
        % represents our expected value for the measurement of radiance by our
        % spectrometer. But the spectrometer has measurement error, which we
        % will assume to be gaussian. There at each wavelength there is a 1D
        % pdf that represents the spread of possible measurements
        radiance = transmittance * (radiance * (1/1e-6) * (1/1e4)); % micro-watts/cm^2/nm/sr
        
        bayes_inputs.measurement.mean = radiance;
        bayes_inputs.measurement.varaince = bayes_inputs.measurement.uncertainty * radiance; % Radiometric variance is some percentage of the mean
        
        
        % create the covaraince matrix of the model parameters
        % if the covariance matrix is diagonal, then we are assuming each
        % measurement (spectral channel) is independent of one another
        bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
        
        % lets create our measurement pdf at each wavelength. We do this by
        % sampling a gaussian pdf with a  mean and variance as defined above
        
        
        
    end
    
    
elseif strcmp(bayes_inputs.data_type,'modis') == true
    
    % This defines the number of spectral channels used to create the
    % measurement covariance
    bayes_inputs.spectral_bins = size(data_struct,3);

    
    if strcmp(measurement_prior,'custom')
        
        % lets create the variance and mean for each model parameter
        bayes_inputs.measurement.variance = [4,10]; % variance for the effective radius (microns squared) and optical thickness respectively
        bayes_inputs.measurement.mean = [10,15]; % expected values for the effective radius (microns) and the optical depth
        bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
        
        error('mvnpdf is not properly set up yet!')
        bayes_inputs.measurement.prior = mvnpdf(bayes_inputs.measurement.mean,bayes_inputs.measurement.covariance);
        
    elseif strcmp(measurement_prior,'standard_normal')
        % if the standard normal option is chosen, the code assumes the model
        % parameters take on the standard normal gaussian, where the the
        % expected value for each parameter is 0, and the variance is 1
        
        bayes_inputs.measurement.varaince = ones(1,bayes_inputs.spectral_bins);
        bayes_inputs.measurement.mean = zeros(1,bayes_inputs.spectral_bins);
        
        % create the covaraince matrix of the model parameters
        bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
        
        error('standard normal option is not properly set up yet!')
        
    elseif strcmp(measurement_prior,'gaussian')
        % if the 'gaussian' option is chosen, a multivariate gaussian pdf with
        % custom variance and mean will be implemented
        
        % we use the modis measurement uncertainty. Radiance is
        % measured in units of watts per square meter per
        % micron per steradian.
        
        % ---**--- Important Quantity ---**---
        bayes_inputs.measurement.uncertainty = 0.02; % percentage of measurement uncertainty for reflectance according To "VALIDATION OF MODIS-DERIVED TOP-OF-ATMOSPHERE SPECTRAL RADIANCES BY MEANS OF VICARIOUS CALIBRATION"
        
        % we will calculate the plank function to estimate the measurement
        % expected values across the MODIS spectrum. The MODIS data is
        % reported in watts/m^2/micron/sr.
        
        wavelengths = [data_struct.EV.m250.bands.center;data_struct.EV.m500.bands.center]; % ---**--- Change this after the modis structure is fixes
        
        
        
        T_sun = 5800; % Kelvin - temperature of the photosphere of the sun
        radiance = plancks_function(wavelengths,T_sun,'microns'); % watts/m^2/micron/sr
        
        % now lets compute the spectral radiance at the top of Earths
        % atmosphere
        r_sun = 696340000; % meters - radius of the sun
        dist_sun2earth = 1.5e11; % meters - average distance from sun to earth
        
        radiance = radiance * (r_sun/dist_sun2earth)^2; % watts/m^2/micron/sr - top of atmophsere spectral radiance
        
        % guess at the transmittance
        transmittance = 0.5; % constant over all wavelengths
        %  The radiance calculation below
        % represents our expected value for the measurement of radiance by our
        % spectrometer. But the spectrometer has measurement error, which we
        % will assume to be gaussian. There at each wavelength there is a 1D
        % pdf that represents the spread of possible measurements
        radiance = transmittance * radiance;
        
        bayes_inputs.measurement.mean = radiance;
        bayes_inputs.measurement.variance = ones(length(wavelengths),1)*bayes_inputs.measurement.uncertainty; % Radiometric variance is some percentage of the mean
        
        
        
        
        
    end
    
    
end


end

