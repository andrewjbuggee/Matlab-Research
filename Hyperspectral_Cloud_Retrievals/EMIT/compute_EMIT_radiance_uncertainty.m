%% Compute the EMIT radiance uncertainty using the noise model defined in the L1B Theoretical Basis document


% By Andrew John Buggee

%%

function [radiance_uncertainty, radiance_uncertainty_percent_perChannel] = compute_EMIT_radiance_uncertainty(emit)


%% Which computer is this?

if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------


    emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/';


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    emitDataPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/';

elseif strcmp(whatComputer, 'curc')==true

    % ----- Folders on the CU super computer -----

    emitDataPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/';


end

%% Read the EMIT noise coefficients

filename = 'emit_noise.txt';
delimiter = ' ';
headerLines = 1;

noise_model = importdata([emitDataPath, filename], delimiter, headerLines);


%%  Compute the radiance uncertainty at each channel

% define the number of pixels for which we are computing the noise
num_pixels = size(emit.radiance.measurements, 2);

% From equation 13 of the EMIT L1B Algorithm: Calibrated Radiance at
% Sensor Theoretical Basis document

eta_1 = repmat(noise_model.data(:,2), 1, num_pixels);
eta_2 = repmat(noise_model.data(:,3), 1, num_pixels);
eta_3 = repmat(noise_model.data(:,4), 1, num_pixels);

% -------------------------------------------------------
% --- Compute the noise-equivelant change in radiance ---
% -------------------------------------------------------

% This is the minimum detectable change in radiance. Any change in the
% signal that is smaller than this value will go undetected. Essentially it
% is the system noise. It has units of radiance

% occasionally the emit radiance measurements are less than 0. How can this
% be? There are no photons in some wavelengths, but the measurement
% dependent uncertainty, like shot noise, causes the recorded signal to be
% negative?

% take the real part of the equation below to ignore the cases where some
% radiance measurements are less than 0
%radiance_uncertainty = real( (eta_1 .* sqrt(eta_2 .* emit.radiance.measurements)) + eta_3 );      % microW/cm^2/nm/sr
%radiance_uncertainty = real( (eta_1 .* sqrt(eta_2 + emit.radiance.measurements)) + eta_3 );      % microW/cm^2/nm/sr


%% ---- The above equations describe a noise model for the instrument ----

% Peter thinks I should set the EMIT radiance at uncertainty at 5% 
% radiance_uncertainty = 0.05 .* emit.radiance.measurements;              % microW/cm^2/nm/sr

% Let's set it to 3%
radiance_uncertainty_percent_perChannel = 0.03;
radiance_uncertainty = 0.03 .* emit.radiance.measurements;              % microW/cm^2/nm/sr



end