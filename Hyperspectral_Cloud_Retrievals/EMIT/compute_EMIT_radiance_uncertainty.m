%% Compute the EMIT radiance uncertainty using the noise model defined in the L1B Theoretical Basis document


% By Andrew John Buggee

%%

function radiance_uncertainty = compute_EMIT_radiance_uncertainty(emit)


%% Which computer is this?

if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------


    emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/';


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    emitDataPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/';


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

radiance_uncertainty = eta_1 .* sqrt(eta_2 .* emit.radiance.measurements) + eta_3;      % microW/cm^2/nm/sr

end