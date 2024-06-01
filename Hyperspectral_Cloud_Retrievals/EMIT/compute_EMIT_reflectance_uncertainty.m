%% Compute the EMIT reflectance uncertainty using the Noise-Equivalent Change in Radiance


% By Andrew John Buggee

%%

function reflectance_uncertainty = compute_EMIT_reflectance_uncertainty(emit, inputs)



%%  Compute the reflectance uncertainty at each spectral channel

% define the solar spectral irradiance
solar_irradiance = inputs.source.flux;          % W/nm/m^2

% convert solar irradiance to the same units as the EMIT measurements
solar_irradiance = solar_irradiance .* 100;      % microW/nm/cm^2

% define the solar zenith angle
sza = emit.obs.solar.zenith;        % deg

% define the noise-equivelant change in radiance
NEdR = emit.radiance.uncertainty;               % % microW/nm/cm^2/sr


% determine the number of pixels
num_pixels = size(emit.radiance.measurements, 2);


reflectance_uncertainty = zeros(size(emit.radiance.measurements));

% We use the same equation to compute spectral reflectance
for nn = 1:num_pixels

    reflectance_uncertainty(:,nn) = ( pi * NEdR(:,nn) )./...
                                ( cosd(sza) * solar_irradiance);                   % 1/sr

end



end