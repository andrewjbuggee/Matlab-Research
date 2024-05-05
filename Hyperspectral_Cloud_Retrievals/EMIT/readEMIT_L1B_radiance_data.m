%% ----- READ IN EMIT L1B DATA -----

% this function will read in L1B EMIT data as .nc files
% will produce all necessary information into a cell array, and the data
% set into a structure


% ---- Description of the different fields within the netCDF File -----

%   (1) radiance - units of microW/cm^2/nm/sr



% By Andrew J. Buggee
%%

function [ds] = readEMIT_L1B_radiance_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retrieve the netCDF info structure
info = ncinfo(fileName);



%% --- Read in the radiance data ---

% read radiance values
radiance0 = ncread(fileName, 'radiance');

% Let's reshape the data so the format is cross-track by down-track by
% spectral-channel
for bb = 1:size(radiance0,1)
    ds.measurements(:,:, bb) = reshape(radiance0(bb,:,:), size(radiance0,2), size(radiance0,3), []);        % - microW/nm/cm^2/sr
end

% read the center wavelength of each band
ds.wavelength = ncread(fileName, 'sensor_band_parameters/wavelengths');     % nm

% read the full-width at half-max - roughly speaking the band width of each
% spectral channel
ds.fwhm = ncread(fileName, 'sensor_band_parameters/fwhm');     % nm

% read the lat and long position of each pixel
ds.lat = ncread(fileName, 'location/lat');      % degrees north
ds.long = ncread(fileName, 'location/lon');    % degrees east






end





