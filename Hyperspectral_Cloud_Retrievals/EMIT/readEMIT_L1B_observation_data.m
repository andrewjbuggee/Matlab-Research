%% ----- READ IN EMIT L1B OBSERVATION DATA -----

% this function will read in L1B EMIT data as .nc files
% will produce all necessary information into a cell array, and the data
% set into a structure


% ---- Description of the different fields within the netCDF File -----

%   (1) radiance - units of microW/cm^2/nm/sr



% By Andrew J. Buggee
%%

function [ds] = readEMIT_L1B_observation_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retrieve the netCDF info structure
info = ncinfo(fileName);



%% --- Read in the observation geometry data ---

% read radiance values
obs_data = ncread(fileName, 'obs');

% in the observation L1B file, the 'bands' are the following variables:
% (1) path length (meters) - dist between sensor and ground
% (2) To-sensor Zenith (degrees) - 0-90 degrees from zenith
% (3) To-sensor Azimuth (degrees) - 0-360 degrees clockwise from N
% (4) To-sun Zenith (degrees) - 0-90 degrees from zenith
% (5) To-sun Azimuth (degrees) - 0-360 degrees clockwise from N
% (6) Phase Angle (degrees) - degrees between to-sensor and to-sun vectors
% in principal plane
% (7) Slope (angle) - local surface slope as derived from DEM in degrees
% (8) Aspect (angle) - local surface aspect 0 to 360 degrees clockwise from
% N
% (9) cosine i (unitless) - apparent local illumination factor based on DEM
% slope and aspect and to sun vector, 0 to 1
% (10) UTC time (fractional hours) - fractional hours since UTC midnight
% (11) Earth-Sun distance (AU's) - Distance between the Earth and Sun

% Let's reshape the data so the format is cross-track by down-track by
% spectral-channel
ds.sensor.zenith = reshape(obs_data(2,:,:), size(obs_data,2), size(obs_data,3), []);
ds.sensor.azimuth = reshape(obs_data(3,:,:), size(obs_data,2), size(obs_data,3), []);
ds.solar.zenith = reshape(obs_data(4,:,:), size(obs_data,2), size(obs_data,3), []);
ds.solar.azimuth = reshape(obs_data(5,:,:), size(obs_data,2), size(obs_data,3), []);
ds.utc_time = reshape(obs_data(10,:,:), size(obs_data,2), size(obs_data,3), []);









end





