%% ----- READ IN EMIT L1B OBSERVATION DATA -----

% this function will read in L1B EMIT data as .nc files
% The observation files contain the solar and viewing geometry among other
% things. The list of the information included is below:
% (1)    "Path length (sensor-to-ground in meters)"
% (2)    "To-sensor azimuth (0 to 360 degrees CW from N)"
% (3)    "To-sensor zenith (0 to 90 degrees from zenith)"
% (4)    "To-sun azimuth (0 to 360 degrees CW from N)"
% (5)    "To-sun zenith (0 to 90 degrees from zenith)"
% (6)    "Solar phase (degrees between to-sensor and to-sun vectors in principal plane)"
% (7)    "Slope (local surface slope as derived from DEM in degrees)"
% (8)    "Aspect (local surface aspect 0 to 360 degrees clockwise from N)"
% (9)    "Cosine(i) (apparent local illumination factor based on DEM slope and aspect and to sun vector)"
% (10)    "UTC Time (decimal hours for mid-line pixels)"
% (11)    "Earth-sun distance (AU)"


% ---- Description of the different fields within the netCDF File -----



% By Andrew J. Buggee
%%

function [ds] = readEMIT_L1B_observation_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retrieve the netCDF info structure
info = ncinfo(fileName);



%% --- Read in the observation geometry data ---

% read radiance values
obs_data = ncread(fileName, 'obs');

% in the observation L1B file, the 'bands' are labeled accoring to the
% netCDF file from 17 Jan 2024 /sensor_band_parameters/observation_bands:

% (1)    "Path length (sensor-to-ground in meters)"
% (2)    "To-sensor azimuth (0 to 360 degrees CW from N)"
% (3)    "To-sensor zenith (0 to 90 degrees from zenith)"
% (4)    "To-sun azimuth (0 to 360 degrees CW from N)"
% (5)    "To-sun zenith (0 to 90 degrees from zenith)"
% (6)    "Solar phase (degrees between to-sensor and to-sun vectors in principal plane)"
% (7)    "Slope (local surface slope as derived from DEM in degrees)"
% (8)    "Aspect (local surface aspect 0 to 360 degrees clockwise from N)"
% (9)    "Cosine(i) (apparent local illumination factor based on DEM slope and aspect and to sun vector)"
% (10)    "UTC Time (decimal hours for mid-line pixels)"
% (11)    "Earth-sun distance (AU)"

% Let's reshape the data so the format is cross-track by down-track by
% spectral-channel
ds.sensor.azimuth = reshape(obs_data(2,:,:), size(obs_data,2), size(obs_data,3), []);
ds.sensor.zenith = reshape(obs_data(3,:,:), size(obs_data,2), size(obs_data,3), []);
ds.solar.azimuth = reshape(obs_data(4,:,:), size(obs_data,2), size(obs_data,3), []);
ds.solar.zenith = reshape(obs_data(5,:,:), size(obs_data,2), size(obs_data,3), []);
ds.utc_time = reshape(obs_data(10,:,:), size(obs_data,2), size(obs_data,3), []);









end





