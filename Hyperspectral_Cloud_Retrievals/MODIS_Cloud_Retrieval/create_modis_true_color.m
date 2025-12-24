%% Create True Color RGB Image from MODIS Data

% This function creates a true color (RGB) composite image from MODIS data
% using bands 1 (red), 4 (green), and 3 (blue)
%
% INPUTS:
%   modis - MODIS data structure containing reflectance or radiance data
%           Must have modis.EV1km.reflectance or modis.EV1km.radiance
%   options - structure with optional parameters:
%       .use_radiance - true to use radiance, false for reflectance (default: false)
%       .gamma - gamma correction factor (default: 2.2)
%       .brightness - brightness adjustment factor (default: 1.0)
%       .contrast - contrast stretch percentiles [low, high] (default: [2, 98])
%
% OUTPUTS:
%   rgb_image - MxNx3 array containing the RGB image (values 0-1)
%   lat_grid - MxN array of latitude coordinates
%   lon_grid - MxN array of longitude coordinates
%
% MODIS Band Information:
%   Band 1 (620-670 nm) - Red
%   Band 4 (545-565 nm) - Green  
%   Band 3 (459-479 nm) - Blue
%
% Example usage:
%   [rgb, lat, lon] = create_modis_true_color(modis);
%   [rgb, lat, lon] = create_modis_true_color(modis, struct('gamma', 1.8, 'brightness', 1.2));
%
% By Andrew John Buggee

%%

function [rgb_image, lat_grid, lon_grid] = create_modis_true_color(modis, options)

% Set default options
if nargin < 2
    options = struct();
end

if ~isfield(options, 'use_radiance')
    options.use_radiance = false;
end

if ~isfield(options, 'gamma')
    options.gamma = 2.2;
end

if ~isfield(options, 'brightness')
    options.brightness = 1.0;
end

if ~isfield(options, 'contrast')
    options.contrast = [2, 98];  % percentile stretch
end


%% Extract RGB bands from MODIS data

% Check if data has reflectance or radiance
if isfield(modis, 'EV1km')==false
    % This is MYD06/MOD06 cloud product - no radiance/reflectance bands
    error('MODIS cloud product (MYD06/MOD06) does not contain radiance/reflectance data. Use MYD021KM/MOD021KM or MYD02HKM/MOD02HKM instead.');
end

% Try to find reflectance or radiance data
% Different MODIS products have different structures
if options.use_radiance
    % Look for radiance data
    if isfield(modis.EV1km, 'radiance')==true
        data_source = modis.EV_1KM_Emissive.radiance;

    else
        error('Could not find radiance data in MODIS structure. Check field names.');
    end
else
    % Look for reflectance data (preferred)
    if isfield(modis.EV1km, 'reflectance')
        data_source = modis.EV1km.reflectance;

    else
        warning('Could not find reflectance data. Trying radiance instead...');
        options.use_radiance = true;
        if isfield(modis, 'EV_1KM_RefSB') && isfield(modis.EV_1KM_RefSB, 'radiance')
            data_source = modis.EV_1KM_RefSB.radiance;
        elseif isfield(modis, 'radiance')
            data_source = modis.radiance;
        else
            error('Could not find reflectance or radiance data in MODIS structure.');
        end
    end
end



% Extract the three bands
% MODIS 1km bands: Band 1 (Red), Band 3 (Blue), Band 4 (Green)
% In the data cube, these are typically indices 1, 3, 4 but may vary

% Get data dimensions
data_size = size(data_source);

% Assume band dimension is the third dimension
if length(data_size) == 3
    % Data format: [rows x cols x bands]
    
    % For MODIS L1B products, bands are ordered as:
    % Band 1 (index 1) - Red (620-670 nm)
    % Band 2 (index 2) - NIR (841-876 nm) 
    % Band 3 (index 3) - Blue (459-479 nm)
    % Band 4 (index 4) - Green (545-565 nm)
    
    if data_size(3) >= 4
        red_band = double(data_source(:, :, 1));
        green_band = double(data_source(:, :, 4));
        blue_band = double(data_source(:, :, 3));
    else
        error('MODIS data does not have enough bands for RGB composite. Need at least 4 bands.');
    end
else
    error('Unexpected MODIS data structure. Expected 3D array [rows x cols x bands].');
end


%% Process each band

% Remove fill values and invalid data
fill_value = -9999;
red_band(red_band <= fill_value) = NaN;
green_band(green_band <= fill_value) = NaN;
blue_band(blue_band <= fill_value) = NaN;

% For radiance, convert to reflectance approximation if needed
if options.use_radiance
    % Simple scaling - proper conversion would need solar irradiance
    % This is a rough approximation for visualization
    red_band = red_band / max(red_band(:), [], 'omitnan') * 0.6;
    green_band = green_band / max(green_band(:), [], 'omitnan') * 0.6;
    blue_band = blue_band / max(blue_band(:), [], 'omitnan') * 0.6;
else
    % Reflectance should already be in reasonable units (0-1 or 0-100%)
    % Check range and scale if needed
    if max(red_band(:), [], 'omitnan') > 2
        % Likely scaled by 10000 (common in MODIS products)
        red_band = red_band / 10000;
        green_band = green_band / 10000;
        blue_band = blue_band / 10000;
    end
end

% Clip to valid range
red_band = max(0, min(1, red_band));
green_band = max(0, min(1, green_band));
blue_band = max(0, min(1, blue_band));


%% Apply contrast stretch

% Remove NaNs for percentile calculation
red_valid = red_band(~isnan(red_band));
green_valid = green_band(~isnan(green_band));
blue_valid = blue_band(~isnan(blue_band));

if ~isempty(red_valid)
    red_min = prctile(red_valid, options.contrast(1));
    red_max = prctile(red_valid, options.contrast(2));
    red_band = (red_band - red_min) / (red_max - red_min);
end

if ~isempty(green_valid)
    green_min = prctile(green_valid, options.contrast(1));
    green_max = prctile(green_valid, options.contrast(2));
    green_band = (green_band - green_min) / (green_max - green_min);
end

if ~isempty(blue_valid)
    blue_min = prctile(blue_valid, options.contrast(1));
    blue_max = prctile(blue_valid, options.contrast(2));
    blue_band = (blue_band - blue_min) / (blue_max - blue_min);
end

% Clip again after stretch
red_band = max(0, min(1, red_band));
green_band = max(0, min(1, green_band));
blue_band = max(0, min(1, blue_band));


%% Apply brightness adjustment
red_band = red_band * options.brightness;
green_band = green_band * options.brightness;
blue_band = blue_band * options.brightness;

% Clip after brightness
red_band = max(0, min(1, red_band));
green_band = max(0, min(1, green_band));
blue_band = max(0, min(1, blue_band));


%% Apply gamma correction
red_band = red_band .^ (1/options.gamma);
green_band = green_band .^ (1/options.gamma);
blue_band = blue_band .^ (1/options.gamma);


%% Create RGB composite
rgb_image = cat(3, red_band, green_band, blue_band);


%% Get geographic coordinates
lat_grid = double(modis.geo.lat);
lon_grid = double(modis.geo.long);

% Handle NaN in coordinates
lat_grid(lat_grid <= -9999) = NaN;
lon_grid(lon_grid <= -9999) = NaN;


%% Display summary
fprintf('True color RGB image created successfully.\n');
fprintf('Image size: %d x %d pixels\n', size(rgb_image, 1), size(rgb_image, 2));
fprintf('Latitude range: [%.2f, %.2f]\n', min(lat_grid(:), [], 'omitnan'), max(lat_grid(:), [], 'omitnan'));
fprintf('Longitude range: [%.2f, %.2f]\n', min(lon_grid(:), [], 'omitnan'), max(lon_grid(:), [], 'omitnan'));

% Check coordinate orientation
fprintf('Top-left corner: lat=%.2f, lon=%.2f\n', lat_grid(1,1), lon_grid(1,1));
fprintf('Top-right corner: lat=%.2f, lon=%.2f\n', lat_grid(1,end), lon_grid(1,end));
fprintf('Bottom-left corner: lat=%.2f, lon=%.2f\n', lat_grid(end,1), lon_grid(end,1));
fprintf('Bottom-right corner: lat=%.2f, lon=%.2f\n', lat_grid(end,end), lon_grid(end,end));

if lat_grid(1,1) < lat_grid(end,1)
    fprintf('Note: Latitude increases from top to bottom (may need flipping for display)\n');
else
    fprintf('Note: Latitude decreases from top to bottom (standard orientation)\n');
end

end
