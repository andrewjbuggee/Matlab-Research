%% Create True Color RGB Image from EMIT Data

% This function creates a true color (RGB) composite image from EMIT 
% hyperspectral data by selecting appropriate wavelength bands for red, 
% green, and blue channels
%
% INPUTS:
%   emit - EMIT data structure containing radiance data
%          Must have emit.radiance.measurements, emit.radiance.wavelength
%   options - structure with optional parameters:
%       .red_wl - target wavelength for red channel in nm (default: 650)
%       .green_wl - target wavelength for green channel in nm (default: 560)
%       .blue_wl - target wavelength for blue channel in nm (default: 470)
%       .gamma - gamma correction factor (default: 2.2)
%       .brightness - brightness adjustment factor (default: 1.0)
%       .contrast - contrast stretch percentiles [low, high] (default: [2, 98])
%       .convert_to_reflectance - true to convert radiance to reflectance (default: true)
%
% OUTPUTS:
%   rgb_image - MxNx3 array containing the RGB image (values 0-1)
%   lat_grid - MxN array of latitude coordinates
%   lon_grid - MxN array of longitude coordinates
%   band_indices - structure with .red, .green, .blue band indices used
%
% EMIT Wavelength Selection:
%   Red channel: ~650 nm (similar to MODIS Band 1: 620-670 nm)
%   Green channel: ~560 nm (similar to MODIS Band 4: 545-565 nm)
%   Blue channel: ~470 nm (similar to MODIS Band 3: 459-479 nm)
%
% Example usage:
%   [rgb, lat, lon] = create_emit_true_color(emit);
%   [rgb, lat, lon, bands] = create_emit_true_color(emit, struct('gamma', 1.8));
%
% By Andrew John Buggee

%%

function [rgb_image, lat_grid, lon_grid, band_indices] = create_emit_true_color(emit, options)

% Set default options
if nargin < 2
    options = struct();
end

if ~isfield(options, 'red_wl')
    options.red_wl = 650;  % nm
end

if ~isfield(options, 'green_wl')
    options.green_wl = 560;  % nm
end

if ~isfield(options, 'blue_wl')
    options.blue_wl = 470;  % nm
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

if ~isfield(options, 'convert_to_reflectance')
    options.convert_to_reflectance = true;
end


%% Check data structure

if ~isfield(emit, 'radiance')
    error('EMIT structure must contain measurements field.');
end


%% Find closest wavelength bands for RGB

% Find index of band closest to target red wavelength
[~, red_idx] = min(abs(emit.radiance.wavelength - options.red_wl));
% Find index of band closest to target green wavelength
[~, green_idx] = min(abs(emit.radiance.wavelength - options.green_wl));
% Find index of band closest to target blue wavelength
[~, blue_idx] = min(abs(emit.radiance.wavelength - options.blue_wl));

% Store band information
band_indices.red = red_idx;
band_indices.green = green_idx;
band_indices.blue = blue_idx;
band_indices.red_wl = emit.radiance.wavelength(red_idx);
band_indices.green_wl = emit.radiance.wavelength(green_idx);
band_indices.blue_wl = emit.radiance.wavelength(blue_idx);

disp(newline)
fprintf('Selected EMIT bands for RGB composite:\n');
fprintf('  Red:   Band %d (%.2f nm, target: %.2f nm)\n', red_idx, band_indices.red_wl, options.red_wl);
fprintf('  Green: Band %d (%.2f nm, target: %.2f nm)\n', green_idx, band_indices.green_wl, options.green_wl);
fprintf('  Blue:  Band %d (%.2f nm, target: %.2f nm)\n', blue_idx, band_indices.blue_wl, options.blue_wl);
disp(newline)

%% Extract RGB bands from EMIT data

% Extract the three bands
red_band = double(emit.radiance.measurements(:, :, red_idx));
green_band = double(emit.radiance.measurements(:, :, green_idx));
blue_band = double(emit.radiance.measurements(:, :, blue_idx));


%% Process each band

% Remove fill values and invalid data
% EMIT uses -9999 as fill value
fill_value = -9999;
red_band(red_band <= fill_value) = NaN;
green_band(green_band <= fill_value) = NaN;
blue_band(blue_band <= fill_value) = NaN;

% Also remove any negative values (shouldn't exist in valid radiance)
red_band(red_band < 0) = NaN;
green_band(green_band < 0) = NaN;
blue_band(blue_band < 0) = NaN;


%% Convert radiance to reflectance if requested

if options.convert_to_reflectance
    % EMIT radiance is in microW/cm^2/nm/sr
    % To convert to reflectance, we need solar irradiance and geometry
    
    % Check if we have solar zenith angle for proper conversion
    if isfield(emit, 'obs') && isfield(emit.obs, 'solar')
        % Use proper conversion with solar zenith angle
        sza = emit.obs.sensor.zenith;  % degrees
        mu0 = cosd(sza);  % cosine of solar zenith angle
        
        % Approximate solar irradiance at top of atmosphere (W/m^2/nm)
        % These are approximate values for the visible spectrum
        F0_red = 1.5;    % W/m^2/nm at ~650 nm
        F0_green = 1.8;  % W/m^2/nm at ~560 nm  
        F0_blue = 1.9;   % W/m^2/nm at ~470 nm
        
        % Convert EMIT radiance from microW/cm^2/nm/sr to W/m^2/nm/sr
        % 1 microW/cm^2 = 0.01 W/m^2
        red_band = red_band * 0.01;
        green_band = green_band * 0.01;
        blue_band = blue_band * 0.01;
        
        % Convert to reflectance: rho = pi * L / (F0 * mu0)
        red_band = pi * red_band ./ (F0_red * mu0);
        green_band = pi * green_band ./ (F0_green * mu0);
        blue_band = pi * blue_band ./ (F0_blue * mu0);
        
        fprintf('Converted radiance to reflectance using solar geometry.\n');
    else
        % Simple normalization approach (less accurate)
        warning('Solar geometry not available. Using simple normalization.');
        red_band = red_band / prctile(red_band(:), 99, 'omitnan');
        green_band = green_band / prctile(green_band(:), 99, 'omitnan');
        blue_band = blue_band / prctile(blue_band(:), 99, 'omitnan');
    end
else
    % Just normalize radiance to 0-1 range for visualization
    red_band = red_band / max(red_band(:), [], 'omitnan');
    green_band = green_band / max(green_band(:), [], 'omitnan');
    blue_band = blue_band / max(blue_band(:), [], 'omitnan');
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

if ~isempty(red_valid) && numel(red_valid) > 1
    red_min = prctile(red_valid, options.contrast(1));
    red_max = prctile(red_valid, options.contrast(2));
    if red_max > red_min
        red_band = (red_band - red_min) / (red_max - red_min);
    end
end

if ~isempty(green_valid) && numel(green_valid) > 1
    green_min = prctile(green_valid, options.contrast(1));
    green_max = prctile(green_valid, options.contrast(2));
    if green_max > green_min
        green_band = (green_band - green_min) / (green_max - green_min);
    end
end

if ~isempty(blue_valid) && numel(blue_valid) > 1
    blue_min = prctile(blue_valid, options.contrast(1));
    blue_max = prctile(blue_valid, options.contrast(2));
    if blue_max > blue_min
        blue_band = (blue_band - blue_min) / (blue_max - blue_min);
    end
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
if isfield(emit.radiance, 'geo') && isfield(emit.radiance.geo, 'lat') && isfield(emit.radiance.geo, 'long')

    lat_grid = double(emit.radiance.geo.lat);
    lon_grid = double(emit.radiance.geo.long);
    
    % Handle fill values in coordinates
    lat_grid(lat_grid <= -9999) = NaN;
    lon_grid(lon_grid <= -9999) = NaN;
    
    % Handle EMIT's lon field name (could be 'long' or 'lon')
    if sum(isnan(lon_grid(:))) == numel(lon_grid) && isfield(emit.radiance.geo, 'long')
        lon_grid = double(emit.radiance.geo.long);
        lon_grid(lon_grid <= -9999) = NaN;
    end
else
    warning('Geographic coordinates not found in EMIT structure.');
    lat_grid = [];
    lon_grid = [];
end


%% Display summary
disp(newline)
fprintf('\nTrue color RGB image created successfully.\n');
fprintf('Image size: %d x %d pixels\n', size(rgb_image, 1), size(rgb_image, 2));

if ~isempty(lat_grid) && ~isempty(lon_grid)
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
disp(newline)


end
