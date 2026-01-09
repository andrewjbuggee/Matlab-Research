%% Plot Instrument Footprints on Geographic Map

% This function plots the ground footprints of different satellite instruments
% to visualize their spatial coverage and overlap
%
% INPUTS:
%   modis - MODIS data structure (can be [] to skip)
%   emit  - EMIT data structure (can be [] to skip)
%   airs  - AIRS data structure (can be [] to skip)
%   amsr  - AMSR data structure (can be [] to skip)
%   pixels - structure specifying which pixels to plot for each instrument:
%       .modis - structure with .row and .col for MODIS pixel(s) 
%                OR [] if data is already filtered to single/few pixels
%       .emit  - structure with .row and .col for EMIT pixel(s)
%                OR [] if data is already filtered to single/few pixels
%       .airs  - structure with .row and .col for AIRS pixel(s)
%                OR [] if data is already filtered to single/few pixels
%       .amsr  - structure with .row and .col for AMSR pixel(s)
%                OR [] if data is already filtered to single/few pixels
%   options - structure with optional plotting parameters:
%       .show_centers - true/false to show pixel centers (default: true)
%       .show_labels - true/false to show instrument labels (default: true)
%       .alpha - transparency for polygons (default: 0.3)
%       .colors - structure with RGB colors for each instrument
%       .show_rgb - true/false to show MODIS RGB image backdrop (default: false)
%       .rgb_image - Pre-computed RGB image (MxNx3) (optional)
%       .rgb_lat - Latitude grid for RGB image (MxN) (optional)
%       .rgb_lon - Longitude grid for RGB image (MxN) (optional)
%       .rgb_alpha - Transparency for RGB image (default: 0.8)
%       .tick_font_size - Font size for lat/lon axis labels (default: 16)
%       .latlim - [min_lat, max_lat] to limit display region (optional)
%       .lonlim - [min_lon, max_lon] to limit display region (optional)
%
% OUTPUT:
%   fig_handle - handle to the generated figure
%
% Example usage with full arrays:
%   pixels.modis.row = 1000;
%   pixels.modis.col = 500;
%   pixels.emit.row = 512;
%   pixels.emit.col = 256;
%   pixels.airs.row = 45;
%   pixels.airs.col = 15;
%   pixels.amsr.row = 100;
%   pixels.amsr.col = 50;
%   
%   fig_handle = plot_instrument_footprints(modis, emit, airs, AMSR, pixels);
%
% Example usage with filtered data:
%   % After using remove_unwanted_*_data functions
%   pixels.modis = [];  % Data already filtered
%   pixels.emit = [];   % Data already filtered
%   pixels.airs = [];   % Data already filtered
%   pixels.amsr = [];   % Data already filtered
%   
%   fig_handle = plot_instrument_footprints(modis_filtered, emit_filtered, airs_filtered, amsr_filtered, pixels);
%
% Example usage with RGB backdrop and geographic limits:
%   % First create the RGB image
%   [rgb_img, rgb_lat, rgb_lon] = create_modis_true_color(modis_full);
%   
%   % Then plot with RGB backdrop, limits, and custom font size
%   options.show_rgb = true;
%   options.rgb_image = rgb_img;
%   options.rgb_lat = rgb_lat;
%   options.rgb_lon = rgb_lon;
%   options.tick_font_size = 18;  % Larger font for axis labels
%   options.latlim = [-30, -20];  % Only show between -30 and -20 degrees latitude
%   options.lonlim = [-75, -65];  % Only show between -75 and -65 degrees longitude
%   
%   fig_handle = plot_instrument_footprints(modis, emit, airs, AMSR, pixels, options);
%
% By Andrew John Buggee

%%

function [fig_handle, gx] = plot_instrument_footprints_3(modis, emit, airs, amsr, pixels, options)

% Set default options
if nargin < 6
    options = struct();
end

if nargin < 5
    pixels = struct();
end

if ~isfield(options, 'show_centers')
    options.show_centers = true;
end

if ~isfield(options, 'show_labels')
    options.show_labels = true;
end

if ~isfield(options, 'alpha')
    options.alpha = 0.3;
end

if ~isfield(options, 'colors')
    % Default colors: MODIS (blue), EMIT (red), AIRS (orange), AMSR (green)
    options.colors.modis = [0.2, 0.4, 0.8];
    options.colors.emit = [0.8, 0.2, 0.2];
    options.colors.airs = [1.0, 0.6, 0.0];  % Orange for AIRS
    options.colors.amsr = [0.2, 0.8, 0.4];
end

if ~isfield(options, 'show_rgb')
    options.show_rgb = false;
end

if ~isfield(options, 'rgb_alpha')
    options.rgb_alpha = 0.8;
end

if ~isfield(options, 'tick_font_size')
    options.tick_font_size = 18;  % Default font size for lat/lon labels
end

% Create figure
fig_handle = figure('Position', [100, 100, 1200, 800]);

% Create geographic axes first (before hold on)
gx = geoaxes(fig_handle);
hold(gx, 'on');
% Set font size for latitude/longitude labels
gx.FontSize = options.tick_font_size;

%% Plot MODIS or EMIT RGB image backdrop if requested
if options.show_rgb && isfield(options, 'rgb_image') && ...
        isfield(options, 'rgb_lat') && isfield(options, 'rgb_lon')
    
    fprintf('Plotting RGB image backdrop...\n');
    
    % WORKAROUND: GeographicAxes doesn't support RGB image overlays directly
    % Solution: Use map axes (axesm) instead, which supports geoshow with images
    
    try
        % Get the geographic bounds from the data
        data_latlim = [min(options.rgb_lat(:), [], 'omitnan'), ...
                       max(options.rgb_lat(:), [], 'omitnan')];
        data_lonlim = [min(options.rgb_lon(:), [], 'omitnan'), ...
                       max(options.rgb_lon(:), [], 'omitnan')];
        
        % Apply user-specified limits if provided
        if isfield(options, 'latlim') && ~isempty(options.latlim)
            latlim = options.latlim;
            fprintf('Using user-specified latitude limits: [%.2f, %.2f]\n', latlim(1), latlim(2));
        else
            latlim = data_latlim;
        end
        
        if isfield(options, 'lonlim') && ~isempty(options.lonlim)
            lonlim = options.lonlim;
            fprintf('Using user-specified longitude limits: [%.2f, %.2f]\n', lonlim(1), lonlim(2));
        else
            lonlim = data_lonlim;
        end
        
        % Close current figure and recreate with map axes instead
        close(fig_handle);
        fig_handle = figure('Position', [100, 100, 1200, 800]);
        
        % Calculate appropriate tick spacing to get approximately 5 ticks on each axis
        lat_range = diff(latlim);
        lon_range = diff(lonlim);
        
        % Target approximately 5-6 tick marks per axis
        % Use nice round numbers: 0.1, 0.2, 0.5, 1, 2, 5, 10, etc.
        nice_intervals = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50];
        
        % Find the interval that gives closest to 5 ticks for latitude
        lat_intervals_per_nice = lat_range ./ nice_intervals;
        [~, idx_lat] = min(abs(lat_intervals_per_nice - 5));
        lat_label_spacing = nice_intervals(idx_lat);
        
        % Find the interval that gives closest to 5 ticks for longitude
        lon_intervals_per_nice = lon_range ./ nice_intervals;
        [~, idx_lon] = min(abs(lon_intervals_per_nice - 5));
        lon_label_spacing = nice_intervals(idx_lon);
        
        fprintf('Calculated label spacing: Lat = %.2f deg, Lon = %.2f deg\n', ...
                lat_label_spacing, lon_label_spacing);
        
        % Create a map axes using axesm
        % Use a geographic (unprojected) coordinate system for simplicity
        ax = axesm('MapProjection', 'eqdcylin', ...  % Equidistant cylindrical (Plate Carrée)
                   'MapLatLimit', latlim, ...
                   'MapLonLimit', lonlim, ...
                   'Frame', 'off', ...               % Remove black border
                   'Grid', 'off', ...
                   'MeridianLabel', 'on', ...
                   'ParallelLabel', 'on', ...
                   'MLabelLocation', lon_label_spacing, ...  % Auto-calculated longitude spacing
                   'PLabelLocation', lat_label_spacing, ...  % Auto-calculated latitude spacing
                   'LabelFormat', 'compass', ...     % Use compass format (N, S, E, W)
                   'LabelUnits', 'dm', ...           % Display in degrees-minutes
                   'FontSize', options.tick_font_size);  % Set font size for labels
        
        % Set label positions
        setm(ax, 'MLabelParallel', 'south');  % Put longitude labels on south edge
        setm(ax, 'PLabelMeridian', 'west');   % Put latitude labels on west edge

        % Explicitly set font size for meridian and parallel labels
        setm(ax, 'FontSize', options.tick_font_size);

        % Also set font size directly on the axes object
        set(ax, 'FontSize', options.tick_font_size);

        % Turn off box
        ax.BoxStyle = "full";
        ax.Box = "off";

        % Create a spatial referencing object for the image
        % IMPORTANT: georasterref expects [nrows, ncols] and the image to be
        % oriented with NORTH at TOP (latitude decreasing as rows increase)
        
        % Work on copies to avoid modifying the original data
        rgb_img = options.rgb_image;
        rgb_lat = options.rgb_lat;
        rgb_lon = options.rgb_lon;
        
        rasterSize = [size(rgb_img, 1), size(rgb_img, 2)];
        
        % Check if we need to flip the image
        % MODIS typically has latitude INCREASING from first row to last row
        % But display needs latitude DECREASING from first row to last row (north at top)
        fprintf('Before flip: First row lat = %.2f, Last row lat = %.2f\n', ...
                rgb_lat(1,1), rgb_lat(end,1));
        
        if rgb_lat(1,1) < rgb_lat(end,1) & strcmp(options.rgb_image_type, 'modis')==true
            % Latitude increases from top to bottom - need to flip
            fprintf('Flipping RGB image vertically (MODIS has south at top, need north at top)...\n');
            rgb_img = flipud(rgb_img);
            rgb_lat = flipud(rgb_lat);
            rgb_lon = flipud(rgb_lon);

            % Now verify orientation is correct
            fprintf('After flip: First row lat = %.2f, Last row lat = %.2f\n', ...
                rgb_lat(1,1), rgb_lat(end,1));

        elseif rgb_lat(1,1) < rgb_lat(end,1) & strcmp(options.rgb_image_type, 'emit')==true
            % Don't flip
            fprintf('Dont need to flip EMIT data.\n');
            rgb_img = flipud(rgb_img);
            rgb_lat = flipud(rgb_lat);
            rgb_lon = flipud(rgb_lon);

            % Now verify orientation is correct
            fprintf('After flip: First row lat = %.2f, Last row lat = %.2f\n', ...
                rgb_lat(1,1), rgb_lat(end,1));
            
        end

        
        
        % Crop RGB data to the specified limits (if provided)
        if isfield(options, 'latlim') || isfield(options, 'lonlim')
            fprintf('Cropping RGB image to specified geographic limits...\n');
            
            % Find pixels within the specified bounds
            in_bounds = (rgb_lat >= latlim(1)) & (rgb_lat <= latlim(2)) & ...
                       (rgb_lon >= lonlim(1)) & (rgb_lon <= lonlim(2));
            
            % Find the bounding box of pixels to keep
            [rows_with_data, cols_with_data] = find(in_bounds);
            
            if ~isempty(rows_with_data)
                row_min = min(rows_with_data);
                row_max = max(rows_with_data);
                col_min = min(cols_with_data);
                col_max = max(cols_with_data);
                
                % Crop to bounding box (with small buffer to avoid edge artifacts)
                rgb_img = rgb_img(row_min:row_max, col_min:col_max, :);
                rgb_lat = rgb_lat(row_min:row_max, col_min:col_max);
                rgb_lon = rgb_lon(row_min:row_max, col_min:col_max);
                
                fprintf('Cropped RGB from original size to [%dx%d] pixels\n', ...
                       size(rgb_img, 1), size(rgb_img, 2));
            else
                warning('No RGB data within specified limits. Using full RGB image.');
            end
        end
        
        % The issue: surfm/meshm don't handle RGB images well
        % Better approach: use geoshow but with lat/lon data directly
        
        fprintf('Displaying RGB using geoshow with lat/lon grids...\n');
        
        % geoshow can accept [lat, lon, RGB] format for textured images
        % This should work on map axes when given explicit coordinates
        h_img = geoshow(rgb_lat, rgb_lon, rgb_img, 'DisplayType', 'texturemap');
        
        % Set transparency
        if ~isempty(h_img)
            try
                set(h_img, 'AlphaData', options.rgb_alpha);
            catch
                % Alpha might not be supported in this mode
            end
        end
        
        fprintf('RGB image displayed on map.\n');
        
        % Set transparency if possible
        h_img = findobj(gca, 'Type', 'image');
        if ~isempty(h_img)
            set(h_img, 'AlphaData', options.rgb_alpha);
        end
        
        hold on;
        
        % Flag that we're using map axes instead of geoaxes
        using_map_axes = true;
        gx = ax;  % Use map axes as the plotting target
        
        fprintf('RGB image displayed successfully on map axes.\n');
        fprintf('Lat limits: [%.2f, %.2f]\n', latlim(1), latlim(2));
        fprintf('Lon limits: [%.2f, %.2f]\n', lonlim(1), lonlim(2));
        
    catch ME
        warning('Could not display RGB image: %s', ME.message);
        fprintf('Error details: %s\n', ME.getReport());
        fprintf('Proceeding without RGB backdrop...\n');
        
        % Recreate geoaxes if map axes failed
        close(fig_handle);
        fig_handle = figure('Position', [100, 100, 1200, 800]);
        gx = geoaxes(fig_handle);
        hold(gx, 'on');
        gx.FontSize = options.tick_font_size;  % Set font size for labels
        using_map_axes = false;
    end
    
else
    % No RGB image - use regular geoaxes
    using_map_axes = false;
end

% Store the axes type for later use
if ~exist('using_map_axes', 'var')
    using_map_axes = false;
end

% Initialize storage for legend entries
legend_entries = {};
legend_handles = [];

%% Plot MODIS footprints
if ~isempty(modis)
    
    % Check if pixels are specified or if data is already filtered
    if isfield(pixels, 'modis') && ~isempty(pixels.modis) && ...
            isfield(pixels.modis, 'row') && isfield(pixels.modis, 'col')
        % Pixels specified - data is full array
        modis_pixels = pixels.modis;
        data_is_filtered = false;
    else
        % No pixels specified - assume data is already filtered
        % Plot all pixels in the filtered data
        data_is_filtered = true;
        modis_pixels = struct();
    end
    
    if data_is_filtered
        % Data has been filtered - plot all remaining pixels
        % Determine number of pixels from data structure
        lat_size = size(modis.geo.lat);
        if length(lat_size) == 2 && lat_size(1) > 1
            % Still a 2D array (shouldn't happen after filtering, but check)
            num_modis = lat_size(1) * lat_size(2);
            [rows, cols] = meshgrid(1:lat_size(1), 1:lat_size(2));
            modis_pixels.row = rows(:);
            modis_pixels.col = cols(:);
        else
            % Filtered to vector
            num_modis = length(modis.geo.lat);
            modis_pixels.row = ones(num_modis, 1);  % Use index as row
            modis_pixels.col = (1:num_modis)';       % Use index as col
        end
    else
        % Using specified pixels
        num_modis = length(modis_pixels.row);
    end
    
    for ii = 1:num_modis
        
        row = modis_pixels.row(ii);
        col = modis_pixels.col(ii);
        
        % Get the footprint corners for this pixel
        [lat_corners, lon_corners] = get_modis_footprint(modis, row, col, data_is_filtered);
        
        % Plot the footprint polygon
        if using_map_axes
            % Use plotm for map axes
            h_modis = plotm(lat_corners, lon_corners, '-', ...
                'Color', options.colors.modis, 'LineWidth', 2);
        else
            % Use geoplot for geographic axes
            h_modis = geoplot(lat_corners, lon_corners, '-', ...
                'Color', options.colors.modis, 'LineWidth', 2);
        end
        
        % Plot center point if requested
        if options.show_centers
            if data_is_filtered
                % For filtered data, use direct indexing
                center_lat = modis.geo.lat(ii);
                center_lon = modis.geo.long(ii);
            else
                % For full array, use row/col indexing
                center_lat = modis.geo.lat(row, col);
                center_lon = modis.geo.long(row, col);
            end
            
            if using_map_axes
                plotm(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.modis, ...
                    'MarkerEdgeColor', options.colors.modis, 'MarkerSize', 6);
            else
                geoplot(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.modis, ...
                    'MarkerEdgeColor', options.colors.modis, 'MarkerSize', 6);
            end
        end
        
        % Store handle for legend (only for first pixel)
        if ii == 1
            legend_handles = [legend_handles, h_modis];
            legend_entries{end+1} = 'MODIS';
        end
    end
end

%% Plot EMIT footprints
if ~isempty(emit)
    
    % Check if pixels are specified or if data is already filtered
    if isfield(pixels, 'emit') && ~isempty(pixels.emit) && ...
            isfield(pixels.emit, 'row') && isfield(pixels.emit, 'col')
        % Pixels specified - data is full array
        emit_pixels = pixels.emit;
        data_is_filtered = false;
    else
        % No pixels specified - assume data is already filtered
        data_is_filtered = true;
        emit_pixels = struct();
    end
    
    if data_is_filtered
        % Data has been filtered - plot all remaining pixels
        lat_size = size(emit.radiance.geo.lat);
        if length(lat_size) == 2 && lat_size(1) > 1
            % Still a 2D array
            num_emit = lat_size(1) * lat_size(2);
            [rows, cols] = meshgrid(1:lat_size(1), 1:lat_size(2));
            emit_pixels.row = rows(:);
            emit_pixels.col = cols(:);
        else
            % Filtered to vector
            num_emit = length(emit.radiance.geo.lat);
            emit_pixels.row = ones(num_emit, 1);
            emit_pixels.col = (1:num_emit)';
        end
    else
        % Using specified pixels
        num_emit = length(emit_pixels.row);
    end
    
    for ii = 1:num_emit
        
        row = emit_pixels.row(ii);
        col = emit_pixels.col(ii);
        
        % Get the footprint corners for this pixel
        [lat_corners, lon_corners] = get_emit_footprint(emit, row, col, data_is_filtered);
        
        % Plot the footprint polygon
        if using_map_axes
            h_emit = plotm(lat_corners, lon_corners, '-', ...
                'Color', options.colors.emit, 'LineWidth', 2);
        else
            h_emit = geoplot(lat_corners, lon_corners, '-', ...
                'Color', options.colors.emit, 'LineWidth', 2);
        end
        
        % Plot center point if requested
        if options.show_centers
            if data_is_filtered
                center_lat = emit.radiance.geo.lat(ii);
                center_lon = emit.radiance.geo.long(ii);
            else
                center_lat = emit.radiance.geo.lat(row, col);
                center_lon = emit.radiance.geo.long(row, col);
            end
            
            if using_map_axes
                plotm(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.emit, ...
                    'MarkerEdgeColor', options.colors.emit, 'MarkerSize', 6);
            else
                geoplot(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.emit, ...
                    'MarkerEdgeColor', options.colors.emit, 'MarkerSize', 6);
            end
        end
        
        % Store handle for legend (only for first pixel)
        if ii == 1
            legend_handles = [legend_handles, h_emit];
            legend_entries{end+1} = 'EMIT';
        end
    end
end

%% Plot AIRS footprints
if ~isempty(airs)
    
    % Check if pixels are specified or if data is already filtered
    if isfield(pixels, 'airs') && ~isempty(pixels.airs) && ...
            isfield(pixels.airs, 'row') && isfield(pixels.airs, 'col')
        % Pixels specified - data is full array
        airs_pixels = pixels.airs;
        data_is_filtered = false;
    else
        % No pixels specified - assume data is already filtered
        data_is_filtered = true;
        airs_pixels = struct();
    end
    
    if data_is_filtered
        % Data has been filtered - plot all remaining pixels
        % AIRS data could be in different structures depending on processing level
        % Try to find latitude data
        if isfield(airs, 'geo') && isfield(airs.geo, 'Latitude')
            lat_data = airs.geo.Latitude;
        elseif isfield(airs, 'Latitude')
            lat_data = airs.Latitude;
        else
            error('Cannot find latitude data in AIRS structure');
        end
        
        lat_size = size(lat_data);
        if length(lat_size) == 2 && lat_size(1) > 1
            % Still a 2D array
            num_airs = lat_size(1) * lat_size(2);
            [rows, cols] = meshgrid(1:lat_size(1), 1:lat_size(2));
            airs_pixels.row = rows(:);
            airs_pixels.col = cols(:);
        else
            % Filtered to vector
            num_airs = length(lat_data);
            airs_pixels.row = ones(num_airs, 1);
            airs_pixels.col = (1:num_airs)';
        end
    else
        % Using specified pixels
        num_airs = length(airs_pixels.row);
    end
    
    for ii = 1:num_airs
        
        row = airs_pixels.row(ii);
        col = airs_pixels.col(ii);
        
        % Get the footprint corners for this pixel
        [lat_corners, lon_corners] = get_airs_footprint(airs, row, col, data_is_filtered);
        
        % Plot the footprint polygon
        if using_map_axes
            % Use plotm for map axes
            h_airs = plotm(lat_corners, lon_corners, '-', ...
                'Color', options.colors.airs, 'LineWidth', 2);
        else
            % Use geoplot for geographic axes
            h_airs = geoplot(lat_corners, lon_corners, '-', ...
                'Color', options.colors.airs, 'LineWidth', 2);
        end
        
        % Plot center point if requested
        if options.show_centers
            if data_is_filtered
                % For filtered data, use direct indexing
                if isfield(airs, 'geo') && isfield(airs.geo, 'Latitude')
                    center_lat = airs.geo.Latitude(ii);
                    center_lon = airs.geo.Longitude(ii);
                elseif isfield(airs, 'Latitude')
                    center_lat = airs.Latitude(ii);
                    center_lon = airs.Longitude(ii);
                end
            else
                % For full array, use row/col indexing
                if isfield(airs, 'geo') && isfield(airs.geo, 'Latitude')
                    center_lat = airs.geo.Latitude(row, col);
                    center_lon = airs.geo.Longitude(row, col);
                elseif isfield(airs, 'Latitude')
                    center_lat = airs.Latitude(row, col);
                    center_lon = airs.Longitude(row, col);
                end
            end
            
            if using_map_axes
                plotm(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.airs, ...
                    'MarkerEdgeColor', options.colors.airs, 'MarkerSize', 6);
            else
                geoplot(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.airs, ...
                    'MarkerEdgeColor', options.colors.airs, 'MarkerSize', 6);
            end
        end
        
        % Store handle for legend (only for first pixel)
        if ii == 1
            legend_handles = [legend_handles, h_airs];
            legend_entries{end+1} = 'AIRS';
        end
    end
end

%% Plot AMSR footprints
if ~isempty(amsr)
    
    % Check if pixels are specified or if data is already filtered
    if isfield(pixels, 'amsr') && ~isempty(pixels.amsr) && ...
            isfield(pixels.amsr, 'row') && isfield(pixels.amsr, 'col')
        % Pixels specified - data is full array
        amsr_pixels = pixels.amsr;
        data_is_filtered = false;
    else
        % No pixels specified - assume data is already filtered
        data_is_filtered = true;
        amsr_pixels = struct();
    end
    
    if data_is_filtered
        % Data has been filtered - plot all remaining pixels
        lat_size = size(amsr.geo.Latitude);
        if length(lat_size) == 2 && lat_size(1) > 1
            % Still a 2D array
            num_amsr = lat_size(1) * lat_size(2);
            [rows, cols] = meshgrid(1:lat_size(1), 1:lat_size(2));
            amsr_pixels.row = rows(:);
            amsr_pixels.col = cols(:);
        else
            % Filtered to vector
            num_amsr = length(amsr.geo.Latitude);
            amsr_pixels.row = ones(num_amsr, 1);
            amsr_pixels.col = (1:num_amsr)';
        end
    else
        % Using specified pixels
        num_amsr = length(amsr_pixels.row);
    end
    
    for ii = 1:num_amsr
        
        row = amsr_pixels.row(ii);
        col = amsr_pixels.col(ii);
        
        % Get the footprint corners for this pixel
        [lat_corners, lon_corners] = get_amsr_footprint(amsr, row, col, data_is_filtered);
        
        % Plot the footprint polygon
        if using_map_axes
            h_amsr = plotm(lat_corners, lon_corners, '-', ...
                'Color', options.colors.amsr, 'LineWidth', 2);
        else
            h_amsr = geoplot(lat_corners, lon_corners, '-', ...
                'Color', options.colors.amsr, 'LineWidth', 2);
        end
        
        % Plot center point if requested
        if options.show_centers
            if data_is_filtered
                center_lat = amsr.geo.Latitude(ii);
                center_lon = amsr.geo.Longitude(ii);
            else
                center_lat = amsr.geo.Latitude(row, col);
                center_lon = amsr.geo.Longitude(row, col);
            end
            
            if using_map_axes
                plotm(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.amsr, ...
                    'MarkerEdgeColor', options.colors.amsr, 'MarkerSize', 6);
            else
                geoplot(center_lat, center_lon, 'o', ...
                    'MarkerFaceColor', options.colors.amsr, ...
                    'MarkerEdgeColor', options.colors.amsr, 'MarkerSize', 6);
            end
        end
        
        % Store handle for legend (only for first pixel)
        if ii == 1
            legend_handles = [legend_handles, h_amsr];
            legend_entries{end+1} = 'AMSR';
        end
    end
end

%% Finalize plot
if using_map_axes
    % Map axes - already has grid and labels
    % No basemap needed since we have RGB image
    title(gx, 'Instrument Footprint Comparison', 'FontSize', 14, 'FontWeight', 'bold');
else
    % Geographic axes - only show basemap if RGB is not displayed
    if ~options.show_rgb
        geobasemap(gx, 'topographic');  % or 'satellite', 'streets', 'landcover', etc.
    end
    title(gx, 'Instrument Footprint Comparison', 'FontSize', 14, 'FontWeight', 'bold');
    grid(gx, 'on');
end

% Add legend if requested
if options.show_labels && ~isempty(legend_handles)
    legend(gx, legend_handles, legend_entries, 'Location', 'best', ...
        'FontSize', 18, 'Interpreter', 'latex');
end

% Update tick label interpreter and font size
if using_map_axes
    % Map axes (when RGB is shown)
    gx.FontSize = options.tick_font_size;
    % Set interpreter for axis labels
    if isfield(gx, 'Title') && ~isempty(gx.Title)
        gx.Title.Interpreter = 'latex';
        gx.Title.FontSize = options.tick_font_size + 4;
    end
else
    % Geographic axes (when RGB is not shown)
    gx.FontSize = options.tick_font_size;
    
    % Set LaTeX interpreter for latitude axis
    % gx.LatitudeAxis.TickLabelInterpreter = 'latex';
    gx.LatitudeAxis.FontSize = options.tick_font_size;
    if ~isempty(gx.LatitudeLabel)
        gx.LatitudeLabel.FontSize = options.tick_font_size + 10;
        gx.LatitudeLabel.Interpreter = 'latex';
    end
    
    % Set LaTeX interpreter for longitude axis
    % gx.LongitudeAxis.TickLabelInterpreter = 'latex';
    gx.LongitudeAxis.FontSize = options.tick_font_size;
    if ~isempty(gx.LongitudeLabel)
        gx.LongitudeLabel.FontSize = options.tick_font_size + 10;
        gx.LongitudeLabel.Interpreter = 'latex';
    end
    
    % Set LaTeX interpreter for title
    if ~isempty(gx.Title)
        gx.Title.Interpreter = 'latex';
        gx.Title.FontSize = options.tick_font_size + 10;
    end
end




hold(gx, 'off');

end


%% Helper function: Get MODIS footprint corners
function [lat_corners, lon_corners] = get_modis_footprint(modis, row, col, data_is_filtered)
    % MODIS pixel size grows with scan angle
    % At nadir: ~1 km x 1 km
    % At scan edge (55 deg): ~2 km (along-track) x 5 km (along-scan)
    
    % Get pixel center location
    if data_is_filtered
        % Data is filtered - use direct indexing
        center_lat = modis.geo.lat(row);
        center_lon = modis.geo.long(row);
        
        % Get sensor zenith angle if available
        if isfield(modis, 'sensor') && isfield(modis.sensor, 'zenith')
            sensor_zenith = modis.sensor.zenith(row);
        else
            % Estimate from neighboring pixels if possible
            sensor_zenith = 25;  % Assume mid-range value
        end
    else
        % Data is full array - use row/col indexing
        center_lat = modis.geo.lat(row, col);
        center_lon = modis.geo.long(row, col);
        
        % Get sensor zenith angle if available (to estimate pixel distortion)
        if isfield(modis, 'sensor') && isfield(modis.sensor, 'zenith')
            sensor_zenith = modis.sensor.zenith(row, col);  % degrees
        else
            % If no sensor geometry, estimate from scan position
            [~, ncols] = size(modis.geo.lat);
            % Assume center column is nadir, edges are ~55 degrees
            scan_fraction = abs(col - ncols/2) / (ncols/2);
            sensor_zenith = scan_fraction * 55;  % degrees
        end
    end
    
    % Calculate pixel size based on zenith angle
    % At nadir (0 deg): 1 km x 1 km
    % At edge (55 deg): ~2 km x 5 km
    % Empirical relationship from MODIS ATBD
    along_track_size = 1.0 * (1 + 0.02 * (sensor_zenith/10)^2);  % km
    along_scan_size = 1.0 * (1 / cosd(sensor_zenith));  % km
    
    % Convert km to approximate degrees (rough approximation)
    lat_per_km = 1.0 / 111.0;  % degrees latitude per km
    lon_per_km = 1.0 / (111.0 * cosd(center_lat));  % degrees longitude per km
    
    half_along_track = (along_track_size / 2) * lat_per_km;
    half_along_scan = (along_scan_size / 2) * lon_per_km;
    
    % Create rectangular footprint (simplified - actual MODIS pixels bow-tie)
    lat_corners = zeros(5, 1);
    lon_corners = zeros(5, 1);
    
    % Top-left corner
    lat_corners(1) = center_lat + half_along_track;
    lon_corners(1) = center_lon - half_along_scan;
    
    % Top-right corner
    lat_corners(2) = center_lat + half_along_track;
    lon_corners(2) = center_lon + half_along_scan;
    
    % Bottom-right corner
    lat_corners(3) = center_lat - half_along_track;
    lon_corners(3) = center_lon + half_along_scan;
    
    % Bottom-left corner
    lat_corners(4) = center_lat - half_along_track;
    lon_corners(4) = center_lon - half_along_scan;
    
    % Close the polygon
    lat_corners(5) = lat_corners(1);
    lon_corners(5) = lon_corners(1);
end


%% Helper function: Get EMIT footprint corners
function [lat_corners, lon_corners] = get_emit_footprint(emit, row, col, data_is_filtered)
    % EMIT is a pushbroom imager with ~60m pixel size
    % Footprint is relatively consistent across the swath
    
    if data_is_filtered
        % Data is filtered - use direct indexing
        center_lat = emit.radiance.geo.lat(row);
        center_lon = emit.radiance.geo.long(row);
    else
        % Data is full array - use row/col indexing
        center_lat = emit.radiance.geo.lat(row, col);
        center_lon = emit.radiance.geo.long(row, col);
    end
    
    % EMIT pixel size: approximately 60m x 60m
    pixel_size_km = 0.060;  % km
    
    % Convert km to approximate degrees
    lat_per_km = 1.0 / 111.0;  % degrees latitude per km
    lon_per_km = 1.0 / (111.0 * cosd(center_lat));  % degrees longitude per km
    
    half_size_lat = (pixel_size_km / 2) * lat_per_km;
    half_size_lon = (pixel_size_km / 2) * lon_per_km;
    
    lat_corners = zeros(5, 1);
    lon_corners = zeros(5, 1);
    
    % Top-left corner
    lat_corners(1) = center_lat + half_size_lat;
    lon_corners(1) = center_lon - half_size_lon;
    
    % Top-right corner
    lat_corners(2) = center_lat + half_size_lat;
    lon_corners(2) = center_lon + half_size_lon;
    
    % Bottom-right corner
    lat_corners(3) = center_lat - half_size_lat;
    lon_corners(3) = center_lon + half_size_lon;
    
    % Bottom-left corner
    lat_corners(4) = center_lat - half_size_lat;
    lon_corners(4) = center_lon - half_size_lon;
    
    % Close the polygon
    lat_corners(5) = lat_corners(1);
    lon_corners(5) = lon_corners(1);
end


%% Helper function: Get AIRS footprint corners
function [lat_corners, lon_corners] = get_airs_footprint(airs, row, col, data_is_filtered)
    % AIRS L2 retrieval footprint
    % The L2 standard retrieval uses a 50 km diameter footprint (circular)
    % This is larger than the L1B infrared footprint (13.5 km) because
    % the retrieval combines multiple L1B footprints
    % Reference: AIRS Version 7 L2 Product User Guide
    
    % Get pixel center location
    if data_is_filtered
        % Data is filtered - use direct indexing
        % Try different possible structure formats
        if isfield(airs, 'geo') && isfield(airs.geo, 'Latitude')
            center_lat = airs.geo.Latitude(row);
            center_lon = airs.geo.Longitude(row);
        elseif isfield(airs, 'Latitude')
            center_lat = airs.Latitude(row);
            center_lon = airs.Longitude(row);
        else
            error('Cannot find latitude/longitude in AIRS structure');
        end
        
        % Get sensor zenith angle if available
        if isfield(airs, 'sensor') && isfield(airs.sensor, 'zen')
            sensor_zenith = airs.sensor.zen(row);
        elseif isfield(airs, 'scanang')
            sensor_zenith = airs.scanang(row);
        elseif isfield(airs, 'satzen')
            sensor_zenith = airs.satzen(row);
        else
            % Estimate from position if not available
            sensor_zenith = 25;  % Assume mid-range value
        end
    else
        % Data is full array - use row/col indexing
        if isfield(airs, 'geo') && isfield(airs.geo, 'Latitude')
            center_lat = airs.geo.Latitude(row, col);
            center_lon = airs.geo.Longitude(row, col);
        elseif isfield(airs, 'Latitude')
            center_lat = airs.Latitude(row, col);
            center_lon = airs.Longitude(row, col);
        else
            error('Cannot find latitude/longitude in AIRS structure');
        end
        
        % Get sensor zenith angle if available
        if isfield(airs, 'sensor') && isfield(airs.sensor, 'zen')
            sensor_zenith = airs.sensor.zen(row, col);
        elseif isfield(airs, 'scanang')
            sensor_zenith = airs.scanang(row, col);
        elseif isfield(airs, 'satzen')
            sensor_zenith = airs.satzen(row, col);
        else
            % If no sensor geometry, estimate from scan position
            % AIRS has 90 footprints per scan line
            [~, ncols] = size(airs.geo.Latitude);
            scan_fraction = abs(col - ncols/2) / (ncols/2);
            sensor_zenith = scan_fraction * 49.5;  % AIRS max scan angle is ~49.5 degrees
        end
    end
    
    % Calculate footprint size based on zenith angle
    % At nadir (0 deg): 50 km diameter → 25 km radius (L2 retrieval footprint)
    % The footprint grows with zenith angle due to slant path through atmosphere
    % Using a simplified model similar to MODIS:
    % The along-track dimension grows modestly with zenith
    % The along-scan dimension grows as 1/cos(zenith)
    
    nadir_radius_km = 50 / 2;  % 25 km at nadir (L2 retrieval footprint)
    
    % Along-track stretching (modest growth)
    along_track_size = 2 * nadir_radius_km * (1 + 0.01 * (sensor_zenith/10)^2);  % km
    
    % Along-scan stretching (more pronounced growth)
    along_scan_size = 2 * nadir_radius_km * (1 / cosd(sensor_zenith));  % km
    
    % Convert km to approximate degrees
    lat_per_km = 1.0 / 111.0;  % degrees latitude per km
    lon_per_km = 1.0 / (111.0 * cosd(center_lat));  % degrees longitude per km
    
    half_along_track = (along_track_size / 2) * lat_per_km;
    half_along_scan = (along_scan_size / 2) * lon_per_km;
    
    % Create rectangular footprint (simplified representation of circular footprint)
    lat_corners = zeros(5, 1);
    lon_corners = zeros(5, 1);
    
    % Top-left corner
    lat_corners(1) = center_lat + half_along_track;
    lon_corners(1) = center_lon - half_along_scan;
    
    % Top-right corner
    lat_corners(2) = center_lat + half_along_track;
    lon_corners(2) = center_lon + half_along_scan;
    
    % Bottom-right corner
    lat_corners(3) = center_lat - half_along_track;
    lon_corners(3) = center_lon + half_along_scan;
    
    % Bottom-left corner
    lat_corners(4) = center_lat - half_along_track;
    lon_corners(4) = center_lon - half_along_scan;
    
    % Close the polygon
    lat_corners(5) = lat_corners(1);
    lon_corners(5) = lon_corners(1);
end


%% Helper function: Get AMSR footprint corners
function [lat_corners, lon_corners] = get_amsr_footprint(amsr, row, col, data_is_filtered)
    % AMSR has much larger footprint (~10-50 km depending on channel)
    % Footprint size varies by frequency band
    % We'll estimate from the spacing between pixels
    
    if data_is_filtered
        % Data is filtered - use direct indexing
        center_lat = amsr.geo.Latitude(row);
        center_lon = amsr.geo.Longitude(row);
        
        % For filtered data with few pixels, use a default footprint size
        % Assume ~10 km footprint (conservative estimate)
        half_dlat = 10 / 111.0 / 2;  % degrees
        half_dlon = 10 / (111.0 * cosd(center_lat)) / 2;  % degrees
        
    else
        % Data is full array - estimate from neighboring pixels
        center_lat = amsr.geo.Latitude(row, col);
        center_lon = amsr.geo.Longitude(row, col);
        
        % Estimate footprint size from neighboring pixels
        [nrows, ncols] = size(amsr.geo.Latitude);
        
        % Calculate approximate pixel spacing
        dlat = 0;
        dlon = 0;
        count = 0;
        
        if row > 1
            dlat = dlat + abs(amsr.geo.Latitude(row, col) - amsr.geo.Latitude(row-1, col));
            dlon = dlon + abs(amsr.geo.Longitude(row, col) - amsr.geo.Longitude(row-1, col));
            count = count + 1;
        end
        
        if row < nrows
            dlat = dlat + abs(amsr.geo.Latitude(row, col) - amsr.geo.Latitude(row+1, col));
            dlon = dlon + abs(amsr.geo.Longitude(row, col) - amsr.geo.Longitude(row+1, col));
            count = count + 1;
        end
        
        if col > 1
            dlat = dlat + abs(amsr.geo.Latitude(row, col) - amsr.geo.Latitude(row, col-1));
            dlon = dlon + abs(amsr.geo.Longitude(row, col) - amsr.geo.Longitude(row, col-1));
            count = count + 1;
        end
        
        if col < ncols
            dlat = dlat + abs(amsr.geo.Latitude(row, col) - amsr.geo.Latitude(row, col+1));
            dlon = dlon + abs(amsr.geo.Longitude(row, col) - amsr.geo.Longitude(row, col+1));
            count = count + 1;
        end
        
        if count > 0
            half_dlat = (dlat / count) / 2;
            half_dlon = (dlon / count) / 2;
        else
            % Fallback: assume ~10 km footprint
            half_dlat = 10 / 111.0 / 2;  % degrees
            half_dlon = 10 / (111.0 * cosd(center_lat)) / 2;  % degrees
        end
    end
    
    lat_corners = zeros(5, 1);
    lon_corners = zeros(5, 1);
    
    % Top-left corner
    lat_corners(1) = center_lat + half_dlat;
    lon_corners(1) = center_lon - half_dlon;
    
    % Top-right corner
    lat_corners(2) = center_lat + half_dlat;
    lon_corners(2) = center_lon + half_dlon;
    
    % Bottom-right corner
    lat_corners(3) = center_lat - half_dlat;
    lon_corners(3) = center_lon + half_dlon;
    
    % Bottom-left corner
    lat_corners(4) = center_lat - half_dlat;
    lon_corners(4) = center_lon - half_dlon;
    
    % Close the polygon
    lat_corners(5) = lat_corners(1);
    lon_corners(5) = lon_corners(1);
end
