%% Finding overlapping data between EMIT and MODIS
%
% This function identifies overlapping pixels between EMIT, Aqua/MODIS, and AIRS
% datasets based on spatial proximity and user-defined cloud property criteria.
% For each MODIS pixel meeting the criteria, it finds the closest EMIT pixel,
% and for each EMIT pixel, it finds the closest AIRS pixel.
%
% INPUTS:
%   folder_paths - Structure containing paths to data files with fields:
%       .coincident_dataPath   - Base path to the coincident data directory
%       .coincident_dataFolder - Folder name containing EMIT, MODIS, and AIRS files
%
%   criteria - Structure defining cloud filtering criteria with fields:
%       .cld_phase    - Cloud phase to filter ('water' for liquid water clouds)
%       .cld_tau_min  - Minimum cloud optical thickness (tau)
%       .cld_tau_max  - Maximum cloud optical thickness (tau)
%       .H            - Maximum horizontal inhomogeneity index threshold
%                       (H_860 < criteria.H will be selected)
%
%   plot_data - (Optional) Flag or structure for plotting control
%               (currently not used in function but reserved for future plotting)
%
% OUTPUTS:
%   overlap - Structure containing matched pixel indices and distances with fields:
%       .modis - MODIS pixel information:
%           .row           - Row indices in MODIS array (n_matches x 1)
%           .col           - Column indices in MODIS array (n_matches x 1)
%           .linear_idx    - Linear indices for MODIS array access (n_matches x 1)
%       .emit - EMIT pixel information:
%           .row           - Row indices in EMIT array (n_matches x 1)
%           .col           - Column indices in EMIT array (n_matches x 1)
%           .linear_idx    - Linear indices for EMIT array access (n_matches x 1)
%       .airs - AIRS pixel information:
%           .row           - Row indices in AIRS array (n_matches x 1)
%           .col           - Column indices in AIRS array (n_matches x 1)
%           .linear_idx    - Linear indices for AIRS array access (n_matches x 1)
%       .distance_modis_emit_km - Distance between MODIS and EMIT pixels (km)
%       .distance_emit_airs_km  - Distance between EMIT and AIRS pixels (km)
%
%   emit - Structure containing EMIT L1B radiance data
%       (returned from retrieveEMIT_data function)
%
%   modis - Structure containing Aqua/MODIS L1B and L2 data
%       (returned from retrieveMODIS_data function)
%
%   airs - Structure containing AIRS L2 atmospheric profile data
%       (returned from readAIRS_L2_data function)
%
%   folder_paths - Updated structure with added fields:
%       .L1B_fileName_emit  - Name of the EMIT L1B file processed
%       .L1B_fileName_modis - Name of the MODIS L1B file processed
%
% FILTERING CRITERIA:
%   The function applies the following filters to MODIS pixels:
%   1. Within EMIT footprint (spatial overlap using polygon boundary)
%   2. Over ocean (land fraction < threshold)
%   3. Liquid water clouds (phase == 2)
%   4. Cloud optical thickness within specified range [tau_min, tau_max]
%   5. Effective radius retrieval uncertainty < 10%
%   6. Horizontal inhomogeneity index < H threshold
%
% MATCHING LOGIC:
%   - For each MODIS pixel meeting criteria → find closest EMIT pixel
%   - For each matched EMIT pixel → find closest AIRS pixel
%   - Distances calculated using WGS84 ellipsoid for geodetic accuracy
%
% EXAMPLE USAGE:
%   folder_paths.coincident_dataPath = '/path/to/data/';
%   folder_paths.coincident_dataFolder = 'coincident_obs_2024/';
%   
%   criteria.cld_phase = 'water';
%   criteria.cld_tau_min = 5;
%   criteria.cld_tau_max = 30;
%   criteria.H = 0.3;
%   
%   [overlap, emit, modis, airs, folder_paths] = ...
%       findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths, criteria, true);
%
% NOTES:
%   - EMIT footprint detection uses convex hull boundary with inpolygon()
%   - All distances computed using WGS84 geodetic calculations
%   - MODIS effective radius uncertainty threshold is hard-coded at 10%
%   - Ocean detection uses land_or_ocean() function with 20 km tolerance
%
% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - AIRS pixel acquisition time (available: airs.Time_UTC)
% - Then compute: 
%   overlap.time_difference_modis_emit_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))
%   overlap.time_difference_emit_airs_seconds(nn) = abs(emit_time(idx_emit) - airs_time(idx_airs))
%
% By Andrew John Buggee
%%

function [overlap, emit, modis, airs, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths, criteria, plot_data)



%% Load EMIT Data
[emit, folder_paths.L1B_fileName_emit] = retrieveEMIT_data([folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder]);

%% Load Aqua/MODIS Data
[modis, folder_paths.L1B_fileName_modis] = retrieveMODIS_data([folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder]);

%% Load AIRS data

% Find which files are AIRS data
airs = readAIRS_L2_data([folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder]);

%% Get MODIS coordinates
modis_lat = double(modis.geo.lat);
modis_long = double(modis.geo.long);

%% Apply criteria filters first (more efficient than checking spatial overlap first)

% Which cloud phase should we look for?
if strcmp(criteria.cld_phase, 'water')==true
    is_liquidWater = modis.cloud.phase==2;           % 2 is the value designated for liquid water
else
    error('What is the phase??')
end

% Set the optical depth limits
is_tauGreater_and_Less = modis.cloud.optThickness17 >= criteria.cld_tau_min &...
    modis.cloud.optThickness17 <= criteria.cld_tau_max;

% Set the criteria for the horizontal homogeneity index
is_H_lessThan = modis.cloud.inhomogeneity_index(:,:,2) < criteria.H;

% Find pixels where the effective radius retrieval uncertainty is less than 10 percent
is_reUncertLessThan10Percent = modis.cloud.effRad_uncert_17<=10;

% Check if over ocean
isOcean = reshape(land_or_ocean(modis_lat(:), modis_long(:), 20, false), size(modis_lat));

% Combine all criteria
idx_criteria = logical(is_liquidWater .* is_tauGreater_and_Less .*...
                       is_reUncertLessThan10Percent .* is_H_lessThan .* isOcean);

%% Find MODIS pixels within EMIT footprint using convex hull

% Get EMIT footprint boundary
emit_lat = emit.radiance.geo.lat;
emit_long = emit.radiance.geo.long;

% Create boundary of EMIT footprint (trace the perimeter)
% Get the four edges
top_edge = [emit_lat(1,:)', emit_long(1,:)'];
bottom_edge = [emit_lat(end,:)', emit_long(end,:)'];
left_edge = [emit_lat(:,1), emit_long(:,1)];
right_edge = [emit_lat(:,end), emit_long(:,end)];

% Combine to form boundary (counter-clockwise)
boundary_pts = [top_edge; 
                flipud(right_edge(2:end-1,:)); 
                flipud(bottom_edge); 
                left_edge(2:end-1,:)];

% Find MODIS pixels inside the EMIT boundary
idx_in_footprint = inpolygon(modis_lat(:), modis_long(:), boundary_pts(:,1), boundary_pts(:,2));
idx_in_footprint = reshape(idx_in_footprint, size(modis_lat));

%% Combine footprint constraint with other criteria
idx_combined_master = logical(idx_in_footprint .* idx_criteria);

fprintf('Found %d MODIS pixels meeting all criteria within EMIT footprint\n', sum(idx_combined_master(:)));

%% Find row/column indices for MODIS pixels and match with closest EMIT and AIRS pixels

n_matches = sum(idx_combined_master(:));

% Pre-allocate output arrays for MODIS
overlap.modis.row = zeros(n_matches, 1);
overlap.modis.col = zeros(n_matches, 1);
overlap.modis.linear_idx = zeros(n_matches, 1);    % Linear index for MODIS

% Pre-allocate output arrays for EMIT
overlap.emit.row = zeros(n_matches, 1);
overlap.emit.col = zeros(n_matches, 1);
overlap.emit.linear_idx = zeros(n_matches, 1);     % Linear index for EMIT

% Pre-allocate output arrays for AIRS
overlap.airs.row = zeros(n_matches, 1);
overlap.airs.col = zeros(n_matches, 1);
overlap.airs.linear_idx = zeros(n_matches, 1);     % Linear index for AIRS

% Pre-allocate distance arrays
overlap.distance_modis_emit_km = zeros(n_matches, 1);  % Distance between MODIS and EMIT
overlap.distance_emit_airs_km = zeros(n_matches, 1);   % Distance between EMIT and AIRS

% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - AIRS pixel acquisition time (available: airs.Time_UTC)
% - Then compute: 
%   overlap.time_difference_modis_emit_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))
%   overlap.time_difference_emit_airs_seconds(nn) = abs(emit_time(idx_emit) - airs_time(idx_airs))

% Setup WGS84 ellipsoid for accurate distance calculation
wgs84 = wgs84Ellipsoid('kilometer');

% Flatten EMIT coordinates for vectorized distance calculation
emit_lat_flat = emit_lat(:);
emit_long_flat = emit_long(:);

% Flatten AIRS coordinates for vectorized distance calculation
airs_lat = double(airs.geo.Latitude);
airs_long = double(airs.geo.Longitude);
airs_lat_flat = airs_lat(:);
airs_long_flat = airs_long(:);

% Get linear indices of matching MODIS pixels
modis_linear_idx = find(idx_combined_master);

for nn = 1:n_matches
    
    % ====== MODIS Pixel Information ======
    % Get the linear index for this MODIS pixel
    overlap.modis.linear_idx(nn) = modis_linear_idx(nn);
    
    % Convert linear index to row/col for MODIS
    [overlap.modis.row(nn), overlap.modis.col(nn)] = ind2sub(size(modis_lat), modis_linear_idx(nn));
    
    % Get lat/lon for this MODIS pixel
    modis_lat_current = modis_lat(modis_linear_idx(nn));
    modis_long_current = modis_long(modis_linear_idx(nn));
    
    
    % ====== Find Closest EMIT Pixel ======
    % First do a quick squared-distance filter to reduce candidates
    dist_sq_emit = (emit_lat_flat - modis_lat_current).^2 + (emit_long_flat - modis_long_current).^2;
    [~, idx_emit] = min(dist_sq_emit);
    
    % Store EMIT linear index
    overlap.emit.linear_idx(nn) = idx_emit;
    
    % Convert linear index to row/col for EMIT
    [overlap.emit.row(nn), overlap.emit.col(nn)] = ind2sub(size(emit_lat), idx_emit);
    
    % Get lat/lon for this EMIT pixel
    emit_lat_current = emit_lat_flat(idx_emit);
    emit_long_current = emit_long_flat(idx_emit);
    
    % Calculate accurate distance between MODIS and EMIT using WGS84 ellipsoid
    overlap.distance_modis_emit_km(nn) = distance(modis_lat_current, modis_long_current, ...
        emit_lat_current, emit_long_current, wgs84);
    
    
    % ====== Find Closest AIRS Pixel ======
    % Use the EMIT pixel location (not MODIS) to find the closest AIRS pixel
    % First do a quick squared-distance filter to reduce candidates
    dist_sq_airs = (airs_lat_flat - emit_lat_current).^2 + (airs_long_flat - emit_long_current).^2;
    [~, idx_airs] = min(dist_sq_airs);
    
    % Store AIRS linear index
    overlap.airs.linear_idx(nn) = idx_airs;
    
    % Convert linear index to row/col for AIRS
    [overlap.airs.row(nn), overlap.airs.col(nn)] = ind2sub(size(airs_lat), idx_airs);
    
    % Get lat/lon for this AIRS pixel
    airs_lat_current = airs_lat_flat(idx_airs);
    airs_long_current = airs_long_flat(idx_airs);
    
    % Calculate accurate distance between EMIT and AIRS using WGS84 ellipsoid
    overlap.distance_emit_airs_km(nn) = distance(emit_lat_current, emit_long_current, ...
        airs_lat_current, airs_long_current, wgs84);
    
end

fprintf('Average distance between matched MODIS-EMIT pixels: %.2f km\n', mean(overlap.distance_modis_emit_km));
fprintf('Median distance: %.2f km\n', median(overlap.distance_modis_emit_km));
fprintf('Max distance: %.2f km\n', max(overlap.distance_modis_emit_km));

fprintf('\nAverage distance between matched EMIT-AIRS pixels: %.2f km\n', mean(overlap.distance_emit_airs_km));
fprintf('Median distance: %.2f km\n', median(overlap.distance_emit_airs_km));
fprintf('Max distance: %.2f km\n', max(overlap.distance_emit_airs_km));

end
