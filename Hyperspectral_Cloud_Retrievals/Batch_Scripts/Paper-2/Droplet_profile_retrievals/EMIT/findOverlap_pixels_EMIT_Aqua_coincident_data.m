%% Finding overlapping data between EMIT and MODIS, AIRS, and AMSR on the Aqua satellite
%
% This function identifies overlapping pixels between EMIT, Aqua/MODIS, AIRS, and AMSR
% datasets based on spatial proximity and user-defined cloud property criteria.
% For each MODIS pixel meeting the criteria, it finds ALL EMIT pixels whose
% centers fall completely within the MODIS pixel boundary (defined by the
% grid spacing of the MODIS lat/lon arrays), then randomly selects 10 of
% those EMIT pixels. For each matched MODIS pixel, it also finds the
% closest AIRS and AMSR pixels.
%
% The output arrays are flat vectors of length n_matches * 10, where each
% group of 10 consecutive rows corresponds to one MODIS pixel. The MODIS,
% AIRS, and AMSR indices are repeated 10 times per match.
%
% INPUTS:
%   folder_paths - Structure containing paths to data files with fields:
%       .coincident_dataPath   - Base path to the coincident data directory
%       .coincident_dataFolder - Folder name containing EMIT, MODIS, AIRS, and AMSR files
%
%   criteria - Structure defining cloud filtering criteria with fields:
%       .cld_phase       - Cloud phase to filter ('water' for liquid water clouds)
%       .cld_tau_min     - Minimum cloud optical thickness (tau)
%       .cld_tau_max     - Maximum cloud optical thickness (tau)
%       .findN_smallest_H - (logical) Flag to select filtering mode:
%           false (default): Use criteria.H as a threshold. All MODIS
%                            pixels with H < criteria.H are selected.
%           true:            Ignore criteria.H. Instead, find all MODIS
%                            pixels meeting phase, tau, re uncertainty,
%                            and ocean criteria, then select the
%                            criteria.H_N_smallest pixels with the
%                            smallest inhomogeneity index values.
%       .H               - Maximum horizontal inhomogeneity index threshold
%                           (used when findN_smallest_H is false)
%       .H_N_smallest    - Number of MODIS pixels with smallest H to keep
%                           (used when findN_smallest_H is true)
%
%   plot_data - (Optional) Flag or structure for plotting control
%               (currently not used in function but reserved for future plotting)
%
% OUTPUTS:
%   overlap - Structure containing matched pixel indices and distances with fields:
%       All index/distance arrays below are flat vectors of length
%       n_total = n_matches * 10 (or fewer if a MODIS pixel has < 10 EMIT pixels).
%       Every group of 10 consecutive entries corresponds to one MODIS pixel.
%
%       .modis - MODIS pixel information (repeated 10x per match):
%           .row           - Row indices in MODIS array (n_total x 1)
%           .col           - Column indices in MODIS array (n_total x 1)
%           .linear_idx    - Linear indices for MODIS array access (n_total x 1)
%       .emit - EMIT pixel information (10 randomly selected per MODIS pixel):
%           .row           - Row indices in EMIT array (n_total x 1)
%           .col           - Column indices in EMIT array (n_total x 1)
%           .linear_idx    - Linear indices for EMIT array access (n_total x 1)
%       .airs - AIRS pixel information (repeated 10x per match):
%           .row           - Row indices in AIRS array (n_total x 1)
%           .col           - Column indices in AIRS array (n_total x 1)
%           .linear_idx    - Linear indices for AIRS array access (n_total x 1)
%       .amsr - AMSR pixel information (repeated 10x per match):
%           .row           - Row indices in AMSR array (n_total x 1)
%           .col           - Column indices in AMSR array (n_total x 1)
%           .linear_idx    - Linear indices for AMSR array access (n_total x 1)
%       .num_emit_pixels_per_modis - Number of EMIT pixels found in each
%                                    MODIS pixel before random selection (n_matches x 1)
%       .distance_modis_emit_km - Distance between each MODIS-EMIT pair (n_total x 1)
%       .distance_modis_airs_km  - Distance between MODIS and AIRS pixels (n_total x 1)
%       .distance_modis_amsr_km  - Distance between MODIS and AMSR pixels (n_total x 1)
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
%   amsr - Structure containing AMSR-E/2 L2B ocean data
%       (returned from readAMSR_L2_data function)
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
%   - For each MODIS pixel meeting criteria → find ALL EMIT pixels within
%     the MODIS pixel boundary (defined by grid spacing half-widths)
%   - Randomly select 10 EMIT pixels from those found (without replacement)
%   - For each matched MODIS pixel → find closest AIRS pixel
%   - For each matched MODIS pixel → find closest AMSR pixel
%   - MODIS, AIRS, and AMSR indices are repeated 10 times (once per EMIT pixel)
%   - Remove any EMIT pixels with all-NaN radiance (cloud masked by EMIT pipeline)
%   - If all 10 EMIT pixels for a MODIS pixel are masked, remove that MODIS match
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
%   [overlap, emit, modis, airs, amsr, folder_paths] = ...
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
% - AMSR pixel acquisition time (available: amsr.Time_UTC)
% - Then compute: 
%   overlap.time_difference_modis_emit_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))
%   overlap.time_difference_emit_airs_seconds(nn) = abs(emit_time(idx_emit) - airs_time(idx_airs))
%   overlap.time_difference_emit_amsr_seconds(nn) = abs(emit_time(idx_emit) - amsr_time(idx_amsr))
%
% By Andrew John Buggee
%%

function [overlap, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths, criteria, n_emit_per_modis)



%% Load EMIT Data

[emit, folder_paths.L1B_fileName_emit] = retrieveEMIT_data([folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder]);

%% Load Aqua/MODIS Data

[modis, folder_paths.L1B_fileName_modis] = retrieveMODIS_data([folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder]);

%% Load AIRS data

% Find which files are AIRS data
airs = readAIRS_L2_data([folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder]);

%% Load AMSR-E/2 data

amsr = readAMSR_L2_data([folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder]);

%% Get MODIS coordinates
modis_lat = double(modis.geo.lat);
modis_long = double(modis.geo.long);

%% Apply criteria filters first (more efficient than checking spatial overlap first)

% Set default for findN_smallest_H if not provided
if ~isfield(criteria, 'findN_smallest_H')
    criteria.findN_smallest_H = false;
end

% Which cloud phase should we look for?
if strcmp(criteria.cld_phase, 'water')==true
    is_liquidWater = modis.cloud.phase==2;           % 2 is the value designated for liquid water
else
    error('What is the phase??')
end

% Set the optical depth limits
is_tauGreater_and_Less = modis.cloud.optThickness17 >= criteria.cld_tau_min &...
    modis.cloud.optThickness17 <= criteria.cld_tau_max;

% Find pixels where the effective radius retrieval uncertainty is less than 10 percent
is_reUncertLessThan10Percent = modis.cloud.effRad_uncert_17<=10;

% Check if over ocean
isOcean = reshape(land_or_ocean(modis_lat(:), modis_long(:), 20, false), size(modis_lat));


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


%% Combine footprint constraint with criteria

if criteria.findN_smallest_H == false

    % ---- Option 1: Use H threshold ----
    % Keep all pixels with H < criteria.H
    is_H_lessThan = modis.cloud.inhomogeneity_index(:,:,2) < criteria.H;

    idx_criteria = logical(is_liquidWater .* is_tauGreater_and_Less .*...
                           is_reUncertLessThan10Percent .* is_H_lessThan .* isOcean);

    idx_combined_master = logical(idx_in_footprint .* idx_criteria);

    fprintf('Found %d MODIS pixels meeting all criteria (H < %.2f) within EMIT footprint\n',...
        sum(idx_combined_master(:)), criteria.H);

else

    % ---- Option 2: Find N pixels with smallest H ----
    % First, find all pixels meeting phase, tau, re uncertainty, ocean,
    % and footprint criteria (no H threshold applied)
    idx_criteria_noH = logical(is_liquidWater .* is_tauGreater_and_Less .*...
                               is_reUncertLessThan10Percent .* isOcean);

    idx_candidates = logical(idx_in_footprint .* idx_criteria_noH);

    n_candidates = sum(idx_candidates(:));
    fprintf('Found %d MODIS pixels meeting all criteria (no H filter) within EMIT footprint\n', n_candidates);

    % Get the H values for the candidate pixels
    H_values = modis.cloud.inhomogeneity_index(:,:,2);
    H_candidates = H_values(idx_candidates);

    % Sort by H value (ascending) and keep the N smallest
    N = criteria.H_N_smallest;

    if n_candidates <= N
        % Fewer candidates than requested — keep all
        idx_combined_master = idx_candidates;
        fprintf('Keeping all %d candidates (fewer than requested N = %d)\n', n_candidates, N);

    else
        % Find the N smallest H values
        [~, sort_idx] = sort(H_candidates, 'ascend');
        idx_keep = sort_idx(1:N);

        % Convert back to a logical mask over the full MODIS grid
        candidate_linear_idx = find(idx_candidates);
        idx_combined_master = false(size(modis_lat));
        idx_combined_master(candidate_linear_idx(idx_keep)) = true;

        fprintf('Selected %d MODIS pixels with smallest H values (H range: %.4f to %.4f)\n',...
            N, H_candidates(sort_idx(1)), H_candidates(sort_idx(N)));

    end

end

%% Find row/column indices for MODIS pixels and match with closest EMIT, AIRS, and AMSR pixels

n_matches = sum(idx_combined_master(:));

% Number of EMIT pixels to randomly select per MODIS pixel
% n_emit_per_modis = 30;

% Maximum total entries (may be fewer if some MODIS pixels have < 10 EMIT pixels)
n_total_max = n_matches * n_emit_per_modis;

% Pre-allocate output arrays for MODIS (repeated n_emit_per_modis times per match)
overlap.modis.row = zeros(n_total_max, 1);
overlap.modis.col = zeros(n_total_max, 1);
overlap.modis.linear_idx = zeros(n_total_max, 1);

% Pre-allocate output arrays for EMIT (n_emit_per_modis per MODIS pixel)
overlap.emit.row = zeros(n_total_max, 1);
overlap.emit.col = zeros(n_total_max, 1);
overlap.emit.linear_idx = zeros(n_total_max, 1);

% Pre-allocate output arrays for AIRS (repeated n_emit_per_modis times per match)
overlap.airs.row = zeros(n_total_max, 1);
overlap.airs.col = zeros(n_total_max, 1);
overlap.airs.linear_idx = zeros(n_total_max, 1);

% Pre-allocate output arrays for AMSR (repeated n_emit_per_modis times per match)
overlap.amsr.row = zeros(n_total_max, 1);
overlap.amsr.col = zeros(n_total_max, 1);
overlap.amsr.linear_idx = zeros(n_total_max, 1);

% Pre-allocate distance arrays
overlap.distance_modis_emit_km = zeros(n_total_max, 1);
overlap.distance_modis_airs_km = zeros(n_total_max, 1);
overlap.distance_modis_amsr_km = zeros(n_total_max, 1);

% Track how many EMIT pixels were found in each MODIS pixel (before selection)
overlap.num_emit_pixels_per_modis = zeros(n_matches, 1);

% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - AIRS pixel acquisition time (available: airs.Time_UTC)
% - AMSR pixel acquisition time (available: amsr.Time_UTC)
% - Then compute:
%   overlap.time_difference_modis_emit_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))
%   overlap.time_difference_emit_airs_seconds(nn) = abs(emit_time(idx_emit) - airs_time(idx_airs))
%   overlap.time_difference_emit_amsr_seconds(nn) = abs(emit_time(idx_emit) - amsr_time(idx_amsr))

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

% Flatten AMSR coordinates for vectorized distance calculation
amsr_lat = double(amsr.geo.Latitude);
amsr_long = double(amsr.geo.Longitude);
amsr_lat_flat = amsr_lat(:);
amsr_long_flat = amsr_long(:);

%% Compute MODIS pixel half-widths from grid spacing
% Use the spacing between neighboring pixels to define each MODIS pixel's
% boundary. For interior pixels, use the average of the forward and backward
% differences. For edge pixels, use one-sided differences.

% % Half-width in latitude direction (row spacing)
% dlat = diff(modis_lat, 1, 1);  % difference along rows
% modis_half_dlat = zeros(size(modis_lat));
% modis_half_dlat(1,:) = abs(dlat(1,:)) / 2;
% modis_half_dlat(end,:) = abs(dlat(end,:)) / 2;
% modis_half_dlat(2:end-1,:) = (abs(dlat(1:end-1,:)) + abs(dlat(2:end,:))) / 4;

% Third-width in latitude direction (row spacing)
dlat = diff(modis_lat, 1, 1);  % difference along rows
modis_half_dlat = zeros(size(modis_lat));
modis_half_dlat(1,:) = abs(dlat(1,:)) / 3;
modis_half_dlat(end,:) = abs(dlat(end,:)) / 3;
modis_half_dlat(2:end-1,:) = (abs(dlat(1:end-1,:)) + abs(dlat(2:end,:))) / 6;

% % Half-width in longitude direction (column spacing)
% dlon = diff(modis_long, 1, 2);  % difference along columns
% modis_half_dlon = zeros(size(modis_long));
% modis_half_dlon(:,1) = abs(dlon(:,1)) / 2;
% modis_half_dlon(:,end) = abs(dlon(:,end)) / 2;
% modis_half_dlon(:,2:end-1) = (abs(dlon(:,1:end-1)) + abs(dlon(:,2:end))) / 4;

% Half-width in longitude direction (column spacing)
dlon = diff(modis_long, 1, 2);  % difference along columns
modis_half_dlon = zeros(size(modis_long));
modis_half_dlon(:,1) = abs(dlon(:,1)) / 3;
modis_half_dlon(:,end) = abs(dlon(:,end)) / 3;
modis_half_dlon(:,2:end-1) = (abs(dlon(:,1:end-1)) + abs(dlon(:,2:end))) / 6;



% Get linear indices of matching MODIS pixels
modis_linear_idx = find(idx_combined_master);

% Running index into the flat output arrays
idx_out = 0;

for nn = 1:n_matches

    % ====== MODIS Pixel Information ======
    % Get the linear index for this MODIS pixel
    modis_lidx = modis_linear_idx(nn);

    % Convert linear index to row/col for MODIS
    [modis_row, modis_col] = ind2sub(size(modis_lat), modis_lidx);

    % Get lat/lon for this MODIS pixel
    modis_lat_current = modis_lat(modis_lidx);
    modis_long_current = modis_long(modis_lidx);

    % Get the half-widths for this MODIS pixel
    half_dlat = modis_half_dlat(modis_lidx);
    half_dlon = modis_half_dlon(modis_lidx);


    % ====== Find ALL EMIT Pixels Within This MODIS Pixel ======
    % Define the MODIS pixel boundary
    lat_min = modis_lat_current - half_dlat;
    lat_max = modis_lat_current + half_dlat;
    lon_min = modis_long_current - half_dlon;
    lon_max = modis_long_current + half_dlon;

    % Find all EMIT pixels whose centers fall within the MODIS pixel boundary
    idx_emit_all = find(emit_lat_flat >= lat_min & emit_lat_flat <= lat_max & ...
                        emit_long_flat >= lon_min & emit_long_flat <= lon_max);

    n_found = numel(idx_emit_all);
    overlap.num_emit_pixels_per_modis(nn) = n_found;

    % ====== Randomly Select 10 EMIT Pixels ======
    if n_found >= n_emit_per_modis
        % Randomly select n_emit_per_modis without replacement
        rand_idx = randperm(n_found, n_emit_per_modis);
        idx_emit_selected = idx_emit_all(rand_idx);
    else
        % Fewer than 10 EMIT pixels found -- keep all and warn
        warning(['MODIS pixel %d (row=%d, col=%d) has only %d EMIT pixels ' ...
                 'within its boundary (expected >= %d). Keeping all available.'], ...
                 nn, modis_row, modis_col, n_found, n_emit_per_modis);
        idx_emit_selected = idx_emit_all;
    end

    n_selected = numel(idx_emit_selected);

    % Convert selected EMIT linear indices to row/col
    [rows_emit, cols_emit] = ind2sub(size(emit_lat), idx_emit_selected);

    % Calculate accurate distances between MODIS and each selected EMIT pixel
    dist_km = zeros(n_selected, 1);
    for ee = 1:n_selected
        dist_km(ee) = distance(modis_lat_current, modis_long_current, ...
            emit_lat_flat(idx_emit_selected(ee)), emit_long_flat(idx_emit_selected(ee)), wgs84);
    end


    % ====== Find Closest AIRS Pixel ======
    % Use the MODIS pixel location to find the closest AIRS pixel
    dist_sq_airs = (airs_lat_flat - modis_lat_current).^2 + (airs_long_flat - modis_long_current).^2;
    [~, idx_airs] = min(dist_sq_airs);

    % Convert to row/col for AIRS
    [airs_row, airs_col] = ind2sub(size(airs_lat), idx_airs);

    % Calculate accurate distance between MODIS and AIRS using WGS84 ellipsoid
    dist_modis_airs = distance(modis_lat_current, modis_long_current, ...
        airs_lat_flat(idx_airs), airs_long_flat(idx_airs), wgs84);


    % ====== Find Closest AMSR Pixel ======
    % Use the MODIS pixel location to find the closest AMSR pixel
    dist_sq_amsr = (amsr_lat_flat - modis_lat_current).^2 + (amsr_long_flat - modis_long_current).^2;
    [~, idx_amsr] = min(dist_sq_amsr);

    % Convert to row/col for AMSR
    [amsr_row, amsr_col] = ind2sub(size(amsr_lat), idx_amsr);

    % Calculate accurate distance between MODIS and AMSR using WGS84 ellipsoid
    dist_modis_amsr = distance(modis_lat_current, modis_long_current, ...
        amsr_lat_flat(idx_amsr), amsr_long_flat(idx_amsr), wgs84);


    % ====== Store Results: one row per selected EMIT pixel ======
    out_range = idx_out + (1:n_selected);

    % MODIS info (repeated for each EMIT pixel)
    overlap.modis.row(out_range) = modis_row;
    overlap.modis.col(out_range) = modis_col;
    overlap.modis.linear_idx(out_range) = modis_lidx;

    % EMIT info (unique per row)
    overlap.emit.row(out_range) = rows_emit;
    overlap.emit.col(out_range) = cols_emit;
    overlap.emit.linear_idx(out_range) = idx_emit_selected;

    % AIRS info (repeated for each EMIT pixel)
    overlap.airs.row(out_range) = airs_row;
    overlap.airs.col(out_range) = airs_col;
    overlap.airs.linear_idx(out_range) = idx_airs;

    % AMSR info (repeated for each EMIT pixel)
    overlap.amsr.row(out_range) = amsr_row;
    overlap.amsr.col(out_range) = amsr_col;
    overlap.amsr.linear_idx(out_range) = idx_amsr;

    % Distances
    overlap.distance_modis_emit_km(out_range) = dist_km;
    overlap.distance_modis_airs_km(out_range) = dist_modis_airs;
    overlap.distance_modis_amsr_km(out_range) = dist_modis_amsr;

    idx_out = idx_out + n_selected;

end

% Trim pre-allocated arrays in case some MODIS pixels had < 10 EMIT pixels
n_total = idx_out;
if n_total < n_total_max
    overlap.modis.row = overlap.modis.row(1:n_total);
    overlap.modis.col = overlap.modis.col(1:n_total);
    overlap.modis.linear_idx = overlap.modis.linear_idx(1:n_total);
    overlap.emit.row = overlap.emit.row(1:n_total);
    overlap.emit.col = overlap.emit.col(1:n_total);
    overlap.emit.linear_idx = overlap.emit.linear_idx(1:n_total);
    overlap.airs.row = overlap.airs.row(1:n_total);
    overlap.airs.col = overlap.airs.col(1:n_total);
    overlap.airs.linear_idx = overlap.airs.linear_idx(1:n_total);
    overlap.amsr.row = overlap.amsr.row(1:n_total);
    overlap.amsr.col = overlap.amsr.col(1:n_total);
    overlap.amsr.linear_idx = overlap.amsr.linear_idx(1:n_total);
    overlap.distance_modis_emit_km = overlap.distance_modis_emit_km(1:n_total);
    overlap.distance_modis_airs_km = overlap.distance_modis_airs_km(1:n_total);
    overlap.distance_modis_amsr_km = overlap.distance_modis_amsr_km(1:n_total);
end

%% Remove EMIT pixels that have been masked out (all-NaN radiance)
% The EMIT cloud mask sets radiance to NaN at every wavelength for masked pixels.
% Check each selected EMIT pixel and remove it if all wavelengths are NaN.

% emit.radiance.measurements is (rows x cols x wavelengths)
is_masked = false(n_total, 1);
for pp = 1:n_total
    emit_spectrum = emit.radiance.measurements(overlap.emit.row(pp), overlap.emit.col(pp), :);
    if all(isnan(emit_spectrum))
        is_masked(pp) = true;
    end
end

n_masked = sum(is_masked);
fprintf('\nFound %d EMIT pixels with all-NaN radiance (cloud masked). Removing...\n', n_masked);

% Remove masked EMIT pixels from all overlap arrays
if n_masked > 0
    idx_keep = ~is_masked;

    overlap.modis.row = overlap.modis.row(idx_keep);
    overlap.modis.col = overlap.modis.col(idx_keep);
    overlap.modis.linear_idx = overlap.modis.linear_idx(idx_keep);
    overlap.emit.row = overlap.emit.row(idx_keep);
    overlap.emit.col = overlap.emit.col(idx_keep);
    overlap.emit.linear_idx = overlap.emit.linear_idx(idx_keep);
    overlap.airs.row = overlap.airs.row(idx_keep);
    overlap.airs.col = overlap.airs.col(idx_keep);
    overlap.airs.linear_idx = overlap.airs.linear_idx(idx_keep);
    overlap.amsr.row = overlap.amsr.row(idx_keep);
    overlap.amsr.col = overlap.amsr.col(idx_keep);
    overlap.amsr.linear_idx = overlap.amsr.linear_idx(idx_keep);
    overlap.distance_modis_emit_km = overlap.distance_modis_emit_km(idx_keep);
    overlap.distance_modis_airs_km = overlap.distance_modis_airs_km(idx_keep);
    overlap.distance_modis_amsr_km = overlap.distance_modis_amsr_km(idx_keep);

    % Check if any MODIS pixel lost ALL of its EMIT pixels
    % If so, remove those MODIS (and corresponding AIRS/AMSR) entries entirely
    unique_modis_remaining = unique(overlap.modis.linear_idx);
    unique_modis_original = unique(modis_linear_idx);
    modis_lost = setdiff(unique_modis_original, unique_modis_remaining);

    if ~isempty(modis_lost)
        fprintf('Removed %d MODIS pixels entirely (all EMIT pixels were masked).\n', numel(modis_lost));
    end

    n_total = numel(overlap.emit.linear_idx);
end

% Summary statistics
fprintf('\n--- Overlap Summary ---\n');
fprintf('MODIS pixels meeting all criteria: %d\n', n_matches);
fprintf('EMIT pixels selected per MODIS pixel: %d\n', n_emit_per_modis);
fprintf('EMIT pixels removed (cloud masked): %d\n', n_masked);
fprintf('Total EMIT pixels stored after masking: %d\n', n_total);
fprintf('Unique MODIS pixels remaining: %d\n', numel(unique(overlap.modis.linear_idx)));
fprintf('Average EMIT pixels found per MODIS pixel (before selection): %.1f\n', mean(overlap.num_emit_pixels_per_modis));
fprintf('Min/Max EMIT pixels found per MODIS pixel: %d / %d\n', min(overlap.num_emit_pixels_per_modis), max(overlap.num_emit_pixels_per_modis));

if n_total > 0
    fprintf('\nAverage distance between matched MODIS-EMIT pixels: %.2f km\n', mean(overlap.distance_modis_emit_km));
    fprintf('Median distance: %.2f km\n', median(overlap.distance_modis_emit_km));
    fprintf('Max distance: %.2f km\n', max(overlap.distance_modis_emit_km));

    fprintf('\nAverage distance between matched MODIS-AIRS pixels: %.2f km\n', mean(overlap.distance_modis_airs_km));
    fprintf('Median distance: %.2f km\n', median(overlap.distance_modis_airs_km));
    fprintf('Max distance: %.2f km\n', max(overlap.distance_modis_airs_km));

    fprintf('\nAverage distance between matched MODIS-AMSR pixels: %.2f km\n', mean(overlap.distance_modis_amsr_km));
    fprintf('Median distance: %.2f km\n', median(overlap.distance_modis_amsr_km));
    fprintf('Max distance: %.2f km\n', max(overlap.distance_modis_amsr_km));
else
    fprintf('\nNo EMIT pixels remaining after cloud mask removal.\n');
end

end