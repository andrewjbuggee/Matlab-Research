function [results, low_cloud_5min] = calc_row15_cloud_top(arm_folder, low_cloud_max_m)
% CALC_ROW15_CLOUD_TOP  Sheet row 15: cloud-top height from KAZR-ARSCL.
%
%   [results, low_cloud_5min] = calc_row15_cloud_top(arm_folder, low_cloud_max_m)
%
%   DATA: KAZR-ARSCL at the Scripps Pier (epcarsclkazr1kolliasM1.c1), file
%   cloud_arscl_M1_20230215_20240214.nc. ARSCL merges the Ka-band radar
%   (KAZR, which sets cloud TOPS), the micropulse lidar and the ceilometer
%   into up to 10 hydrometeor layers every 4 s:
%       cloud_layer_top_height(layer, time), cloud_layer_base_height(layer, time)
%   in m above ground; layer 1 = the lowest. Empty layer slots hold -1.
%   ARSCL has no QC field for these variables.
%
%   WHAT IS COMPUTED: the top of the LOWEST layer, kept when
%     * the layer exists (its base and top are real heights, >= 0), and
%     * its top is above ground and <= low_cloud_max_m (3500 m), i.e. a
%       low cloud (the marine stratocumulus deck, not cirrus).
%   Every 4-s sample that passes counts once (no time averaging).
%
%   SECOND OUTPUT, used for the PBL height (row 17):
%   low_cloud_5min = start times (datetime) of every 5-min interval that
%   holds at least one low-cloud sample.

file = fullfile(arm_folder, 'cloud_arscl_M1_20230215_20240214.nc');
t = read_netcdf_time(file);

% Read only layer 1. In MATLAB the dimensions are (layer, time), so
% start = [1 1] and count = [1 Inf] = "layer 1, every time".
top0  = double(ncread(file, 'cloud_layer_top_height',  [1 1], [1 Inf])).';   % m above ground
base0 = double(ncread(file, 'cloud_layer_base_height', [1 1], [1 Inf])).';

layer_exists = isfinite(top0) & top0 >= 0 & isfinite(base0) & base0 >= 0;   % -1 = empty slot
is_low_cloud = layer_exists & top0 > 0 & top0 <= low_cloud_max_m;

cloud_top_m = top0;
cloud_top_m(~is_low_cloud) = NaN;

stats = season_stats(cloud_top_m, sheet_season_masks(t));
results = make_result(15, "Cloud Top Height (low clouds only)", "m", "Pier", ...
    "epcarsclkazr1kolliasM1.c1 (KAZR-ARSCL)", stats, ...
    sprintf("filters: top of the lowest cloud layer <= %g m AGL; cloud_source_flag codes any; ", low_cloud_max_m) + ...
    "averaging window 0 min (0 = every 4-s profile); sheet seasons (UTC)");

% 5-min intervals that contain at least one low-cloud sample (for row 17)
low_cloud_5min = unique(floor_time(t(is_low_cloud), 5));

end
