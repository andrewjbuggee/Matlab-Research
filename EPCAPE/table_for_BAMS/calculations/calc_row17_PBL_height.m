function results = calc_row17_PBL_height(arm_folder, low_cloud_5min, low_cloud_max_m)
% CALC_ROW17_PBL_HEIGHT  Sheet row 17: boundary-layer height, "PBLH-Thermo (Zhang)".
%
%   results = calc_row17_PBL_height(arm_folder, low_cloud_5min, low_cloud_max_m)
%
%   DATA: ARM's PBLHTTHERMO VAP at the Scripps Pier (epcpblhtthermoM1.c1),
%   file pblh_thermo_M1_20230215_20240214.nc. Variable pbl_height (km above
%   ground), the PBL height from combined Raman-lidar and TROPoe
%   potential-temperature profiles (ARM Data Discovery: "Planetary Boundary
%   Layer Height Best Estimate using Machine Learning"). One value every
%   10 min (:00, :10, ..., :50); no QC field; no data 11-31 Oct 2023.
%
%   TWO NEW LINES (m above ground, sheet seasons, each 10-min value counts once):
%   Line 1: only while KAZR-ARSCL sees a low cloud: the 5-min interval
%           nearest the PBL time stamp holds a low-cloud sample (from
%           calc_row15_cloud_top). This reading reproduces the sheet best
%           of the ones tried; Kavin's method was not shared.
%   Line 2: every 10-min value (all times).
%
%   low_cloud_5min  start times of 5-min intervals with a low cloud (row 15)

file = fullfile(arm_folder, 'pblh_thermo_M1_20230215_20240214.nc');
t = read_netcdf_time(file);
pbl_km = double(ncread(file, 'pbl_height'));
units  = strtrim(string(ncreadatt(file, 'pbl_height', 'units')));
if units ~= "km"
    error('calc_row17_PBL_height: expected pbl_height in km, found "%s"', units);
end
pbl_m = 1000 * pbl_km;                                    % km -> m

% Line 1: is there a low cloud in the 5-min interval nearest each PBL time stamp?
% (PBL time stamps fall on whole 10-min marks, so the nearest 5-min interval
%  is the one that STARTS at the time stamp.)
nearest_5min = datetime(round(posixtime(t) / 300) * 300, 'ConvertFrom', 'posixtime');
cloudy = ismember(nearest_5min, low_cloud_5min);
pbl_cloudy = pbl_m;
pbl_cloudy(~cloudy) = NaN;
stats1 = season_stats(pbl_cloudy, sheet_season_masks(t));
results = make_result(17, "Boundary Layer Height", "m", "Pier", ...
    "epcpblhtthermoM1.c1 (PBLHTTHERMO) + KAZR-ARSCL", stats1, ...
    sprintf("filters: pbl_height (10-min values; no QC field) only when KAZR-ARSCL has a cloud top <= %g m ", low_cloud_max_m) + ...
    "within ±5 min; the best match among six readings; sheet seasons (UTC).");

% Line 2: all times
stats2 = season_stats(pbl_m, sheet_season_masks(t));
results(2) = make_result(17, "Boundary Layer Height", "m", "Pier", "epcpblhtthermoM1.c1 (PBLHTTHERMO)", stats2, ...
    "filters: none (pbl_height has no QC field); every 10-min value; sheet seasons (UTC).");

end
