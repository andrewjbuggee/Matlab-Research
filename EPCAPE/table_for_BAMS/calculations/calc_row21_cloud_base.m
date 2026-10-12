function results = calc_row21_cloud_base(arm_folder, low_cloud_max_m)
% CALC_ROW21_CLOUD_BASE  Sheet row 21: cloud-base height from the ceilometer.
%
%   results = calc_row21_cloud_base(arm_folder, low_cloud_max_m)
%
%   DATA: Vaisala ceilometer at the Scripps Pier (epcceilM1.b1), file
%   cbh_ceil_M1_20230215_20240214.nc, one sample every 16 s:
%     first_cbh         lowest cloud base (m above the instrument)
%     qc_first_cbh      0 = passed every QC test
%     detection_status  0 clear; 1, 2, 3 = one, two, three cloud bases found;
%                       4 = obscured with no base (fog, heavy rain)
%
%   WHAT IS COMPUTED: first_cbh where qc_first_cbh = 0, detection_status is
%   1-3 (a cloud base was found) and the base is <= low_cloud_max_m
%   (3500 m, low cloud only). Every 16-s sample counts once; sheet seasons.

file = fullfile(arm_folder, 'cbh_ceil_M1_20230215_20240214.nc');
t = read_netcdf_time(file);
cbh_m  = read_arm_qc0(file, 'first_cbh');                  % m, NaN where QC ~= 0
status = double(ncread(file, 'detection_status'));

keep = ismember(status, [1 2 3]) & isfinite(cbh_m) & cbh_m <= low_cloud_max_m;
cbh_m(~keep) = NaN;

stats = season_stats(cbh_m, sheet_season_masks(t));
results = make_result(21, "Cloud Base Height", "m", "Pier", "epcceilM1.b1 (Vaisala ceilometer)", stats, ...
    sprintf("filters: qc_first_cbh = 0; detection_status 1-3 (a cloud base detected); first_cbh <= %g m ", low_cloud_max_m) + ...
    "above the instrument; averaging window 0 min (0 = every 16-s sample); sheet seasons (UTC).");

end
