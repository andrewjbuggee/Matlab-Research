function results = calc_row20_LCL(arm_folder)
% CALC_ROW20_LCL  Sheet row 20: lifting condensation level (LCL) from radiosondes.
%
%   results = calc_row20_LCL(arm_folder)
%
%   DATA: ARM's SONDEPARAM VAP at the Scripps Pier (epcsondeparamM1.c1), in
%   the folder sondeparam_M1/: ONE FILE PER RADIOSONDE LAUNCH (~3-4 a day).
%   In each file:
%     lcl(parcel_type, time)  LCL in km for three parcels (parcel_type 1 =
%                             surface, 2 = most unstable, 3 = mixed layer);
%                             the same value is repeated at every time step
%                             of the launch; -9999 = missing
%     data_quality(time)      bit mask of input checks on the profile
%                             (0 = every check passed). There is no qc_lcl,
%                             so data_quality = 0 for the whole launch is
%                             used as the "QC = 0" test.
%
%   WHAT IS COMPUTED: one value per launch = the surface-parcel LCL (km -> m),
%   kept only when data_quality = 0 throughout the launch; mean ± std over the
%   launches of each season (sheet seasons; launch time = first time step).
%   (The VAP documentation does not say whether lcl is above ground or sea
%   level; the site is ~7 m above sea level, so it hardly matters.)

files = dir(fullfile(arm_folder, 'sondeparam_M1', 'epcsondeparamM1.c1.*.nc'));
files = sort({files.name});
n_launch = numel(files);

launch_time = NaT(n_launch, 1);
lcl_m       = nan(n_launch, 1);

for i = 1:n_launch
    path = fullfile(arm_folder, 'sondeparam_M1', files{i});
    t = read_netcdf_time(path);
    launch_time(i) = t(1);

    parcel_type = double(ncread(path, 'parcel_type'));          % e.g. [1; 2; 3]
    lcl_km      = double(ncread(path, 'lcl'));                  % (parcel_type, time)
    lcl_km(lcl_km < -9000) = NaN;                               % -9999 = missing
    dq          = double(ncread(path, 'data_quality'));         % (time)

    surface = lcl_km(parcel_type == 1, :);                      % surface parcel, every time step
    surface = surface(isfinite(surface));
    passed_checks = all(isfinite(dq)) && all(dq == 0);
    if ~isempty(surface) && passed_checks
        lcl_m(i) = 1000 * surface(1);                           % km -> m
    end
end

% Only launches inside the campaign (the season masks do this as well)
stats = season_stats(lcl_m, sheet_season_masks(launch_time));
results = make_result(20, "Lifting Condensation Level", "m", "Pier", "epcsondeparamM1.c1 (SONDEPARAM)", stats, ...
    "filters: surface-parcel LCL (parcel type 1), data_quality = 0 (every input check passed); one value " + ...
    "per radiosonde launch; km -> m; sheet seasons (UTC).");

end
