function results = calc_rows22to24_rain(arm_folder, rain_day_min_mm)
% CALC_ROWS22TO24_RAIN  Sheet rows 22-24: "Rain Accumulation" at three disdrometers.
%
%   results = calc_rows22to24_rain(arm_folder, rain_day_min_mm)
%
%   DATA (1-min rain rates, mm/h):
%     row 22  LD M1   Parsivel2 laser disdrometer, Scripps Pier (epcldM1.b1),
%                     precip_rate with qc_precip_rate = 0
%     row 23  LD S2   Parsivel2 at Mt. Soledad (epcldS2.b1; the sheet's "LDS1":
%                     there is no S1 facility), precip_rate with QC = 0
%     row 24  VDIS M1 2D video disdrometer VAP, Scripps Pier
%                     (epcvdisquantsM1.c1), rain_rate. No QC field. The VAP
%                     writes NaN (not 0) when it sees no drops, so NaN is
%                     counted as 0 (no rain).
%
%   WHAT IS COMPUTED (Kavin did not document his method; of seven readings
%   this one matches the sheet best, and its numbers are mm/DAY, not mm/hr):
%     1. daily (UTC) rain accumulation = sum of rate x (sample length in h);
%        a day counts only if >= 80 % of its 1-min samples are valid, so an
%        outage does not look like a dry day;
%     2. rain days = days with >= rain_day_min_mm (0.1 mm);
%     3. the mean over the rain days of each calendar month;
%     4. value = mean ± std ACROSS those monthly values in each season. Each
%        month is placed in the season of its first rain day (sheet seasons).

specs = { ...   % sheet row, file, variable, location, product text, QC text
    22, 'rain_ld_M1_20230215_20240214.nc',   'precip_rate', "Pier", ...
        "epcldM1.b1 (Parsivel LD)", "QC = 0 on precip_rate", "Rain Accumulation - LDM1"; ...
    23, 'rain_ld_S2_20230215_20240214.nc',   'precip_rate', "Soledad", ...
        "epcldS2.b1 (Parsivel LD at Mt. Soledad; the sheet's 'LDS1')", "QC = 0 on precip_rate", "Rain Accumulation - LDS1"; ...
    24, 'rain_vdis_M1_20230215_20240214.nc', 'rain_rate',   "Pier", ...
        "epcvdisquantsM1.c1 (VDISQUANTS)", "rain_rate has no QC field; NaN (no drops) counted as 0", "Rain Accumulation - VDIS"};

results = [];
for i = 1:size(specs, 1)
    [row, fname, var, location, product, qc_text, quantity] = specs{i, :};
    file = fullfile(arm_folder, fname);
    t    = read_netcdf_time(file);
    rate = read_arm_qc0(file, var);                        % mm/h (QC = 0 where a qc_ field exists)
    if row == 24
        rate(isnan(rate)) = 0;                             % VDIS: NaN = no drops seen = no rain
    end

    % --- 1. daily accumulation (mm) and how complete each day is ------------------------
    dt_h = median(diff(posixtime(t))) / 3600;              % sample length in hours (1 min)
    day  = dateshift(t, 'start', 'day');                   % UTC day of each sample
    [days_u, ~, which_day] = unique(day);
    valid = isfinite(rate);
    total_mm  = accumarray(which_day(valid), rate(valid) * dt_h, [numel(days_u) 1]);
    n_valid   = accumarray(which_day(valid), 1,                  [numel(days_u) 1]);
    coverage  = n_valid / (24 / dt_h);                     % fraction of the day's samples that are valid

    % --- 2. rain days ---------------------------------------------------------------
    is_rain_day = coverage >= 0.8 & n_valid > 0 & total_mm >= rain_day_min_mm;
    rain_days   = days_u(is_rain_day);
    rain_mm     = total_mm(is_rain_day);

    % --- 3. monthly mean of the rain days; each month dated by its first rain day ----------
    month_id = year(rain_days) * 100 + month(rain_days);   % e.g. 202307
    [~, first_day_idx, which_month] = unique(month_id, 'first');
    monthly_mean_mm = accumarray(which_month, rain_mm, [], @mean);
    month_time      = rain_days(first_day_idx);

    % --- 4. mean ± std across months, per season ----------------------------------------
    stats = season_stats(monthly_mean_mm, sheet_season_masks(month_time));
    r = make_result(row, quantity, "mm/day", location, product, stats, ...
        "filters: " + qc_text + sprintf("; daily accumulation on rain days (>= %g mm), averaged over each ", rain_day_min_mm) + ...
        "calendar month; value = mean ± std across months; sheet seasons (UTC); the best of seven readings.");
    results = [results, r]; %#ok<AGROW>
end

end
