function results = calc_row08_LW_down(arm_folder)
% CALC_ROW08_LW_DOWN  Sheet row 8: downwelling longwave (LW) irradiance.
%
%   results = calc_row08_LW_down(arm_folder)
%
%   DATA: RADFLUX at the Scripps Pier (ARM datastream epcradflux1longM1.c2),
%   file radflux_M1_20230215_20240214.nc. Variable downwelling_longwave
%   (W m-2): 1-min means of the 1-s SKYRAD pyrgeometer samples, time-stamped
%   at the END of each minute. Only samples with qc_downwelling_longwave = 0.
%   (The sheet says "From Lubin, Soledad", but Lubin's LW tables are this
%   Pier record: his monthly means agree with it within 0.4 W m-2.)
%
%   TWO NEW LINES:
%   Line 1, Kavin's statistic applied to the full record:
%     1. shift the time stamps 1 min back (to the START of each minute),
%        keep 2023-02-15 .. 2024-02-14;
%     2. hourly means, with hours in Pacific Standard Time (UTC - 8 h), as in
%        Lubin's tables;
%     3. for each calendar month, the MEDIAN of each hour of day across the
%        month's days -> a 24-hour "median diurnal cycle" per month;
%     4. season = mean over its 3 months of those cycles (24 values);
%        the value = mean of the 24 values, ± = their std (the spread of the
%        diurnal cycle, NOT day-to-day variability). EPCAPE = mean ± std of
%        all 12 x 24 monthly medians. Calendar quarters.
%   Line 2: every QC = 0 1-min value, day and night; mean ± std; sheet seasons.

file = fullfile(arm_folder, 'radflux_M1_20230215_20240214.nc');
t  = read_netcdf_time(file);                          % end of each minute (UTC)
lw = read_arm_qc0(file, 'downwelling_longwave');      % W m-2, NaN where QC ~= 0

% ======================= Line 1: Kavin's statistic ===================================
t_start  = t - minutes(1);                                      % start of each minute
keep     = isfinite(lw) & t_start >= datetime(2023, 2, 15) & t_start < datetime(2024, 2, 15);
t_local  = t_start(keep) - hours(8);                            % Pacific Standard Time
lw_keep  = lw(keep);

% 2. hourly means: every sample goes to the hour that contains it
hour_start = dateshift(t_local, 'start', 'hour');
[hours_u, ~, which_hour] = unique(hour_start);
hourly_mean = accumarray(which_hour, lw_keep, [], @mean);       % one mean per hour

% 3. median of each (hour of day, month) over the month's days
hod = hour(hours_u) + 1;                                        % hour of day 1..24 (MATLAB index)
mon = month(hours_u);                                           % calendar month 1..12
medians = accumarray([hod mon], hourly_mean, [24 12], @median, NaN);  % 24 x 12 matrix

% 4. season statistics
season_months = {1:12, [4 5 6], [7 8 9], [10 11 12], [1 2 3]};
stats1.mean = nan(1, 5); stats1.std = nan(1, 5); stats1.n = zeros(1, 5);
for k = 1:5
    if k == 1
        v = medians(:);                                         % all 12 x 24 medians
    else
        v = mean(medians(:, season_months{k}), 2, 'omitnan');   % the season's 24-hour cycle
    end
    v = v(isfinite(v));
    stats1.mean(k) = mean(v);
    stats1.std(k)  = std(v);
    stats1.n(k)    = numel(v);
end
results = make_result(8, "Longwave (LW) Flux", "W/m2", "Pier", ...
    "epcradflux1longM1.c2 (RADFLUX; the record Lubin's tables are made from)", stats1, ...
    "Kavin's method on the full record. filters: QC = 0 on downwelling_longwave; hourly means (PST); " + ...
    "each month's median of each hour; season = mean of its 3 months' median diurnal cycles over 24 hours, " + ...
    "± = spread of that cycle; calendar quarters.");

% ======================= Line 2: every QC = 0 sample ==================================
stats2 = season_stats(lw, sheet_season_masks(t));
results(2) = make_result(8, "Longwave (LW) Flux", "W/m2", "Pier", "epcradflux1longM1.c2 (RADFLUX)", stats2, ...
    "filters: QC = 0 on downwelling_longwave; every QC = 0 sample (no further averaging); day and night; " + ...
    "± = std of all values; sheet seasons (UTC).");

end
