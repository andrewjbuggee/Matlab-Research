function results = calc_rows31_32_FM120_hours(fm, vis_files, settings)
% CALC_ROWS31_32_FM120_HOURS  Rows 31-32: hours of on-site cloud and haze sampling.
%
%   results = calc_rows31_32_FM120_hours(fm, vis_files, settings)
%     fm         read_fm120_daily_files(...)   FM-120 5-min rows (time, lwc_gm3, nd_cm3)
%     vis_files  {vis1, vis2}: the two Mt. Soledad visibility files, each from
%                read_visibility_file(...), in time order
%     settings   thresholds (see make_table_for_BAMS.m)
%
%   Each FM-120 row is a 5-min interval. It is sorted into "cloud" or "haze"
%   using the Mt. Soledad visibility and the FM-120 liquid water content
%   (LWC); the hours are (number of rows) x 5 min.
%
%   TWO NEW LINES PER ROW:
%   Line 1, Kavin's code re-run (not truncated; truncating each season to
%   whole hours gives the sheet's numbers):
%     * visibility of row i = MEAN of the raw ~1-s samples from row i's time
%       to row i+1's time, both ends included (so across an FM-120 gap the
%       window spans the whole gap); zeros and spikes are NOT removed;
%     * rows before the 2nd visibility file starts use file 1; the row at the
%       start of file 2 (2023-07-19 04:00) gets 0 m; later rows use file 2;
%       the last row keeps 0 m (these reproduce his code exactly);
%     * cloud = visibility < 1 km and LWC > 0.01 g m-3;
%       haze  = 1 km < visibility < 5 km, NO LWC condition;
%     * calendar quarters (UTC); EPCAPE = every row.
%   Line 2, the BAMS-paper definitions (Lynn Russell's EPCAPE BAMS paper):
%     * visibility of a row = MEDIAN of the screened 1-s samples inside its
%       5-min interval [t, t + 5 min); screening removes samples <= 0 m and
%       > 48,270 m (the sensor's 30-mile maximum); an interval needs at least
%       settings.vis_min_samples (125) valid samples, of ~250;
%     * cloud = visibility < 1 km and LWC > 0.01 g m-3;
%       haze  = visibility < 5 km and LWC < 0.01 g m-3 (and not cloud);
%     * the FM-120 must have reported N and LWC; sheet seasons (UTC).

t_fm = fm.time;
lwc  = fm.lwc_gm3;
n_rows = numel(t_fm);

% =================== Line 1: Kavin's code =================================================
avgvis = zeros(n_rows, 1);                                  % his avgvis = zeros(...)
t2_start = floor_time(vis_files{2}.time(1), 5);             % first 5-min mark of file 2
split = find(t_fm >= t2_start, 1, 'first');                 % first FM-120 row on or after it
rows_1 = 1 : min(split - 1, n_rows - 1);                    % rows that read file 1
rows_2 = split + 1 : n_rows - 1;                            % rows that read file 2
avgvis(rows_1) = kavin_window_mean(vis_files{1}, t_fm(rows_1), t_fm(rows_1 + 1));
avgvis(rows_2) = kavin_window_mean(vis_files{2}, t_fm(rows_2), t_fm(rows_2 + 1));
avgvis(split)  = 0;                                         % his avgvis(36952) = 0
avgvis(n_rows) = 0;                                         % never assigned in his loop

kavin_cloud = avgvis < 1000 & lwc > 0.01;
kavin_haze  = avgvis > 1000 & avgvis < 5000;
q_masks = calendar_quarter_masks(t_fm);
row_h = 5 / 60;                                             % hours per FM-120 row

% =================== Line 2: BAMS-paper definitions ===========================================
% one visibility record: both files joined, each time stamp once (the files overlap on
% 2023-07-19 04:00-13:57 with identical values), screened
t_all = [vis_files{1}.time; vis_files{2}.time];
v_all = [vis_files{1}.m;    vis_files{2}.m];
[t_all, first] = unique(t_all, 'first');
v_all = v_all(first);
ok = v_all > 0 & v_all <= settings.vis_max_valid_m;          % remove zeros and spikes
t_all = t_all(ok);
v_all = v_all(ok);

% median and count of the samples in each FM-120 row's interval [t, t + 5 min)
[has_row, which_row] = ismember(floor_time(t_all, 5), t_fm);
vis_median = accumarray(which_row(has_row), v_all(has_row), [n_rows 1], @median, NaN);
vis_count  = accumarray(which_row(has_row), 1,              [n_rows 1]);
vis_median(vis_count < settings.vis_min_samples) = NaN;      % too few samples

operating = isfinite(fm.nd_cm3) & isfinite(lwc);
bams_cloud = operating & vis_median < settings.cloud_vis_max_m & lwc > settings.cloud_lwc_min_gm3;
bams_haze  = operating & vis_median > 0 & vis_median < settings.haze_vis_max_m & ...
             lwc < settings.haze_lwc_max_gm3 & ~bams_cloud;
s_masks = sheet_season_masks(t_fm);

% =================== Results ================================================================
results = [ ...
    make_result(31, "On-Site Cloud Sampling by Fog Drop Monitor", "hr", "Soledad", "FM120 & msd_visibility", ...
        season_sums(row_h * kavin_cloud, q_masks), ...
        "Kavin's code re-run on the same files, not truncated (truncating each season to whole hours gives the sheet)."), ...
    make_result(31, "On-Site Cloud Sampling by Fog Drop Monitor", "hr", "Soledad", "FM120 & msd_visibility", ...
        season_sums(row_h * bams_cloud, s_masks), ...
        sprintf("filters: BAMS-paper definition: median visibility < %g km and LWC > %g g/m3; ", ...
            settings.cloud_vis_max_m / 1000, settings.cloud_lwc_min_gm3) + vis_note(settings)), ...
    make_result(32, "On-Site Haze Sampling by Fog Drop Monitor", "hr", "Soledad", "FM120 & msd_visibility", ...
        season_sums(row_h * kavin_haze, q_masks), ...
        "Kavin's code re-run on the same files, not truncated (truncating each season to whole hours gives the sheet)."), ...
    make_result(32, "On-Site Haze Sampling by Fog Drop Monitor", "hr", "Soledad", "FM120 & msd_visibility", ...
        season_sums(row_h * bams_haze, s_masks), ...
        sprintf("filters: BAMS-paper definition: median visibility < %g km and LWC < %g g/m3; ", ...
            settings.haze_vis_max_m / 1000, settings.haze_lwc_max_gm3) + vis_note(settings))];

end


function avg = kavin_window_mean(vis, t_lo, t_hi)
% Mean of the visibility samples with t_lo(i) <= time <= t_hi(i), for every i (NaN if none).
%
% Same result as Kavin's loop
%     avgvis(i) = mean(vis(time <= t5min(i+1) & time >= t5min(i)), 'omitmissing')
% but without a loop over ~90,000 rows: with the samples in time order, the
% samples in a window are numbers count(< t_lo) + 1 to count(<= t_hi), and
% their sum is a difference of the running sum ("cumulative sum").
ts = posixtime(vis.time);                  % sample times (s), sorted
x  = vis.m;
good = isfinite(x);                        % 'omitmissing': NaN samples do not count
running_sum   = [0; cumsum(x .* good, 'omitnan')];
running_count = [0; cumsum(good)];
first = count_sorted(ts, posixtime(t_lo), 'below');        % samples before the window
last  = count_sorted(ts, posixtime(t_hi), 'at_or_below');  % samples up to the window's end
n_in  = running_count(last + 1) - running_count(first + 1);
avg   = (running_sum(last + 1) - running_sum(first + 1)) ./ n_in;
avg(n_in == 0) = NaN;                      % the mean of no samples
end


function txt = vis_note(settings)
txt = sprintf("interval visibility = median of the screened 1-s samples in each 5-min FM-120 interval " + ...
    "(0 < vis <= %s m; >= %d samples); hours = rows x 5 min; sheet seasons (UTC).", ...
    "48,270", settings.vis_min_samples);
end
