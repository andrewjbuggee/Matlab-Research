function results = calc_rows25to29_FM120_5min(fm, F)
% CALC_ROWS25TO29_FM120_5MIN  Rows 25, 26, 28, 29 from the FM-120 5-min files.
%
%   results = calc_rows25to29_FM120_5min(fm, F)
%     fm = read_fm120_daily_files(...)   FM-120 5-min rows
%     F  = read_15min_file(...)          Abbey Williams' 15-min classes
%
%   Same cloud/haze classes as Kavin (Abbey Williams' ClrCldHaz_class), but
%   with the FM-120's own 5-min values instead of the 15-min averages:
%     * each 5-min FM-120 row takes the class of the 15-min interval that
%       holds it (its time rounded down to 15 min);
%     * number concentration = the FM-120 "Number Conc" column (2-50 um, cm-3);
%     * effective diameter  = the FM-120 "ED" column (um), set to NaN when
%       no droplet was counted (N = 0);
%     * mean ± std over the rows of the class; each row counts once;
%       sheet seasons (UTC).

% class of the 15-min interval holding each 5-min row (NaN if there is none)
[found, idx] = ismember(floor_time(fm.time, 15), F.time);
row_class = nan(size(fm.time));
row_class(found) = F.class(idx(found));

ed_um = fm.ed_um;
ed_um(~(fm.nd_cm3 > 0)) = NaN;                 % no droplets -> no diameter

masks = sheet_season_masks(fm.time);
src = "FM120 daily 5-min files + Abbey's 15-min classes";
note_nd = @(code, kind) sprintf("filters: FM-120 5-min rows, each taking the ClrCldHaz_class of the 15-min interval " + ...
    "holding it (= %d, %s); nd_cm3 (Number Conc, 2-50 um); mean ± std over rows; sheet seasons (UTC).", code, kind);
note_ed = @(code, kind) sprintf("filters: FM-120 5-min rows, each taking the ClrCldHaz_class of the 15-min interval " + ...
    "holding it (= %d, %s); ed_um (the FM-120 ED column, the effective diameter), NaN when N = 0; mean ± std " + ...
    "over rows; sheet seasons (UTC).", code, kind);

only = @(values, code) values .* nan_unless(row_class == code);

results = [ ...
    make_result(25, "Cloud Droplet Number Concentration", "cm-3", "Soledad", src, ...
        season_stats(only(fm.nd_cm3, 1), masks), note_nd(1, "cloud")), ...
    make_result(26, "Cloud Droplet Effective Diameter", "um", "Soledad", src, ...
        season_stats(only(ed_um, 1), masks), note_ed(1, "cloud")), ...
    make_result(28, "Haze Droplet Number Concentration", "cm-3", "Soledad", src, ...
        season_stats(only(fm.nd_cm3, 2), masks), note_nd(2, "haze")), ...
    make_result(29, "Haze Droplet Effective Diameter", "um", "Soledad", src, ...
        season_stats(only(ed_um, 2), masks), note_ed(2, "haze"))];

end
