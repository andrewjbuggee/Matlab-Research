function results = calc_row07_SW_transmittance(kavin_folder)
% CALC_ROW07_SW_TRANSMITTANCE  Sheet row 7: shortwave (SW) transmittance.
%
%   results = calc_row07_SW_transmittance(kavin_folder)
%
%   DATA: Dan Lubin's hourly SW transmittance tables, <MON>_SW_Boxplot.csv,
%   in Kavin's folder (the same files Kavin used). Each file is one calendar
%   month:
%     * line 1 = header "0,1,...,23,Means"
%     * then one line per day of the month (day 1 first)
%     * columns 1-24 = hour of day 0-23 in Pacific STANDARD time (UTC-8);
%       column 25 ("Means") is not used
%     * each value = measured SW down / a clear-sky model (sheet note:
%       Atwater & Ball 1981); night hours are empty
%     * day 1 holds -0.05 in the night hours: a plotting placeholder, not a
%       measurement
%
%   WHAT IS COMPUTED (one new line for the table):
%     every (day, hour) cell of the season's three months, with the -0.05
%     placeholders removed; mean ± std over all those cells. Seasons are
%     calendar quarters (Kavin's), so winter = Jan-Mar.
%
%   HOW THIS DIFFERS FROM KAVIN'S CODE: his MATLAB readtable read the header
%   line as a data row, so his rows 2:31 are days 1-30 (day 31 dropped) and
%   his Var2:Var24 are hours 1-23 (hour 0 dropped); his means also include
%   the -0.05 placeholders. Here every day and hour is used and the
%   placeholders are removed.
%
%   MISSING FILES: JAN, FEB, MAR and JUN are not in the folder, so only
%   summer and fall can be computed; the other columns are "n/a".

month_names = ["JAN" "FEB" "MAR" "APR" "MAY" "JUN" "JUL" "AUG" "SEP" "OCT" "NOV" "DEC"];

% --- 1. Read every month that exists -------------------------------------------------
cells = cell(1, 12);          % cells{m} = matrix (days x 24 hours) for month m, [] if missing
for m = 1:12
    path = fullfile(kavin_folder, month_names(m) + "_SW_Boxplot.csv");
    if ~isfile(path)
        continue
    end
    M = readmatrix(path, 'NumHeaderLines', 1, 'Delimiter', ',');  % empty cells -> NaN
    M = M(:, 1:24);           % hours 0-23 (drop the "Means" column)
    M(M < 0) = NaN;           % remove the -0.05 placeholders (no real value is negative)
    cells{m} = M;
end

% --- 2. Mean ± std per season (calendar quarters) --------------------------------------
season_months = {1:12, [4 5 6], [7 8 9], [10 11 12], [1 2 3]};   % EPCAPE, Spring, Summer, Fall, Winter
stats.mean = nan(1, 5);
stats.std  = nan(1, 5);
stats.n    = zeros(1, 5);
for k = 1:5
    months = season_months{k};
    if any(cellfun(@isempty, cells(months)))
        continue                       % a month file is missing: leave this season as n/a
    end
    v = vertcat(cells{months});        % stack the months' tables
    v = v(isfinite(v));                % every non-empty (day, hour) cell
    stats.mean(k) = mean(v);
    stats.std(k)  = std(v);
    stats.n(k)    = numel(v);
end

results = make_result(7, "Shortwave (SW) Transmittance", "fraction", "Pier", ...
    "Lubin's hourly SW transmittance tables (Kavin's files)", stats, ...
    "Same data as Kavin (Lubin's <MON>_SW_Boxplot.csv). filters: every day and hour of each month; " + ...
    "-0.05 placeholders removed; mean ± std over all (day, hour) cells; calendar quarters (winter = Jan-Mar). " + ...
    "JAN, FEB, MAR and JUN tables not provided, so EPCAPE, Spring and Winter are n/a.");

end
