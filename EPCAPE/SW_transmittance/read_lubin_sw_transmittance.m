%% Read Dan Lubin's hourly SW transmittance tables (<MON>_SW_Boxplot.csv) as one hourly UTC series
%
% DESCRIPTION:
%   Dan Lubin computed an hourly shortwave transmittance at the Scripps Pier
%   for EPCAPE: measured SW down divided by a clear-sky model, which Kavin's
%   table note names as the analytical model of Atwater & Ball (1981). He
%   shared one CSV per calendar month (Kavin used them for the BAMS paper's
%   Fig. 5 and the table's row 7). This function reads them and returns one
%   hourly series in UTC. It combines read_lubin_sw and lubin_sw_series of
%   the Python module comparisons/seasonal_averages/kavin.py.
%
%   FILE FORMAT (<MON>_SW_Boxplot.csv, MON = JAN ... DEC)
%     * line 1: header "0,1,...,23,Means" (the files start with a byte-order
%       mark, which readmatrix skips with the header line)
%     * then one line per day of the month
%     * columns 1-24 = hour of the day 0-23 in Pacific STANDARD time (UTC-8,
%       no daylight saving: daylight is centred near the same column all
%       year, and the Python notebook's Section 1b found the best match with
%       RADFLUX for an 8-h shift); column 25 ("Means") is not used
%     * night hours are empty, except on the first row, whose night hours
%       hold -0.05: a plotting placeholder, not a measurement. No real value
%       is negative, so every value < 0 is set to NaN here.
%   Each value is taken to describe the hour [h, h + 1 h) PST, i.e. the
%   column label is the hour's START. compare_sw_transmittance_with_lubin.m
%   tests this by shifting the series and checking the correlation.
%
%   WHICH DAY IS EACH ROW
%   * MAR-DEC are 2023 and JAN is 2024; row r is day r of the month.
%   * FEB holds days from both Februaries of the campaign (15 Feb 2023 -
%     14 Feb 2024), 28 rows. Two readings exist:
%       "chronological" (default): rows 1-14 = 15-28 Feb 2023 and rows 15-28
%           = 1-14 Feb 2024. This is what matching each row against the
%           RADFLUX SW / clear-sky ratio day by day showed (EPCAPE BAMS
%           Fig. 4 work, Oct 2026: rows 3-14 matched 17-28 Feb 2023 and rows
%           15-23 matched 1-9 Feb 2024; the rest were mostly clear days that
%           matched several dates equally well).
%       "dayOfMonth": row r = day r, with days 1-14 in 2024 and 15-28 in
%           2023. This is the Python lubin_sw_series assumption.
%     The comparison script reports the correlation with RADFLUX for both.
%
%   WHERE THE FILES ARE
%   Pass one folder or several. Each month is read from the FIRST folder
%   that has its file. (Kavin_data_for_BAMS_paper holds APR, MAY, JUL-DEC;
%   Kavin_code_for_BAMS/Annual Echo Plot holds all twelve; the APR files in
%   the two folders are byte-identical.)
%
% INPUTS:
%   folders - folder path (string/char) or a string array of folders, in
%             order of preference
%   Optional name-value pairs:
%   'FebMapping'  "chronological" (default) | "dayOfMonth"
%   'UtcOffset_h' hours to ADD to the PST hour to get UTC (default 8)
%
% OUTPUTS:
%   T       - timetable, one row per non-empty (day, hour) cell, sorted:
%               time (row times)  hour START in UTC
%               T_lubin           Dan's transmittance (unitless)
%               month, row        calendar month of the file and its row
%               date_local        the PST date the row was assigned to
%               hour_pst          column hour 0-23 (PST)
%   sources - table: month, file (full path) or "" if no folder had it
%
% REFERENCES:
%   Atwater, M. A., and J. T. Ball, 1981 (clear-sky model named in Kavin's
%     table note; the reference was not checked here).
%   Python original: VS_CODE/Python-Research/EPCAPE/comparisons/
%     seasonal_averages/kavin.py (read_lubin_sw, corrected_sw_cells,
%     lubin_sw_series, LUBIN_UTC_OFFSET_H).
%   FEB row mapping: EPCAPE/Kavin_code_for_BAMS/Fig 4/EPCAPE_BAMS_Fig4_AJB_v4_2.m,
%     header assumption b.
%
% By Andrew John Buggee
%%

function [T, sources] = read_lubin_sw_transmittance(folders, varargin)


% -------------------------------------------------------------------
% ------------------------ Parse the inputs -------------------------
% -------------------------------------------------------------------

p = inputParser;
addParameter(p, 'FebMapping', "chronological", @(x) any(strcmpi(string(x), ["chronological" "dayOfMonth"])));
addParameter(p, 'UtcOffset_h', 8, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});
opt = p.Results;

folders = string(folders);
month_names = ["JAN" "FEB" "MAR" "APR" "MAY" "JUN" "JUL" "AUG" "SEP" "OCT" "NOV" "DEC"];



% -------------------------------------------------------------------
% -------------- Read each month from the first folder --------------
% -------------------------------------------------------------------

file = strings(12, 1);
parts = cell(12, 1);

for m = 1:12

    % the first folder that has this month's file
    for f = 1:numel(folders)
        candidate = fullfile(folders(f), month_names(m) + "_SW_Boxplot.csv");
        if isfile(candidate)
            file(m) = candidate;
            break
        end
    end
    if file(m) == ""
        continue
    end

    % one row per day, columns = PST hours 0-23 (empty cells -> NaN)
    M = readmatrix(file(m), 'NumHeaderLines', 1, 'Delimiter', ',');
    M = M(:, 1:24);                % drop the "Means" column
    M(M < 0) = NaN;                % the -0.05 night placeholders (no real value is negative)
    n_rows = size(M, 1);
    rows = (1:n_rows)';

    % the PST date of each row
    if m == 2
        switch lower(string(opt.FebMapping))
            case "chronological"
                % rows 1-14 = 15-28 Feb 2023, rows 15-28 = 1-14 Feb 2024
                date_local = NaT(n_rows, 1);
                first = rows <= 14;
                date_local(first)  = datetime(2023, 2, 14 + rows(first));
                date_local(~first) = datetime(2024, 2, rows(~first) - 14);
            case "dayofmonth"
                % row r = day r; days 1-14 in 2024, 15-28 in 2023
                yr = 2023 * ones(n_rows, 1);
                yr(rows <= 14) = 2024;
                date_local = datetime(yr, 2, rows);
        end
    else
        yr = 2023;
        if m == 1
            yr = 2024;
        end
        date_local = datetime(yr, m, rows);
        % a row beyond the month's length (e.g. a 31st row in a 30-day
        % month) has no date; datetime would roll it into the next month
        bad = month(date_local) ~= m;
        date_local(bad) = NaT;
    end

    % one entry per (row, hour) cell, at the hour's START in UTC
    [hh, rr] = meshgrid(0:23, rows);                  % rr = row, hh = PST hour
    vals = M(:);
    d = date_local(rr(:));
    keep = isfinite(vals) & ~isnat(d);
    time_utc = d(keep) + hours(hh(keep)) + hours(opt.UtcOffset_h);

    parts{m} = table(time_utc, vals(keep), m * ones(nnz(keep), 1), rr(keep), d(keep), hh(keep), ...
        'VariableNames', {'time', 'T_lubin', 'month', 'row', 'date_local', 'hour_pst'});
end



% -------------------------------------------------------------------
% ------------------------- Assemble outputs ------------------------
% -------------------------------------------------------------------

T = vertcat(parts{:});
T = sortrows(table2timetable(T, 'RowTimes', 'time'));

sources = table((1:12)', file, 'VariableNames', {'month', 'file'});


end
