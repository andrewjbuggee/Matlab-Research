function T = read_posted_product(file_path)
% READ_POSTED_PRODUCT  Read one EPCAPE data file (Jeramy Dedrick's posted products, or Sourita Saha's CSV).
%
%   T = read_posted_product(file_path)
%
%   Two file layouts are recognised from the header line:
%
%   1. Jeramy Dedrick's posted products, from the UC San Diego Library item
%      "Aerosol Quantification Products at Mt. Soledad and Scripps Pier during
%      EPCAPE 2023-24" (https://doi.org/10.6075/J0WM1DSC). Tab-separated .txt;
%      the first column, "time", is a MATLAB datenum (days since year 0, UTC),
%      and the other columns are numbers (missing = NaN).
%
%   2. Sourita Saha's Gaussian-corrected CCN file at the 0.1 % setpoint
%      (ss01_ccn_gauss_corrected_multiple_peak_fitting.csv, shared by email on
%      20 Oct 2026). Comma-separated; the time is split over the columns
%      year, month, day, hour, minute, second (1-s samples). The CCN in it is
%      straight from ARM and NOT multiplied by Jeramy's 1.57 leak factor
%      (Sourita Saha, email of 20 Oct 2026). The file does not state its time
%      zone. It is UTC like the ARM data it comes from (my inference, supported
%      by the check in Part 9 of the live script: read as UTC, every posted
%      0.1 % row of the UTC CCN file contains one of its 10-min stages, and
%      1.57 x the stage mean correlates with the posted value at r = 0.996).
%
%   Output: a timetable with the file's other columns as variables. The column
%   names are kept exactly as written (e.g. '0.1%' in the CCN file), so use
%   T{:, k} or T.("0.1%") to read them.

% The delimiter must be given: with a tab delimiter forced on a CSV, the whole
% header line becomes one variable name, which is truncated to namelengthmax and
% the 'time' column is never found (the error that prompted this layout check).
[~, ~, ext] = fileparts(file_path);
if strcmpi(ext, '.csv')
    delimiter = ',';
else
    delimiter = '\t';
end
opts = detectImportOptions(file_path, 'FileType', 'text', 'Delimiter', delimiter, ...
                           'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'double');   % every column is numeric (NaN where missing)
raw = readtable(file_path, opts);

names = raw.Properties.VariableNames;
date_parts = {'year', 'month', 'day', 'hour', 'minute', 'second'};

if ismember('time', names)
    % layout 1: datenum -> datetime, rounded to the nearest second (the datenums carry
    % float round-off, e.g. 738932.041666667 = 01:00:00 minus a few microseconds)
    t = datetime(raw.time, 'ConvertFrom', 'datenum');
    t = dateshift(t, 'start', 'second', 'nearest');
    raw.time = [];

elseif all(ismember(date_parts, names))
    % layout 2: calendar columns -> datetime (whole seconds already, no rounding needed)
    t = datetime(raw.year, raw.month, raw.day, raw.hour, raw.minute, raw.second);
    raw(:, date_parts) = [];

else
    error('read_posted_product:noTime', ...
          'No time column in %s: expected "time" (datenum) or year, month, day, hour, minute, second.', file_path);
end

T = table2timetable(raw, 'RowTimes', t);
T.Properties.DimensionNames{1} = 'time';

end
