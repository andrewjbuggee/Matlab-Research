function T = read_posted_product(file_path)
% READ_POSTED_PRODUCT  Read one tab-separated file of Jeramy Dedrick's posted EPCAPE products.
%
%   T = read_posted_product(file_path)
%
%   The files come from the UC San Diego Library item "Aerosol Quantification
%   Products at Mt. Soledad and Scripps Pier during EPCAPE 2023-24"
%   (https://doi.org/10.6075/J0WM1DSC). Each is a tab-separated text file whose
%   first column, "time", is a MATLAB datenum (days since year 0, UTC) and whose
%   other columns are numbers (missing = NaN).
%
%   Output: a timetable with the file's other columns as variables. The column
%   names are kept exactly as written (e.g. '0.1%' in the CCN file), so use
%   T{:, k} or T.("0.1%") to read them.

opts = detectImportOptions(file_path, 'FileType', 'text', 'Delimiter', '\t', ...
                           'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'double');   % every column is numeric (NaN where missing)
raw = readtable(file_path, opts);

% datenum -> datetime, rounded to the nearest second (the datenums carry float round-off,
% e.g. 738932.041666667 = 01:00:00 minus a few microseconds)
t = datetime(raw.time, 'ConvertFrom', 'datenum');
t = dateshift(t, 'start', 'second', 'nearest');

raw.time = [];
T = table2timetable(raw, 'RowTimes', t);
T.Properties.DimensionNames{1} = 'time';

end
