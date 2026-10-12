function fm = read_fm120_daily_files(folder)
% READ_FM120_DAILY_FILES  Read the FM-120 fog monitor's daily 5-min files.
%
%   fm = read_fm120_daily_files(folder)
%
%   The DMT FM-120 fog monitor ran at Mt. Soledad, next to the GCVI inlet,
%   and measured droplets 2-50 um across. Lauren Robinson and Rachel Chang
%   (Dalhousie) provided one CSV per day, named YYYY-MM-DD.csv, holding
%   5-min averages of the 1-s data (see readme_FM120.md in the same folder).
%   Column 1 is the START of each 5-min interval (UTC). There is no QC field.
%
%   File quirks handled here:
%   * Some files end with a row for 00:00 of the NEXT day. The README says
%     to ignore it (that interval is in the next day's file), so only rows
%     dated on the file's own day are kept.
%   * Bin column names sometimes end with a space; only the columns below
%     are used, and they are found by name.
%   * A blank cell (one LWC value on 2023-11-14) is read as NaN.
%   * 2024-01-02.csv has a header and no data; it is skipped.
%
%   Output: struct with column vectors (one element per 5-min interval)
%     fm.time     interval start (datetime, UTC)
%     fm.nd_cm3   droplet number concentration, 2-50 um (# cm-3) ["Number Conc (#/cm^3)"]
%     fm.lwc_gm3  liquid water content (g m-3)                    ["LWC (g/m^3)"]
%     fm.ed_um    "ED" = equivalent (effective) diameter (um)      ["ED (um)"]

files = dir(fullfile(folder, '*.csv'));
files = sort({files.name});
files = files(~cellfun(@isempty, regexp(files, '^\d{4}-\d{2}-\d{2}\.csv$')));   % keep YYYY-MM-DD.csv only
if isempty(files)
    error('No FM-120 daily files (YYYY-MM-DD.csv) in %s', folder);
end

wanted = ["Number Conc (#/cm^3)", "LWC (g/m^3)", "ED (um)"];
time = {}; nd = {}; lwc = {}; ed = {};

for i = 1:numel(files)
    path = fullfile(folder, files{i});
    file_day = datetime(files{i}(1:10), 'InputFormat', 'yyyy-MM-dd');

    % Read the header line to find the columns we need (names compared without spaces at the ends)
    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve', 'Delimiter', ',');
    names = strtrim(string(opts.VariableNames));
    cols = zeros(1, numel(wanted));
    for k = 1:numel(wanted)
        hit = find(names == wanted(k), 1, 'first');   % 'first': a repeated column keeps its first copy
        if isempty(hit)
            cols = [];   % e.g. 2024-01-02.csv, whose header differs and which has no data
            break
        end
        cols(k) = hit;
    end
    if isempty(cols)
        continue
    end

    % Column 1 = time stamp text, the wanted columns = numbers
    opts.SelectedVariableNames = opts.VariableNames([1 cols]);
    opts = setvartype(opts, opts.VariableNames{1}, 'datetime');
    opts = setvaropts(opts, opts.VariableNames{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    opts = setvartype(opts, opts.VariableNames(cols), 'double');
    T = readtable(path, opts);
    if height(T) == 0
        continue
    end

    % Keep only rows dated on the file's own day
    own_day = T{:, 1} >= file_day & T{:, 1} < file_day + days(1);
    time{end + 1} = T{own_day, 1};   %#ok<AGROW>
    nd{end + 1}   = T{own_day, 2};   %#ok<AGROW>
    lwc{end + 1}  = T{own_day, 3};   %#ok<AGROW>
    ed{end + 1}   = T{own_day, 4};   %#ok<AGROW>
end

fm.time    = vertcat(time{:});
fm.nd_cm3  = vertcat(nd{:});
fm.lwc_gm3 = vertcat(lwc{:});
fm.ed_um   = vertcat(ed{:});

% Sort by time and drop any repeated time stamp (keeping its first copy)
[fm.time, order] = unique(fm.time, 'first');
fm.nd_cm3  = fm.nd_cm3(order);
fm.lwc_gm3 = fm.lwc_gm3(order);
fm.ed_um   = fm.ed_um(order);

if any(isnat(fm.time))
    error('read_fm120_daily_files: some time stamps could not be read');
end

end
