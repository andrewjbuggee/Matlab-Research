%% make_table_for_BAMS.m
% Re-derives the rows of the EPCAPE seasonal-averages table
% (Table-EPCAPE_SeasonalAverages.xlsx, sheet "Cloud&MetQuantities") that are
% attributed to Kavin, and writes the new lines ("Andrew Buggee" rows) to an
% Excel file and a CSV file.
%
% HOW TO RUN
%   1. Edit the three folder paths in Section 1 below.
%   2. Press Run. It takes about a minute; the slowest steps are the 22 million
%      MWRRET v2 samples and the 25 million visibility samples.
%   3. The results are printed in the Command Window and saved in
%      output_folder (see Section 4).
%
% WHAT EACH LINE IS
%   Every result line says, in its Notes, which data file was used, which
%   samples were kept (filters) and how the average was taken. There are two
%   kinds of lines:
%     * "Kavin's method re-run": his own calculation on his own input file,
%       written out step by step (sometimes with a bug fixed; the Notes say so);
%     * an independent calculation from the original data.
%
% SEASONS (one column each, as in the sheet)
%   EPCAPE = 2023-02-15 .. 2024-02-14, Spring = Apr-Jun 2023, Summer = Jul-Sep
%   2023, Fall = Oct-Dec 2023, Winter = 15 Feb-Mar 2023 + Jan-14 Feb 2024 (UTC).
%   Some of Kavin's calculations use calendar quarters instead (winter =
%   Jan-Mar); the Notes say "calendar quarters" where that is the case.
%   Every value is the mean ± (sample) standard deviation over all samples
%   in the season, so the campaign value is NOT the mean of the 4 seasons.
%
% QUALITY CONTROL
%   Every ARM variable that has a qc_<variable> field keeps only samples
%   with QC = 0 (every test passed). Variables without a QC field are used
%   as distributed (the Notes say so).
%
% FOLDERS
%   calculations/  one function per sheet row (or group of rows)
%   readers/       functions that read the data files
%   helpers/       seasons, statistics, number formatting, time handling
%
% Written for MATLAB R2024b and later (tested in R2024b); base MATLAB only,
% no toolboxes needed.

clear; clc;

%% 1. Settings: EDIT THESE FOLDERS ------------------------------------------------------
% Folder with the data Andrew shared (subfolders ARM/, FM120_daily_5min/, visibility/)
data_folder  = "/Users/andrewbuggee/Documents/Scripps - UCSD/EPCAPE/table_for_BAMS_data";
% Folder with Kavin's own input files (Lubin's *_SW_Boxplot.csv, EPCAPE_15mindata.nc,
% ET_Dcrit_nativetime.csv, epcmwrret_turnM1_hourly_averages.nc)
kavin_folder = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper";
% Where the results are written (created if it does not exist)
output_folder = fullfile(fileparts(mfilename('fullpath')), 'output');

% Thresholds and choices (each one changes a number in the table)
settings.low_cloud_max_m   = 3500;   % "low cloud": cloud top (KAZR) or base (ceilometer) <= 3.5 km
settings.rain_window_min   = 30;     % LWP "rain-free": no Pier rain within ±30 min
settings.rain_day_min_mm   = 0.1;    % a "rain day" has >= 0.1 mm
settings.cloud_vis_max_m   = 1000;   % FM-120 cloud: visibility < 1 km ...
settings.cloud_lwc_min_gm3 = 0.01;   % ... and LWC > 0.01 g m-3   (BAMS-paper definition)
settings.haze_vis_max_m    = 5000;   % FM-120 haze: visibility < 5 km ...
settings.haze_lwc_max_gm3  = 0.01;   % ... and LWC < 0.01 g m-3   (BAMS-paper definition)
settings.vis_max_valid_m   = 48270;  % visibility sensor maximum (30 statute miles); above = spike
settings.vis_min_samples   = 125;    % a 5-min interval needs >= 125 of its ~250 visibility samples

% Make the functions in the subfolders visible to MATLAB
code_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(code_folder, 'calculations'), fullfile(code_folder, 'readers'), fullfile(code_folder, 'helpers'));

arm_folder = fullfile(data_folder, 'ARM');
tic_all = tic;

%% 2. Rows 7-24: radiation, LWP, cloud boundaries, PBL height, LCL, rain ----------------------
fprintf('Row 7  SW transmittance ...\n');           r07 = calc_row07_SW_transmittance(kavin_folder);
fprintf('Row 8  LW down ...\n');                    r08 = calc_row08_LW_down(arm_folder);
fprintf('Row 10 LWP (slowest step) ...\n');         r10 = calc_row10_LWP(arm_folder, kavin_folder, settings.rain_window_min);
fprintf('Row 15 cloud top ...\n');                  [r15, low_cloud_5min] = calc_row15_cloud_top(arm_folder, settings.low_cloud_max_m);
fprintf('Row 17 PBL height ...\n');                 r17 = calc_row17_PBL_height(arm_folder, low_cloud_5min, settings.low_cloud_max_m);
fprintf('Row 20 LCL (one file per launch) ...\n');  r20 = calc_row20_LCL(arm_folder);
fprintf('Row 21 cloud base ...\n');                 r21 = calc_row21_cloud_base(arm_folder, settings.low_cloud_max_m);
fprintf('Rows 22-24 rain ...\n');                   r22 = calc_rows22to24_rain(arm_folder, settings.rain_day_min_mm);

%% 3. Rows 25-55: FM-120, activated fraction, GCVI-AMS (Mt. Soledad) ----------------------------
fprintf('Reading the FM-120 daily files ...\n');
fm = read_fm120_daily_files(fullfile(data_folder, 'FM120_daily_5min'));
fprintf('Reading the 15-min file ...\n');
F15 = read_15min_file(fullfile(kavin_folder, 'EPCAPE_15mindata.nc'));
fprintf('Reading the visibility files ...\n');
vis_files = {read_visibility_file(fullfile(data_folder, 'visibility', 'msd_visibility_20230224-20230719.csv')), ...
             read_visibility_file(fullfile(data_folder, 'visibility', 'msd_visibility_20230719-20240214.csv'))};

fprintf('Rows 25-55 ...\n');
r15min = calc_rows_15min_file(F15);
r25    = calc_rows25to29_FM120_5min(fm, F15);
r27    = calc_rows27_30_activated_fraction(kavin_folder);
r31    = calc_rows31_32_FM120_hours(fm, vis_files, settings);

%% 4. Collect, order like the sheet, print and save -----------------------------------------
% Listed so that, within each sheet row, Kavin's-method lines come first (r15min before r25)
results = [r07, r08, r10, r15, r17, r20, r21, r22, r15min, r25, r27, r31];
% Order by sheet row. MATLAB's sort is "stable": lines with the same sheet row keep the order above.
[~, order] = sort([results.sheet_row]);
results = results(order);

% One table: sheet row, columns A-D of the sheet, the five season cells, author, notes
n = numel(results);
T = table('Size', [n 12], ...
    'VariableTypes', ["double", repmat("string", 1, 11)], ...
    'VariableNames', ["SheetRow", "Quantity", "Units", "Location", "Instrument_Product_Method", ...
                      "EPCAPE", "Spring", "Summer", "Fall", "Winter", "Author", "Notes"]);
for i = 1:n
    r = results(i);
    T(i, :) = {r.sheet_row, r.quantity, r.units, r.location, r.product, ...
               r.text(1), r.text(2), r.text(3), r.text(4), r.text(5), "Andrew Buggee", r.notes};
end

% The unrounded numbers, one line per (result, season), for checking
seasons = ["EPCAPE", "Spring", "Summer", "Fall", "Winter"];
N = table();
for i = 1:n
    r = results(i);
    for k = 1:5
        N = [N; table(r.sheet_row, i, r.product, seasons(k), r.mean(k), r.std(k), r.n(k), r.text(k), ...
            'VariableNames', ["SheetRow", "Line", "Product", "Season", "Mean", "Std", "N", "Text"])]; %#ok<AGROW>
    end
end

disp(T(:, ["SheetRow", "Quantity", "EPCAPE", "Spring", "Summer", "Fall", "Winter"]));

if ~isfolder(output_folder)
    mkdir(output_folder);
end
xlsx_file = fullfile(output_folder, 'table_for_BAMS_MATLAB.xlsx');
if isfile(xlsx_file)
    delete(xlsx_file);                    % start from an empty workbook each run
end
writetable(T, xlsx_file, 'Sheet', 'Table');
writetable(N, xlsx_file, 'Sheet', 'Numbers');
writetable(N, fullfile(output_folder, 'table_for_BAMS_MATLAB_numbers.csv'));
fprintf('\nWrote %s\n(sheet "Table": the new lines; sheet "Numbers": unrounded means, stds and counts)\n', xlsx_file);
fprintf('Total time: %.1f min\n', toc(tic_all) / 60);
