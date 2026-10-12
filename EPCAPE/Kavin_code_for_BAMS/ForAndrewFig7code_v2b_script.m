%% EPCAPE Figure 7 (version 2b): seasonal cloud vs. haze droplets, CCN activation, aerosol composition and cloud-base height
% Version 2b of Kavin's ForAndrewFig7code.mlx (v1), updated 2026-10-08. Same
% figure as version 2a (Kavin's numbers) except for the numbers and the D_eff
% panel listed here, which come from Andrew's verification table,
% "Table-EPCAPE_SeasonalAverages - AJB verification - ver 2.xlsx"
% (sheet Cloud&MetQuantities), produced with check_seasonal_averages.ipynb.
%
% *From the verification table*
%
% Row 66 / row 70, FM-120 cloud / haze sampling hours with the BAMS-paper
% definitions (cloud = median visibility < 1 km and LWC > 0.01 g m^-3; haze =
% median visibility < 5 km and LWC < 0.01 g m^-3; visibility = median of the
% screened 1-s samples in each 5-min FM-120 interval): the blue "Cloud: N Hours"
% and red "Haze: N Hours" labels, rounded to the nearest hour.
%
% Rows 51 / 60, cloud / haze droplet effective diameter: the D_eff panel. Each
% daily FM-120 5-min row takes the ClrCldHaz_class of the 15-min interval of
% EPCAPE_15mindata.nc that holds it; D_eff = the FM-120 ED column (the
% instrument's effective diameter, which uses bin mid-points), NaN when N = 0;
% sheet seasons. Kavin's D_eff (version 2a) instead takes each 15-min
% FM120_DSD with every bin's UPPER edge as its diameter. The notebook finds the
% ED column within 1 % (median) of sum(N D^3)/sum(N D^2) recomputed from the
% 5-min bins with mid-points.
%
% Rows 53 / 62, cloud / haze activated fraction: the "% Activated" labels
% (Kavin's count re-run; the same values as before, 94/100, 85/56, 32).
%
% The script reproduces rows 51 and 60 and prints them next to the table.
% Everything else (size distributions, D_crit, pies, cloud-base heights) is
% unchanged from version 2a: Kavin's method on the 15-min file, calendar-quarter
% seasons, and Lynn's edits and the bug fixes below.
%
% *Lynn Russell's review comments (Google Doc, Sep 24 - Oct 3)*
%
% L1. Extend the CBH axis to 0 m and move the "CBH < 250" arrow to clarify.
% -> y-axis now 0-500 m with a new 0-50 m bin; the arrow now points DOWN from
% the 250 m Soledad line, i.e. it marks the whole region below it.
%
% L2. Make the organics green brighter, like the example.
% -> organics wedge = pure green [0 1 0] (#00FF00), drawn opaque.
%
% L3. NO3 -> Nitrate, SO4 -> Sulfate, NH4 -> Ammonium, Cl -> Chloride.
% -> legend names are read from the AMS_species variable of the NetCDF, which
% are exactly these names (this also removes the incorrect ion charges "SO4+"
% and "Cl+" of the v1 legend).
%
% L4. Close the boxes in the top row; less space between panels; more space
% below the "Droplet Diameter" x-axis label.
% -> manual (point-based) layout replaces tiledlayout.
%
% *Corrections to v1 code* (controlled by applyV1BugFixes below; false
% reproduces v1 behaviour so the two versions can be compared)
%
% B1. DSD shading. shadedErrorBar(x,y,errBar) expects errBar = OFFSETS from y
% (upper = y + errBar(1,:), lower = y - errBar(2,:); see shadedErrorBar.m).
% v1 passed the absolute 75th/25th percentiles, so the band was
% [median - p25, median + p75] instead of the interquartile range [p25, p75].
%
% B2. AMS pies. v1 drew the HAZE composition directly under the
% "Cloud: N Hours" label and the CLOUD composition under "Haze: N Hours".
% v2 puts the cloud pie under the cloud label and the haze pie under the haze
% label. (Check: Spring haze nitrate = 36 %, cloud nitrate = 28 %.)
%
% B3. Spring D_crit % activated labels. v1 passed the labels as [cloud, haze]
% for a [haze, cloud] data matrix (Summer used the correct [haze, cloud]
% order), and boxplot sorts numeric group labels, so Spring showed
% cloud = 100 %, haze = 94 %. Actual: cloud 94 %, haze 100 % (table rows 27, 30).
% Box-plot outliers are also drawn in the box colour (v1: always red).
%
% Always fixed: v1 used nativetime without defining it after clear all.
%
% *Not changed, but flagged for the authors* (see Fig7_open_questions.txt)
%
% (i) The y-axis is labelled dN/dlog(Dp), but dlogDp is computed with MATLAB's
% log (= natural log), i.e. the plotted quantity is dN/dlnDp. The usual
% aerosol convention dN/dlog10(Dp) is ln(10) = 2.303x larger
% (Seinfeld & Pandis, 2016, Atmospheric Chemistry and Physics, 3rd ed., Ch. 8).
%
% (ii) The "Cloud/Haze: N Hours" labels (rows 66, 70) come from a 5-min
% visibility + FM-120 LWC classification, whereas the size distributions, pies,
% D_crit and cloud-base panels use the 15-min ClrCldHaz_class (e.g. Spring:
% 137 h vs. 106 h of cloud).
%
% (iii) ARM_CBH_clean has no values below 50 m (raw ARM_CBH has 54 values in
% 40-50 m, all "other low clouds"), so the new 0-50 m bin is empty.
%
% (iv) Dta_critical exists only in Jun-Sep and Dec 2023, so the Spring % activated
% is June alone (941 cloud, 19 haze samples) and the Fall value is December alone.
%
% (v) The size-distribution row plots each bin at FM120_diam, the bin's UPPER
% edge (3, 4, ..., 50 um), as in v1, while the D_eff row now uses mid-points.
%
% *Inputs*
%
% EPCAPE_15mindata.nc: 15-min merged EPCAPE data set (Kavin). Its AMS_mass and
% AMS_massfrac equal AMS_mass_15min_EF and AMS_massfrac_15min_EF from the
% EPCAPE_15min_data.mat file used in v1 (checked against the means stored in
% the v1 output).
%
% ET_Dcrit_15min.csv: 15-min wet critical diameter (same 15-min grid).
%
% ET_Dcrit_nativetime.csv: native-time Dta, Dta_critical and cloud/haze class.
%
% Daily Averaged FM-120 FIles: 5-min FM-120 files (UTC), Mt. Soledad.
%
% Table-EPCAPE_SeasonalAverages - AJB verification - ver 2.xlsx: Andrew's
% verification of the seasonal-averages table (rows 51, 53, 60, 62, 66, 70).
%
% shadedErrorBar.m: in the Kavin_code_for_BAMS folder.

%% User settings
clear

dataDir = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper';
procDir = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed';
codeDir = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS';

ncFile        = fullfile(dataDir, 'EPCAPE_15mindata.nc');
dcrit15File   = fullfile(dataDir, 'ET_Dcrit_15min.csv');
dcritNatFile  = fullfile(dataDir, 'ET_Dcrit_nativetime.csv');
fm120Dir      = fullfile(procDir, 'Daily Averaged FM-120 FIles');
% Andrew's verification of the seasonal-averages table (rows read in section 4b)
verTableFile  = fullfile(procDir, 'derived', 'seasonal_averages_check', ...
                         'Table-EPCAPE_SeasonalAverages - AJB verification - ver 2.xlsx');
verTableSheet = 'Cloud&MetQuantities';

applyV1BugFixes = true;      % B1-B3 above; false = v1 behaviour (with Lynn's edits)
saveFigure      = true;      % write <figName>.png (300 dpi) and .pdf to outDir
outDir          = codeDir;
figName         = 'Fig7_v2b';
% 'off' in the live script (ForAndrewFig7code_v2b.mlx): in R2024b the Live Editor takes
% more than 10 min to export a visible copy of this figure (the same code run as a plain
% script takes ~40 s); the exported files are the same. 'on' in ForAndrewFig7code_v2b_script.m.
figureVisible   = 'on';

% shadedErrorBar.m lives next to this script
if isempty(which('shadedErrorBar')), addpath(codeDir); end

% Seasons as defined in v1 (calendar quarters, UTC months). The D_eff panel uses the
% sheet's seasons instead (local function inSheetSeason), as rows 51 and 60 do.
seasonName   = {'Spring', 'Summer', 'Fall', 'Winter'};
seasonMonths = [4 6; 7 9; 10 12; 1 3];          % [first last] month of each season

%% 1. 15-min data set
% time: nanoseconds since 2023-02-15 00:00 UTC on a 15-min grid (35041 steps)
time_ns_15min = ncread(ncFile, 'time');
t_15min = datetime(time_ns_15min, 'ConvertFrom', 'epochtime', ...
    'Epoch', datetime(2023,2,15,0,0,0), 'TicksPerSecond', 1e9);

% Cloud/haze class at Mt. Soledad: 0 = neither, 1 = cloud, 2 = haze
clrCldHaz_class = ncread(ncFile, 'ClrCldHaz_class');

% FM-120 droplet size distribution: dN per bin [cm^-3], 30 bins x time
FM120_DSD_cm3 = ncread(ncFile, 'FM120_DSD');
FM120_diam_um = ncread(ncFile, 'FM120_diam');          % bin UPPER edge [um]: 3, 4, ..., 50 (the 2-3 um bin is 3)
FM120_Nd_cm3  = ncread(ncFile, 'FM120_Nd');            % total droplet number [cm^-3]
FM120_LWC_gm3 = ncread(ncFile, 'FM120_LWC');           % liquid water content [g m^-3]

% Ceilometer cloud-base height (Scripps Pier, M1), cleaned [m]
ARM_CBH_clean_m = ncread(ncFile, 'ARM_CBH_clean');

% AMS non-refractory mass, CVI enhancement-factor corrected [ug m^-3], 5 species x time
AMS_mass_ugm3    = ncread(ncFile, 'AMS_mass');
AMS_massfrac     = ncread(ncFile, 'AMS_massfrac');      % mass fraction [-]
AMS_species_name = cellstr(ncread(ncFile, 'AMS_species'));   % Nitrate, Sulfate, Organics, Ammonium, Chloride

% FM-120 bin edges [um] (readme_FM120.md: 1-um bins from 2 to 14 um, 2-um bins to 50 um)
FM120_edges_um = [2:14, 16:2:50];
% Bin width in log space. NOTE: natural log, as in v1 (see header).
dlogDp = log(FM120_edges_um(2:end) ./ FM120_edges_um(1:end-1))';     % 30 x 1

%% 2. Critical diameters
% 15-min mean wet critical diameter on the same 15-min grid [nm -> um]
dcrit15 = readtable(dcrit15File);
Dcrit_wet_um_15min = dcrit15.Dta_critical_wet_fKDd_mean / 1e3;
assert(height(dcrit15) == numel(t_15min), 'ET_Dcrit_15min.csv is not on the 15-min grid of %s', ncFile)

% Native-time particle diameter (Dta) and critical diameter (Dta_critical), both [nm]
dcritNat = readtable(dcritNatFile);
nativetime    = dcritNat.time;               % undefined in v1 after `clear all`
Dta_nm        = dcritNat.Dta;
Dta_crit_nm   = dcritNat.Dta_critical;
classNative   = dcritNat.ClrCldHaz_class;    % 1 = cloud, 2 = haze

%% 3. FM-120 5-min rows for the D_eff panel (verification-table rows 51, 60)
% Daily FM-120 files: 5-min averages of 1-s data, each row labelled by its interval
% START (UTC). Only rows dated on the file's own day are kept (some files end with
% 00:00 of the next day, which readme_FM120.md says to ignore).
[t_5min, Nd_5min_cm3, ED_5min_um] = readFM120rows(fm120Dir);
% The FM-120 ED column is the instrument's effective diameter (bin mid-points); it is
% undefined when no droplet was counted (the files write 0), so N = 0 -> NaN.
ED_5min_um(~(Nd_5min_cm3 > 0)) = NaN;

% Class of each 5-min row = ClrCldHaz_class of the 15-min interval [t, t + 15 min)
% holding it (15-min intervals are labelled by their START; the file's time stamps
% carry a few ns of float jitter, removed by rounding to whole seconds).
t_15min_s = datetime(2023,2,15,0,0,0) + seconds(round(double(time_ns_15min) / 1e9));
t_5min_floor15 = dateshift(t_5min, 'start', 'hour') + minutes(15 * floor(minute(t_5min) / 15));
[found, loc] = ismember(t_5min_floor15, t_15min_s);
class_5min = nan(size(t_5min));
class_5min(found) = clrCldHaz_class(loc(found));

%% 4. Per-season quantities
% Everything the figure needs is computed here and stored in the struct array S,
% so the plotting block below only draws.
mo_15min  = month(t_15min);
mo_native = month(nativetime);
inSeason  = @(mo, s) mo >= seasonMonths(s,1) & mo <= seasonMonths(s,2);

for s = 1:4
    % 15-min cloud and haze samples in this season
    iCld = find(clrCldHaz_class == 1 & inSeason(mo_15min, s));
    iHaz = find(clrCldHaz_class == 2 & inSeason(mo_15min, s));
    S(s).iCld = iCld;                                        %#ok<SAGROW>
    S(s).iHaz = iHaz;

    % Droplet size distribution: median and 25th/75th percentiles of dN/dlogDp
    [S(s).dsdCld_med, S(s).dsdCld_p25, S(s).dsdCld_p75] = dsdPercentiles(FM120_DSD_cm3(:, iCld), dlogDp);
    [S(s).dsdHaz_med, S(s).dsdHaz_p25, S(s).dsdHaz_p75] = dsdPercentiles(FM120_DSD_cm3(:, iHaz), dlogDp);

    % Effective diameter for the D_eff panel (rows 51, 60): FM-120 ED of the 5-min rows
    % of each class, sheet seasons. ED = sum(n D^3) / sum(n D^2) with bin mid-points
    % (Hansen & Travis, 1974, Space Sci. Rev. 16, 527-610, eq. 2.53, for r_eff).
    S(s).DeffCld_um = ED_5min_um(class_5min == 1 & inSheetSeason(t_5min, seasonName{s}));
    S(s).DeffHaz_um = ED_5min_um(class_5min == 2 & inSheetSeason(t_5min, seasonName{s}));
    % Kavin's D_eff (version 2a), kept for the printed comparison only: each 15-min
    % FM120_DSD with every bin's UPPER edge as its diameter
    S(s).DeffCld_upper_um = effectiveDiameter(FM120_DSD_cm3(:, iCld), FM120_diam_um);
    S(s).DeffHaz_upper_um = effectiveDiameter(FM120_DSD_cm3(:, iHaz), FM120_diam_um);

    % 15-min wet critical diameter
    S(s).DcritCld_um = Dcrit_wet_um_15min(iCld);
    S(s).DcritHaz_um = Dcrit_wet_um_15min(iHaz);

    % Percent of native-time particles with Dta > Dta_critical (NaN when no data)
    S(s).pctActCld = percentActivated(Dta_nm, Dta_crit_nm, classNative == 1 & inSeason(mo_native, s));
    S(s).pctActHaz = percentActivated(Dta_nm, Dta_crit_nm, classNative == 2 & inSeason(mo_native, s));

    % AMS composition for the pies: mean concentration of each species over the samples
    % [ug m^-3] = table rows 35-42 (+ chloride). v1 used sums, see header.
    S(s).amsCld_ugm3 = mean(AMS_mass_ugm3(:, iCld), 2, 'omitmissing');
    S(s).amsHaz_ugm3 = mean(AMS_mass_ugm3(:, iHaz), 2, 'omitmissing');

    % Cloud-base height histograms [hours per 50-m bin], CBH < 2 km only
    % columns: other low clouds (class 0), cloud (1), haze (2)
    isLow = ARM_CBH_clean_m < 2000 & inSeason(mo_15min, s);
    S(s).cbhHours = [cbhHistogramHours(ARM_CBH_clean_m(isLow & clrCldHaz_class == 0)), ...
                     cbhHistogramHours(ARM_CBH_clean_m(isLow & clrCldHaz_class == 1)), ...
                     cbhHistogramHours(ARM_CBH_clean_m(isLow & clrCldHaz_class == 2))];
end

%% 4b. Numbers printed on the figure, from the verification table
% Rows of "Table-EPCAPE_SeasonalAverages - AJB verification - ver 2.xlsx"; each row's
% quantity, author and method are checked, so an inserted or moved row stops the script.
%   66 / 70  FM-120 cloud / haze sampling hours (BAMS-paper definitions) -> hour labels
%   53 / 62  cloud / haze activated fraction (Kavin's count re-run)       -> % Activated
%   51 / 60  cloud / haze effective diameter (5-min rows, ED column)      -> checked below
rowCldHours = readVerificationRow(verTableFile, verTableSheet, 66, ...
    "On-Site Cloud Sampling by Fog Drop Monitor", "filters: BAMS-paper definition");
rowHazHours = readVerificationRow(verTableFile, verTableSheet, 70, ...
    "On-Site Haze Sampling by Fog Drop Monitor", "filters: BAMS-paper definition");
rowCldAct = readVerificationRow(verTableFile, verTableSheet, 53, ...
    "Cloud Droplet Activated Fraction", "Kavin's count re-run");
rowHazAct = readVerificationRow(verTableFile, verTableSheet, 62, ...
    "Haze Droplet Activated Fraction", "Kavin's count re-run");
rowCldDeff = readVerificationRow(verTableFile, verTableSheet, 51, ...
    "Cloud Droplet Effective Diameter", "filters: FM-120 5-min rows");
rowHazDeff = readVerificationRow(verTableFile, verTableSheet, 60, ...
    "Haze Droplet Effective Diameter", "filters: FM-120 5-min rows");
for s = 1:4
    S(s).cloudHours_table = rowCldHours.mean.(seasonName{s});     % [h]
    S(s).hazeHours_table  = rowHazHours.mean.(seasonName{s});     % [h]
    S(s).fracActCld_table = rowCldAct.mean.(seasonName{s});       % fraction
    S(s).fracActHaz_table = rowHazAct.mean.(seasonName{s});       % fraction
end

fprintf('\nNumbers printed on the figure, from %s\n', verTableFile)
fprintf('%-7s %-30s %10s %8s\n', 'Season', 'Quantity (table row)', 'table', 'printed')
for s = 1:4
    fprintf('%-7s %-30s %10.4g %8d\n', seasonName{s}, 'cloud hours (66)', S(s).cloudHours_table, round(S(s).cloudHours_table))
    fprintf('%-7s %-30s %10.4g %8d\n', seasonName{s}, 'haze hours (70)',  S(s).hazeHours_table,  round(S(s).hazeHours_table))
    fprintf('%-7s %-30s %10.4g %8d   (this script: %g)\n', seasonName{s}, 'cloud % activated (53)', ...
        100*S(s).fracActCld_table, round(100*S(s).fracActCld_table), S(s).pctActCld)
    fprintf('%-7s %-30s %10.4g %8d   (this script: %g)\n', seasonName{s}, 'haze % activated (62)', ...
        100*S(s).fracActHaz_table, round(100*S(s).fracActHaz_table), S(s).pctActHaz)
end

% D_eff panel data next to rows 51 and 60 (mean +- std over the 5-min rows)
fprintf('\nD_eff panel [um]: this script vs. verification-table rows 51 (cloud) and 60 (haze)\n')
fprintf('%-7s %-6s %18s %16s %8s\n', 'Season', 'Class', 'this script', 'table', 'n rows')
allSeasons = ['EPCAPE', seasonName];
className = {'cloud', 'haze'};
rowDeff = [rowCldDeff, rowHazDeff];
for s = 1:numel(allSeasons)
    for c = 1:2
        x = ED_5min_um(class_5min == c & inSheetSeason(t_5min, allSeasons{s}));
        tbl = rowDeff;
        fprintf('%-7s %-6s %9.2f +- %5.2f %8.1f +- %4.1f %8d\n', allSeasons{s}, className{c}, ...
            mean(x, 'omitmissing'), std(x, 'omitmissing'), tbl(c).mean.(allSeasons{s}), tbl(c).std.(allSeasons{s}), sum(~isnan(x)))
    end
end

%% 5. Figure layout (all sizes in points)
% Panels keep the v1 sizes; only the gaps change (Lynn: less space between
% panels, more space below the "Droplet Diameter" label).
L.colW     = 177.5;   % panel width (v1)
L.colGap   = 16;      % between seasons (v1: 24); must fit the 2-digit % activated labels
L.dsdH     = 168;     % size-distribution panel height (v1)
L.boxH     = 36;      % D_crit and D_eff panel height (v1)
L.rowGap   = 6;       % DSD -> D_crit -> D_eff (v1: 31)
L.labelGap = 46;      % D_eff axis -> CBH panel, holds tick labels + x label + ~16 pt of space (v1: 30)
L.cbhH     = 168;     % cloud-base-height panel height (v1)
L.left     = 52;      % room for y tick labels + y labels
L.bottom   = 36;      % room for "Hours" tick labels + label
L.top      = 22;      % room for season titles
L.legendW  = 100;     % room right of Winter for the composition legend

L.W = L.left + 4*L.colW + 3*L.colGap + L.legendW;
L.H = L.bottom + L.cbhH + L.labelGap + 2*L.boxH + 2*L.rowGap + L.dsdH + L.top;

y.cbh   = L.bottom;
y.deff  = y.cbh   + L.cbhH + L.labelGap;
y.dcrit = y.deff  + L.boxH + L.rowGap;
y.dsd   = y.dcrit + L.boxH + L.rowGap;
x0      = L.left + (0:3) * (L.colW + L.colGap);       % left edge of each season column

% points -> normalized figure units (normalized keeps the layout intact when the
% Live Editor rescales the figure)
toNorm = @(p) [p(1)/L.W, p(2)/L.H, p(3)/L.W, p(4)/L.H];

% AMS pie boxes, relative to the lower-left corner of the size-distribution panel.
% Same size and approximately the same place as in v1 (v1 axes [x y 0.07 0.07] of a 975 x 613 pt figure).
pie.dx    = 88;                    % [pt] from the panel's left edge
pie.wh    = [68.25 42.9];          % [pt] box size; the pie diameter is the smaller side
pie.dyTop = 85.5;                  % [pt] upper pie (under the "Cloud" label at y = 1180)
pie.dyBot = 24.2;                  % [pt] lower pie (under the "Haze" label at y = 630)

% Colours
blue = [0 0 1];  red = [1 0 0];
cbhColors = [0.369 0.529 0.388;    % other low clouds
             0.643 0.663 0.941;    % cloud
             0.941 0.451 0.482];   % haze
% AMS standard species colours (nitrate, sulfate, organics, ammonium, chloride). v1 drew
% them with piechart's default FaceAlpha = 0.6, i.e. blended 60:40 with white. To make
% ONLY the organics brighter (L2), wedges are drawn opaque, the other four keep their
% v1 appearance by pre-blending with white, and organics is pure green.
amsBaseColors = [0 0 1; 1 0 0; 0 1 0; 1 0.647 0; 1 0 0.82];
v1PieAlpha = 0.6;
pieColors  = 1 - v1PieAlpha * (1 - amsBaseColors);   % = alpha*c + (1-alpha)*white
pieColors(3,:) = [0 1 0];                             % organics, like Lynn's example (#00FF00)
pieEdgeColor = [0.15 0.15 0.15];                      % piechart default edge colour

% Diameter axis shared by the three upper rows
Dlim_um   = [3 50];
Dticks_um = [3 10 20 30 40 50];

%% 6. Draw the figure
fig = figure('Color', 'w', 'Units', 'points', 'Position', [20 20 L.W L.H], 'Visible', figureVisible);
% Figure-level defaults = MATLAB factory values, so user startup.m settings
% (e.g. dark themes, large fonts) do not change this figure.
set(fig, 'DefaultAxesFontSize', 10, 'DefaultAxesFontWeight', 'normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1.1, 'DefaultAxesLabelFontSizeMultiplier', 1.1, ...
    'DefaultAxesTitleFontWeight', 'bold', 'DefaultAxesColor', 'w', ...
    'DefaultAxesXColor', [0.15 0.15 0.15], 'DefaultAxesYColor', [0.15 0.15 0.15], ...
    'DefaultAxesZColor', [0.15 0.15 0.15], 'DefaultAxesGridColor', [0.15 0.15 0.15], ...
    'DefaultTextColor', 'k', 'DefaultTextFontSize', 10, ...
    'DefaultLineLineWidth', 0.5, 'DefaultLineMarkerSize', 6);

for s = 1:4
    isFirstCol = (s == 1);

    % ---------- Row 1: droplet size distributions + AMS pies ----------
    ax = axes(fig, 'Units', 'normalized', 'Position', toNorm([x0(s) y.dsd L.colW L.dsdH]));
    shadedErrorBar(FM120_diam_um, S(s).dsdCld_med, ...
        bandErrors(S(s).dsdCld_med, S(s).dsdCld_p25, S(s).dsdCld_p75, applyV1BugFixes), ...
        'lineprops', '-b', 'plotAxes', ax);
    hold(ax, 'on')
    shadedErrorBar(FM120_diam_um, S(s).dsdHaz_med, ...
        bandErrors(S(s).dsdHaz_med, S(s).dsdHaz_p25, S(s).dsdHaz_p75, applyV1BugFixes), ...
        'lineprops', '-r', 'plotAxes', ax);
    hold(ax, 'off')
    set(ax, 'XScale', 'log', 'XLim', Dlim_um, 'XTick', Dticks_um, 'XTickLabel', [], ...
        'YLim', [0 1400], 'Box', 'on', 'Layer', 'top')      % L4: closed box
    if isFirstCol
        ylabel(ax, 'dN/dlog(D_p) [cm^{-3}]')
    else
        ax.YTickLabel = [];
    end
    title(ax, seasonName{s})
    % hours from verification-table rows 66 / 70 (section 4b), rounded to whole hours
    text(ax, 10, 1180, "Cloud: " + round(S(s).cloudHours_table) + " Hours", 'Color', blue, 'FontSize', 7)
    text(ax, 10,  630, "Haze: "  + round(S(s).hazeHours_table)  + " Hours", 'Color', red,  'FontSize', 7)

    % Pies (B2: v1 had the haze pie on top, under the "Cloud" label)
    if applyV1BugFixes
        dyCld = pie.dyTop;  dyHaz = pie.dyBot;
    else
        dyCld = pie.dyBot;  dyHaz = pie.dyTop;
    end
    drawPie(fig, toNorm([x0(s)+pie.dx, y.dsd+dyCld, pie.wh]), S(s).amsCld_ugm3, pieColors, pieEdgeColor);
    drawPie(fig, toNorm([x0(s)+pie.dx, y.dsd+dyHaz, pie.wh]), S(s).amsHaz_ugm3, pieColors, pieEdgeColor);

    % ---------- Row 2: wet critical diameter, labelled with % activated ----------
    posDcrit = toNorm([x0(s) y.dcrit L.colW L.boxH]);
    ax = axes(fig, 'Units', 'normalized', 'Position', posDcrit);
    plotDcritBoxes(ax, S(s), s, applyV1BugFixes);
    % boxplot resizes the axes and switches it to 'outerposition', which would let the
    % tick labels shrink it further; pin the planned inner position again
    set(ax, 'PositionConstraint', 'innerposition', 'Position', posDcrit, ...
        'XScale', 'log', 'XLim', Dlim_um, 'XTick', Dticks_um, 'XTickLabel', [], ...
        'YLim', [0 100], 'Box', 'on')
    text(ax, 3, 85, 'D_{\it crit}')
    if isFirstCol, ylabel(ax, '% Activated'), end

    % ---------- Row 3: effective diameter ----------
    posDeff = toNorm([x0(s) y.deff L.colW L.boxH]);
    ax = axes(fig, 'Units', 'normalized', 'Position', posDeff);
    % FM-120 5-min rows (verification-table rows 51, 60); haze at y = 1, cloud at y = 2
    boxplot(ax, padColumns(S(s).DeffHaz_um, S(s).DeffCld_um), 'Orientation', 'horizontal', ...
        'Colors', 'rb', 'Symbol', outlierSymbol(applyV1BugFixes));
    set(ax, 'PositionConstraint', 'innerposition', 'Position', posDeff, ...   % see D_crit row
        'XScale', 'log', 'XLim', Dlim_um, 'XTick', Dticks_um, 'XTickLabelRotation', 0, ...
        'YTick', [], 'Box', 'on')
    % trailing space: the Live Editor clips the hook of a final italic 'f' ("eff" -> "efl")
    text(ax, 3, 2.20, 'D_{\it eff }')
    xlabel(ax, 'Droplet Diameter [\mum]')

    % ---------- Row 4: cloud-base height histogram ----------
    ax = axes(fig, 'Units', 'normalized', 'Position', toNorm([x0(s) y.cbh L.colW L.cbhH]));
    cbhCenters_m = 25:50:475;                        % L1: bins now start at 0 m
    bh = barh(ax, cbhCenters_m, S(s).cbhHours, 'stacked');
    for k = 1:3, bh(k).FaceColor = cbhColors(k,:); end
    hold(ax, 'on')
    plot(ax, [0 120], [250 250], 'k-')               % Mt. Soledad height (~250 m)
    hold(ax, 'off')
    set(ax, 'XLim', [0 120], 'YLim', [0 500], 'Box', 'on', ...
        'YTick', 0:50:500, 'YTickLabel', {'0','','100','','200','','300','','400','','500'})
    xlabel(ax, 'Hours')
    if isFirstCol
        ylabel(ax, 'Cloud Base Height [m]')
    else
        ax.YTickLabel = [];
    end

    if s == 4
        legend(bh, {'Other Low Clouds', 'Cloud', 'Haze'}, 'Location', 'northeast', ...
            'Color', 'w', 'TextColor', 'k')
        text(ax, 65, 265, 'Soledad Height', 'FontSize', 7)
        % L1: arrow points down from the 250 m line into the region it describes
        arrowX_h = 55;  arrowY_m = [243 135];
        [ax1, ay1] = data2norm(ax, arrowX_h, arrowY_m(1));
        [ax2, ay2] = data2norm(ax, arrowX_h, arrowY_m(2));
        an = annotation(fig, 'arrow', [ax1 ax2], [ay1 ay2]);
        an.LineWidth = 1;  an.HeadWidth = 7;  an.HeadLength = 7;
        text(ax, arrowX_h + 5, mean(arrowY_m), {'CBH < 250 m is', 'below Soledad peak'}, ...
            'FontSize', 7, 'VerticalAlignment', 'middle')
    end
end

% Composition legend (L3), right of the Winter size-distribution panel. Drawn from
% patches with the pie colours because piechart's built-in legend cannot be placed.
axLeg = axes(fig, 'Units', 'normalized', 'Position', toNorm([x0(4)+L.colW+6, y.dsd, 10, 10]), ...
    'Visible', 'off');
hp = gobjects(1, 5);
for k = 1:5
    hp(k) = patch(axLeg, NaN, NaN, pieColors(k,:), 'EdgeColor', pieEdgeColor);
end
lgd = legend(hp, AMS_species_name, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.15 0.15 0.15]);
lgd.Units = 'normalized';
lgd.Position(1:2) = [(x0(4)+L.colW+6)/L.W, (y.dsd+L.dsdH)/L.H - lgd.Position(4)];   % top-aligned with the panel

%% 7. Export
if saveFigure
    exportgraphics(fig, fullfile(outDir, [figName '.png']), 'Resolution', 300)
    exportgraphics(fig, fullfile(outDir, [figName '.pdf']), 'ContentType', 'vector')
    fprintf('Saved %s.png / .pdf in %s\n', figName, outDir)
end

%% 8. Seasonal statistics (the numbers printed throughout v1)
% Hrs(15m) = 15-min class hours (Kavin); Hrs(label) = verification-table rows 66/70;
% Deff(panel) = 5-min ED rows (rows 51/60); Deff(2a) = Kavin's 15-min upper-edge D_eff
fprintf('\n%-7s %-6s %9s %10s %18s %18s %16s %16s %7s\n', 'Season', 'Class', 'Hrs(15m)', 'Hrs(label)', ...
    'Nd [cm^-3]', 'LWC [g m^-3]', 'Deff(panel) [um]', 'Deff(2a) [um]', '%Activ')
for s = 1:4
    cls = {'cloud', 'haze'};
    idx = {S(s).iCld, S(s).iHaz};
    hrsLabel = [S(s).cloudHours_table, S(s).hazeHours_table];
    deff = {S(s).DeffCld_um, S(s).DeffHaz_um};
    deffUpper = {S(s).DeffCld_upper_um, S(s).DeffHaz_upper_um};
    act  = [S(s).pctActCld, S(s).pctActHaz];
    for c = 1:2
        i = idx{c};
        fprintf('%-7s %-6s %9.2f %10.4g %8.2f +- %6.2f %8.4f +- %6.4f %7.2f +- %5.2f %7.2f +- %5.2f %7.0f\n', ...
            seasonName{s}, cls{c}, numel(i)/4, hrsLabel(c), ...
            mean(FM120_Nd_cm3(i), 'omitmissing'),  std(FM120_Nd_cm3(i), 'omitmissing'), ...
            mean(FM120_LWC_gm3(i), 'omitmissing'), std(FM120_LWC_gm3(i), 'omitmissing'), ...
            mean(deff{c}, 'omitmissing'), std(deff{c}, 'omitmissing'), ...
            mean(deffUpper{c}, 'omitmissing'), std(deffUpper{c}, 'omitmissing'), act(c))
    end
end
fprintf('\nAMS mass fraction, mean +- std over 15-min samples (%s)\n', strjoin(AMS_species_name, ', '))
for s = 1:4
    fprintf('%-7s cloud: %s +- %s\n', seasonName{s}, mat2str(mean(AMS_massfrac(:, S(s).iCld), 2, 'omitmissing')', 3), ...
        mat2str(std(AMS_massfrac(:, S(s).iCld), [], 2, 'omitmissing')', 3))
    fprintf('%-7s haze : %s +- %s\n', seasonName{s}, mat2str(mean(AMS_massfrac(:, S(s).iHaz), 2, 'omitmissing')', 3), ...
        mat2str(std(AMS_massfrac(:, S(s).iHaz), [], 2, 'omitmissing')', 3))
end

%% Local functions

function [t_5min, Nd_cm3, ED_um] = readFM120rows(fm120Dir)
% Daily FM-120 5-min rows (UTC, interval START), as check_seasonal_averages.ipynb
% reads them (sources.load_fm120): only rows dated on the file's own day are kept
% (some files end with 00:00 of the next day, which readme_FM120.md says to ignore),
% empty files are skipped, and a repeated time keeps its first row.
%   Nd_cm3  "Number Conc (#/cm^3)", droplets 2-50 um [cm^-3]   (column 4)
%   ED_um   "ED (um)", the instrument's effective diameter [um] (column 7)
% Import options come from the largest (full-day) file and are reused for every file,
% because readtable's auto-detection misreads files with only 1-2 data rows. The
% bin-column headers differ between files, but columns 1-7 must match.
files = dir(fullfile(fm120Dir, '20*.csv'));
[~, iRef] = max([files.bytes]);
refFile = fullfile(files(iRef).folder, files(iRef).name);
opts = detectImportOptions(refFile);
opts.SelectedVariableNames = opts.VariableNames([1 4 7]);
refHeader = leadingFields(refFile, 7);
t_5min = cell(numel(files), 1);  Nd_cm3 = cell(numel(files), 1);  ED_um = cell(numel(files), 1);
for j = 1:numel(files)
    f = fullfile(files(j).folder, files(j).name);
    tbl = readtable(f, opts);
    if height(tbl) == 0, continue, end              % e.g. 2024-01-02.csv: header only
    if ~isequal(leadingFields(f, 7), refHeader)
        error('FM-120 file %s has a different column layout than %s', files(j).name, files(iRef).name)
    end
    day = datetime(files(j).name(1:10), 'InputFormat', 'yyyy-MM-dd');
    ownDay = tbl{:, 1} >= day & tbl{:, 1} < day + days(1);
    t_5min{j} = tbl{ownDay, 1};
    Nd_cm3{j} = tbl{ownDay, 2};
    ED_um{j}  = tbl{ownDay, 3};
end
t_5min = vertcat(t_5min{:});  Nd_cm3 = vertcat(Nd_cm3{:});  ED_um = vertcat(ED_um{:});
[t_5min, iFirst] = unique(t_5min, 'first');         % sorted; a repeated time keeps its first row
Nd_cm3 = Nd_cm3(iFirst);
ED_um  = ED_um(iFirst);
end

function fields = leadingFields(fileName, n)
% First n comma-separated fields of the header line of a csv file (trailing spaces removed).
fid = fopen(fileName, 'r');
header = fgetl(fid);
fclose(fid);
fields = strtrim(strsplit(header, ','));
fields = fields(1:n);
end

function tf = inSheetSeason(t, season)
% Seasons of the seasonal-averages sheet (its column headers), whole UTC days:
% EPCAPE 2023-02-15..2024-02-14, Spring Apr-Jun 2023, Summer Jul-Sep 2023,
% Fall Oct-Dec 2023, Winter 2023-02-15..03-31 + 2024-01-01..02-14.
switch season
    case 'EPCAPE', r = [datetime(2023,2,15), datetime(2024,2,15)];
    case 'Spring', r = [datetime(2023,4,1),  datetime(2023,7,1)];
    case 'Summer', r = [datetime(2023,7,1),  datetime(2023,10,1)];
    case 'Fall',   r = [datetime(2023,10,1), datetime(2024,1,1)];
    case 'Winter', r = [datetime(2023,2,15), datetime(2023,4,1); datetime(2024,1,1), datetime(2024,2,15)];
    otherwise, error('Unknown season %s', season)
end
tf = false(size(t));
for k = 1:size(r, 1)
    tf = tf | (t >= r(k, 1) & t < r(k, 2));
end
end

function row = readVerificationRow(file, sheet, rowNumber, quantity, notesStart)
% One row of Andrew's verification table: quantity (column A), author (J) and notes (K)
% are checked, then columns E-I (EPCAPE, Spring, Summer, Fall, Winter) are parsed from
% text such as "12.7 +- 4.3", "137" or "n/a" into row.mean.<season> and row.std.<season>.
c = readcell(file, 'Sheet', sheet, 'Range', sprintf('A%d:K%d', rowNumber, rowNumber));
c(cellfun(@(x) isa(x, 'missing'), c)) = {''};
if ~strcmp(strtrim(string(c{1})), quantity) || ~strcmp(strtrim(string(c{10})), "Andrew Buggee") ...
        || ~startsWith(strtrim(string(c{11})), notesStart)
    error('Row %d of %s is not "%s" (Andrew Buggee, notes "%s..."); found "%s" / "%s".', ...
        rowNumber, file, quantity, notesStart, string(c{1}), string(c{10}))
end
seasons = {'EPCAPE', 'Spring', 'Summer', 'Fall', 'Winter'};
for k = 1:5
    [row.mean.(seasons{k}), row.std.(seasons{k})] = parseMeanStd(c{4 + k});
end
end

function [m, sd] = parseMeanStd(cellValue)
% "12.7 +- 4.3" -> 12.7, 4.3 ; "137" or 137 -> 137, NaN ; "n/a" or empty -> NaN, NaN
m = NaN;  sd = NaN;
if isnumeric(cellValue)
    m = double(cellValue);
    return
end
parts = strtrim(split(string(cellValue), char(177)));    % char(177) = plus-minus sign
m = str2double(parts(1));
if numel(parts) > 1, sd = str2double(parts(2)); end
end

function [med, p25, p75] = dsdPercentiles(dN_cm3, dlogDp)
% Median, 25th and 75th percentile of dN/dlogDp in each bin across samples [cm^-3].
% dN_cm3: 30 bins x samples, dlogDp: 30 x 1
x = dN_cm3 ./ dlogDp;
med = median(x, 2, 'omitmissing');
p25 = prctile(x, 25, 2);
p75 = prctile(x, 75, 2);
end

function err = bandErrors(med, p25, p75, fixBand)
% Error input for shadedErrorBar, which draws [y - err(2,:), y + err(1,:)].
if fixBand
    err = [p75 - med, med - p25]';     % band = [p25, p75]  (B1)
else
    err = [p75, p25]';                 % v1: band = [med - p25, med + p75]
end
end

function Deff_um = effectiveDiameter(dN_cm3, D_um)
% Effective diameter of each DSD (column) [um]: sum(n D^3) / sum(n D^2).
% NaN when the DSD is empty or missing.
Deff_um = sum(dN_cm3 .* D_um.^3, 1) ./ sum(dN_cm3 .* D_um.^2, 1);
Deff_um = Deff_um(:);
end

function pct = percentActivated(Dta_nm, Dcrit_nm, keep)
% Percent of samples with particle diameter above the critical diameter, rounded
% to an integer as in v1 (round(piechart Proportions * 100)). NaN when no samples.
k = keep & ~isnan(Dta_nm) & ~isnan(Dcrit_nm);
pct = round(100 * sum(Dta_nm(k) > Dcrit_nm(k)) / sum(k));
end

function hours = cbhHistogramHours(cbh_m)
% Hours of 15-min samples in 50-m cloud-base-height bins centred on 25:50:475 m.
% The 50-500 m bins reproduce v1 exactly (v1: histcounts(250 - CBH, -250:50:200),
% i.e. (a, b] bins plus the closed 50-100 m bin); the new 0-50 m bin is [0, 50).
nV1 = fliplr(histcounts(250 - cbh_m, -250:50:200));   % ascending CBH: [50,100], (100,150], ..., (450,500]
n0  = sum(cbh_m >= 0 & cbh_m < 50);
hours = [n0, nV1]' / 4;                                % 15-min samples -> hours, 10 x 1
end

function M = padColumns(a, b)
% Two vectors as the columns of one matrix, padding the shorter one with NaN (for boxplot).
n = max(numel(a), numel(b));
M = nan(n, 2);
M(1:numel(a), 1) = a(:);
M(1:numel(b), 2) = b(:);
end

function sym = outlierSymbol(fix)
% '+' alone draws outliers in the box colour; v1 used the default red '+'.
if fix, sym = '+'; else, sym = 'r+'; end
end

function plotDcritBoxes(ax, Ss, s, fixLabels)
% Wet critical diameter box plots. Haze at y = 33 (red), cloud at y = 66 (blue),
% each labelled with its % activated from the summary table (section 4b). A box is
% drawn only when its % activated exists (v1: Fall shows the cloud box only, at
% y = 56; Winter is empty).
pctCld = round(100 * Ss.fracActCld_table);
pctHaz = round(100 * Ss.fracActHaz_table);
hasCld = ~isnan(pctCld);
hasHaz = ~isnan(pctHaz);
if hasCld && hasHaz
    labels = {num2str(pctHaz), num2str(pctCld)};
    if ~fixLabels && s == 1
        labels = fliplr(labels);       % B3: v1 Spring labels were swapped
    end
    boxplot(ax, padColumns(Ss.DcritHaz_um, Ss.DcritCld_um), 'Orientation', 'horizontal', ...
        'Positions', [33 66], 'Widths', 15, 'Colors', 'rb', 'Labels', labels, ...
        'Symbol', outlierSymbol(fixLabels));
elseif hasCld
    boxplot(ax, Ss.DcritCld_um, 'Orientation', 'horizontal', 'Positions', 56, 'Widths', 15, ...
        'Colors', 'b', 'Labels', {num2str(pctCld)}, 'Symbol', outlierSymbol(fixLabels));
elseif hasHaz
    boxplot(ax, Ss.DcritHaz_um, 'Orientation', 'horizontal', 'Positions', 33, 'Widths', 15, ...
        'Colors', 'r', 'Labels', {num2str(pctHaz)}, 'Symbol', outlierSymbol(fixLabels));
else
    ax.YTick = [];
end
end

function p = drawPie(fig, posNorm, mass_ugm3, colors, edgeColor)
% AMS composition pie. As in v1, an axes is created at the target position and
% piechart takes over the current axes' position. The data are normalised to
% mass fractions first: like pie, piechart draws only a PARTIAL pie when the
% values sum to less than 1, which mean concentrations [ug m^-3] often do.
axes(fig, 'Units', 'normalized', 'Position', posNorm);
massFraction = mass_ugm3 / sum(mass_ugm3, 'omitmissing');
p = piechart(massFraction, 'LabelStyle', 'none', 'ColorOrder', colors, ...
    'FaceAlpha', 1, 'EdgeColor', edgeColor);
end

function [xn, yn] = data2norm(ax, x, y)
% Data coordinates of a linear axes -> normalized figure coordinates (for annotation).
p = ax.Position;      % normalized
xn = p(1) + (x - ax.XLim(1)) / diff(ax.XLim) * p(3);
yn = p(2) + (y - ax.YLim(1)) / diff(ax.YLim) * p(4);
end
