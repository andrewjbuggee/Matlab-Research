%% EPCAPE_BAMS_Fig4_AJB_v1  -- reproduces BAMS Figure 4 (Kavin's slide 4)
% Monthly time series at EPCAPE, Feb 2023 - Feb 2024 (13 campaign months; Feb
% 2023 = 15-28 Feb 2023 and Feb 2024 = 1-14 Feb 2024):
%   A  cloud-top, boundary-layer and cloud-base heights (notched boxes, one
%      colour per instrument/product; no whiskers)
%   B  ARSCL echo (hydrometeor) fraction profile, monthly mean (colour) with
%      ceilometer boundary-layer cloud base height (orange boxes)
%   C  MWR liquid water path (blue boxes)
%   D  upwind cloud vertical fraction, UCVF (blue boxes, left axis), with the
%      monthly mean relative humidity and "temperature as % of 20 C" (orange
%      lines, right axis)
%   E  SW cloud transmittance (blue boxes, right axis) and downwelling LW flux
%      (dark-red line, left axis)
%   F  rainfall from the two laser disdrometers (LD) and the video disdrometer
%      (VDIS). Default (rainPanel = "monthlyTotal"): monthly accumulation [mm]
%      as grouped bars, the same numbers as Die Wang's Fig. A3 (slide 12).
%      rainPanel = "slideMeanStd" gives the slide's version instead: boxes =
%      mean +/- 1 std of the DAILY totals on rain days (not quartiles).
% Output: EPCAPE_BAMS_Fig4_AJB_v1.png (300 dpi) and .pdf, and a table of
% every plotted variable, month by month, with the file and variable it
% comes from and the numbers drawn (section 9). The table is printed and
% saved as EPCAPE_BAMS_Fig4_AJB_v1_plotted_data.csv.
%
% WHY THIS FILE EXISTS
%   The code that made the slide's Figure 4 is NOT in Kavin's folder. The last
%   block of EPCAPEDanLWSWtransmittance.mlx ('echoyearly') is only a
%   calendar-month (Jan-Dec) precursor with no UCVF/RH row. This script
%   rebuilds the slide figure from the same input files, written from scratch.
%   The layout, colours, axis limits and labels follow the slide image
%   (EPCAPEoverviewFigs261005.pptx, slide 4, image10.png).
%
% DATA SOURCES (every number in the figure comes from one of these files)
%   A  EPCAPE_CTH_CBH_PBLH_SChellappan_LaRC.csv  GOES18_CTH, KAZR_CTH,
%        PBLH_THERMO, PBLH_MPL, CEIL_CBH (Chellappan, LaRC; 1 Mar 2023 -
%        10 Feb 2024, so there is NO Feb 2023 data for these five; Lynn's slide
%        comment "Missing Feb 23 data from Damao").
%      BL_Depth.csv (Depth), EPCAPE_sonde_LCL.csv (LCL)
%   B  arscl_echomask_hourly.nc (echo_fraction, height x time)
%      epcceil_hourly_cbh.nc (bl_cbh = hourly mean BL cloud base, 0-3000 m)
%   C  epcmwrret_turnM1_hourly_averages.nc (phys_lwp, g m^-2)
%   D  EPCAPE_15mindata.nc (UCVF, MET_rh, MET_temp)
%   E  <MON>_SW_Boxplot.csv (Dan's hourly SW transmittance, day x hour, PST)
%      radflux_M1_20230215_20240214.nc (ARM RADFLUX1LONG downwelling_longwave)
%   F  monthly_rain_accumulation_ver_3.csv (Die Wang's disdrometer code,
%        disdro_epc_ver_3.ipynb; ARM epcldquantsM1/S2 and epcvdisquantsM1).
%        Monthly total = sum of daily totals (rain_rate/60 per UTC day).
%      Slide version only: MEAN-EPCAPE Russel Rainfall.xlsx and
%        STD-Russel EPCAPE Rainfall.xlsx. Their 'mean' = monthly total /
%        number of rain days in the ver_3 CSV (all 36 values agree to 0.001
%        mm), i.e. mean daily rain on days with rain [mm/day]; 'std' = the
%        spread of those daily totals.
%
% DIFFERENCES FROM KAVIN'S SLIDE VERSION
%   Plain-language version for anyone: EPCAPE_BAMS_Fig4_changes_and_differences.txt
%   (same folder). Items 1-4 are data-handling problems in the slide version;
%   5-6 are requested changes; 7-10 are small choices.
%
%   1. ZERO PADDING -> boxes pulled toward 0 (panels B and C). Biggest visible
%      difference: on the slide the Feb 2023 and Feb 2024 CBH boxes reach the
%      surface (0 m). Kavin's code (precursor block of EPCAPEDanLWSWtransmittance.mlx)
%      builds each panel's monthly matrix as
%          A = nan(500,12);  for each month: A(1:n,month) = values;
%      Each month's column holds EVERY hour of the month, NaN hours included
%      (up to 744). When a column longer than 500 is written, MATLAB grows
%      the matrix and fills the new rows of all OTHER columns with 0, not
%      NaN. boxplot skips NaN but counts the 0s as real values. Every month
%      shorter than the longest month therefore gets fake zeros: 24 for each
%      30-day month and 244 for each half-February. A simulation of that code
%      on the 13 campaign months (AJB, 2026-10-14) reproduces the slide boxes:
%        CBH  Feb 2023: 176 real hours + 244 zeros -> Q1 = median = 0 m,
%             Q3 = 1219 m, whisker 2919 m (slide: 0 / ~1220 / ~2930).
%             Real data: Q1/median/Q3 = 1012/1280/1550 m.
%        CBH  Feb 2024: 219 real + 244 zeros -> 0 / 0 / 1214 m, whisker 2843 m.
%             Real: 939/1246/1650 m.
%        CBH  Nov 2023 (23 zeros among 187 real values): Q1 356 vs real 486 m.
%             Apr, Jun, Sep are lowered slightly.
%        LWP  same effect. Slide whisker tops Feb 2023 ~25, Oct ~33, Feb 2024
%             ~48 g m^-2 match the padded values 26.1, 33.3, 48.8. The real
%             values are 98.7, 72.2, 168.6 g m^-2.
%      Here each month's group holds only its real values (groupByMonth), so
%      no box is padded. The bl_cbh QC keeps 0 < CBH <= 3000 m, so a real 0 m
%      cloud base cannot occur.
%   2. UCVF boxes (panel D) are drawn at the month ticks. On the slide they sit
%      on a different x scale (13 boxes spread from just right of 'Feb 23' to
%      just right of 'Feb 24'), so they do not line up with the other panels.
%      The box values themselves match the 15-min UCVF monthly quartiles.
%   3. LW line (panel E). The slide line follows Dan's monthly LW numbers
%      (mean over the day of the hourly medians) for Mar 2023 - Jan 2024, but
%      Dan's FEB file pools Feb 2023 and Feb 2024, so it cannot give two
%      February points. The slide's Feb 2023 point (~350 W m^-2) matches no
%      data (real 296.7). Here the same statistic is computed for all 13
%      months from the 1-min ARM RADFLUX data (qc == 0). That recipe matches
%      Dan's files exactly for 8 months and within 0.3 W m^-2 for the rest
%      (EPCAPE_BAMS_Fig4_LW_validation_AJB.m). The LW axis starts at 290
%      instead of 300 so Feb 2023 shows.
%   4. SW boxes for the two Februaries (panel E). Dan's FEB file pools both
%      years; here it is split by row (assumption b). The slide's two Feb
%      boxes (medians ~0.5 and ~0.55) cannot be traced to a method; here the
%      medians are 0.70 (Feb 2023) and 0.68 (Feb 2024).
%   5. Rain (panel F). Slide: boxes of mean +/- 1 std of DAILY rain on rain
%      days (Russel spreadsheets). Now: monthly accumulation in mm (Lynn
%      asked for the amount per month, 2026-10-08), the same numbers as Fig.
%      A3, with 'x' = no data and '*' = more than half the days missing. The
%      old panel is kept: rainPanel = "slideMeanStd".
%   6. Lynn's slide comments (panel letters A-F, 'Hydrometeor Fraction') are
%      available but off by default: lynnEdits.
%   7. Panel A, Feb 2023: product/month boxes need at least
%      minSamplesPerBox = 10 samples. This drops 5-sample boxes of
%      GOES/KAZR/PBLH/CEIL that would otherwise appear on 28 Feb PST
%      (assumption a). The slide has no Feb 2023 boxes for these products.
%   8. Months are assigned in fixed PST (UTC-8). Kavin used DST-aware
%      America/Los_Angeles time. Only samples within 7-8 h of a month boundary
%      can change month.
%   9. RH line (panel D): the slide's Nov 2023 point (~63%) is not the monthly
%      mean (48.7%) or median (43.8%) of MET_rh. Here: monthly mean. The
%      slide's other RH points and the "% of 20 C" line match the means.
%  10. Legend typo 'KAZRCTH' -> 'KAZR CTH'. Whiskers and boxes are drawn by
%      the local function drawBoxes (same Tukey 1.5 IQR whiskers and prctile
%      quartiles as boxplot) instead of boxplot, so boxes can share axes with
%      the echo image and the lines.
%
% ASSUMPTIONS (check with Kavin)
%   a. Months are assigned in PST (UTC-8, no daylight saving), the same time
%      base as BAMS Fig. 5 v2. ARM netCDF times are UTC. The CSV times and the
%      15-min netCDF times carry no time zone; they are taken as UTC when
%      naiveTimesAreUTC = true (below). This only moves samples within 8 h of
%      a month boundary.
%   b. Dan's FEB_SW_Boxplot.csv has 28 rows. Matching each row against the
%      RADFLUX SW/clear-sky ratio day by day shows rows 3-14 = 17-28 Feb 2023
%      and rows 15-23 = 1-9 Feb 2024 (rows 1-2 and 24-28 are mostly clear days
%      and match several dates equally well). Rows 1-14 are therefore used as
%      Feb 2023 and rows 15-28 as Feb 2024.
%   c. (Resolved 2026-10-14) Rain units: the ver_3 CSV is in mm per month;
%      the Russel spreadsheets are mm per rain day (see DATA SOURCES, F).
%   d. The rain files call the second disdrometer 'LD/S1'; it is ARM facility
%      S2 (Mt. Soledad, epcldquantsS2; no S1 exists). Label kept as 'LD/S1'
%      to match Fig. A3.
%
% NOTE ON THE CAPTION (slide notes) -- the caption does not describe this figure:
%   * It says whiskers are the 10th/90th percentiles. Kavin's boxplot default
%     (and the slide; e.g. the Mar CBH whisker ends at 2680 m = the largest
%     value <= Q3 + 1.5 IQR) is Tukey: 1.5 x IQR. Set whiskerDef below.
%   * The slide's rain boxes were mean +/- 1 std, not 25th/75th percentiles.
%     With the monthly-total bars (default) the panel has no boxes; the
%     caption's "(e) rain accumulation" then fits.
%   * Panel letters/order in the caption ((f) RH and temperature, (g) UCVF)
%     do not match the slide's six panels.
%   * '% of 20 C' is temperature in deg C divided by 20 C. A ratio of Celsius
%     temperatures has no physical meaning (0 C is arbitrary); a separate
%     temperature axis would be clearer.
%
% Andrew J. Buggee, Scripps Institution of Oceanography, 2026-10-14

clear; close all

%% 1. Settings
fprintf('[%s] 1. Settings\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

% --- Input folders ---------------------------------------------------------
% Most inputs sit in Kavin's 'Annual Echo Plot' folder. The 15-min merged file
% and the RADFLUX file live in Andrew's Python project.
echoDir_andrew  = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS/Annual Echo Plot';
fifteenDir_andrew = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper';
radfluxDir_andrew = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed';
rainDir_andrew    = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed/Monthly Precip - 2DVD - data from Wang Die';
if isfolder(echoDir_andrew)
    echoDir    = echoDir_andrew;      % Andrew's machine
    fifteenDir = fifteenDir_andrew;
    radfluxDir = radfluxDir_andrew;
    rainDir    = rainDir_andrew;
else
    echoDir    = pwd;                 % anyone else: put every input file in
    fifteenDir = pwd;                 % the current folder
    radfluxDir = pwd;
    rainDir    = pwd;
end
% Outputs (PNG/PDF) go next to this script, in the 'Fig 4' folder; the input
% data stay in 'Annual Echo Plot' (Kavin's Fig. 5 code reads them too).
outDir = fileparts(mfilename('fullpath'));
figTag = 'EPCAPE_BAMS_Fig4_AJB_v1';      % base name of every output file
if isempty(outDir), outDir = pwd; end   % e.g. code pasted into the Command Window

% The figure is drawn and exported while INVISIBLE, then shown at the end
% (showFigureAtEnd). Drawing on a visible figure hung R2024b in -batch runs
% on 2026-10-14 (3 of 3 runs, no PNG written); the same code with a hidden
% figure finished in 13 s. Exports are identical either way.
showFigureAtEnd = ~batchStartupOptionUsed;   % no window in matlab -batch
exportFigs    = true;     % false = draw only, write no files

% --- Time handling ----------------------------------------------------------
% PST = fixed UTC-8 (no daylight saving), matching Dan's SW CSV hours and the
% BAMS Fig. 5 v2 time base. Used to assign every sample to a month.
tz = '-08:00';
% true  = the CSV times and the 15-min netCDF times are UTC (converted to PST)
% false = they are already local time (used as they are, like Kavin's code)
naiveTimesAreUTC = true;

% Campaign window [start, end) in PST. 13 months: Feb 2023 ... Feb 2024.
tStart = datetime(2023,2,15,0,0,0,'TimeZone',tz);
tEnd   = datetime(2024,2,15,0,0,0,'TimeZone',tz);
nMonth = 13;
monthFirstDay = datetime(2023,2,1) + calmonths(0:nMonth-1);   % for labels
monthLabels   = cellstr(string(monthFirstDay,'MMM yy'));      % 'Feb 23', ...

% --- Box style ---------------------------------------------------------------
% "tukey"  : whiskers reach the furthest point within 1.5 x IQR of the box
%            (MATLAB boxplot default; what the slide shows)
% "p10p90" : whiskers at the 10th and 90th percentiles (what the caption says)
whiskerDef = "tukey";

% Panel F. "monthlyTotal" = monthly rain accumulation [mm], grouped bars,
% matching Fig. A3 (Lynn asked for the amount per month, 2026-10-08).
% "slideMeanStd" = the slide's mean +/- std of daily rain on rain days.
rainPanel = "monthlyTotal";

% Lynn Russell's slide comments: panel letters A-F and 'Hydrometeor Fraction'
% in place of 'Echo Fraction'. false = the slide as Kavin made it.
lynnEdits = false;

% --- Colours (read off the slide / Kavin's precursor code) -----------------
orange   = [0.8500 0.3250 0.0980];   % CBH boxes, RH and T lines
blue     = [0      0.4470 0.7410];   % LWP, UCVF and SW boxes; LD/M1
yellow   = [0.9290 0.6940 0.1250];   % VDIS/M1
lwColor  = [0.6350 0.0780 0.1840];   % LW line (dark red)
heightColors = [0.129 0.408 0.910;   % GOES18 CTH   (Kavin's precursor colours)
                0.129 0.635 0.910;   % KAZR CTH
                0.129 0.910 0.400;   % PBLH THERMO
                0.098 0.702 0.224;   % PBLH MPL
                0.871 0.855 0.118;   % BL DEPTH
                0.702 0.098 0.400;   % SONDE LCL
                0.871 0.118 0.180];  % CEIL CBH
boxWidth = 0.5;                       % box width in month units (panels B-E)

% Panel A: a product/month box is drawn only if it has at least this many
% samples. With naiveTimesAreUTC = true, the first five LaRC samples (1 Mar
% 2023 01:30-07:30 UTC) fall on 28 Feb PST and would otherwise draw 5-sample
% Feb 2023 boxes for GOES/KAZR/PBLH/CEIL, which the slide (and Lynn's
% "Missing Feb 23 data" comment) do not have. Set to 1 to draw every box.
minSamplesPerBox = 10;

%% 2. Panel A inputs: heights from LaRC, BL depth and sonde LCL
fprintf('[%s] 2. Panel A inputs: heights from LaRC, BL depth and sonde LCL\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

larc   = readtable(fullfile(echoDir,'EPCAPE_CTH_CBH_PBLH_SChellappan_LaRC.csv'));
blDep  = readtable(fullfile(echoDir,'BL_Depth.csv'));
sonde  = readtable(fullfile(echoDir,'EPCAPE_sonde_LCL.csv'));

% Month index (1 = Feb 2023 ... 13 = Feb 2024) of every sample; NaN outside
% the campaign window.
mLarc  = monthIndex(naiveToZone(larc.Time,  tz, naiveTimesAreUTC), tStart, tEnd);
mBl    = monthIndex(naiveToZone(blDep.Time, tz, naiveTimesAreUTC), tStart, tEnd);
mSonde = monthIndex(naiveToZone(sonde.time, tz, naiveTimesAreUTC), tStart, tEnd);

% Stack all seven height products into one long vector with a group label so
% that boxchart can dodge them by colour within each month (Kavin's approach).
heightNames = ["GOES18 CTH","KAZR CTH","PBLH THERMO","PBLH MPL", ...
               "BL DEPTH","SONDE LCL","CEIL CBH"];
hVals  = [larc.GOES18_CTH; larc.KAZR_CTH; larc.PBLH_THERMO; larc.PBLH_MPL; ...
          blDep.Depth;     sonde.LCL;     larc.CEIL_CBH];                 % [m]
hMonth = [mLarc; mLarc; mLarc; mLarc; mBl; mSonde; mLarc];
hGroup = [repmat(heightNames(1),height(larc),1); repmat(heightNames(2),height(larc),1);
          repmat(heightNames(3),height(larc),1); repmat(heightNames(4),height(larc),1);
          repmat(heightNames(5),height(blDep),1); repmat(heightNames(6),height(sonde),1);
          repmat(heightNames(7),height(larc),1)];
keep   = ~isnan(hMonth) & ~isnan(hVals);
hVals = hVals(keep);  hMonth = hMonth(keep);
% Ordered categorical: fixes the colour/legend order (replaces Kavin's
% '\color{white}1' legend trick, which existed only to sort the groups).
hGroup = categorical(hGroup(keep), heightNames, 'Ordinal', true);
% Drop product/month combinations with fewer than minSamplesPerBox samples.
[gid, ~, ~] = findgroups(hGroup, hMonth);
nPerGroup = accumarray(gid, 1);
tooFew = nPerGroup(gid) < minSamplesPerBox;
if any(tooFew)
    fprintf('Panel A: dropped %d samples in product/month boxes with < %d samples\n', ...
            nnz(tooFew), minSamplesPerBox);
end
hVals = hVals(~tooFew);  hMonth = hMonth(~tooFew);  hGroup = hGroup(~tooFew);

%% 3. Panel B inputs: echo fraction profiles and ceilometer CBH
fprintf('[%s] 3. Panel B inputs: echo fraction profiles and ceilometer CBH\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

echoFile = fullfile(echoDir,'arscl_echomask_hourly.nc');
tEcho    = readNcTime(echoFile, tz);
zEcho_m  = double(ncread(echoFile,'height'));        % [m], 160-18010 m, 30 m step
echoFrac = double(ncread(echoFile,'echo_fraction')); % [height x time], 0-1
if size(echoFrac,1) ~= numel(zEcho_m)                % guard against orientation
    echoFrac = echoFrac.';
end
mEcho = monthIndex(tEcho, tStart, tEnd);
echoProfile = nan(numel(zEcho_m), nMonth);           % monthly mean profile
for k = 1:nMonth
    echoProfile(:,k) = mean(echoFrac(:, mEcho == k), 2, 'omitmissing');
end
clear echoFrac                                       % ~40 MB, no longer needed

ceilFile = fullfile(echoDir,'epcceil_hourly_cbh.nc');
mCeil    = monthIndex(readNcTime(ceilFile, tz), tStart, tEnd);
cbh_m    = double(ncread(ceilFile,'bl_cbh'));        % hourly mean BL CBH [m]

%% 4. Panel C inputs: MWR liquid water path
fprintf('[%s] 4. Panel C inputs: MWR liquid water path\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

mwrFile = fullfile(echoDir,'epcmwrret_turnM1_hourly_averages.nc');
mMwr    = monthIndex(readNcTime(mwrFile, tz), tStart, tEnd);
lwp_gm2 = double(ncread(mwrFile,'phys_lwp'));        % hourly mean LWP [g m^-2]

%% 5. Panel D inputs: UCVF, RH and temperature (15-min merged file)
fprintf('[%s] 5. Panel D inputs: UCVF, RH and temperature (15-min merged file)\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

f15  = fullfile(fifteenDir,'EPCAPE_15mindata.nc');
t15  = readNcTime(f15, tz, naiveTimesAreUTC);
m15  = monthIndex(t15, tStart, tEnd);
ucvf = double(ncread(f15,'UCVF'));                   % 0-1, 24-h trajectory mean
rh_pct = double(ncread(f15,'MET_rh'));               % relative humidity [%]
T_degC = double(ncread(f15,'MET_temp'));             % air temperature [deg C]
ucvf(ucvf < 0 | ucvf > 1) = NaN;                     % valid range 0-1 (fill -9999)

rhMean_pct  = nan(1,nMonth);
TMean_degC  = nan(1,nMonth);
for k = 1:nMonth
    rhMean_pct(k) = mean(rh_pct(m15 == k), 'omitmissing');
    TMean_degC(k) = mean(T_degC(m15 == k), 'omitmissing');
end
% Slide quantity "% of 20 C": monthly mean temperature in deg C / 20 deg C x 100.
TPct20_pct = TMean_degC / 20 * 100;

%% 6. Panel E inputs: SW transmittance (Dan's CSVs) and LW flux (RADFLUX)
fprintf('[%s] 6. Panel E inputs: SW transmittance (Dan''s CSVs) and LW flux (RADFLUX)\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

% Dan's <MON>_SW_Boxplot.csv: one row per day, columns 1-24 = PST hours 0-23,
% column 25 = 'Means' (not used). -0.05 is a night-time fill value.
% Campaign month index of each calendar-month file (Feb is split, see below).
swMonths = {'MAR',2; 'APR',3; 'MAY',4; 'JUN',5; 'JUL',6; 'AUG',7; ...
            'SEP',8; 'OCT',9; 'NOV',10; 'DEC',11; 'JAN',12};
swGroups = cell(1,nMonth);                           % all hourly values per month
for i = 1:size(swMonths,1)
    swGroups{swMonths{i,2}} = readSW(fullfile(echoDir,[swMonths{i,1} '_SW_Boxplot.csv']));
end
% FEB file pools both Februaries: rows 1-14 = 15-28 Feb 2023, rows 15-28 =
% 1-14 Feb 2024 (assumption b in the header).
[~, febSW] = readSW(fullfile(echoDir,'FEB_SW_Boxplot.csv'));
swGroups{1}      = reshape(febSW(1:14,:),[],1);
swGroups{nMonth} = reshape(febSW(15:end,:),[],1);

% LW: 1-min ARM RADFLUX downwelling LW (qc == 0 only) -> hourly means (PST)
% -> for each month, the median of each hour of the day -> mean of those 24
% medians. This is the statistic of Dan's LW_Stepplot files ('Median' column)
% that Kavin averaged for the precursor's LW line.
radFile = fullfile(radfluxDir,'radflux_M1_20230215_20240214.nc');
tRad    = readNcTime(radFile, tz);
lwRaw   = double(ncread(radFile,'downwelling_longwave'));       % [W m^-2]
qcLw    = ncread(radFile,'qc_downwelling_longwave');
lwRaw(qcLw ~= 0) = NaN;
lwHourly = retime(timetable(tRad, lwRaw), 'hourly', 'mean');    % hour-start stamps
mLw  = monthIndex(lwHourly.tRad, tStart, tEnd);
hLw  = hour(lwHourly.tRad);
lwMonthly_Wm2 = nan(1,nMonth);
for k = 1:nMonth
    hourlyMedian = nan(24,1);
    for h = 0:23
        hourlyMedian(h+1) = median(lwHourly.lwRaw(mLw == k & hLw == h), 'omitmissing');
    end
    lwMonthly_Wm2(k) = mean(hourlyMedian, 'omitmissing');
end

% For the printed comparison: Dan's monthly LW (mean of his hourly medians).
danLW_Wm2 = nan(1,nMonth);
for i = 1:size(swMonths,1)
    lwTbl = readmatrix(fullfile(echoDir,[swMonths{i,1} '_LW_Stepplot.csv']), 'NumHeaderLines',1);
    danLW_Wm2(swMonths{i,2}) = mean(lwTbl(:,4), 'omitmissing');  % column 4 = Median
end

%% 7. Panel F inputs: rainfall per instrument
fprintf('[%s] 7. Panel F inputs: rainfall per instrument\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

rainNames  = ["LD/M1","LD/S1","VDIS/M1"];

% Monthly totals (Fig. A3 numbers). One row per campaign month, 2023-02 ...
% 2024-02; columns <inst>_mm (total), <inst>_rain_days, <inst>_days_recorded,
% and days_in_period (14 for the two half Februaries).
accTbl = readtable(fullfile(rainDir,'monthly_rain_accumulation_ver_3.csv'), ...
                   'TextType','string');
assert(isequal(accTbl.month(:), string(monthFirstDay,'yyyy-MM').'), ...
       'monthly_rain_accumulation_ver_3.csv months do not match Feb 2023 - Feb 2024');
rainCols = ["LD_M1","LD_S1","VDIS_M1"];          % same order as rainNames
rainTotal_mm  = nan(3,nMonth);                   % [instrument x month]
rainDaysRec   = nan(3,nMonth);                   % days with data recorded
for g = 1:3
    rainTotal_mm(g,:) = accTbl.(rainCols(g) + "_mm").';
    rainDaysRec(g,:)  = accTbl.(rainCols(g) + "_days_recorded").';
end
rainDaysInPeriod = accTbl.days_in_period.';
% Coverage flags, as in Fig. A3:
%   noData   : no day recorded -> no bar, 'x' (VDIS Aug, Oct, Nov 2023; Die
%              Wang's screen removed VDIS Aug and Oct, the VDIS was off in Nov)
%   halfData : < 50% of the period's days recorded -> '*', total is a lower
%              bound (LD/M1 Oct 2023, 15/31 days)
rainNoData   = rainDaysRec == 0;
rainHalfData = ~rainNoData & rainDaysRec < 0.5*rainDaysInPeriod;
rainTotal_mm(rainNoData) = NaN;                  % the CSV writes 0.00 for these

% Slide version: mean and std of the daily rain on rain days.
rain = readRain(fullfile(echoDir,'MEAN-EPCAPE Russel Rainfall.xlsx'), ...
                fullfile(echoDir,'STD-Russel EPCAPE Rainfall.xlsx'), monthFirstDay);
rainColors = [blue; orange; yellow];

%% 8. Figure
fprintf('[%s] 8. Figure\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)

fig = figure('Visible','off','Units','pixels','Position',[50 50 640 703]);   % slide image aspect (879 x 966)
resetToFactoryDefaults(fig);   % undo Andrew's startup.m dark theme for this figure
fig.Color = 'w';
% Font sizes chosen to match the slide image at this figure size
% (factory 10 pt looked ~25% smaller than the slide).
set(fig,'DefaultAxesFontSize',12,'DefaultLegendFontSize',11, ...
        'DefaultColorbarFontSize',12,'DefaultTextFontSize',12);

% Axes positions measured from the slide image (879 x 966 px), normalized.
L = 95/879;  W = (658-95)/879;          % common left edge and width
pos = @(yTop,yBot) [L, 1-yBot/966, W, (yBot-yTop)/966];
legX  = 675/879;                        % legend/colorbar left edge, panels A, B, F
legX2 = 718/879;                        % panels D, E (room for the right y label)
xl  = [0.5 nMonth+0.5];

% ---- A: heights -------------------------------------------------------------
axA = axes(fig,'Position',pos(3,255));
bA = boxchart(axA, hMonth, hVals, 'GroupByColor', hGroup, ...
              'MarkerStyle','none','Notch','on','WhiskerLineStyle','none');
for g = 1:numel(bA)
    bA(g).BoxFaceColor = heightColors(g,:);
end
ylim(axA,[0 1750]); xlim(axA,xl)
set(axA,'XTick',1:nMonth,'XTickLabel',[],'Box','on')
ylabel(axA,'Altitude [m]')
lgA = legend(axA, bA, cellstr(heightNames), 'Units','normalized');
lgA.Position(1:2) = [legX, axA.Position(2)+axA.Position(4)-lgA.Position(4)-25/966];

% ---- B: echo fraction + CBH -------------------------------------------------
axB = axes(fig,'Position',pos(272,532));
hold(axB,'on')
imagesc(axB, 1:nMonth, zEcho_m, echoProfile, 'AlphaData', ~isnan(echoProfile));
set(axB,'YDir','normal')
colormap(axB, sky(256)); clim(axB,[0 0.5])
cb = colorbar(axB,'Units','normalized');
cb.Position = [legX, axB.Position(2), 32/879, axB.Position(4)];
if lynnEdits, ylabel(cb,'Hydrometeor Fraction'), else, ylabel(cb,'Echo Fraction'), end
% Real values only (no zero padding): this is why the Feb boxes do not reach
% 0 m as they do on Kavin's slide (header, difference 1).
drawBoxes(axB, 1:nMonth, groupByMonth(cbh_m, mCeil, nMonth), orange, boxWidth, whiskerDef);
ylim(axB,[0 3000]); xlim(axB,xl)
set(axB,'XTick',1:nMonth,'XTickLabel',[],'Box','on','Layer','top')
ylabel(axB,'CBH [m]')

% ---- C: LWP -----------------------------------------------------------------
axC = axes(fig,'Position',pos(550,624));
% Same as panel B: real values only (header, difference 1).
drawBoxes(axC, 1:nMonth, groupByMonth(lwp_gm2, mMwr, nMonth), blue, boxWidth, whiskerDef);
ylim(axC,[0 500]); xlim(axC,xl)
set(axC,'XTick',1:nMonth,'XTickLabel',[],'Box','on')
ylabel(axC,'LWP')

% ---- D: UCVF + RH + T -------------------------------------------------------
axD = axes(fig,'Position',pos(643,713));
yyaxis(axD,'left')
hU = drawBoxes(axD, 1:nMonth, groupByMonth(ucvf, m15, nMonth), blue, boxWidth, whiskerDef);
ylim(axD,[0 1]); ylabel(axD,'UCVF')
yyaxis(axD,'right')
hRH = plot(axD, 1:nMonth, rhMean_pct, '-',  'Color',orange,'LineWidth',0.75);
hold(axD,'on')
hT  = plot(axD, 1:nMonth, TPct20_pct, '--', 'Color',orange,'LineWidth',0.75);
ylim(axD,[50 100]); set(axD,'YTick',[60 80 100]); ylabel(axD,'%')
axD.YAxis(1).Color = 'k'; axD.YAxis(2).Color = 'k';
xlim(axD,xl); set(axD,'XTick',1:nMonth,'XTickLabel',[],'Box','on')
lgD = legend(axD,[hU hRH hT],{'UCVF','RH','% of 20^{\circ}C'},'Units','normalized');
lgD.Position(1:2) = [legX2, axD.Position(2)+axD.Position(4)-lgD.Position(4)];

% ---- E: LW line + SW boxes --------------------------------------------------
axE = axes(fig,'Position',pos(733,807));
yyaxis(axE,'right')
hSW = drawBoxes(axE, 1:nMonth, swGroups, blue, boxWidth, whiskerDef);
ylim(axE,[-0.15 1.1]); set(axE,'YTick',[0 0.5 1]); ylabel(axE,'SW')
yyaxis(axE,'left')
hLW = plot(axE, 1:nMonth, lwMonthly_Wm2, '-', 'Color',lwColor,'LineWidth',0.75);
% Slide: [300 410]. Lowered to 290 so the Feb 2023 value (297 W m^-2) shows.
ylim(axE,[290 410]); set(axE,'YTick',[300 350 400]); ylabel(axE,'LW')
axE.YAxis(1).Color = 'k'; axE.YAxis(2).Color = 'k';
xlim(axE,xl); set(axE,'XTick',1:nMonth,'XTickLabel',[],'Box','on')
lgE = legend(axE,[hLW hSW],{'LW','SW'},'Units','normalized');
lgE.Position(1:2) = [legX2, axE.Position(2)+axE.Position(4)-lgE.Position(4)];

% ---- F: rainfall ------------------------------------------------------------
axF = axes(fig,'Position',pos(824,893));
hold(axF,'on')
dx = [-0.25 0 0.25];                    % dodge offsets for the three instruments
hR = gobjects(1,3);
if rainPanel == "monthlyTotal"
    % Grouped bars of the monthly accumulation (Fig. A3 numbers).
    yTop = 145;                                          % mm; max is 135.5 (LD/S1 Mar 2023)
    for g = 1:3
        hR(g) = bar(axF, (1:nMonth)+dx(g), rainTotal_mm(g,:), 0.22, ...
                    'FaceColor',rainColors(g,:), 'EdgeColor','k', 'LineWidth',0.5);
        % coverage markers just above the bar (or the axis, if no bar)
        yMark = max(rainTotal_mm(g,:), 0, 'includemissing');
        yMark(isnan(yMark)) = 0;
        k = find(rainNoData(g,:));
        text(axF, k+dx(g), yMark(k) + 0.02*yTop, 'x', 'FontSize',9, ...
             'HorizontalAlignment','center','VerticalAlignment','bottom','Color',[0.4 0.4 0.4])
        k = find(rainHalfData(g,:));
        text(axF, k+dx(g), yMark(k) + 0.02*yTop, '*', 'FontSize',11, ...
             'HorizontalAlignment','center','VerticalAlignment','bottom','Color',[0.4 0.4 0.4])
    end
    ylim(axF,[0 yTop])
    ylabel(axF,'Rain [mm]')
    text(axF, 0.5, 0.95, 'x  no data     *  > 50% of days missing', 'Units','normalized', ...
         'FontSize',9,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center','VerticalAlignment','top')
else
    for g = 1:3
        hR(g) = drawMeanStdBoxes(axF, (1:nMonth)+dx(g), rain.mean(g,:), rain.std(g,:), ...
                                 rainColors(g,:), 0.14);   % narrow boxes, as on the slide
    end
    ylim(axF,[0 20])
    ylabel(axF,'Rainfall')
end
xlim(axF,xl)
set(axF,'XTick',1:nMonth,'XTickLabel',monthLabels,'XTickLabelRotation',30,'Box','on')
xlabel(axF,'Month')
lgF = legend(axF,hR,cellstr(rainNames),'Units','normalized');
lgF.Position(1:2) = [legX, axF.Position(2)+axF.Position(4)-lgF.Position(4)];

% ---- Panel letters (Lynn's comment) -----------------------------------------
if lynnEdits
    allAx = [axA axB axC axD axE axF];
    for k = 1:numel(allAx)
        text(allAx(k), 0.01, 0.97, char('A'+k-1), 'Units','normalized', ...
             'FontSize',12,'FontWeight','bold','VerticalAlignment','top', ...
             'BackgroundColor','w','Margin',1);
    end
end

if exportFigs
    fprintf('[%s] exporting PNG/PDF\n', string(datetime('now','Format','HH:mm:ss')));
    exportgraphics(fig, fullfile(outDir,[figTag '.png']), 'Resolution',300);
    exportgraphics(fig, fullfile(outDir,[figTag '.pdf']), 'ContentType','vector');
end

fprintf('Figure has %d graphics objects\n', numel(findall(fig)));
if showFigureAtEnd
    fig.Visible = 'on';    % show only after the export (see Settings)
end

%% 9. Table of every plotted variable and the data behind it
fprintf('[%s] 9. Table of every plotted variable\n', string(datetime('now','Format','HH:mm:ss')));  % progress (find hangs)
% One row per plotted variable per month. Columns:
%   Panel, Variable (as in the legend/axis), SourceFile, FileVariable (name
%   inside the file), Units, Shown (how it is drawn), Month, N (samples used
%   that month), Value (for lines/bars/shading), Q1/Median/Q3 and
%   WhiskerLow/WhiskerHigh (for boxes; same numbers drawBoxes uses).
% Boxes that are not drawn (N below minSamplesPerBox in panel A, or N < 2)
% show N with NaN statistics.

srcLaRC = "EPCAPE_CTH_CBH_PBLH_SChellappan_LaRC.csv";
plotted = table();

% ---- A: heights (boxchart: box and median only, no whiskers) -------------
fileA = [srcLaRC srcLaRC srcLaRC srcLaRC "BL_Depth.csv" "EPCAPE_sonde_LCL.csv" srcLaRC];
varA  = ["GOES18_CTH","KAZR_CTH","PBLH_THERMO","PBLH_MPL","Depth","LCL","CEIL_CBH"];
for g = 1:numel(heightNames)
    inG = hGroup == heightNames(g);
    plotted = [plotted; boxRows("A", heightNames(g), fileA(g), varA(g), "m", ...
        "notched box, no whiskers", groupByMonth(hVals(inG), hMonth(inG), nMonth), ...
        monthLabels, whiskerDef)]; %#ok<AGROW>
end
% boxchart draws no whiskers in panel A, so blank them in the table.
plotted.WhiskerLow(plotted.Panel == "A")  = NaN;
plotted.WhiskerHigh(plotted.Panel == "A") = NaN;

% ---- B: echo fraction (shading) and CBH (boxes) --------------------------
% Shading is the monthly mean profile (596 heights); the table gives its mean
% over the plotted height range, 160-3000 m, as one number per month.
inRange = zEcho_m <= 3000;
plotted = [plotted; valueRows("B", "Echo fraction", "arscl_echomask_hourly.nc", ...
    "echo_fraction", "1", "colour: monthly mean profile (table: mean over 160-3000 m)", ...
    mean(echoProfile(inRange,:),1,'omitmissing'), accumarray(mEcho(~isnan(mEcho)),1,[nMonth 1]).', monthLabels)];
plotted = [plotted; boxRows("B", "CBH", "epcceil_hourly_cbh.nc", "bl_cbh", "m", ...
    "box + whiskers", groupByMonth(cbh_m, mCeil, nMonth), monthLabels, whiskerDef)];

% ---- C: LWP ---------------------------------------------------------------
plotted = [plotted; boxRows("C", "LWP", "epcmwrret_turnM1_hourly_averages.nc", ...
    "phys_lwp", "g m^-2", "box + whiskers", groupByMonth(lwp_gm2, mMwr, nMonth), ...
    monthLabels, whiskerDef)];

% ---- D: UCVF (boxes), RH and T (lines) ------------------------------------
plotted = [plotted; boxRows("D", "UCVF", "EPCAPE_15mindata.nc", "UCVF", "1", ...
    "box + whiskers (left axis)", groupByMonth(ucvf, m15, nMonth), monthLabels, whiskerDef)];
nRH = cellfun(@numel, groupByMonth(rh_pct, m15, nMonth));
nT  = cellfun(@numel, groupByMonth(T_degC, m15, nMonth));
plotted = [plotted; valueRows("D", "RH", "EPCAPE_15mindata.nc", "MET_rh", "%", ...
    "line: monthly mean (right axis)", rhMean_pct, nRH, monthLabels)];
plotted = [plotted; valueRows("D", "% of 20 C", "EPCAPE_15mindata.nc", "MET_temp", "%", ...
    "dashed line: 100 x monthly mean T [C] / 20 C (right axis)", TPct20_pct, nT, monthLabels)];
plotted = [plotted; valueRows("D", "T (for reference, not plotted)", "EPCAPE_15mindata.nc", ...
    "MET_temp", "deg C", "not drawn: monthly mean behind '% of 20 C'", TMean_degC, nT, monthLabels)];

% ---- E: SW (boxes) and LW (line) ------------------------------------------
plotted = [plotted; boxRows("E", "SW", "<MON>_SW_Boxplot.csv (Dan); FEB split by row", ...
    "hourly SW transmittance", "1", "box + whiskers (right axis)", swGroups, monthLabels, whiskerDef)];
nLw = cellfun(@numel, groupByMonth(lwHourly.lwRaw, mLw, nMonth));
plotted = [plotted; valueRows("E", "LW", "radflux_M1_20230215_20240214.nc", ...
    "downwelling_longwave (qc == 0)", "W m^-2", ...
    "line: mean over 24 h of hourly medians (N = hourly means)", lwMonthly_Wm2, nLw, monthLabels)];
plotted = [plotted; valueRows("E", "LW, Dan (for reference, not plotted)", ...
    "<MON>_LW_Stepplot.csv (Dan)", "Median column", "W m^-2", ...
    "not drawn: Dan's value of the same statistic (Feb pools 2023+2024)", ...
    danLW_Wm2, nan(1,nMonth), monthLabels)];

% ---- F: rain ----------------------------------------------------------------
for g = 1:3
    if rainPanel == "monthlyTotal"
        plotted = [plotted; valueRows("F", rainNames(g), "monthly_rain_accumulation_ver_3.csv", ...
            rainCols(g) + "_mm", "mm", "bar: monthly total (N = days recorded)", ...
            rainTotal_mm(g,:), rainDaysRec(g,:), monthLabels)]; %#ok<AGROW>
    else
        plotted = [plotted; valueRows("F", rainNames(g), "MEAN-/STD-...Russel Rainfall.xlsx", ...
            "mean", "mm per rain day", "box: mean (line) +/- 1 std (Value = mean)", ...
            rain.mean(g,:), nan(1,nMonth), monthLabels)]; %#ok<AGROW>
    end
end

% Round for display (the CSV keeps full precision).
show = plotted;
numCols = ["Value","Q1","Median","Q3","WhiskerLow","WhiskerHigh"];
for c = numCols
    show.(c) = round(show.(c), 3, 'significant');
end
disp(show)
writetable(plotted, fullfile(outDir, [figTag '_plotted_data.csv']));
fprintf('Plotted-data table (%d rows) written to %s\n', height(plotted), ...
        fullfile(outDir, [figTag '_plotted_data.csv']));

%% Local functions

function t = naiveToZone(t, tz, isUTC)
% naiveToZone  Give a zone-less datetime vector the time zone tz.
%   isUTC = true : the stamps are UTC and are converted to tz.
%   isUTC = false: the stamps are already local time in tz.
if isUTC
    t.TimeZone = 'UTC';
    t.TimeZone = tz;
else
    t.TimeZone = tz;
end
end

function t = readNcTime(file, tz, isUTC)
% readNcTime  Read a netCDF 'time' variable with CF units '<unit> since <date>'
% and return a datetime in time zone tz. ARM files say UTC in the units or
% long_name; for files that do not (EPCAPE_15mindata.nc), isUTC decides.
if nargin < 3, isUTC = true; end
raw   = double(ncread(file,'time'));
units = ncreadatt(file,'time','units');
tok   = regexp(units,'^\s*(\w+)\s+since\s+(\d{4}-\d{2}-\d{2})[ T]?(\d{2}:\d{2}:\d{2})?','tokens','once');
if isempty(tok{3}), tok{3} = '00:00:00'; end           % 'since 2023-02-15'
epoch = datetime([tok{2} ' ' tok{3}], 'InputFormat','yyyy-MM-dd HH:mm:ss');
switch lower(tok{1})
    case 'nanoseconds', dt = seconds(raw/1e9);
    case 'seconds',     dt = seconds(raw);
    case 'minutes',     dt = minutes(raw);
    case 'hours',       dt = hours(raw);
    case 'days',        dt = days(raw);
    otherwise, error('readNcTime: unknown time unit "%s" in %s', tok{1}, file);
end
t = naiveToZone(epoch + dt, tz, isUTC);
end

function m = monthIndex(t, tStart, tEnd)
% monthIndex  Campaign-month number of each time: 1 = Feb 2023 ... 13 = Feb
% 2024 (month of t in the zone of tStart). NaN outside [tStart, tEnd).
t.TimeZone = tStart.TimeZone;
m = (year(t) - 2023)*12 + month(t) - 1;
m(t < tStart | t >= tEnd) = NaN;
m = m(:);
end

function c = groupByMonth(x, m, nMonth)
% groupByMonth  Split x into one column vector per campaign month (NaNs kept
% out). Avoids the fixed-size, zero-padded arrays of the precursor code.
x = x(:);
c = cell(1,nMonth);
for k = 1:nMonth
    v = x(m == k);
    c{k} = v(~isnan(v));
end
end

function [vals, M] = readSW(file)
% readSW  Read one of Dan's <MON>_SW_Boxplot.csv files.
%   M    : [days x 24] SW transmittance, columns = PST hours 0-23
%   vals : M as one column vector
% The header row (hours 0..23, 'Means') is skipped; column 25 ('Means') is
% dropped; the -0.05 night fill becomes NaN.
M = readmatrix(file, 'NumHeaderLines', 1);
M = M(:,1:24);
M(abs(M + 0.05) < 1e-9) = NaN;
vals = M(:);
vals = vals(~isnan(vals));
end

function rain = readRain(meanFile, stdFile, monthFirstDay)
% readRain  Read the monthly rain mean and std spreadsheets.
% Each has rows of {'yyyy/MM', value, 'MEAN-<inst>' or 'SD-<inst>'}; extra
% columns in the MEAN file (unlabelled numbers to the right) are ignored.
% Output rain.mean / rain.std are [3 instruments x 13 months] (NaN = no data),
% instruments in the order LD/M1, LD/S1, VDIS/M1.
inst = ["LD/M1","LD/S1","VDIS/M1"];
monthKey = string(monthFirstDay,'yyyy/MM');
rain.mean = nan(3,numel(monthFirstDay));
rain.std  = nan(3,numel(monthFirstDay));
files = {meanFile, stdFile};  fields = {'mean','std'};
for f = 1:2
    C = readcell(files{f});
    for r = 2:size(C,1)                                   % row 1 = header
        key = C{r,1};
        if isdatetime(key), key = string(key,'yyyy/MM'); end  % if Excel stored a date
        if ~ischar(key) && ~isstring(key), continue, end
        k = find(monthKey == string(key));
        g = find(inst == regexprep(string(C{r,3}),'^(MEAN|SD)-',''));
        v = C{r,2};
        if isempty(k) || isempty(g), continue, end
        if isnumeric(v), rain.(fields{f})(g,k) = v; end   % 'missing' -> NaN
    end
end
end

function h = drawBoxes(ax, x, groups, boxColor, w, whiskerDef)
% drawBoxes  Box-and-whisker plot drawn with lines, so that boxes can share
% an axes (or a yyaxis side) with images and lines at the same x positions.
%   box    : 25th-75th percentile (prctile, the same method as boxplot)
%   line   : median
%   whisker: "tukey"  -> furthest point within 1.5 IQR of the box
%            "p10p90" -> 10th and 90th percentiles
% Boxes in boxColor; whiskers and caps black (as on the slide). Returns the
% box patch handle for the legend.
% All boxes of one call are drawn as just THREE graphics objects (one
% NaN-separated whisker line, one multi-face box patch, one NaN-separated
% median line) rather than ~6 per box: Andrew has seen MATLAB hang on figures
% with many graphics objects.
hold(ax,'on')
xw = [];  yw = [];          % whiskers + caps (black), NaN-separated
vb = [];  fb = [];          % box patch vertices and faces
xm = [];  ym = [];          % median segments, NaN-separated
for k = 1:numel(x)
    v = groups{k};
    v = v(~isnan(v));
    if numel(v) < 2, continue, end
    [p, lo, hi] = boxStats(v, whiskerDef);
    c = x(k);  cap = w/4*[-1 1];
    xw = [xw, c c NaN, c c NaN, c+cap NaN, c+cap NaN];            %#ok<AGROW>
    yw = [yw, hi p(3) NaN, p(1) lo NaN, hi hi NaN, lo lo NaN];   %#ok<AGROW>
    fb = [fb; size(vb,1) + (1:4)];                                 %#ok<AGROW>
    vb = [vb; c+w/2*[-1;1;1;-1], [p(1);p(1);p(3);p(3)]];           %#ok<AGROW>
    xm = [xm, c+w/2*[-1 1], NaN];                                  %#ok<AGROW>
    ym = [ym, p(2) p(2), NaN];                                     %#ok<AGROW>
end
h = gobjects(1);
if isempty(fb), return, end
plot(ax, xw, yw, 'k-', 'LineWidth',1);
% box as an unfilled patch so the legend icon is a box (as on the slide)
h = patch(ax, 'Vertices',vb, 'Faces',fb, 'FaceColor',boxColor, 'FaceAlpha',0, ...
          'EdgeColor',boxColor, 'LineWidth',1);
plot(ax, xm, ym, '-', 'Color',boxColor, 'LineWidth',1);
end

function [p, lo, hi] = boxStats(v, whiskerDef)
% boxStats  Quartiles p = [Q1 median Q3] (prctile, as boxplot) and whisker
% ends lo/hi: "tukey" = furthest values within 1.5 IQR of the box,
% "p10p90" = 10th/90th percentiles. Used by drawBoxes and by the table in
% section 9, so the table holds exactly the drawn numbers.
v = v(~isnan(v));
p = prctile(v,[25 50 75]);
if whiskerDef == "tukey"
    iqrV = p(3) - p(1);
    lo = min(v(v >= p(1) - 1.5*iqrV));
    hi = max(v(v <= p(3) + 1.5*iqrV));
else
    lo = prctile(v,10);
    hi = prctile(v,90);
end
end

function T = boxRows(panel, var, file, fileVar, units, shown, groups, months, whiskerDef)
% boxRows  Table rows (one per month) for a variable drawn as boxes.
n = numel(groups);
stats = nan(n,5);  N = zeros(n,1);
for k = 1:n
    v = groups{k};  v = v(~isnan(v));  N(k) = numel(v);
    if N(k) >= 2
        [p, lo, hi] = boxStats(v, whiskerDef);
        stats(k,:) = [p lo hi];
    end
end
T = rowsTable(panel, var, file, fileVar, units, shown, months, N, nan(n,1), stats);
end

function T = valueRows(panel, var, file, fileVar, units, shown, values, N, months)
% valueRows  Table rows (one per month) for a variable drawn as a line, bar
% or shading (one number per month).
n = numel(values);
T = rowsTable(panel, var, file, fileVar, units, shown, months, N(:), values(:), nan(n,5));
end

function T = rowsTable(panel, var, file, fileVar, units, shown, months, N, value, stats)
% rowsTable  Common layout of the section-9 table.
n = numel(months);
T = table(repmat(string(panel),n,1), repmat(string(var),n,1), repmat(string(file),n,1), ...
    repmat(string(fileVar),n,1), repmat(string(units),n,1), repmat(string(shown),n,1), ...
    string(months(:)), N(:), value(:), stats(:,1), stats(:,2), stats(:,3), stats(:,4), stats(:,5), ...
    'VariableNames',{'Panel','Variable','SourceFile','FileVariable','Units','Shown', ...
    'Month','N','Value','Q1','Median','Q3','WhiskerLow','WhiskerHigh'});
end

function h = drawMeanStdBoxes(ax, x, mu, sd, c, w)
% drawMeanStdBoxes  Rain panel: a box from mean-std to mean+std with a line at
% the mean (what Kavin's boxchart of [mean-2s mean-s mean mean mean+s
% mean+2s] draws, without the 6-point boxchart detour). Months with no data
% (NaN) are skipped. Light fill + coloured edge, like boxchart.
h = gobjects(1);
for k = 1:numel(x)
    if isnan(mu(k)), continue, end
    lo = mu(k) - sd(k);  hi = mu(k) + sd(k);
    if isnan(sd(k)), lo = mu(k); hi = mu(k); end
    hp = patch(ax, x(k)+w/2*[-1 1 1 -1], [lo lo hi hi], c, ...
               'FaceAlpha',0.2,'EdgeColor',c,'LineWidth',1);
    plot(ax, x(k)+w/2*[-1 1], [mu(k) mu(k)], '-', 'Color',c,'LineWidth',1);
    if ~isgraphics(h), h = hp; end
end
end

function resetToFactoryDefaults(fig)
% resetToFactoryDefaults  Give one figure MATLAB's factory graphics defaults.
% Andrew's startup.m sets ROOT defaults (black background, white text, 18 pt
% bold fonts, 4 pt lines). Setting each changed default back to its factory
% value at the FIGURE level fixes this figure only. defaultFigure* properties
% are skipped (they only apply to figures created later).
rootDefaults = fieldnames(get(groot,'default'));
for k = 1:numel(rootDefaults)
    name = rootDefaults{k};
    if startsWith(name,'defaultFigure','IgnoreCase',true), continue, end
    set(fig, name, get(groot, regexprep(name,'^default','factory','ignorecase')));
end
fig.Color = get(groot,'factoryFigureColor');
end
