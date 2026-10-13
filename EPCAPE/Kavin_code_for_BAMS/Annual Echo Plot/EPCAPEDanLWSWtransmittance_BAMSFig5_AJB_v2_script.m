%% EPCAPEDanLWSWtransmittance_BAMSFig5_AJB_v2  -- makes BAMS Figure 5
% Seasonal diurnal cycle at the EPCAPE M1 site (Scripps Pier), 15 Feb 2023 -
% 17 Feb 2024, in four seasonal columns (Spring = Apr-Jun, Summer = Jul-Sep,
% Fall = Oct-Dec, Winter = Jan-Mar):
%   row A: ARSCL hydrometeor fraction profiles (colour) with ceilometer
%          boundary-layer cloud base height (orange boxes)
%   row B: MWR liquid water path (boxes)
%   row C: SW cloud transmittance (blue boxes, right axis) and mean
%          downwelling LW flux (dark-red line, left axis)
% Output: EPCAPEDanLWSWtransmittance_BAMSFig5_AJB_v2.png (300 dpi) and .pdf.
%
% LINEAGE
%   EPCAPEDanLWSWtransmittance.mlx         Kavin's original (R2025b), untouched.
%   EPCAPEDanLWSWtransmittance_AJB.mlx     v1: runs on Andrew's machine and
%                                          reproduces Kavin's figures exactly.
%   ..._BAMSFig5_AJB_v2.mlx (this file)    v2: bug fixes + Lynn Russell's comments;
%                                          Fig. 4 code deleted.
%   Plain-language change log for Kavin/Lynn:
%   EPCAPE_BAMS_Fig5_changes_and_open_questions.txt (same folder).
%
% EDIT TAGS (search for these to find every change):
%   % AJB edits:       v1 change (file paths, figure handling, exports, ';').
%   % AJB note:        v1 comment only; flags an issue in Kavin's code.
%   % AJB edits (v2):  v2 change. Kavin's replaced lines are kept, commented
%                      out, next to the new code wherever a block was rewritten.
%
% v2 CHANGES TO THE DATA (each changes what Fig. 5 shows)
%   1. SW hour offset fixed. readtable reads the SW CSV header row of hours as
%      data, and Kavin's Var2:Var24 put hour k in column k, while every other
%      panel puts hour h in column h+1. SW was drawn 1 h early.
%   2. One time base for all panels: PST (UTC-8, no daylight saving). The SW
%      CSV hours are PST: valid daytime hours are centred on ~11.5 in every
%      month with no jump at the DST changes, and June data start in hour 4
%      (sunrise 4:42 PST = 5:42 PDT). The netCDF panels were in DST-aware
%      local time, so during PDT (12 Mar - 5 Nov 2023) SW was a further 1 h
%      early. The netCDF times are now converted to PST (setting ncTimeZone).
%      ASSUMPTION: the LW CSVs use the same PST hours as the SW CSVs.
%   3. Day 31 of Jan/Mar/May/Jul/Aug/Oct/Dec is now included in the SW boxes
%      (Kavin's rows 2:31 kept only 30 days).
%   4. The SW fill value -0.05 (night hours only) is set to NaN. No visible
%      change in Fig. 5; it does change the printed annual/seasonal SW means.
%   5. d2 (ceilometer time) now uses the ceilometer's own time vector. No
%      change in values: both files have identical time stamps.
%
% v2 CHANGES FOR LYNN RUSSELL'S COMMENTS
%   "Typo in Transmittance"          -> 'SW Transmittance' (was 'Trasmittance').
%   "Add panel labels -- A, B, etc." -> A, B, C down the first column (one per
%                                       row); columns are labelled by season.
%   "Add units to LW Flux and Height if they fit"
%                                    -> 'LW Flux [W m^{-2}]'; the left height
%                                       axis is now 'Cloud Base Height [m]'
%                                       (orange, Andrew's edit).
%   "why is right axis label needed?" / "remove label for Cloud Base Height on
%    right"                          -> removed. The right axis has no tick
%                                       labels and the same 0-3000 m range as
%                                       the left axis. (An orange 'CBH' key was
%                                       tried in panel D; Andrew switched it off
%                                       in favour of the orange left label.)
%   "Change echo to Hydrometeor"     -> colorbar 'Hydrometeor Fraction'.
%
% v2 CHANGES TO MATCH THE SLIDE VERSION OF FIG. 5 (Kavin's later revision)
%   Black axis lines/ticks/labels; no '23' tick label; 'SW' (blue) and 'LW'
%   (red) keys in the Winter bottom panel. Also: LWP unit written g m^{-2}, and
%   the x label is 'Local Hour [PST]' (set by timeAxisLabel below).
%
% v2 CODE DELETED (still in v1 and in Kavin's original)
%   - The Fig. 4 precursor (calendar-month version; the slide's Fig. 4 is made by
%     other code) and the inputs only it used (BL depth, sonde LCL).
%   - Kavin's Spring+Summer test figure (not in BAMS; errors once day 31 is kept).
%   - The LWP mean/std lines from the Fig. 4 block are REPLACED by a corrected
%     table, lwpStats (zero-padding and one-month column offset fixed).
%   The test figures before Fig. 5 (LW lines, SW boxcharts) still run in the
%   same figure window; each is wiped by the next 'clf'.
%
% EDITS BY ANDREW (2026-10-13, in this file): x label 'Local Hour [PST]'; font
%   sizes ax_fnt/cb_fnt; left y label 'Cloud Base Height [m]' in orange; the
%   orange 'CBH' key commented out.
%
% Andrew J. Buggee, Scripps Institution of Oceanography, 2026-10-13

% AJB edits: settings block (v1). Kavin's original read every file from
% MATLAB's current folder; these variables make the paths explicit.
dataDir_andrew = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS/Annual Echo Plot';
if isfolder(dataDir_andrew)
    dataDir = dataDir_andrew;      % Andrew's machine
else
    dataDir = pwd;                 % anyone else: run from the data folder (Kavin's original behaviour)
end
outDir = dataDir;                  % where the exported figures are written

% AJB edits: figure visibility. 'off' in the .mlx (a visible figure can hang
% exportgraphics in the R2024b Live Editor); 'on' in the plain _script.m.
figureVisible = 'on';

% AJB edits: set false to skip writing the PNG/PDF files.
exportFigs = true;

% AJB edits (v2): time zone for the hourly netCDF files (ARSCL, ceilometer,
% MWR), whose time stamps are hour-start UTC. '-08:00' = PST, a fixed UTC-8
% offset with no daylight saving, to match the hours in the SW (and, assumed,
% LW) CSVs. Set to 'America/Los_Angeles' to get Kavin's DST-aware local time
% back (the SW panels would then be 1 h off from the others during PDT).
ncTimeZone = '-08:00';
if strcmp(ncTimeZone,'-08:00')
    timeAxisLabel = 'Local Hour [PST]';
else
    timeAxisLabel = 'Local Time [Hour]';
end

%%
APR_LW_Stepplot = readtable(fullfile(dataDir,"APR_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
APR_SW_Boxplot = readtable(fullfile(dataDir,"APR_SW_Boxplot.csv"))  % AJB edits: read from dataDir (was a bare file name)
AUG_LW_Stepplot = readtable(fullfile(dataDir,"AUG_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
AUG_SW_Boxplot = readtable(fullfile(dataDir,"AUG_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
DEC_LW_Stepplot = readtable(fullfile(dataDir,"DEC_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
DEC_SW_Boxplot = readtable(fullfile(dataDir,"DEC_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
DJF_LW_Stepplot = readtable(fullfile(dataDir,"DJF_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
FEB_LW_Stepplot = readtable(fullfile(dataDir,"FEB_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
FEB_SW_Boxplot = readtable(fullfile(dataDir,"FEB_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
JAN_LW_Stepplot = readtable(fullfile(dataDir,"JAN_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
JAN_SW_Boxplot = readtable(fullfile(dataDir,"JAN_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
JJA_LW_Stepplot = readtable(fullfile(dataDir,"JJA_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
JUL_LW_Stepplot = readtable(fullfile(dataDir,"JUL_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
JUL_SW_Boxplot = readtable(fullfile(dataDir,"JUL_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
JUN_LW_Stepplot = readtable(fullfile(dataDir,"JUN_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
JUN_SW_Boxplot = readtable(fullfile(dataDir,"JUN_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
MAM_LW_Stepplot = readtable(fullfile(dataDir,"MAM_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
MAR_LW_Stepplot =readtable(fullfile(dataDir,"MAR_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
MAR_SW_Boxplot = readtable(fullfile(dataDir,"MAR_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
MAY_LW_Stepplot = readtable(fullfile(dataDir,"MAY_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
MAY_SW_Boxplot = readtable(fullfile(dataDir,"MAY_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
NOV_LW_Stepplot = readtable(fullfile(dataDir,"NOV_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
NOV_SW_Boxplot = readtable(fullfile(dataDir,"NOV_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
OCT_LW_Stepplot = readtable(fullfile(dataDir,"OCT_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
OCT_SW_Boxplot = readtable(fullfile(dataDir,"OCT_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
SEP_LW_Stepplot = readtable(fullfile(dataDir,"SEP_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
SEP_SW_Boxplot = readtable(fullfile(dataDir,"SEP_SW_Boxplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
SON_LW_Stepplot = readtable(fullfile(dataDir,"SON_LW_Stepplot.csv"));  % AJB edits: read from dataDir (was a bare file name)
%%
decLW = DJF_LW_Stepplot{:,"DEC"};
janLW = DJF_LW_Stepplot{:,"JAN"};
febLW = table2array(DJF_LW_Stepplot(:,"FEB"));
djftime = table2array(DJF_LW_Stepplot(:,"Hour"));

marLW = table2array(MAM_LW_Stepplot(:,"MAR"));
aprLW = table2array(MAM_LW_Stepplot(:,"APR"));
mayLW = table2array(MAM_LW_Stepplot(:,"MAY"));
mamtime = table2array(MAM_LW_Stepplot(:,"Hour"));

junLW = table2array(JJA_LW_Stepplot(:,"JUN"));
julLW = table2array(JJA_LW_Stepplot(:,"JUL"));
augLW = table2array(JJA_LW_Stepplot(:,"AUG"));
jjatime = table2array(JJA_LW_Stepplot(:,"Hour"));

sepLW = table2array(SON_LW_Stepplot(:,"SEP"));
octLW = table2array(SON_LW_Stepplot(:,"OCT"));
novLW = table2array(SON_LW_Stepplot(:,"NOV"));
sontime = table2array(SON_LW_Stepplot(:,"Hour"));


% AJB edits: create the working figure explicitly (Kavin relied on an
% implicit figure). 'Visible' comes from the figureVisible setting at the top
% ('off' in the .mlx to avoid the R2024b Live Editor exportgraphics hang).
% resetToFactoryDefaults (local function at the end of this file) undoes the
% root graphics defaults from Andrew's startup.m for this figure only.
fig_work = figure('Visible',figureVisible);
resetToFactoryDefaults(fig_work);
t = tiledlayout(2,2,'TileSpacing','compact');
nexttile
hold on
plot(djftime,decLW)
plot(djftime,janLW)
plot(djftime,febLW)
hold off
legend("DEC","JAN","FEB")

nexttile
hold on
plot(mamtime, marLW)
plot(mamtime,aprLW)
plot(mamtime,mayLW)
hold off
legend("MAR","APR","MAY")

nexttile
hold on
plot(jjatime, junLW)
plot(jjatime,julLW)
plot(jjatime,augLW)
hold off
legend("JUN","JUL","AUG")

nexttile
hold on
plot(sontime,sepLW)
plot(sontime,octLW)
plot(sontime,novLW)
hold off
legend("SEP","OCT","NOV")
%%
% AJB edits (v2): SW transmittance block rewritten. Kavin's original lines are kept,
% commented out, just below this new code.
%
% (1) HOUR OFFSET FIX. readtable reads the CSV header row (0,1,...,23,Means)
%     as data: row 1 holds the hour labels and Var1 = hour 0 ... Var24 = hour 23.
%     Kavin took Var2:Var24 (hours 1-23) into columns 1-23 and padded column 24
%     with NaN, so column k held hour k. Every other diurnal panel in Fig. 5 puts
%     hour h in column h+1 (boxplot position h+1, tick label h), so SW was drawn
%     1 h early. Now column k = hour k-1, like the other panels. The assert
%     below stops the script if a CSV ever has a different hour row.
% (2) ALL DAYS. Kavin kept rows 2:31 (30 days), which dropped day 31 of Jan,
%     Mar, May, Jul, Aug, Oct and Dec. Rows 2:end keep every day (Feb has 28).
% (3) FILL VALUE. -0.05 marks night hours with no SW retrieval (it never occurs
%     in a daytime hour). It is set to NaN so it can't enter any mean or box.
%     This doesn't change the Fig. 5 boxes (night boxes sit below ylim 0) but it
%     does change the annual/seasonal SW means printed later in this script.
% NOTE on time base: these CSV hours are PST (UTC-8, no daylight saving); see
% the ncTimeZone setting at the top, which puts the netCDF panels on PST too.
swHourVars = "Var" + (1:24);    % Var1 = hour 0 ... Var24 = hour 23 (PST)
sepSW = SEP_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
octSW = OCT_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
novSW = NOV_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
decSW = DEC_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
janSW = JAN_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
febSW = FEB_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
marSW = MAR_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
aprSW = APR_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
maySW = MAY_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
junSW = JUN_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
julSW = JUL_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance
augSW = AUG_SW_Boxplot{2:end, swHourVars};   % [days x 24 hours], unitless transmittance

% Check that row 1 of every SW CSV really is the hour row 0..23.
swTables = {SEP_SW_Boxplot, OCT_SW_Boxplot, NOV_SW_Boxplot, DEC_SW_Boxplot, JAN_SW_Boxplot, FEB_SW_Boxplot, MAR_SW_Boxplot, APR_SW_Boxplot, MAY_SW_Boxplot, JUN_SW_Boxplot, JUL_SW_Boxplot, AUG_SW_Boxplot};
for k = 1:numel(swTables)
    assert(isequal(swTables{k}{1,swHourVars}, 0:23), ...
        'SW CSV %d: row 1 is not the hour row 0..23; check the file layout.', k)
end

% Replace the -0.05 night fill with NaN (exact match is safe: readtable parses
% the text '-0.05' to the same double as the literal -0.05).
sepSW(sepSW == -0.05) = NaN;
octSW(octSW == -0.05) = NaN;
novSW(novSW == -0.05) = NaN;
decSW(decSW == -0.05) = NaN;
janSW(janSW == -0.05) = NaN;
febSW(febSW == -0.05) = NaN;
marSW(marSW == -0.05) = NaN;
aprSW(aprSW == -0.05) = NaN;
maySW(maySW == -0.05) = NaN;
junSW(junSW == -0.05) = NaN;
julSW(julSW == -0.05) = NaN;
augSW(augSW == -0.05) = NaN;

% --- Kavin's original SW block (kept for reference, not run) ---
% % AJB note (SW hour offset): readtable reads the CSV header row (0,1,...,23,
% % Means) as data, so row 1 = hour labels and Var1 = hour 0 ... Var24 = hour 23
% % (checked in R2024b; Kavin's 'timecheck' line in a later block shows the same).
% % Var2:Var24 therefore puts hours 1-23 in columns 1-23, and column 24 is NaN.
% % Every other diurnal panel puts hour h in column h+1, so in the boxplots below
% % the SW boxes sit ONE HOUR EARLIER than the other panels' hours. On top of that
% % the SW CSV hours are PST (UTC-8, no daylight saving): valid daytime hours are
% % centred on ~11.5 every month with no jump at the DST changes, and June data
% % start in hour 4 (sunrise 4:42 PST = 5:42 PDT). The netCDF panels use
% % DST-aware local time, so during PDT (12 Mar - 5 Nov) SW is 2 h early in total.
% % Left as-is here so this version reproduces Kavin's figure; fixed in v2.
% % AJB note (rows): rows 2:31 keep 30 days, so day 31 of Jan/Mar/May/Jul/Aug/
% % Oct/Dec is dropped. Values of -0.05 are a fill value; they occur only in
% % night hours, so they don't change the daytime boxes, but they do enter the
% % annual/seasonal SW means printed further down.
% sepSW = SEP_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% sepSW(:,24) = nan(30,1);
% octSW = OCT_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% octSW(:,24) = nan(30,1);
% novSW = NOV_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% novSW(:,24) = nan(30,1);
% decSW = DEC_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% decSW(:,24) = nan(30,1);
% janSW = JAN_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% janSW(:,24) = nan(30,1);
% febSW = FEB_SW_Boxplot{2:29,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% febSW(:,24) = nan(28,1);
% marSW = MAR_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% marSW(:,24) = nan(30,1);
% aprSW = APR_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% aprSW(:,24) = nan(30,1);
% maySW = MAY_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% maySW(:,24) = nan(30,1);
% junSW = JUN_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% junSW(:,24) = nan(30,1);
% julSW = JUL_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% julSW(:,24) = nan(30,1);
% augSW = AUG_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
% augSW(:,24) = nan(30,1);

clf
t2 = tiledlayout(2,2,"TileSpacing","compact")
nexttile
hold on
boxchart(sepSW)
boxchart(octSW)
boxchart(novSW)
ylim([0 1])
hold off

nexttile
hold on
boxchart(decSW)
boxchart(janSW)
boxchart(febSW)
ylim([0 1])
hold off

nexttile
hold on
boxchart(marSW)
boxchart(aprSW)
boxchart(maySW)
ylim([0 1])
hold off

nexttile
hold on
boxchart(junSW)
boxchart(julSW)
boxchart(augSW)
ylim([0 1])
hold off
%%
clf
% AJB note: seasons here are Winter = Jan-Mar, Spring = Apr-Jun, Summer =
% Jul-Sep, Fall = Oct-Dec (not DJF/MAM/JJA/SON). The figure caption should say
% so. The SW 'avg' lines take mean(...,3) of a 2-D matrix, which returns the
% input unchanged, so the SW 'averages' are the stacked daily values of the
% three months (boxes span all days), while the LW curves are means of the
% three monthly mean diurnal cycles.
winLWavg = mean([janLW, febLW, marLW],2,'omitmissing');
winSWavg = mean([janSW; febSW(1:28,:); marSW],3,'omitmissing');
sprLWavg = mean([aprLW, mayLW, junLW],2,'omitmissing');
sprSWavg = mean([aprSW; maySW; junSW],3,'omitmissing');
sumLWavg = mean([julLW, augLW, sepLW], 2,'omitmissing');
sumSWavg = mean([julSW; augSW; sepSW],3,'omitmissing');
falLWavg = mean([octLW, novLW, decLW], 2,'omitmissing');
falSWavg = mean([octSW; novSW; decSW],3,'omitmissing');

AnnualLWmean = mean([janLW,febLW,marLW,aprLW,mayLW,junLW,julLW,augLW,sepLW,octLW,novLW,decLW],"all",'omitmissing')
AnnualLWstd = std([janLW,febLW,marLW,aprLW,mayLW,junLW,julLW,augLW,sepLW,octLW,novLW,decLW],[],"all",'omitmissing')

AnnualSWmean = mean([janSW;febSW; marSW; aprSW; maySW; junSW; julSW; augSW; sepSW; octSW; novSW; decSW],'all','omitmissing')
AnnualSWstd = std([janSW;febSW; marSW; aprSW; maySW; junSW; julSW; augSW; sepSW; octSW; novSW; decSW],[],'all','omitmissing')


winLWavgmean = mean(winLWavg,'all','omitmissing')
winLWavgstd = std(winLWavg,[],'all','omitmissing')
winSWavgmean = mean(winSWavg,'all','omitmissing')
winSWavgstd = std(winSWavg,[],'all','omitmissing')

sprLWavgmean = mean(sprLWavg,'all','omitmissing')
sprLWavgstd = std(sprLWavg,[],'all','omitmissing')
sprSWavgmean = mean(sprSWavg,'all','omitmissing')
sprSWavgstd = std(sprSWavg,[],'all','omitmissing')

sumLWavgmean = mean(sumLWavg,'all','omitmissing')
sumLWavgstd = std(sumLWavg,[],'all','omitmissing')
sumSWavgmean = mean(sumSWavg,'all','omitmissing')
sumSWavgstd = std(sumSWavg,[],'all','omitmissing')

falLWavgmean = mean(falLWavg,'all','omitmissing')
falLWavgstd = std(falLWavg,[],'all','omitmissing')
falSWavgmean = mean(falSWavg,'all','omitmissing')
falSWavgstd = std(falSWavg,[],'all','omitmissing')

% AJB edits (v2): no second figure. Kavin's 'f = figure' opened Figure 2 here after the
% 'clf' above had already emptied Figure 1 (fig_work), so Figure 1 was left
% as a blank window. All later blocks, including Fig. 5, now reuse fig_work
% (f is kept as a name for it in case anything refers to f).
% f = figure('Visible',figureVisible);     % Kavin's line (v1 had resetToFactoryDefaults(f) too)
f = fig_work;
t4 = tiledlayout(2,2,'TileSpacing',"tight");

t11 = tiledlayout(t4,1,1);
t11.Layout.Tile = 1
ax11 = axes(t11);
boxchart(ax11,winSWavg(:,1:23))
text(17, .9, 'WINTER')
ax11.XColor ='b';
ax11.YColor ='b';
ax11.Box = 'off';
ax11.XTick = ["1" "5" "10" "15" "20"];
ax11.YLim = [0 1];
ax11.XTickLabel = [];
ylabel(ax11,'SW Transmittance','Color','b')

xtickangle(ax11,0)
ax12 = axes(t11);
plot(ax12,djftime,winLWavg,'r','LineWidth',1)
ax12.XAxisLocation = 'Top';
ax12.XColor = 'r';
ax12.YAxisLocation = 'Right';
ax12.YColor = 'r';
ax12.YLim = [300 410];
ax12.XLim = [0 24];
ax12.XTick = [1 5 10 15 20];
ax12.Color = 'none';
ax12.Box = 'off';
ax12.YTickLabel = [];
xlabel(ax12,'Local Time [Hour]','Color','r')



t12 = tiledlayout(t4,1,1);
t12.Layout.Tile = 2;
ax21 = axes(t12);
boxchart(ax21,sprSWavg(:,1:23))
text(17, .9, 'SPRING')
ax21.XColor = 'b';
ax21.YColor = 'b';
ax21.Box = 'off';
ax21.YLim = [0 1];
ax21.XTick = ["1" "5" "10" "15" "20"];
ax21.XTickLabel = [];
ax21.YTickLabel = [];
xtickangle(ax21,0)
ax22 = axes(t12);
plot(ax22,mamtime,sprLWavg,'r','LineWidth',1)
ax22.XAxisLocation = 'Top';
ax22.XLim = [0 24];
ax22.XTick = [1 5 10 15 20];
ax22.XColor = 'r';
ax22.YAxisLocation = 'Right';
ax22.YColor = 'r';
ax22.YLim = [300 410];
ax22.Color = 'none';
ax22.Box = 'off';
xlabel(ax22,'Local Time [Hour]','Color','r')
ylabel(ax22,'LW Flux [W/m^2]')


t13 = tiledlayout(t4,1,1);
t13.Layout.Tile = 3;
ax31 = axes(t13);
boxchart(ax31,sumSWavg(:,1:23))
text(17, .9, 'SUMMER')
ax31.XColor = 'b';
ax31.YColor = 'b';
ax31.Box = 'off';
ax31.XTick = ["1" "5" "10" "15" "20"];
ax31.YLim = [0 1];
xtickangle(ax31,0)
xlabel(ax31,'Local Time [Hour]','Color','b')
ylabel(ax31,'SW Transmittance','Color','b')

ax32 = axes(t13);
plot(ax32,jjatime,sumLWavg,'r','LineWidth',1)
ax32.XAxisLocation = 'Top';
ax32.XLim = [0 24];
ax32.XTick = [1 5 10 15 20];
ax32.XColor = 'r';
ax32.YAxisLocation = 'Right';
ax32.YColor = 'r';
ax32.Color = 'none';
ax32.Box = 'off';
ax32.YLim = [300 410];
ax32.XTickLabel = [];
ax32.YTickLabel = [];


t14 = tiledlayout(t4,1,1);
t14.Layout.Tile = 4;
ax41 = axes(t14);
boxchart(ax41,falSWavg(:,1:23))
text(17, .9, 'FALL')
ax41.XColor = 'b';
ax41.YColor = 'b';
ax41.Box = 'off';
ax41.XTick = ["1" "5" "10" "15" "20"];
ax41.YLim = [0 1];
ax41.YTickLabel = [];
xtickangle(ax41,0)
xlabel(ax41,'Local Time [Hour]','Color','b')
ax42 = axes(t14);
hold on
plot(ax42,sontime,falLWavg,'r','LineWidth',1)
hold off
ax42.XAxisLocation = 'Top';
ax42.XLim = [0 24];
ax42.XColor = 'r';
ax42.XTick = [1 5 10 15 20];
ax42.XTickLabel = [];
ax42.YAxisLocation = 'Right';
ax42.YColor = 'r';
ax42.YLim = [300 410];
ax42.Color = 'none';
ax42.Box = 'off';
ylabel(ax42,'LW Flux [W/m^2]',Color='r')
%%
clf
ta = tiledlayout(2,2,'TileSpacing','compact')
nexttile
yyaxis left
boxchart(winSWavg)
xtickangle(0)
ylim([0 1])
ylabel('SW Cloud Transmittance')
yyaxis right
plot(winLWavg)
ylim([300 400])
ylabel('LW Downwelling Flux (W/m^-2)')
xlabel('Local Time (Hour)')
text(20,390,'Winter')

nexttile
yyaxis left
boxchart(sprSWavg)
xtickangle(0)
ylim([0 1])
ylabel('SW Cloud Transmittance')
yyaxis right
plot(sprLWavg)
ylim([300 400])
ylabel('LW Downwelling Flux (W/m^-2)')
xlabel('Local Time (Hour)')
text(20,390,'Spring')


nexttile
yyaxis left
boxchart(sumSWavg)
xtickangle(0)
ylim([0 1])
ylabel('SW Cloud Transmittance')
yyaxis right
plot(sumLWavg)
ylim([300 400])
ylabel('LW Downwelling Flux (W/m^-2)')
xlabel('Local Time (Hour)')
text(20,390,'Summer')

nexttile
yyaxis left
boxchart(falSWavg)
ylim([0 1])
xtickangle(0)
ylabel('SW Cloud Transmittance')
yyaxis right
plot(falLWavg)
ylim([300 400])
ylabel('LW Downwelling Flux (W/m^-2)')
xlabel('Local Time (Hour)')
text(20,390,'Fall')
%%

%%

%%
timecheck = SEP_SW_Boxplot{1,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]},

clf
t5 = tiledlayout(4,2,"TileSpacing","tight");
nexttile
hold on
boxchart(janSW,'BoxFaceColor','[0.8500, 0.3250, 0.0980 ]','MarkerColor','[0.8500, 0.3250, 0.0980 ]')
boxchart(febSW,'BoxFaceColor','#66CCFF','MarkerColor','#66CCFF')
boxchart(marSW,'BoxFaceColor','#0066FF','MarkerColor','#0066FF')
hold off
xtickangle(0)
title('Shortwave Cloud Transmittance')
legend('JAN','FEB','MAR')
ylim([0 1])
ax = gca;
ax.YColor = 'k';
xticks([])
ylabel('Transmittance')

nexttile
hold on
plot(djftime,janLW,'Color','[0.8500, 0.3250, 0.0980 ]','LineWidth',1.5)
plot(djftime,febLW,'Color','#66CCFF','LineWidth',1.5)
plot(djftime,marLW,'Color','#0066FF','LineWidth',1.5)
hold off
ylim([280 420])
title('Longwave Downwelling Flux')
yticks([])
yyaxis right
ylim([280 420])
ax = gca;
ax.YColor = 'k';
xlim([0 24])
xticks([1 5 10 15 20])
xticks([])
ylabel('Flux [W/m^2]')


nexttile
hold on
boxchart(aprSW,'BoxFaceColor','[0.8500, 0.3250, 0.0980 ]','MarkerColor','[0.8500, 0.3250, 0.0980 ]')
boxchart(maySW,'BoxFaceColor','#66CCFF','MarkerColor','#66CCFF')
boxchart(junSW,'BoxFaceColor','#0066FF','MarkerColor','#0066FF')
hold off
xtickangle(0)
legend("APR","MAY","JUN")
ylim([0 1])
ax = gca;
ax.YColor = 'k';
xticks([])
ylabel('Transmittance')


nexttile
hold on
plot(mamtime, aprLW,'Color','[0.8500, 0.3250, 0.0980 ]','LineWidth',1.5)
plot(mamtime, mayLW,'Color','#66CCFF','LineWidth',1.5)
plot(mamtime,junLW,'Color','#0066FF','LineWidth',1.5)
hold off
ylim([280 420])
yticks([])
yyaxis right
ylim([280 420])
ax = gca;
ax.YColor = 'k';
xlim([0 24])
xticks([1 5 10 15 20])
xticks([])
ylabel('Flux [W/m^2]')

nexttile
hold on
boxchart(julSW,'BoxFaceColor','[0.8500, 0.3250, 0.0980 ]','MarkerColor','[0.8500, 0.3250, 0.0980 ]')
boxchart(augSW,'BoxFaceColor','#66CCFF','MarkerColor','#66CCFF')
boxchart(sepSW,'BoxFaceColor','#0066FF','MarkerColor','#0066FF')
hold off
xtickangle(0)
legend('JUL','AUG','SEP')
ylim([0 1])
ax = gca;
ax.YColor = 'k';
xticks([])
ylabel('Transmittance')

nexttile
hold on
plot(jjatime,julLW,'Color','[0.8500, 0.3250, 0.0980 ]','LineWidth',1.5)
plot(jjatime,augLW,'Color','#66CCFF','LineWidth',1.5)
plot(jjatime,sepLW,'Color','#0066FF','LineWidth',1.5)
hold off
yticks([])
yyaxis right
ylim([280 420])
ax = gca;
ax.YColor = 'k';
xlim([0 24])
xticks([1 5 10 15 20])
xticks([])
ylabel('Flux [W/m^2]')


nexttile
hold on
boxchart(octSW,'BoxFaceColor','[0.8500, 0.3250, 0.0980 ]','MarkerColor','[0.8500, 0.3250, 0.0980 ]')
boxchart(novSW,'BoxFaceColor','#66CCFF','MarkerColor','#66CCFF')
boxchart(decSW,'BoxFaceColor','#0066FF','MarkerColor','#0066FF')
hold off
xtickangle(0)
xticks(["1" "5" "10" "15" "20" ])
legend('OCT','NOV','DEC')
ylim([0 1])
ax = gca;
ax.YColor = 'k';
ylabel('Transmittance')

nexttile
hold on
plot(sontime,octLW,'Color','[0.8500, 0.3250, 0.0980 ]','LineWidth',1.5)
plot(sontime,novLW,'Color','#66CCFF','LineWidth',1.5)
plot(sontime,decLW,'Color','#0066FF','LineWidth',1.5)
hold off
yticks([])
yyaxis right
ylim([280 420])
ax = gca;
ax.YColor = 'k';
xlim([0 24])
xticks([1 5 10 15 20])
ylabel('Flux [W/m^2]')


xlabel(t5,'Local Time (Hour)','FontSize',13)
%%
clf
ncdisp(fullfile(dataDir,"arscl_echomask_hourly.nc"));  % AJB edits: read from dataDir (was a bare file name)
time = ncread(fullfile(dataDir,"arscl_echomask_hourly.nc"),"time");  % AJB edits: read from dataDir (was a bare file name)
d = datetime(time, 'ConvertFrom',"epochtime", "Epoch", datetime(1970,1,1,00,00,00),"TicksPerSecond",1, "Format","default",'TimeZone','UTC');
% AJB edits (v2): time zone from the ncTimeZone setting (was 'America/Los_Angeles').
% PST (fixed UTC-8) puts the netCDF hours on the same clock as the SW CSVs.
d = datetime(d,'TimeZone',ncTimeZone);  % AJB edits: added ';' to suppress a very long printout (no change to values)
% AJB note: the result of dateshift below is not assigned to anything, so it
% has no effect on d (the times stay at the hour stamp of the hourly file).
dateshift(d,"end","hour");  % AJB edits: added ';' to suppress a very long printout (no change to values)
height = ncread(fullfile(dataDir,"arscl_echomask_hourly.nc"),"height");  % AJB edits: read from dataDir (was a bare file name)
echo_fraction = ncread(fullfile(dataDir,"arscl_echomask_hourly.nc"), "echo_fraction");  % AJB edits: read from dataDir (was a bare file name)

ncdisp(fullfile(dataDir,"epcceil_hourly_cbh.nc"));  % AJB edits: read from dataDir (was a bare file name)
time2 = ncread(fullfile(dataDir,"epcceil_hourly_cbh.nc"),"time");  % AJB edits: read from dataDir (was a bare file name)
% AJB note: d2 is built from 'time' (ARSCL) instead of 'time2' (ceilometer),
% and the next line converts d, not d2. Harmless here: the two files have
% identical time vectors (8809 hourly stamps, 15 Feb 2023 - 17 Feb 2024).
% AJB edits (v2): d2 now comes from the ceilometer's own time vector 'time2' and is
% converted from d2 (Kavin used 'time' and 'd', the ARSCL times). The values
% are identical today (same 8809 hourly stamps), so the figure doesn't change,
% but this keeps CBH tied to its own clock if either file is ever regenerated.
% AJB edits (v2): time zone from the ncTimeZone setting (was 'America/Los_Angeles').
% PST (fixed UTC-8) puts the netCDF hours on the same clock as the SW CSVs.
d2 = datetime(time2, 'ConvertFrom',"epochtime", "Epoch", datetime(1970,1,1,00,00,00),"TicksPerSecond",1, "Format","default",'TimeZone','UTC');
d2 = datetime(d2,'TimeZone',ncTimeZone);  % AJB edits: added ';' to suppress a very long printout (no change to values)
bl_cbh = ncread(fullfile(dataDir,"epcceil_hourly_cbh.nc"),'bl_cbh');  % AJB edits: read from dataDir (was a bare file name)

sprcbheight = nan(368,24);
for i = 0:24
    heldsize = size(bl_cbh(hour(d2) == i & month(d2) >= 4 & month(d2) <= 6),1);
    sprcbheight(1:heldsize,i+1) = bl_cbh(hour(d2) == i & month(d2) >= 4 & month(d2) <= 6);
end


ncdisp(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"))  % AJB edits: read from dataDir (was a bare file name)
time7913 = ncread(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"),"time");  % AJB edits: read from dataDir (was a bare file name)
d7913 = datetime(time7913, 'ConvertFrom',"epochtime", "Epoch", datetime(1970,1,1,00,00,00),"TicksPerSecond",1, "Format","default",'TimeZone','UTC');
% AJB edits (v2): time zone from the ncTimeZone setting (was 'America/Los_Angeles').
% PST (fixed UTC-8) puts the netCDF hours on the same clock as the SW CSVs.
d7913 = datetime(d7913,'TimeZone',ncTimeZone);

lwp7913 = ncread(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"),"phys_lwp");  % AJB edits: read from dataDir (was a bare file name)
sprlwpbox = nan(330,24);
for i = 0:24
    heldsize = size(lwp7913(hour(d7913) == i & month(d7913) >= 4 & month(d7913) <= 6),1);
    sprlwpbox(1:heldsize,i+1) = lwp7913(hour(d7913) == i & month(d7913) >= 4 & month(d7913) <= 6);
end

tbl = readtable(fullfile(dataDir,"EPCAPE_CTH_CBH_PBLH_SChellappan_LaRC.csv"));  % AJB edits: read from dataDir (was a bare file name)

monthOrder = {'1','2','3','4','5','6','7','8','9','10','11','12'};


timesche= tbl{:,'Time'};
tbl.("month") = month(timesche);  % AJB edits: added ';' to suppress a very long printout (no change to values)
%+ (year(timesche)-2023).*12

% AJB edits (v2): Kavin's Fig. 4-only inputs (BL depth, sonde LCL, combined instrument
% table) were here. Deleted (Andrew, 2026-10-13); they are still in v1 and in
% Kavin's original. tbl (above) stays: it feeds the CBH statistics below.
%%
% AJB note: 'var' (and 'sum' just below, 'min'/'max' in the rainfall block,
% 'std2' at the end) shadow built-in MATLAB functions. Nothing later calls
% those functions, so results are unaffected, but avoid reusing the names.
var = tbl.CEIL_CBH;  % AJB edits: added ';' to suppress a very long printout (no change to values)

win = tbl.month >= 1 & tbl.month <= 3;
spr = tbl.month >= 4 & tbl.month <= 6;
sum = tbl.month >= 7 & tbl.month <= 9;
fal = tbl.month >= 10 & tbl.month <= 12;


mean(var,'omitmissing')
std(var,'omitmissing')

mean(var(spr),'omitmissing')
std(var(spr),'omitmissing')

mean(var(sum),'omitmissing')
std(var(sum),'omitmissing')

mean(var(fal),'omitmissing')
std(var(fal),'omitmissing')

mean(var(win),'omitmissing')
std(var(win),'omitmissing')

%%

% Set up constants used in figure

% Colour keys. 'SW'/'LW' copy the keys on the slide version (Winter, bottom
% row = panel L). 'CBH' replaces the removed right-axis label (panel D).
% Colours are the same RGB triplets Kavin used for the data.
swColor  = [0, 0.4470, 0.7410];        % SW boxes (MATLAB default blue)
lwColor  = [0.6350, 0.0780, 0.1840];   % LW line (dark red)
cbhColor = [0.8500, 0.3250, 0.0980];   % CBH boxes (orange)

ax_fnt = 20;
cb_fnt = 16;

% AJB note: THIS BLOCK PRODUCES BAMS FIGURE 5 (seasonal diurnal cycle of
% echo fraction + ceilometer CBH, MWR LWP, SW transmittance + LW flux).
% The slide version differs only cosmetically (black 'Height [m]' label, no
% '23' tick label, 'SW'/'LW' text in the Winter bottom panel), so the final
% slide was made from a slightly later revision of this block.
% AJB edits (v2): all panels are now on one time base, PST (UTC-8), set by ncTimeZone.
% v1 mixed DST-aware local time (netCDF) with PST (SW CSV), so SW was 1 h early
% Nov-Mar and 2 h early in PDT months. The LW CSV time base can't be checked
% from the data; it is assumed to match the SW CSVs (same source). Ask Kavin.
clf
% AJB edits: fix the figure size so the export matches the slide's aspect
% ratio (slide image is 1365 x 822 px). The Live Editor otherwise uses its
% own small output size, which squeezes the tick labels together.
set(gcf,'Units','pixels','Position',[50 50 1365 822]);
sprechomap = zeros(596,24);
for i = 0:24
    sprechomap(:,i+1) = mean(echo_fraction(:,hour(d) == i & month(d) <= 6 & month(d) >= 4),2,'omitmissing');

end

echospr = tiledlayout(4,4,'TileSpacing','tight');

nexttile(echospr,1,[2 1])
yyaxis right
boxplot(sprcbheight,'Symbol','','Colors',[0.8500, 0.3250, 0.0980],'Widths',0.6)
xlim([0.5 24.5])
ylim([0 3000])
xticklabels([])
yticklabels([])

yyaxis left
pcolor(0.5:24.5,height,sprechomap,'EdgeColor','none')
colormap('sky')
ylim([0 3000])
title('Spring')
% ylabel('Height [m]')
ylabel('Cloud Base Height [m]', 'FontSize', ax_fnt, 'Color', cbhColor)


nexttile(echospr,9,[1 1])
boxplot(sprlwpbox(:,1:24),'Symbol','')
xticklabels([])
% AJB edits (v2): unit written as g m^{-2} to match the W m^{-2} style (was g/m^2).
ylabel('Liquid Water Path [g m^{-2}]', 'FontSize', ax_fnt-4)
ylim([0 500])

nexttile(echospr,13,[1,1])
yyaxis right
boxplot(sprSWavg,'Symbol','','Colors',[0, 0.4470, 0.7410])
yticklabels([])
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" ""])  % AJB edits (v2): '23' label dropped, as on the slide version
xtickangle(0)
ylim([0 1])

yyaxis left
plot(mamtime+1,sprLWavg,Color=[0.6350, 0.0780, 0.1840])
% AJB edits (v2): Lynn: 'Add units to LW Flux and Height if they fit'. LW Flux gets
% [W m^{-2}] (downwelling LW irradiance at the surface); 'Height [m]' already had
% units. Black text as on the slide; the red 'LW' key identifies the line.
ylabel('LW_{\downarrow} Flux [W m^{-2}]','Color',lwColor, 'FontSize', ax_fnt-4)
ax = gca;
ax.YColor = 'k';
ylim([300 410])

xlabel(timeAxisLabel, 'FontSize', ax_fnt)  % AJB edits (v2): label set from ncTimeZone (was 'Local Time [Hour]')
sumcbheight = nan(368,24);
for i = 0:24
    heldsize = size(bl_cbh(hour(d2) == i & month(d2) >= 7 & month(d2) <= 9),1);
    sumcbheight(1:heldsize,i+1) = bl_cbh(hour(d2) == i & month(d2) >= 7 & month(d2) <= 9);
end




sumlwpbox = nan(330,24);
for i = 0:24
    heldsize = size(lwp7913(hour(d7913) == i & month(d7913) >= 7 & month(d7913) <= 9),1);
    sumlwpbox(1:heldsize,i+1) = lwp7913(hour(d7913) == i & month(d7913) >= 7 & month(d7913) <= 9);
end


sumechomap = zeros(596,24);
for i = 0:24
    sumechomap(:,i+1) = mean(echo_fraction(:,hour(d) == i & month(d) >= 7 & month(d) <= 9),2,'omitmissing');

end

%echosum = tiledlayout(4,2,'TileSpacing','compact');

nexttile(echospr,2,[2 1])
yyaxis right
boxplot(sumcbheight,'Symbol','','Colors',[0.8500, 0.3250, 0.0980],'Widths',0.3)
ylim([0 3000])
xlim([0.5 24.5])
xticklabels([])
yticklabels([])

yyaxis left
pcolor(0.5:24.5,height,sumechomap,'EdgeColor','none')
ylim([0 3000])
colormap('sky')
title('Summer')
yticklabels([])

nexttile(echospr,10,[1 1])
boxplot(sumlwpbox(:,1:24),'Symbol','')
xticklabels([])
ylim([0 500])
yticklabels([])

nexttile(echospr,14,[1,1])
yyaxis right
boxplot(sumSWavg,'Symbol','','Colors',[0, 0.4470, 0.7410])
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" ""])  % AJB edits (v2): '23' label dropped, as on the slide version
xtickangle(0)
yticklabels([])
ylim([0 1])

yyaxis left
plot(mamtime+1,sumLWavg,Color=[0.6350, 0.0780, 0.1840])
yticklabels([])
ylim([300 410])

xlabel(timeAxisLabel, 'FontSize', ax_fnt)  % AJB edits (v2): label set from ncTimeZone (was 'Local Time [Hour]'); FontSize added to match the other 3 columns


falcbheight = nan(368,24);
for i = 0:24
    heldsize = size(bl_cbh(hour(d2) == i & month(d2) >= 10 & month(d2) <= 12),1);
    falcbheight(1:heldsize,i+1) = bl_cbh(hour(d2) == i & month(d2) >= 10 & month(d2) <= 12);
end



lwp7913 = ncread(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"),"phys_lwp");  % AJB edits: read from dataDir (was a bare file name)
fallwpbox = nan(330,24);
for i = 0:24
    heldsize = size(lwp7913(hour(d7913) == i & month(d7913) >= 10 & month(d7913) <= 12),1);
    fallwpbox(1:heldsize,i+1) = lwp7913(hour(d7913) == i & month(d7913) >= 10 & month(d7913) <= 12);
end


falechomap = zeros(596,24);
for i = 0:24
    falechomap(:,i+1) = mean(echo_fraction(:,hour(d) == i & month(d) >= 10 & month(d) <= 12),2,'omitmissing');

end


nexttile(echospr,3,[2 1])
yyaxis right
boxplot(falcbheight,'Symbol','','Colors',[0.8500, 0.3250, 0.0980],'Widths',0.3)
ylim([0 3000])
xlim([0.5 24.5])
xticklabels([])
yticklabels([])

yyaxis left
pcolor(0.5:24.5,height,falechomap,'EdgeColor','none')
ylim([0 3000])
colormap('sky')
title('Fall')
yticklabels([])

nexttile(echospr,11,[1 1])
boxplot(fallwpbox(:,1:24),'Symbol','','Colors',[0, 0.4470, 0.7410])
xticklabels([])
ylim([0 500])
yticklabels([])


nexttile(echospr,15,[1,1])
yyaxis right
boxplot(falSWavg,'Symbol','','Colors',[0, 0.4470, 0.7410])
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" ""])  % AJB edits (v2): '23' label dropped, as on the slide version
xtickangle(0)
yticklabels([])
ylim([0 1])


yyaxis left
plot(mamtime+1,falLWavg, Color=[0.6350, 0.0780, 0.1840])
ylim([300 410])
yticklabels([])

xlabel(timeAxisLabel, 'FontSize', ax_fnt)  % AJB edits (v2): label set from ncTimeZone (was 'Local Time [Hour]')

wincbheight = nan(368,24);
for i = 0:24
    heldsize = size(bl_cbh(hour(d2) == i & month(d2) >= 1 & month(d2) <= 3),1);
    wincbheight(1:heldsize,i+1) = bl_cbh(hour(d2) == i & month(d2) >= 1 & month(d2) <= 3);
end


lwp7913 = ncread(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"),"phys_lwp");  % AJB edits: read from dataDir (was a bare file name)
winlwpbox = nan(330,24);
for i = 0:24
    heldsize = size(lwp7913(hour(d7913) == i & month(d7913) >= 1 & month(d7913) <= 3),1);
    winlwpbox(1:heldsize,i+1) = lwp7913(hour(d7913) == i & month(d7913) >= 1 & month(d7913) <= 3);
end


winechomap = zeros(596,24);
for i = 0:24
    winechomap(:,i+1) = mean(echo_fraction(:,hour(d) == i & month(d) >= 1 & month(d) <= 3),2,'omitmissing');

end


nexttile(echospr,4,[2 1])
yyaxis right
boxplot(wincbheight,'Symbol','','Colors',[0.8500, 0.3250, 0.0980],'Widths',0.3)
% AJB edits (v2): Lynn: 'remove label for Cloud Base Height on right' (and 'why is right
% axis label needed?'). It isn't: the right axis has no tick labels and the same
% 0-3000 m range as the left 'Height [m]' axis, so the CBH boxes are read off the
% left axis. An orange 'CBH' key is added in panel D instead (see below).
% ylabel('Cloud Base Height')
xlim([0.5 24.5])
xticklabels([])
yticklabels([])
ylim([0 3000])

yyaxis left
pcolor(0.5:24.5,height,winechomap,'EdgeColor','none')
ylim([0 3000])
colormap('sky')
c = colorbar;
clim([0 0.5])
ylim([0 3000])
title('Winter')
yticklabels([])

% AJB edits (v2): Lynn: 'Change echo to Hydrometeor'. ARSCL echo_fraction = fraction of
% valid samples in the hour with a radar and/or lidar return (cloud or
% precipitation), per the file's description attribute.
ylabel(c, "Hydrometeor Fraction", "FontSize", cb_fnt)

nexttile(echospr,12,[1 1])
boxplot(winlwpbox(:,1:24),'Symbol','')
xticklabels([])
yticklabels([])
ylim([0 500])

nexttile(echospr,16,[1,1])
yyaxis right
boxplot(winSWavg,'Symbol','','Colors',[0, 0.4470, 0.7410])
% AJB edits (v2): Lynn: 'Typo in Transmittance' (was 'SW Trasmittance'). Text is black, as
% on the slide version; the blue 'SW' key added below identifies the boxes.
ylabel('SW Transmittance', 'FontSize', ax_fnt, 'Color', swColor)
ax = gca;
ax.YColor = 'k';
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" ""])  % AJB edits (v2): '23' label dropped, as on the slide version
xtickangle(0)
ylim([0 1])

yyaxis left
plot(mamtime+1,winLWavg,Color=[0.6350, 0.0780, 0.1840])
yticklabels([])
ylim([300 410])

xlabel(timeAxisLabel, 'FontSize', ax_fnt)  % AJB edits (v2): label set from ncTimeZone (was 'Local Time [Hour]')

% AJB edits (v2): finishing touches applied once all 12 panels exist.
% Collect the 12 axes of the tiled layout and sort them by tile number:
% tiles 1-4 = echo/CBH row, 9-12 = LWP row, 13-16 = SW/LW row. (The colorbar
% is a 'colorbar' object, so findobj(...,'Type','axes') skips it.)
fig5Axes = findobj(echospr,'Type','axes');
tileNum  = arrayfun(@(a) a.Layout.Tile, fig5Axes);
[~,order] = sort(tileNum);
fig5Axes = fig5Axes(order);
assert(numel(fig5Axes) == 12, 'Expected 12 panels in Fig. 5, found %d.', numel(fig5Axes))

% Black axes: the slide version of Fig. 5 has black axis lines, ticks and
% labels. Colour now identifies the data through the in-panel keys below.
for k = 1:numel(fig5Axes)
    for iy = 1:numel(fig5Axes(k).YAxis)        % 2 y axes on yyaxis panels
        fig5Axes(k).YAxis(iy).Color = 'k';
    end
end
% Setting an axis colour also recolours its label, which turned Andrew's orange
% 'Cloud Base Height [m]' label black. Re-apply the orange to that one label
% (Spring panel, left y axis).
fig5Axes(1).YAxis(1).Label.Color = cbhColor;

% Lynn: 'Add panel labels -- A, B, etc.' One letter per ROW, in the first
% (Spring) column only (Andrew, 2026-10-13): A = hydrometeor fraction + CBH,
% B = LWP, C = SW/LW. The columns are already labelled by season, so e.g.
% "panel B, Winter" is unambiguous. (The first version labelled all 12 panels
% A-L.) fig5Axes is sorted by tile, so elements 1, 5 and 9 are tiles 1, 9 and
% 13: the Spring panel of each row. Plain 14 pt, not bold, top-left corner (as
% Lynn asked for in BAMS Fig. 8). A white box behind each letter keeps it
% readable where whiskers reach the top of the panel (row B, hours 0-1).
panelLetters = 'ABC';
firstColumnAxes = fig5Axes([1 5 9]);
for k = 1:numel(firstColumnAxes)
    text(firstColumnAxes(k), 0.02, 0.97, panelLetters(k), 'Units','normalized', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'FontSize',14, 'Color','k', 'BackgroundColor','w', 'Margin',1);
end



% text(fig5Axes(12), 0.97, 0.90, 'SW', 'Units','normalized', 'HorizontalAlignment','right', ...
%     'VerticalAlignment','top', 'FontSize',14, 'Color',swColor);
% text(fig5Axes(12), 0.97, 0.72, 'LW', 'Units','normalized', 'HorizontalAlignment','right', ...
%     'VerticalAlignment','top', 'FontSize',14, 'Color',lwColor);
% text(fig5Axes(4), 0.97, 0.97, 'CBH', 'Units','normalized', 'HorizontalAlignment','right', ...
%     'VerticalAlignment','top', 'FontSize',14, 'Color',cbhColor);

% AJB edits: export Fig. 5 (Kavin copied figures out of the Live Editor by
% hand). PNG at 300 dpi for slides; PDF as vector for the manuscript.
% (v2: new file names, matching this script's name.)
if exportFigs
    exportgraphics(gcf, fullfile(outDir,'EPCAPEDanLWSWtransmittance_BAMSFig5_AJB_v2.png'),'Resolution',300);
    exportgraphics(gcf, fullfile(outDir,'EPCAPEDanLWSWtransmittance_BAMSFig5_AJB_v2.pdf'),'ContentType','vector');
end
%%

%%
% AJB edits (v2): Kavin's combined Spring+Summer test figure was here. Deleted: not a
% BAMS figure, and it errors once day 31 is kept (still in v1 / the original).
%%

% AJB edits (v2): corrected LWP statistics (replaces Kavin's mean/std lines that were
% inside the old Fig. 4 block, now deleted). Two problems in the old lines:
%   (a) yearlylwpbox was pre-allocated too short (500 rows < 744 hours in a
%       month), so MATLAB grew it with ZEROS (1023 fake 0 g m^-2 values);
%   (b) column i+1 held month i, so the 'season' columns 4:6, 7:9, 10:12 and
%       1:3 were Mar-May, Jun-Aug, Sep-Nov and (empty, Jan, Feb).
% Here each season is selected directly from the hourly MWR record by month,
% with the same season definitions as Fig. 5. Units: g m^-2. N = number of
% hourly averages with a valid (non-NaN) LWP.
seasonNames  = ["All data"; "Spring (Apr-Jun)"; "Summer (Jul-Sep)"; "Fall (Oct-Dec)"; "Winter (Jan-Mar)"];
seasonMonths = {1:12, 4:6, 7:9, 10:12, 1:3};
lwpMean_gm2 = nan(numel(seasonNames),1);
lwpStd_gm2  = nan(numel(seasonNames),1);
lwpN_hours  = nan(numel(seasonNames),1);
for k = 1:numel(seasonNames)
    inSeason = ismember(month(d7913), seasonMonths{k}) & ~isnan(lwp7913);
    lwpMean_gm2(k) = mean(lwp7913(inSeason));
    lwpStd_gm2(k)  = std(lwp7913(inSeason));
    lwpN_hours(k)  = nnz(inSeason);
end
lwpStats = table(seasonNames, lwpMean_gm2, lwpStd_gm2, lwpN_hours)
%%
% AJB edits (v2): Kavin's old calendar-month Fig. 4 code was here. Deleted
% (Andrew, 2026-10-13): the slide's Fig. 4 is made by other code. It is still in
% v1 and in Kavin's original; its bugs are listed in the change-log text file.

%% Local functions
% AJB edits: new local function (Kavin's original had none).

function resetToFactoryDefaults(fig)
% resetToFactoryDefaults  Give one figure MATLAB's factory graphics defaults.
%
% Andrew's startup.m sets ROOT defaults (black background, white text,
% 18 pt bold axes fonts, 4 pt lines, 15 pt markers). Kavin's code assumes the
% factory defaults, so on Andrew's machine the figures would come out with
% white-on-white text and oversized fonts. Setting each changed default back
% to its factory value AT THE FIGURE LEVEL fixes this figure only; the root
% defaults (and any other figure) are left alone.
%
% defaultFigure* properties are skipped because they cannot be set on a figure
% that already exists (they only apply to figures created later).
rootDefaults = fieldnames(get(groot,'default'));
for k = 1:numel(rootDefaults)
    name = rootDefaults{k};                         % e.g. 'defaultAxesFontSize'
    if startsWith(name,'defaultFigure','IgnoreCase',true)
        continue
    end
    factoryName = regexprep(name,'^default','factory','ignorecase');
    set(fig, name, get(groot, factoryName));        % figure-level override
end
% The figure's own colour is a defaultFigure* property, so set it directly.
fig.Color = get(groot,'factoryFigureColor');
end
