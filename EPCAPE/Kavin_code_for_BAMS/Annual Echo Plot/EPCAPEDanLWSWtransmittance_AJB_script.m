%% EPCAPEDanLWSWtransmittance_AJB
% AJB version of Kavin's live script EPCAPEDanLWSWtransmittance.mlx (saved in
% R2025b, folder "Annual Echo Plot"). Kavin's original is kept untouched.
%
% WHAT THIS CODE MAKES (checked against EPCAPEoverviewFigs261005.pptx):
%   - BAMS Figure 5 (slide 5): the "Spring / Summer / Fall / Winter" block
%     (tiledlayout 'echospr'). Same data and layout as the slide; the slide
%     differs only cosmetically, so it came from a slightly later revision.
%   - A PRECURSOR of BAMS Figure 4 (slide 4): the last block ('echoyearly'),
%     binned by calendar month Jan-Dec with 5 rows. The slide's Figure 4 runs
%     Feb 23 - Feb 24 and has a UCVF/RH row; its code is not in this folder.
%     (The 3-row "Annual Echo Plot" image on slide 4 is parked off the slide
%     and is an older version again.)
%   - Several exploratory figures, each wiped by the next 'clf'.
%
% HOW TO RUN: set dataDir/outDir below (Andrew's path is used when it exists;
% otherwise the current folder, as in Kavin's original) and run the script.
% Figure 5 and the Fig. 4 precursor are exported to outDir as
% EPCAPEDanLWSWtransmittance_AJB_Fig5.png/.pdf and _Fig4precursor.png.
%
% EDIT TAGS (search for these to see every change):
%   % AJB edits: ...  code was changed or added (paths, figure handling,
%                     exports, suppressed printouts). No data values change.
%   % AJB note:  ...  comment only; Kavin's code below it is unchanged. These
%                     flag issues for Kavin. The biggest ones for the paper:
%                     (1) Fig. 5 SW boxes are plotted 1 h early (Nov-Mar)
%                         to 2 h early (PDT months) relative to the other panels;
%                     (2) Fig. 4 precursor: zero padding biases CBH/LWP boxes
%                         and the printed LWP means, and the echo map is
%                         shifted one month from the CBH boxes;
%                     (3) seasons are JFM/AMJ/JAS/OND, not DJF/MAM/JJA/SON.
%
% Andrew J. Buggee, Scripps Institution of Oceanography, 2026-10-13

% AJB edits: settings block (new). Kavin's original read every file from
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
% AJB note (SW hour offset): readtable reads the CSV header row (0,1,...,23,
% Means) as data, so row 1 = hour labels and Var1 = hour 0 ... Var24 = hour 23
% (checked in R2024b; Kavin's 'timecheck' line in a later block shows the same).
% Var2:Var24 therefore puts hours 1-23 in columns 1-23, and column 24 is NaN.
% Every other diurnal panel puts hour h in column h+1, so in the boxplots below
% the SW boxes sit ONE HOUR EARLIER than the other panels' hours. On top of that
% the SW CSV hours are PST (UTC-8, no daylight saving): valid daytime hours are
% centred on ~11.5 every month with no jump at the DST changes, and June data
% start in hour 4 (sunrise 4:42 PST = 5:42 PDT). The netCDF panels use
% DST-aware local time, so during PDT (12 Mar - 5 Nov) SW is 2 h early in total.
% Left as-is here so this version reproduces Kavin's figure; fixed in v2.
% AJB note (rows): rows 2:31 keep 30 days, so day 31 of Jan/Mar/May/Jul/Aug/
% Oct/Dec is dropped. Values of -0.05 are a fill value; they occur only in
% night hours, so they don't change the daytime boxes, but they do enter the
% annual/seasonal SW means printed further down.
sepSW = SEP_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
sepSW(:,24) = nan(30,1);
octSW = OCT_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
octSW(:,24) = nan(30,1);
novSW = NOV_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
novSW(:,24) = nan(30,1);
decSW = DEC_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
decSW(:,24) = nan(30,1);
janSW = JAN_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
janSW(:,24) = nan(30,1);
febSW = FEB_SW_Boxplot{2:29,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
febSW(:,24) = nan(28,1);
marSW = MAR_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
marSW(:,24) = nan(30,1);
aprSW = APR_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
aprSW(:,24) = nan(30,1);
maySW = MAY_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
maySW(:,24) = nan(30,1);
junSW = JUN_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
junSW(:,24) = nan(30,1);
julSW = JUL_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
julSW(:,24) = nan(30,1);
augSW = AUG_SW_Boxplot{2:31,["Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22","Var23","Var24"]};
augSW(:,24) = nan(30,1);

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

% AJB edits: same treatment for Kavin's second figure: visibility from the
% settings block and factory graphics defaults. All later blocks (including
% Fig. 5 and the Fig. 4 precursor) draw into this figure after clf.
f = figure('Visible',figureVisible);
resetToFactoryDefaults(f);
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
d = datetime(d,'TimeZone','America/Los_Angeles');  % AJB edits: added ';' to suppress a very long printout (no change to values)
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
d2 = datetime(time, 'ConvertFrom',"epochtime", "Epoch", datetime(1970,1,1,00,00,00),"TicksPerSecond",1, "Format","default",'TimeZone','UTC');
d2 = datetime(d,'TimeZone','America/Los_Angeles');  % AJB edits: added ';' to suppress a very long printout (no change to values)
bl_cbh = ncread(fullfile(dataDir,"epcceil_hourly_cbh.nc"),'bl_cbh');  % AJB edits: read from dataDir (was a bare file name)

sprcbheight = nan(368,24);
for i = 0:24
    heldsize = size(bl_cbh(hour(d2) == i & month(d2) >= 4 & month(d2) <= 6),1);
    sprcbheight(1:heldsize,i+1) = bl_cbh(hour(d2) == i & month(d2) >= 4 & month(d2) <= 6);
end


ncdisp(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"))  % AJB edits: read from dataDir (was a bare file name)
time7913 = ncread(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"),"time");  % AJB edits: read from dataDir (was a bare file name)
d7913 = datetime(time7913, 'ConvertFrom',"epochtime", "Epoch", datetime(1970,1,1,00,00,00),"TicksPerSecond",1, "Format","default",'TimeZone','UTC');
d7913 = datetime(d7913,'TimeZone','America/Los_Angeles');

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

tbl2 = readtable(fullfile(dataDir,"BL_Depth.csv"));  % AJB edits: read from dataDir (was a bare file name)
timesche2= tbl2{:,'Time'};
tbl2.("month") = month(timesche2);  % AJB edits: added ';' to suppress a very long printout (no change to values)
%+ (year(timesche2)-2023).*12;

tbl3 = readtable(fullfile(dataDir,"EPCAPE_sonde_LCL.csv"));  % AJB edits: read from dataDir (was a bare file name)
timesche3 = tbl3{:,'time'};
tbl3.("month") = month(timesche3);  % AJB edits: added ';' to suppress a very long printout (no change to values)
%+ (year(timesche3)-2023).*12;
tbl.month;
tbl2.month;
tbl3.month;

newmonth = cat(1,tbl.month,tbl.month,tbl.month,tbl.month, tbl2.month,tbl3.month, tbl.month);
newdata = cat(1,tbl.GOES18_CTH,tbl.KAZR_CTH,tbl.PBLH_THERMO,tbl.PBLH_MPL,tbl2.Depth,tbl3.LCL,tbl.CEIL_CBH);
index = strings(14685+1098+35041,1);
% AJB note: consecutive ranges below share their end/start index (e.g.
% 2937+2937 is in both the KAZR and PBLH THERMO ranges), so the last point of
% each group is labelled as the next instrument (1 of 2937 points per group).
index(1:2937) = "1. GOES18 CTH";
index(2938:2937+2937) = "2. KAZR CTH";
index(2937+2937:2937+2937+2937) = "3. PBLH THERMO";
index(2937+2937+2937:2937+2937+2937+2937) = "4. PBLH MPL";
index(2937+2937+2937+2937:2937+2937+2937+2937+1098)  = "5. BL DEPTH";
index(2937+2937+2937+2937+1098:2937+2937+2937+2937+1098+35041) = "6. LCL";
index(2937+2937+2937+2937+1098+35041:2937+2937+2937+2937+2937+1098+35041) = "7. CEIL CBH";

newmonth = categorical(string(newmonth),monthOrder);  % AJB edits: added ';' to suppress a very long printout (no change to values)
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

% AJB note: THIS BLOCK PRODUCES BAMS FIGURE 5 (seasonal diurnal cycle of
% echo fraction + ceilometer CBH, MWR LWP, SW transmittance + LW flux).
% The slide version differs only cosmetically (black 'Height [m]' label, no
% '23' tick label, 'SW'/'LW' text in the Winter bottom panel), so the final
% slide was made from a slightly later revision of this block.
% AJB note: hours are local (America/Los_Angeles, DST-aware) for the netCDF
% panels, but the SW CSV hours are PST (see the SW hour-offset note above).
% The LW CSV time base can't be checked from the data (assumed same as SW).
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
boxplot(sprcbheight,'Symbol','','Colors',[0.8500, 0.3250, 0.0980],'Widths',0.3)
xlim([0.5 24.5])
ylim([0 3000])
xticklabels([])
yticklabels([])

yyaxis left
pcolor(0.5:24.5,height,sprechomap,'EdgeColor','none')
colormap('sky')
ylim([0 3000])
title('Spring')
ylabel('Height [m]')

nexttile(echospr,9,[1 1])
boxplot(sprlwpbox(:,1:24),'Symbol','')
xticklabels([])
ylabel('Liquid Water Path [g/m^2]')
ylim([0 500])

nexttile(echospr,13,[1,1])
yyaxis right
boxplot(sprSWavg,'Symbol','','Colors',[0, 0.4470, 0.7410])
yticklabels([])
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" "23"])
xtickangle(0)
ylim([0 1])

yyaxis left
plot(mamtime+1,sprLWavg,Color=[0.6350, 0.0780, 0.1840])
ylabel('LW Flux',Color=[0.6350, 0.0780, 0.1840])
ax = gca;
ax.YColor = [0.6350, 0.0780, 0.1840];
ylim([300 410])

xlabel('Local Time [Hour]')
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
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" "23"])
xtickangle(0)
yticklabels([])
ylim([0 1])

yyaxis left
plot(mamtime+1,sumLWavg,Color=[0.6350, 0.0780, 0.1840])
yticklabels([])
ylim([300 410])

xlabel('Local Time [Hour]')


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
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" "23"])
xtickangle(0)
yticklabels([])
ylim([0 1])


yyaxis left
plot(mamtime+1,falLWavg, Color=[0.6350, 0.0780, 0.1840])
ylim([300 410])
yticklabels([])

xlabel('Local Time [Hour]')

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
ylabel('Cloud Base Height')
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

ylabel(c, "Echo Fraction")

nexttile(echospr,12,[1 1])
boxplot(winlwpbox(:,1:24),'Symbol','')
xticklabels([])
yticklabels([])
ylim([0 500])

nexttile(echospr,16,[1,1])
yyaxis right
boxplot(winSWavg,'Symbol','','Colors',[0, 0.4470, 0.7410])
ylabel('SW Trasmittance','Color',[0, 0.4470, 0.7410])
ax = gca;
ax.YColor = '[0, 0.4470, 0.7410]'
xticklabels(["0" "" "" "" "" "5" "" "" "" "" "10" "" "" "" "" "15" "" "" "" "" "20" "" "" "23"])
xtickangle(0)
ylim([0 1])

yyaxis left
plot(mamtime+1,winLWavg,Color=[0.6350, 0.0780, 0.1840])
yticklabels([])
ylim([300 410])

xlabel('Local Time [Hour]')

% AJB edits: export Fig. 5 (Kavin copied figures out of the Live Editor by
% hand). PNG at 300 dpi for slides; PDF as vector for the manuscript.
if exportFigs
    exportgraphics(gcf, fullfile(outDir,'EPCAPEDanLWSWtransmittance_AJB_Fig5.png'),'Resolution',300);
    exportgraphics(gcf, fullfile(outDir,'EPCAPEDanLWSWtransmittance_AJB_Fig5.pdf'),'ContentType','vector');
end
%%

%%
sprsumcbheight = nan(368,24);
for i = 0:24
    heldsize = size(bl_cbh(hour(d2) == i & month(d2) >= 4 & month(d2) <= 9),1);
    sprsumcbheight(1:heldsize,i+1) = bl_cbh(hour(d2) == i & month(d2) >= 4 & month(d2) <= 9);
end



lwp7913 = ncread(fullfile(dataDir,"epcmwrret_turnM1_hourly_averages.nc"),"phys_lwp");  % AJB edits: read from dataDir (was a bare file name)
sprsumlwpbox = nan(330,24);
for i = 0:24
    heldsize = size(lwp7913(hour(d7913) == i & month(d7913) >= 4 & month(d7913) <= 9),1);
    sprsumlwpbox(1:heldsize,i+1) = lwp7913(hour(d7913) == i & month(d7913) >= 4 & month(d7913) <= 9);
end

clf
sprsumechomap = zeros(596,24);
for i = 0:24
    sprsumechomap(:,i+1) = mean(echo_fraction(:,hour(d) == i & month(d) >= 4 & month(d) <= 9),2,'omitmissing');

end

echosprsum = tiledlayout(4,1,'TileSpacing','Compact');

nexttile(echosprsum,[2 1])
yyaxis right
boxplot(sprsumcbheight,'Symbol','','Colors',[0.8500, 0.3250, 0.0980])
ylabel('Cloud Base Height')
xlim([0.5 24.5])
xticklabels([])
yticklabels([])

yyaxis left
pcolor(0.5:24.5,height,sprsumechomap,'EdgeColor','none')
colormap('sky')
c = colorbar;
clim([0 0.5])
ylim([0 3000])
title('Spring Summer Echo Plot')
ylabel('Height [m]')
ylabel(c, "Echo Fraction")

nexttile(echosprsum,[1 1])
boxplot(sprsumlwpbox(:,1:24),'Symbol','')
xticklabels([])
ylabel('Liquid Water Path [g/m^2]')
ylim([0 500])

nexttile(echosprsum,[1,1])
yyaxis right
boxplot(mean(cat(3,sumSWavg, sprSWavg),3),'Symbol','')
ylabel('SW Trasmittance')
xticklabels(["0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23"])

yyaxis left
plot(mamtime+1,mean([sumLWavg sprLWavg],2))
ylabel('LW Flux')
ylim([300 420])

xlabel('Local Time [Hour]')

%%

% AJB note: THIS BLOCK (to the end of the file) MAKES A PRECURSOR OF BAMS
% FIGURE 4, not the slide version. The slide's Figure 4 has 6 rows (adds UCVF/
% RH), runs chronologically Feb 23 - Feb 24 (13 months), and has no title;
% this block bins by calendar month (Jan-Dec, so Feb 2023 and Feb 2024 are
% pooled) and has 5 rows. The code for the final Figure 4 is not in this folder.
% AJB note (zero padding): nan(500,12) is shorter than a month of hourly data
% (up to 744 values), and the loop writes column i+1, so the array grows to
% 744 x 13 and MATLAB pads the new rows/column with ZEROS, not NaN. Months with
% fewer than 744 samples get extra zeros (0 m CBH, 0 g/m^2 LWP) in the boxes
% and in the LWP mean/std printed below. Column 1 is empty (all NaN + padded
% zeros), which is the flat line at the left of the CBH panel.
% Same for yearlylwpbox.
% AJB edits: figure size for the Fig. 4 precursor export (tall 7-row layout).
set(gcf,'Units','pixels','Position',[50 50 900 1000]);
yearlycbheight = nan(500,12);
for i = 1:12
    heldsize = size(bl_cbh(month(d2) == i ),1);
    yearlycbheight(1:heldsize,i+1) = bl_cbh(month(d2) == i );
end
%for i = 1:2
%    heldsize = size(bl_cbh(month(d2) == i),1);
%    yearlycbheight(1:heldsize,i+11) = bl_cbh(month(d2) == i );
%end

yearlylwpbox = nan(500,12);
for i = 1:12
    heldsize = size(lwp7913(month(d7913) == i),1);
    yearlylwpbox(1:heldsize,i+1) = lwp7913(month(d7913) == i );
end

% AJB note (LWP statistics below): column i+1 of yearlylwpbox holds month i,
% so (:,4:6) is Mar-May, (:,7:9) Jun-Aug, (:,10:12) Sep-Nov and (:,1:3) is the
% empty column + Jan + Feb, i.e. one month off from the Apr-Jun/Jul-Sep/
% Oct-Dec/Jan-Mar seasons used elsewhere. With the zero padding (1023 extra
% zeros) the printed means differ from the raw data (R2024b, AJB check):
%   all data 56.4 vs 64.3 g/m^2; 'AMJ' 92.2 vs 89.9; 'JAS' 58.5 vs 60.4;
%   'OND' 42.4 vs 34.0; 'JFM' 33.0 vs 65.1. Check before quoting these in the paper.
mean(yearlylwpbox,'all','omitmissing')
std(yearlylwpbox,[],'all','omitmissing')
mean(yearlylwpbox(:,4:6),'all','omitmissing')
std(yearlylwpbox(:,4:6),[],'all','omitmissing')

mean(yearlylwpbox(:,7:9),'all','omitmissing')
std(yearlylwpbox(:,7:9),[],'all','omitmissing')

mean(yearlylwpbox(:,10:12),'all','omitmissing')
std(yearlylwpbox(:,10:12),[],'all','omitmissing')

mean(yearlylwpbox(:,1:3),'all','omitmissing')
std(yearlylwpbox(:,1:3),[],'all','omitmissing')

%for i = 1:2
%    heldsize = size(lwp7913(month(d7913) == i ),1);
%    yearlylwpbox(1:heldsize,i+11) = lwp7913(month(d7913) == i);
%end

yearlyechomap = zeros(596,12);
clear mean
for i = 1:12
        yearlyechomap(:,i+1) = mean(echo_fraction(:,month(d) == i),2,'omitmissing');
end

%for i = 1:2
%    yearlyechomap(:,i+11) = mean(echo_fraction(:,month(d) == i & year(d) == 2024),2,'omitmissing');
%end

echoyearly = tiledlayout(7,1,'TileSpacing','Tight');

nexttile(echoyearly,1,[2 1])
b = boxchart(newmonth,newdata,'GroupByColor',index,'MarkerStyle','none','Notch','on','WhiskerLineStyle','none');
b(1).BoxFaceColor = [0.129, 0.408, 0.91];
b(2).BoxFaceColor = [0.129, 0.635, 0.91];
b(3).BoxFaceColor = [0.129, 0.91, 0.4];
b(4).BoxFaceColor = [0.098, 0.702, 0.224];
b(5).BoxFaceColor = [0.871, 0.855, 0.118];
b(6).BoxFaceColor = [0.702, 0.098, 0.4];
b(7).BoxFaceColor = [0.871, 0.118, 0.18];

lgd = legend("\color{white}1\color{black}GOES18 CTH","\color{white}2\color{black}KAZRCTH","\color{white}3\color{black}PBLH THERMO","\color{white}4\color{black}PBLH MPL","\color{white}5\color{black}BL DEPTH","\color{white}6\color{black}SONDE LCL","\color{white}7\color{black}CEIL CBH")
lgd.NumColumns = 3;

ylim([0 1750])
xticklabels([])
ylabel("H [m]")

nexttile(echoyearly,3,[2 1])
yyaxis left
% AJB note (alignment): pcolor with 13 x-vertices (1.5:13.5) draws 12 cells
% coloured by columns 1-12, so cell [1.5 2.5] = empty column 1, cell [2.5 3.5]
% = Jan, ..., and Dec (column 13) is not drawn. The CBH boxplot below puts Jan
% at x = 2, so the echo map is shifted one month right of the CBH boxes. The
% CBH boxes also sit on the right y axis with no ylim and no tick labels, so
% their scale is not the 0-3000 m of the left axis.
pcolor(1.5:13.5,height,yearlyechomap,'EdgeColor','none')
colormap('sky')
c = colorbar;
clim([0 0.5])
ylim([0 3000])
ylabel('H [m]')
ylabel(c, "Echo Fraction")

yyaxis right
boxplot(yearlycbheight,'Symbol','','Colors',[0.8500, 0.3250, 0.0980])
ylabel('CBH [m]')
xticklabels([])
yticklabels([])

nexttile(echoyearly,5,[1 1])
boxplot(yearlylwpbox(:,2:13),'Symbol','')
xticklabels([])
ylabel('LWP [g/m^2]')
ylim([0 500])


nexttile(echoyearly,6,[1 1])
SWyearlong = nan(30*24,12);
SW = cat(3,janSW,zeros(30,24), marSW, aprSW, maySW, junSW, julSW, augSW, sepSW, octSW, novSW, decSW);

for i = 1:12
    if i == 2
        SWyearlong(1:28*24,i) = reshape(febSW,28*24,1);
    else
        SWyearlong(:,i) = reshape(SW(:,:,i),30*24,1);
    end
end
yyaxis right
boxplot(SWyearlong,'Colors',[0.8500, 0.3250, 0.0980])
xticklabels([])
ylim([0 1])
ylabel('SW %')

yyaxis left
clear mean
yearlyLW  = mean([janLW febLW marLW aprLW mayLW junLW julLW augLW sepLW octLW novLW decLW],1);
plot(1:12, yearlyLW)
ylim([300 410])
xticklabels([])
ylabel('LW Flux [W/m^2]')

nexttile(echoyearly,7,[1 1])
stdRainfall = readtable(fullfile(dataDir,"STD-Russel EPCAPE Rainfall.xlsx"));  % AJB edits: read from dataDir (was a bare file name)
std2 = stdRainfall.data;
meanRainfall = readtable(fullfile(dataDir,"MEAN-EPCAPE Russel Rainfall.xlsx"));  % AJB edits: read from dataDir (was a bare file name)
mean2 = meanRainfall.Var2;
index2 = meanRainfall.Var3;
time = month(datetime(meanRainfall.Var1,"Format","uuuu/MM"));



q1 = mean2-std2;
q3 = mean2+std2;
min = mean2 - 2.*std2;
max = mean2+ 2.*std2;
bigdata = [min q1 mean2 mean2 q3 max];

datagraph = zeros(size(bigdata,1)*size(bigdata,2),1);
indexgraph = strings(size(bigdata,1)*size(bigdata,2),1);
monthgraph = zeros(size(bigdata,1)*size(bigdata,2),1);
for i = 1:6:(size(bigdata,1)*size(bigdata,2))
    datagraph(i:i+5) = bigdata((i+5)/6,:)';
    indexgraph(i:i+5) = index2((i+5)/6);
    monthgraph(i:i+5) = time((i+5)/6);
end

monthgraph = categorical(string(monthgraph),monthOrder);
boxchart(monthgraph,datagraph,'GroupByColor',indexgraph,'MarkerStyle','none','WhiskerLineStyle','none')
legend("LD/M1","LD/S1","VDIS/M1")
ylabel('Rain [mm/day]')
ylim([0 20])
xlabel('Month')
xticklabels(["Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec"])
clear mean

% AJB edits: export the Fig. 4 precursor so it can be compared with the slide.
if exportFigs
    exportgraphics(gcf, fullfile(outDir,'EPCAPEDanLWSWtransmittance_AJB_Fig4precursor.png'),'Resolution',300);
end
%%
std2
mean2
index2 
time

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
