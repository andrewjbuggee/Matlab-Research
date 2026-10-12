%% EPCAPE: UCVF correlation figure, v1 (reproduces Kavin's figure)
% AJB edits: plain-script copy of EPCAPEUCVFplots_kp_AJB_v1.mlx, made on 2026-10-08
% because the Live Editor hung at the export step (exportgraphics). Same code as the
% .mlx except that the figure is visible here (the .mlx hides it to avoid the hang).
% Run this file (F5, or run('EPCAPEUCVFplots_kp_AJB_v1_script.m')).
% It saves the same figure files as the .mlx: EPCAPEUCVFplots_kp_AJB_v1.png and .pdf.
%
% AJB edits (Andrew Buggee, 2026-10-08): version 1 is Kavin's EPCAPEUCVFplots_kp.mlx
% (saved in MATLAB R2025b), updated so that it reads its inputs from the shared data
% folder and reproduces Kavin's original figure, plus an x-axis label. Kavin's
% statistics and plotting commands are unchanged, so the known issues listed below are
% still in this version. Version 2 (EPCAPEUCVFplots_kp_AJB_v2.mlx) fixes them.
%
% How to find the changes: every line of code that AJB added, removed or changed is
% marked with a comment that starts with "AJB edits:" and says what changed. Comments
% that start with "AJB note:" only explain Kavin's code and change nothing. Comments
% without either tag are Kavin's. The text sections (this header, Panels, Known issues)
% are AJB's.
%
% Summary of the AJB edits:
%
% # Every input file is read from dataDir. Kavin's absolute paths and the files he read
%   from the current folder are replaced by fullfile(dataDir, ...) (User settings and
%   section 1).
% # The clf calls in sections 1 and 2 are removed, and the clf in section 3 is replaced
%   by a new figure with a size that reproduces the panel proportions of Kavin's figure
%   and with MATLAB's factory graphics defaults, so that personal defaults set in a
%   startup.m (e.g. a dark theme) do not change its appearance.
% # The legend positions of Kavin's figure are set explicitly (section 3). R2024b puts
%   all three legends in the top-right corner.
% # In MATLAB releases before R2025a, scatter points drawn with the '.' marker ignore
%   alpha(), so the faint points in panels A-C plot opaque. On those releases the '.'
%   markers are swapped for filled squares of the same rendered size, which honour
%   alpha() (end of section 3).
% # An x-axis label, "UCVF", is added below panel C. Kavin's figure has none. The label
%   fits in the layout's existing bottom margin, so the panels keep the size they have in
%   Kavin's figure.
% # The figure is saved as PNG and PDF (section 4).
%
%% Panels
% x axis of every panel: UCVF. Each variable is multiplied and offset so that several
% relations share one panel (hence "Arbitrary Units"); each fit line is a least-squares
% line (fitlm) and gets the same multiplier and offset as its points. The legends give
% Pearson's r between UCVF and the variable (corrcoef, pairwise-complete rows).
%
% * A: AOS and AOS F_thickly (90 nm) from the Petters and Liggio cloud aqueous reaction
%   file; SO4, MO-OOA, SR-OOA and f44 from KAVIN_TABLE.
% * B: HOA and rBC mass (LANL file) and daily mean ultrafine salt.
% * C: OH burst (UCLA overview file, against UCVF24) and July 2023 dD.
%
%% Known issues in Kavin's original figure (reproduced here, not fixed)
% # Panel B legend: the names rBC and HOA are swapped. The blue line and points are
%   HOA (r = -0.23); the orange ones are rBC (r = -0.27).
% # Panel B legend: "NaCl2, r = 0.43" prints string3, which is the MO-OOA r of panel A.
%   The UCVF vs. ultrafine salt correlation (string14) is r = -0.68 (n = 49 days).
% # Panel B: the HOA and rBC points are plotted multiplied by 2 but their fit lines are
%   not, so each line sits at half the height of its points. With the factor 2, ylim
%   [0 0.6] also cuts off 32 % of the HOA points and 13 % of the rBC points.
% # Panel A: the MO-OOA and SR-OOA fit lines and r values use April-June 2023 only (578
%   of the 2230 UCVF-factor pairs), while the points show the whole record (February
%   2023 - February 2024). With all pairs, r = 0.11 (MO-OOA) and r = 0.52 (SR-OOA).
%
%% User settings
% AJB edits: this whole section is new.

% AJB edits: Kavin's "clear all" moved here from the top of section 1, so that it runs
% before the settings below are defined.
clear all

% AJB edits: folder that holds all of Kavin's input files for this figure
dataDir = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper';

% AJB edits: input files (all in dataDir)
ultrafineFile = fullfile(dataDir, 'ucvf vs ultrafine salt.xlsx');      % daily mean UCVF and ultrafine salt (panel B)
lanlFile      = fullfile(dataDir, 'LANL_UCVF_Correlations.xlsx');      % UCVF, HOA and rBC mass in ug m^-3 (panel B)
ohFile        = fullfile(dataDir, 'EPCAPE_overview_UCLAgraph.csv');    % OH burst (uM), RH and UCVF24 per sample (panel C)
kavinFile     = fullfile(dataDir, 'kavin_UCVF.mat');                   % KAVIN_TABLE: SR-OOA, MO-OOA, f44, NO3, SO4, UCVF (panel A)
dDFile        = fullfile(dataDir, 'dD_UCVFreal.xlsx');                  % time, dD and UCVF (panel C)
aqueousFile   = fullfile(dataDir, 'Cloud aqueous reaction data (Petters and Liggio) - April.csv');  % UCVF, AOS, AOS F_thickly 90 nm (panel A)

% AJB edits: figure output
saveFigure = true;          % save <figName>.png (300 dpi) and <figName>.pdf in outDir
outDir     = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS';
figName    = 'EPCAPEUCVFplots_kp_AJB_v1';

% AJB edits: figure size in points, [width height]. Kavin's code does not set a size.
% This one reproduces the axes of his figure as measured on the image: 471 pt wide,
% panel A 272 pt tall, panels B and C 132 pt tall (with the default 10 pt tick labels).
figSize_pt = [811 675];

% AJB edits: true = draw the '.' scatter markers as filled squares of the same rendered
% size, so that alpha() works (see the header). Needed before R2025a.
replacePointMarkers = isMATLABReleaseOlderThan("R2025a");
%% 1. Load the data
% AJB note: lines without a semicolon print the table, as in Kavin's version.

% AJB edits: Kavin's "clear all" moved to User settings, and his "clf" removed (the
% figure is created in section 3).

% AJB edits: path "/Users/2fast2kavin/Desktop/Atmospheric Aerosol Group/EPCAPE/UCVF plot/ucvf vs ultrafine salt.xlsx" -> ultrafineFile
% AJB note: daily mean UCVF and ultrafine salt (DateAndTime, avgUCVF, avgUfineSalt)
ultrafinetbl = readtable(ultrafineFile)
% AJB edits: path "/Users/2fast2kavin/Desktop/Atmospheric Aerosol Group/EPCAPE/UCVF plot/LANL_UCVF_Correlations.xlsx" -> lanlFile
% AJB note: LANL HOA and rBC mass concentrations in ug m^-3 (DateTime_UTC, UCVF, HOAMass_ugM_3_, rBCMass_ugM_3_)
LANLdata_HOA_rBC = readtable(lanlFile)
% AJB edits: "EPCAPE_overview_UCLAgraph.csv" (current folder) -> ohFile
% AJB note: OH burst concentration (uM) and RH per sample, with UCVF24 and UCVF6 (UCVF24 is used)
OHvsUCVF = readtable(ohFile)
% AJB edits: "load kavin_UCVF.mat" (current folder) -> load(kavinFile)
% AJB note: KAVIN_TABLE: Time, SR-OOA, MO-OOA, f44, NO3, SO4, UCVF (15 Jan 2023 - 15 Feb 2024)
load(kavinFile)
KAVIN_TABLE
% AJB note: running month number, 1 = January 2023, 13 = January 2024; dmoso 4-6 = April-June 2023
dmoso = month(KAVIN_TABLE.Time) + (year(KAVIN_TABLE.Time) -2023)*12
% AJB edits: "dD_UCVFreal.xlsx" (current folder) -> dDFile
% AJB note: dD and UCVF time series, 14 Jul 2023 - 6 Feb 2024. In R2024b readtable returns
% the time column as text such as '2023-07-14 22:20:00+00:00'; the first 19 characters
% (without the UTC offset) are parsed to datetime.
NewMexico = readtable(dDFile)
stringythingy = extractBetween(string(NewMexico.time),1,19);
ddD = datetime(stringythingy,"Format","dd-MMM-uuuu HH:mm:ss");
% AJB note: July values only (the record starts on 14 July 2023, so this is July 2023).
% UCVF is read as text and converted to numbers.
julydD = NewMexico.dD(month(ddD) == 7)
julyUCVF = str2double(string(NewMexico.UCVF(month(ddD) == 7)))

% AJB edits: "Cloud aqueous reaction data (Petters and Liggio) - April.csv" (current folder) -> aqueousFile
% AJB note: Petters and Liggio cloud aqueous reaction data, April (UCVF, AOS, AOSFthickly_90nm, ...)
elavarasi = readtable(aqueousFile)
%% 2. Linear fits and correlation coefficients
% AJB note: each block fits y = b + m*UCVF (fitlm), evaluates the fit on 100 points
% between the smallest and largest UCVF (xFitN, yFitN), and builds the legend text
% stringN = "r = " followed by Pearson's r rounded to 2 decimals.

% AJB edits: Kavin's "clf" removed from the top of this section (nothing is plotted here).

% AJB note: panel B, HOA (fit 1) and rBC (fit 2), LANL data
model = fitlm(LANLdata_HOA_rBC.UCVF,LANLdata_HOA_rBC.HOAMass_ugM_3_);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit1 = linspace(min(LANLdata_HOA_rBC.UCVF), max(LANLdata_HOA_rBC.UCVF), 100);
yFit1 = feval(model, xFit1);

[r,p]=corrcoef(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.HOAMass_ugM_3_,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string1=append(str1,str2);

model = fitlm(LANLdata_HOA_rBC.UCVF,LANLdata_HOA_rBC.rBCMass_ugM_3_);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit2 = linspace(min(LANLdata_HOA_rBC.UCVF), max(LANLdata_HOA_rBC.UCVF), 100);
yFit2 = feval(model, xFit2);

[r,p]=corrcoef(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.rBCMass_ugM_3_,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string2=append(str1,str2);

%%%%%%TILE 2
% AJB note: panel A, MO-OOA (fit 3), April-June 2023 only (dmoso 4-6)

model = fitlm(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4),KAVIN_TABLE.("MO-OOA")(dmoso <= 6 & dmoso >= 4));
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit3 = linspace(min(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4)), max(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4)), 100);
yFit3 = feval(model, xFit3);

[r,p]=corrcoef(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4), KAVIN_TABLE.("MO-OOA")(dmoso <= 6 & dmoso >= 4),'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string3=append(str1,str2);

% AJB note: NO3 (fit 4), all rows; computed but not plotted
model = fitlm(KAVIN_TABLE.UCVF,KAVIN_TABLE.NO3);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit4 = linspace(min(KAVIN_TABLE.UCVF), max(KAVIN_TABLE.UCVF), 100);
yFit4 = feval(model, xFit4);

[r,p]=corrcoef(KAVIN_TABLE.UCVF, KAVIN_TABLE.NO3,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string4=append(str1,str2);

% AJB note: panel A, SO4 (fit 5), all rows
model = fitlm(KAVIN_TABLE.UCVF,KAVIN_TABLE.SO4);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit5 = linspace(min(KAVIN_TABLE.UCVF), max(KAVIN_TABLE.UCVF), 100);
yFit5 = feval(model, xFit5);

[r,p]=corrcoef(KAVIN_TABLE.UCVF, KAVIN_TABLE.SO4,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string5=append(str1,str2);

% AJB note: panel A, SR-OOA (fit 6), April-June 2023 only (dmoso 4-6)
model = fitlm(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4),KAVIN_TABLE.("SR-OOA")(dmoso <= 6 & dmoso >= 4));
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit6 = linspace(min(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4)), max(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4)), 100);
yFit6 = feval(model, xFit6);

[r,p]=corrcoef(KAVIN_TABLE.UCVF(dmoso <= 6 & dmoso >= 4), KAVIN_TABLE.("SR-OOA")(dmoso <= 6 & dmoso >= 4),'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string6 =append(str1,str2);

% AJB note: panel A, AOS (fit 7), Petters and Liggio data
model = fitlm(elavarasi.UCVF,elavarasi.AOS);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit7 = linspace(min(elavarasi.UCVF), max(elavarasi.UCVF), 100);
yFit7 = feval(model, xFit7);

[r,p]=corrcoef(elavarasi.UCVF, elavarasi.AOS,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string7=append(str1,str2);

% AJB note: panel A, f44 (fit 8), all rows

model = fitlm(KAVIN_TABLE.UCVF,KAVIN_TABLE.f44);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit8 = linspace(min(KAVIN_TABLE.UCVF), max(KAVIN_TABLE.UCVF), 100);
yFit8 = feval(model, xFit8);

[r,p]=corrcoef(KAVIN_TABLE.UCVF, KAVIN_TABLE.f44,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string8=append(str1,str2);

%%%%TILE 3
% AJB note: July dD (fit 9); computed but not plotted (fit 13 below is the same fit)

model = fitlm(julyUCVF,julydD);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit9 = linspace(min(julyUCVF), max(julyUCVF), 100);
yFit9 = feval(model, xFit9);

[r,p]=corrcoef(julyUCVF, julydD,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string9=append(str1,str2);

% AJB note: panel C, OH burst vs. UCVF24 (fit 10)

model = fitlm(OHvsUCVF.UCVF24,OHvsUCVF.OH_burst_25000_vv_dilution_ratio_uM);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit10 = linspace(min(OHvsUCVF.UCVF24), max(OHvsUCVF.UCVF24), 100);
yFit10 = feval(model, xFit10);

[r,p]=corrcoef(OHvsUCVF.UCVF24, OHvsUCVF.OH_burst_25000_vv_dilution_ratio_uM,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string10=append(str1,str2);

%xL=xlim; xdiff = diff(xlim); xpos = xL(1) + 0.7*xdiff;
%yL=ylim; ydiff = diff(ylim); ypos= yL(1) + 0.9*ydiff;
%text(xpos,ypos,str,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',10)


%title("Isotopes and Oxidants")

% AJB note: panel A, AOS F_thickly 90 nm (fit 11), Petters and Liggio data
model = fitlm(elavarasi.UCVF,elavarasi.AOSFthickly_90nm);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit11 = linspace(min(elavarasi.UCVF), max(elavarasi.UCVF), 100);
yFit11 = feval(model, xFit11);

[r,p]=corrcoef(elavarasi.UCVF,elavarasi.AOSFthickly_90nm,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string11=append(str1,str2);

% AJB note: RH vs. UCVF24 (fit 12); computed but not plotted
model = fitlm(OHvsUCVF.UCVF24,OHvsUCVF.RH);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit12 = linspace(min(OHvsUCVF.UCVF24), max(OHvsUCVF.UCVF24), 100);
yFit12 = feval(model, xFit12);

[r,p]=corrcoef(OHvsUCVF.UCVF24,OHvsUCVF.RH,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string12=append(str1,str2);


% AJB note: panel C, July dD (fit 13)
model = fitlm(julyUCVF,julydD);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit13 = linspace(min(julyUCVF), max(julyUCVF), 100);
yFit13 = feval(model, xFit13);

[r,p]=corrcoef(julyUCVF,julydD,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string13=append(str1,str2);

% AJB note: panel B, daily mean ultrafine salt (fit 14). The panel B legend prints
% string3 (MO-OOA) instead of string14 (known issue 2 in the header).
model = fitlm(ultrafinetbl.avgUCVF,ultrafinetbl.avgUfineSalt);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit14 = linspace(min(ultrafinetbl.avgUCVF), max(ultrafinetbl.avgUCVF), 100);
yFit14 = feval(model, xFit14);

[r,p]=corrcoef(ultrafinetbl.avgUCVF,ultrafinetbl.avgUfineSalt,'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string14=append(str1,str2);


% AJB note: binned means and standard deviations of SO4 (0.05-wide UCVF bins) and of
% MO-OOA, SR-OOA, f44 and NO3 (0.1-wide bins). Not used in the figure. Each array has
% 489 rows (the largest count in a 0.1-wide bin) and is padded with zeros, not NaN, so
% the 'omitmissing' means and standard deviations include the padding zeros.

discso4index = discretize(KAVIN_TABLE.UCVF,0:0.05:1);
discso4 = zeros(489,20);
for i = 1:20
    discsize = size(KAVIN_TABLE.SO4(discso4index == i));
    discso4(1:discsize(1),i) = KAVIN_TABLE.SO4(discso4index == i);
end
meanso4 = mean(discso4,1,'omitmissing');
stdso4 = std(discso4,[],1,'omitmissing');

discso4index = discretize(KAVIN_TABLE.UCVF,0:0.1:1);
discMO = zeros(489,10);
for i = 1:10
    discsize = size(KAVIN_TABLE.("MO-OOA")(discso4index == i));
    discMO(1:discsize(1),i) = KAVIN_TABLE.('MO-OOA')(discso4index == i);
end
meanMO = mean(discMO,1,'omitmissing');
stdMO = std(discMO,[],1,'omitmissing');

discso4index = discretize(KAVIN_TABLE.UCVF,0:0.1:1);
discSR = zeros(489,10);
for i = 1:10
    discsize = size(KAVIN_TABLE.("SR-OOA")(discso4index == i));
    discSR(1:discsize(1),i) = KAVIN_TABLE.('SR-OOA')(discso4index == i);
end
meanSR = mean(discSR,1,'omitmissing');
stdSR = std(discSR,[],1,'omitmissing');

discso4index = discretize(KAVIN_TABLE.UCVF,0:0.1:1);
discf44 = zeros(489,10);
for i = 1:10
    discsize = size(KAVIN_TABLE.f44(discso4index == i));
    discf44(1:discsize(1),i) = KAVIN_TABLE.f44(discso4index == i);
end
meanf44 = mean(discf44,1,'omitmissing');
stdf44 = std(discf44,[],1,'omitmissing');

discso4index = discretize(KAVIN_TABLE.UCVF,0:0.1:1);
discNO3 = zeros(489,10);
for i = 1:10
    discsize = size(KAVIN_TABLE.NO3(discso4index == i));
    discNO3(1:discsize(1),i) = KAVIN_TABLE.NO3(discso4index == i);
end
meanNO3 = mean(discNO3,1,'omitmissing');
stdNO3 = std(discNO3,[],1,'omitmissing');
%% 3. Figure

% AJB edits: Kavin's "clf" replaced by a new figure of size figSize_pt (points, so the
% size is the same on any screen) with a white background.
fig = figure('Units', 'points', 'Position', [40 40 figSize_pt], 'Color', 'w');

% AJB edits: reset, on this figure only, every graphics default set at the root level
% (for example by startup.m) to MATLAB's factory value. Root defaults are named
% 'default<Type><Prop>' and the matching factory values 'factory<Type><Prop>'. Figure
% properties cannot be figure-level defaults, so they are skipped (the figure colour is
% set above). On a MATLAB without custom root defaults this loop changes nothing.
rootDefaults = fieldnames(get(groot, 'default'));
for k = 1:numel(rootDefaults)
    if startsWith(rootDefaults{k}, 'defaultFigure', 'IgnoreCase', true)
        continue
    end
    factoryName = regexprep(rootDefaults{k}, '^default', 'factory', 'ignorecase');
    set(fig, rootDefaults{k}, get(groot, factoryName))
end

% AJB note: 8 x 4 grid of tiles. The panels span columns 1-3 only (column 4 stays empty
% and exportgraphics crops it away): A = rows 1-4, B = rows 5-6, C = rows 7-8.
tiledlayout(8,4,'TileSpacing','compact')

% AJB note: panel A. Fit lines first, so that the legend entries (in plotting order)
% belong to them; then the points, with the same multiplier and offset as their line and
% in the line's default colour (HandleVisibility off keeps them out of the legend).
nexttile(1,[4,3])
plot(xFit7,yFit7)
hold on
plot(xFit5,yFit5.*2 + 3)
plot(xFit3,yFit3.*5 + 2.5)
plot(xFit6,yFit6.*7 + 1.8)
plot(xFit8,yFit8.*25 - 3)
plot(xFit11,yFit11.*15-4)

scatter(elavarasi.UCVF,elavarasi.AOS,[20],[0, 0.4470, 0.7410],HandleVisibility='off',Marker='.')
s = scatter(KAVIN_TABLE.UCVF,KAVIN_TABLE.SO4.*2 + 3,[20],[0.8500, 0.3250, 0.0980],HandleVisibility='off',Marker='.')
alpha(s,0.1)
s2 = scatter(KAVIN_TABLE.UCVF,KAVIN_TABLE.('MO-OOA').*5 + 2.5,[20],[0.9290, 0.6940, 0.1250],HandleVisibility='off',Marker='.')
alpha(s2,0.3)
s3 = scatter(KAVIN_TABLE.UCVF,KAVIN_TABLE.('SR-OOA').*7 + 1.8,[20],[0.4940, 0.1840, 0.5560],HandleVisibility="off",Marker = '.')
alpha(s3,0.1)
s4 = scatter(KAVIN_TABLE.UCVF,KAVIN_TABLE.f44.*25 -3,[20],[0.4660, 0.6740, 0.1880],HandleVisibility="off",Marker = '.')
alpha(s4,0.1)
s6 = scatter(elavarasi.UCVF,elavarasi.AOSFthickly_90nm.*15-4,[20],[0.3010, 0.7450, 0.9330],HandleVisibility="off",Marker = '.')
alpha(s6,0.9)
% AJB note: no3index is not used
no3index = KAVIN_TABLE.NO3.*95 -29 < 6
hold off
ylim([0 7])
ylabel("Arbitrary Units")
yticklabels("")
lgnd = legend("AOS, "+ string7,"SO4, "+string5,"MO-OOA, "+string3,"SR-OOA, " + string6,"f44, "+string8,"FThickly, "+string11)
lgnd.FontSize = 5
% AJB edits: legend in the bottom-right corner, where it is in Kavin's figure
lgnd.Location = 'southeast';
text(0.03, 0.97*7,'A')
xticklabels([""])
xticks([0:0.1:1])

% AJB note: panel B. Fit lines HOA (yFit1), rBC (yFit2), ultrafine salt (yFit14). The
% HOA and rBC points are multiplied by 2 but their lines are not (known issue 3).
nexttile(17,[2,3])
plot(xFit1,yFit1)
hold on
plot(xFit2,yFit2)
plot(xFit14,yFit14)
sone = scatter(LANLdata_HOA_rBC.UCVF,LANLdata_HOA_rBC.HOAMass_ugM_3_.*2,[20],[0, 0.4470, 0.7410],HandleVisibility='off',Marker='.')
alpha(sone,0.4)
stwo = scatter(LANLdata_HOA_rBC.UCVF,LANLdata_HOA_rBC.rBCMass_ugM_3_.*2,[20],[0.8500, 0.3250, 0.0980],HandleVisibility='off',Marker='.')
alpha(stwo,0.5)
sthr = scatter(ultrafinetbl.avgUCVF,ultrafinetbl.avgUfineSalt,[20],[0.9290, 0.6940, 0.1250]	,HandleVisibility='off',Marker='.')
alpha(sthr,0.5)
hold off
ylabel("Arbitrary Units")
yticklabels("")
% AJB note: known issues 1 and 2. Entry 1 is the HOA line (string1 = HOA r) but is
% named rBC, entry 2 is the rBC line named HOA, and entry 3 prints string3 (the MO-OOA
% r) instead of string14.
legend("rBC, "+ string1,"HOA, "+ string2,"NaCl2, "+ string3, FontSize=5)
% AJB edits: legend in the top-right corner, where it is in Kavin's figure
axB = gca;
axB.Legend.Location = 'northeast';
text(0.03, 0.95*0.6,'B')
xticklabels("")
ylim([0 0.6])
xticks([0:0.1:1])

% AJB note: panel C. OH burst (x20) and July dD (+85); lines and points get the same
% transformation.
nexttile(25,[2,3])
plot(xFit10,yFit10.*20)
hold on
plot(xFit13,yFit13+ 85)
suno = scatter(julyUCVF,julydD+85,[20],[0.8500, 0.3250, 0.0980],HandleVisibility='off',Marker='.')
alpha(suno,0.1)
sdos = scatter(OHvsUCVF.UCVF24,OHvsUCVF.OH_burst_25000_vv_dilution_ratio_uM.*20,[20],[0, 0.4470, 0.7410],HandleVisibility='off',Marker='.')
alpha(sdos,1)
hold off
ylabel("Arbitrary Units")
yticklabels("")
legend("OH, "+ string10,"7/2023 dD, "+string13, FontSize=5)
text(0.03, 0.95*40-20,'C')
ylim([-20 20])
xticks([0:0.1:1])
xticklabels([0:0.1:1])
% AJB edits: x-axis label (not in Kavin's figure)
xlabel("UCVF")

% AJB edits: panel C legend at the position it has in Kavin's figure, measured on his
% image: lower-left corner 4.8 % of the axes width right of the left axis and 19.7 % of
% the axes height above the bottom axis. The axes and its legend are both children of
% the tiled layout, so their Position values share one coordinate frame. drawnow first,
% so the layout has its final axes positions.
axC = gca;
drawnow
axC.Legend.Position(1:2) = axC.Position(1:2) + [0.048 0.197].*axC.Position(3:4);

% AJB edits: before R2025a, '.' scatter markers draw opaque whatever their alpha. Swap
% them for filled squares, which honour alpha. A square with SizeData/4 renders with the
% same area as the '.' marker (both about 12 px at 150 dpi in R2024b). alpha() set the
% marker's face and edge transparency to the same value; the square gets it as its face
% transparency. findall (not findobj) because the scatters have HandleVisibility off.
if replacePointMarkers
    pointScatters = findall(fig, 'Type', 'scatter', 'Marker', '.');
    for k = 1:numel(pointScatters)
        sc = pointScatters(k);
        set(sc, 'Marker', 's', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', sc.MarkerEdgeAlpha, 'SizeData', sc.SizeData/4)
    end
end
%% 4. Save the figure
% AJB edits: this whole section is new. exportgraphics crops the figure to its
% contents, which removes the empty fourth tile column.

if saveFigure
    exportgraphics(fig, fullfile(outDir, [figName '.png']), 'Resolution', 300)
    exportgraphics(fig, fullfile(outDir, [figName '.pdf']), 'ContentType', 'vector')
    fprintf('Saved %s.png and %s.pdf in %s\n', figName, figName, outDir)
end
