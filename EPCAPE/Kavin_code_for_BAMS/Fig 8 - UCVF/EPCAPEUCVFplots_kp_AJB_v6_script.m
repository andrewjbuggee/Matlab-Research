%% EPCAPE: UCVF correlation figure, v6 (panel A split into six rows)
% AJB edits: plain-script copy of EPCAPEUCVFplots_kp_AJB_v6.mlx, made on 2026-10-08
% because the Live Editor hung at the export step (exportgraphics). Same code as the
% .mlx except that figureVisible = 'on' here, so the figure is shown in a window.
% Run this file (Run button, or run('EPCAPEUCVFplots_kp_AJB_v6_script.m')).
% It saves the same figure files as the .mlx: EPCAPEUCVFplots_kp_AJB_v6.png and .pdf.
%
% AJB edits (Andrew Buggee, 2026-10-08): version 6 is version 2
% (EPCAPEUCVFplots_kp_AJB_v2.mlx) with panel A split into six thin panels, one per
% species, stacked close together, so that each species and its fit can be seen on its
% own. Panels B and C are as in version 2. See "Changes in version 6" below.
%
% How to find the changes: lines changed for version 6 are marked with a comment that
% starts with "AJB edits (v6):". "AJB edits (v2):" marks the version 2 bug fixes, and
% "AJB edits:" the version 1 changes to Kavin's code. Comments that start with
% "AJB note:" only explain the code and change nothing. Comments without a tag are
% Kavin's. The text sections (this header, How to run, Panels, Fixes, Changes, Open
% questions) are AJB's.
%
% Summary of the version 1 changes to Kavin's code (kept in version 2):
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
%% How to run this file
% # Folders: in "User settings" below, switch to your two folder lines (dataDir and
%   outDir) as described there. All six input files must be in dataDir. Kavin's
%   original code read two of them from ".../EPCAPE/UCVF plot/" and the other four from
%   MATLAB's current folder, so copy those four into that folder if needed.
% # Press Run. The figure is saved in outDir as EPCAPEUCVFplots_kp_AJB_v6.png and EPCAPEUCVFplots_kp_AJB_v6.pdf.
% # In this .mlx file the figure is hidden while it is made (figureVisible = 'off'); open
%   the saved PNG to see it. EPCAPEUCVFplots_kp_AJB_v6_script.m is the same code as a plain script and shows
%   the figure in a window.
% Needs the Statistics and Machine Learning Toolbox (fitlm), as Kavin's original code
% did. Tested in MATLAB R2024b; written to run in R2025a and later as well.
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
%% Fixes in version 2
% # Panel B legend: the HOA and rBC names are swapped back. Entry 1 (blue) is HOA,
%   r = -0.23; entry 2 (orange) is rBC, r = -0.27. In Kavin's figure the names were
%   swapped while the r values and colours belonged to the right lines.
% # Panel B legend: the ultrafine salt entry prints its own correlation, string14
%   (r = -0.68, n = 49 days), instead of string3, the MO-OOA r of panel A (0.43).
% # Panel B: the HOA and rBC points are no longer multiplied by 2, so they match their
%   fit lines. ylim stays [0 0.6]; it now cuts off 12.5 % of the HOA points and 3.7 % of
%   the rBC points (32 % and 13 % in version 1). Multiplying the lines by 2 instead and
%   doubling ylim would give the same picture, because the y axis has no tick values.
% # Panel A: the MO-OOA and SR-OOA points, fit lines and r values now cover the same
%   period, chosen with ooaAprJunOnly (User settings). true (default): Kavin's
%   April-June 2023 subset; the points are restricted to it (578 UCVF-factor pairs), the
%   legend says "(Apr-Jun 2023)", and r is unchanged (MO-OOA 0.43, SR-OOA 0.4). false:
%   points, fits and r use the whole record, February 2023 - February 2024 (2230 pairs);
%   r = 0.11 (MO-OOA) and 0.52 (SR-OOA). With the default, the dense band of SR-OOA
%   points at low UCVF in Kavin's figure disappears: 533 of its 541 points are from
%   outside April-June 2023, mostly cool-season samples (November-March of both
%   winters in the record).
%
%% Changes in version 6
% # Panel A is split into six panels stacked on top of each other, one per species (AOS,
%   SO4, MO-OOA, SR-OOA, f44, FThickly). They share the x axis and sit close together.
% # Each of the six panels has its own y axis, so every species is shown at its real
%   values with numbers on the axis. Kavin's multipliers and offsets are no longer
%   needed in panel A. The species name is the y-axis label; r is written inside the
%   panel.
% # The figure is taller (figSize_pt) to make room for the six panels. Panels B and C
%   are unchanged.
% Everything else, including ooaAprJunOnly = true (April-June 2023 for MO-OOA and
% SR-OOA), is as in version 2.
%
%% Open questions for Kavin and Lynn
% The v2 changes and the open questions are listed in plain language in
% EPCAPEUCVFplots_v2_changes_and_open_questions.txt (same folder).
%
%% User settings
% AJB edits: this whole section is new.

% AJB edits: Kavin's "clear all" moved here from the top of section 1, so that it runs
% before the settings below are defined.
clear all

% AJB edits: folders. Use ONE of the two blocks below and comment out the other one.
% To switch, remove the "% " at the start of Kavin's two lines and add "% " at the
% start of Andrew's two lines.
%   dataDir: folder with all six input files listed below
%   outDir:  folder where the figure (PNG and PDF) is saved
% --- Andrew ---
dataDir = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper';
outDir  = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS/Fig 8 - UCVF';   % AJB edits (2026-10-14): Figure 8 files moved to the subfolder "Fig 8 - UCVF"
% --- Kavin (the folder used in Kavin's original EPCAPEUCVFplots_kp.mlx) ---
% dataDir = '/Users/2fast2kavin/Desktop/Atmospheric Aerosol Group/EPCAPE/UCVF plot';
% outDir  = dataDir;

% AJB edits: input files (all in dataDir)
ultrafineFile = fullfile(dataDir, 'ucvf vs ultrafine salt.xlsx');      % daily mean UCVF and ultrafine salt (panel B)
lanlFile      = fullfile(dataDir, 'LANL_UCVF_Correlations.xlsx');      % UCVF, HOA and rBC mass in ug m^-3 (panel B)
ohFile        = fullfile(dataDir, 'EPCAPE_overview_UCLAgraph.csv');    % OH burst (uM), RH and UCVF24 per sample (panel C)
kavinFile     = fullfile(dataDir, 'kavin_UCVF.mat');                   % KAVIN_TABLE: SR-OOA, MO-OOA, f44, NO3, SO4, UCVF (panel A)
dDFile        = fullfile(dataDir, 'dD_UCVFreal.xlsx');                  % time, dD and UCVF (panel C)
aqueousFile   = fullfile(dataDir, 'Cloud aqueous reaction data (Petters and Liggio) - April.csv');  % UCVF, AOS, AOS F_thickly 90 nm (panel A)

% AJB edits: figure output
saveFigure = true;          % save <figName>.png (300 dpi) and <figName>.pdf in outDir (set above)
figName    = 'EPCAPEUCVFplots_kp_AJB_v6';  % AJB edits (v6): v2 -> v6

% AJB edits: figure size in points, [width height]. Kavin's code does not set a size.
% AJB edits (v6): taller than v2 (675 pt) to make room for the six panels of panel A.
% AJB edits: 'off' = the figure is made hidden; 'on' = it is shown in a window. The
% saved PNG and PDF are the same either way. Keep 'off' in the .mlx file: with a visible
% figure, the Live Editor (MATLAB R2024b) hangs when it saves the figure. The plain-script
% copy (..._script.m) uses 'on'.
figureVisible = 'on';   % AJB edits: 'on' in this plain-script copy

figSize_pt = [811 860];

% AJB edits: true = draw the '.' scatter markers as filled squares of the same rendered
% size, so that alpha() works (see the header). Needed before R2025a.
replacePointMarkers = isMATLABReleaseOlderThan("R2025a");

% AJB edits (v2): period used for MO-OOA and SR-OOA in panel A (points, fit lines and r).
% true = April-June 2023 only, as in Kavin's fits; false = the whole record
% (February 2023 - February 2024). See fix 4 in the header.
ooaAprJunOnly = true;
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
% AJB edits (v2): rows used for MO-OOA and SR-OOA (points, fits and r) and the matching
% legend suffix. Kavin selected April-June 2023 (dmoso 4-6) inside the fits only.
if ooaAprJunOnly
    ooaRows = dmoso <= 6 & dmoso >= 4;
    ooaLabel = " (Apr-Jun 2023)";
else
    ooaRows = true(height(KAVIN_TABLE), 1);
    ooaLabel = "";
end
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
% AJB note: panel A, MO-OOA (fit 3)
% AJB edits (v2): rows ooaRows instead of Kavin's (dmoso <= 6 & dmoso >= 4) in the four
% lines below, so that the period follows ooaAprJunOnly (same rows when it is true).

model = fitlm(KAVIN_TABLE.UCVF(ooaRows),KAVIN_TABLE.("MO-OOA")(ooaRows));
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit3 = linspace(min(KAVIN_TABLE.UCVF(ooaRows)), max(KAVIN_TABLE.UCVF(ooaRows)), 100);
yFit3 = feval(model, xFit3);

[r,p]=corrcoef(KAVIN_TABLE.UCVF(ooaRows), KAVIN_TABLE.("MO-OOA")(ooaRows),'rows','pairwise');
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

% AJB note: panel A, SR-OOA (fit 6)
% AJB edits (v2): rows ooaRows instead of Kavin's (dmoso <= 6 & dmoso >= 4) in the four
% lines below, so that the period follows ooaAprJunOnly (same rows when it is true).
model = fitlm(KAVIN_TABLE.UCVF(ooaRows),KAVIN_TABLE.("SR-OOA")(ooaRows));
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit6 = linspace(min(KAVIN_TABLE.UCVF(ooaRows)), max(KAVIN_TABLE.UCVF(ooaRows)), 100);
yFit6 = feval(model, xFit6);

[r,p]=corrcoef(KAVIN_TABLE.UCVF(ooaRows), KAVIN_TABLE.("SR-OOA")(ooaRows),'rows','pairwise');
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

% AJB note: panel B, daily mean ultrafine salt (fit 14). string14 is printed in the
% panel B legend (fix 2 in the header).
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
% size is the same on any screen) with a white background, hidden or shown according to
% figureVisible (User settings).
fig = figure('Units', 'points', 'Position', [40 40 figSize_pt], 'Color', 'w', 'Visible', figureVisible);

% AJB edits: newer MATLAB releases (R2025a and later) draw figures in the desktop's
% theme, which can be dark. Use the light theme, so that the figure always has a white
% background and dark text. (If a release has no Theme property, this is skipped.)
if isprop(fig, 'Theme')
    fig.Theme = "light";
end

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

% AJB edits (v6): layout. As before, the panels span columns 1-3 only (column 4 stays
% empty and exportgraphics crops it away), but the grid has 10 rows (v2: 8):
% A = rows 1-6, B = rows 7-8, C = rows 9-10. Panel A is itself a small grid of 6 rows
% and 1 column (one row per species) with tight spacing, so the six species sit close
% together.
layoutAll = tiledlayout(10,4,'TileSpacing','compact');
layoutA = tiledlayout(layoutAll, 6, 1, 'TileSpacing', 'tight');
layoutA.Layout.Tile = 1;            % panel A starts in the top-left tile ...
layoutA.Layout.TileSpan = [6 3];    % ... and covers 6 rows and 3 columns

% AJB edits (v6): panel A as six separate panels, one per species, drawn by
% plotSpeciesPanel (end of this file): the species' samples and fit line at their real
% values, with its own y axis. This replaces Kavin's single panel A (six plot() and six
% scatter() calls with multipliers and offsets, one legend). Colours and transparency
% (last two numbers but one) are Kavin's. Arguments: layout, UCVF, values, fit line x,
% fit line y, colour, transparency, y-axis label, text written inside the panel.
axA1 = plotSpeciesPanel(layoutA, elavarasi.UCVF, elavarasi.AOS, xFit7, yFit7, [0, 0.4470, 0.7410], 1, "AOS", string7);
plotSpeciesPanel(layoutA, KAVIN_TABLE.UCVF, KAVIN_TABLE.SO4, xFit5, yFit5, [0.8500, 0.3250, 0.0980], 0.1, "SO4", string5);
plotSpeciesPanel(layoutA, KAVIN_TABLE.UCVF(ooaRows), KAVIN_TABLE.('MO-OOA')(ooaRows), xFit3, yFit3, [0.9290, 0.6940, 0.1250], 0.3, "MO-OOA", string3 + ooaLabel);
plotSpeciesPanel(layoutA, KAVIN_TABLE.UCVF(ooaRows), KAVIN_TABLE.('SR-OOA')(ooaRows), xFit6, yFit6, [0.4940, 0.1840, 0.5560], 0.1, "SR-OOA", string6 + ooaLabel);
plotSpeciesPanel(layoutA, KAVIN_TABLE.UCVF, KAVIN_TABLE.f44, xFit8, yFit8, [0.4660, 0.6740, 0.1880], 0.1, "f44", string8);
plotSpeciesPanel(layoutA, elavarasi.UCVF, elavarasi.AOSFthickly_90nm, xFit11, yFit11, [0.3010, 0.7450, 0.9330], 0.9, "FThickly", string11);
% AJB edits (v6): panel letter on the top panel, in axes units (Kavin: text(0.03, 0.97*7,'A'))
text(axA1, 0.01, 0.80, 'A', 'Units', 'normalized')

% AJB note: panel B. Fit lines HOA (yFit1), rBC (yFit2), ultrafine salt (yFit14).
% AJB edits (v6): panel B starts in tile 25 (row 7) of the 10-row grid (v2: tile 17)
nexttile(layoutAll, 25, [2,3])
plot(xFit1,yFit1)
hold on
plot(xFit2,yFit2)
plot(xFit14,yFit14)
% AJB edits (v2): HOA and rBC points no longer multiplied by 2 (Kavin: .*2), so that
% they match their unscaled fit lines (fix 3 in the header)
sone = scatter(LANLdata_HOA_rBC.UCVF,LANLdata_HOA_rBC.HOAMass_ugM_3_,[20],[0, 0.4470, 0.7410],HandleVisibility='off',Marker='.')
alpha(sone,0.4)
stwo = scatter(LANLdata_HOA_rBC.UCVF,LANLdata_HOA_rBC.rBCMass_ugM_3_,[20],[0.8500, 0.3250, 0.0980],HandleVisibility='off',Marker='.')
alpha(stwo,0.5)
sthr = scatter(ultrafinetbl.avgUCVF,ultrafinetbl.avgUfineSalt,[20],[0.9290, 0.6940, 0.1250]	,HandleVisibility='off',Marker='.')
alpha(sthr,0.5)
hold off
ylabel("Arbitrary Units")
yticklabels("")
% AJB edits (v2): names of entries 1 and 2 swapped back (entry 1 is the HOA line with
% the HOA r, string1; entry 2 the rBC line with string2), and entry 3 prints string14,
% the ultrafine salt r, instead of string3, the MO-OOA r (fixes 1 and 2). Kavin's line:
% legend("rBC, "+ string1,"HOA, "+ string2,"NaCl2, "+ string3, FontSize=5)
legend("HOA, "+ string1,"rBC, "+ string2,"NaCl2, "+ string14, FontSize=5)
% AJB edits: legend in the top-right corner, where it is in Kavin's figure
axB = gca;
axB.Legend.Location = 'northeast';
text(0.03, 0.95*0.6,'B')
xticklabels("")
ylim([0 0.6])
xticks([0:0.1:1])

% AJB note: panel C. OH burst (x20) and July dD (+85); lines and points get the same
% transformation.
% AJB edits (v6): panel C starts in tile 33 (row 9) of the 10-row grid (v2: tile 25)
nexttile(layoutAll, 33, [2,3])
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

%% Local functions
% AJB edits (v6): new.

function ax = plotSpeciesPanel(layout, ucvf, values, xFit, yFit, color, transparency, yLabelText, insideText)
% Draws one species of panel A in its own axes (the next tile of LAYOUT):
% its samples as points, its least-squares fit line (xFit, yFit from section 2), both at
% their real values, a y-axis label and a short text inside the panel (here: r).
% TRANSPARENCY is the points' alpha value (1 = solid).
ax = nexttile(layout);
plot(ax, xFit, yFit, 'Color', color)                                   % fit line
hold(ax, 'on')
points = scatter(ax, ucvf, values, 20, color, '.', 'HandleVisibility', 'off');
alpha(points, transparency)                                           % see-through points
hold(ax, 'off')
xlim(ax, [0 1])
xticks(ax, 0:0.1:1)
xticklabels(ax, "")         % x numbers only on the bottom panel of the figure (panel C)
ylim(ax, 'padded')          % y range from the data, with a small margin
ylabel(ax, yLabelText)
text(ax, 0.05, 0.80, insideText, 'Units', 'normalized', 'FontSize', 8)
end
