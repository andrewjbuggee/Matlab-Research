%% EPCAPE: UCVF correlation figure, v3 (bin means, all-month fits, new colours)
% AJB edits: plain-script copy of EPCAPEUCVFplots_kp_AJB_v3.mlx, made on 2026-10-08
% because the Live Editor hung at the export step (exportgraphics). Same code as the
% .mlx except that the figure is visible here (the .mlx hides it to avoid the hang).
% Run this file (F5, or run('EPCAPEUCVFplots_kp_AJB_v3_script.m')).
% It saves the same figure files as the .mlx: EPCAPEUCVFplots_kp_AJB_v3.png and .pdf.
%
% AJB edits (Andrew Buggee, 2026-10-08): version 3 is version 2
% (EPCAPEUCVFplots_kp_AJB_v2.mlx) redrawn to make it easier to read: in all three panels
% each species is shown by its fit line and its mean in UCVF bins, in colour-blind-friendly
% colours. See "Changes in version 3" below.
%
% How to find the changes: lines changed for version 3 are marked with a comment that
% starts with "AJB edits (v3):". "AJB edits (v2):" marks the version 2 bug fixes, and
% "AJB edits:" the version 1 changes to Kavin's code. Comments that start with
% "AJB note:" only explain the code and change nothing. Comments without a tag are
% Kavin's. The text sections (this header, Panels, Fixes, Changes, Open questions) are
% AJB's.
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
%% Changes in version 3
% # All fit lines and r values use every month. ooaAprJunOnly is false, so MO-OOA and
%   SR-OOA now use February 2023 - February 2024, like SO4 and f44 (r = 0.11 and 0.52
%   instead of 0.43 and 0.4).
% # In all three panels, the individual samples are replaced by the mean of each
%   species in 14 UCVF bins (bin edges linspace(0, 1, 15), each bin 1/14 = 0.071 wide),
%   drawn at the bin centres. The fit lines still use every sample. Some species have
%   few samples (AOS 18, FThickly 21, ultrafine salt 49, OH 71), so their bins hold only
%   a few samples each, and some bins are empty.
% # One helper function, plotFitAndBinMeans (end of this file), draws every species:
%   its fit line and its bin means, with the same multiplier and offset.
% # New colours. Panel A: colour-blind-friendly entries of mySavedColors.m (entries
%   61-73). Panels B and C: five further colour-blind-safe colours from Okabe and Ito
%   (2008) and Tol (2021), each with at least 3.5:1 contrast on white. Fit lines are
%   1 pt wide (was 0.5 pt) and the markers are filled, so the paler colours stay visible.
% # The legends show each species' line and marker together.
% # Panel A: the MO-OOA offset is 0 instead of Kavin's 2.5, and the top of the panel is
%   7.2 instead of 7 (see the comments in section 3).
% # v2's workaround for transparent '.' markers is removed, because v3 draws no
%   individual samples any more.
% Kavin's other multipliers and offsets are unchanged.
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
outDir     = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS/Fig 8 - UCVF';   % AJB edits (2026-10-14): Figure 8 files moved to the subfolder "Fig 8 - UCVF"
figName    = 'EPCAPEUCVFplots_kp_AJB_v3';  % AJB edits (v3): v2 -> v3

% AJB edits: figure size in points, [width height]. Kavin's code does not set a size.
% This one reproduces the axes of his figure as measured on the image: 471 pt wide,
% panel A 272 pt tall, panels B and C 132 pt tall (with the default 10 pt tick labels).
figSize_pt = [811 675];

% AJB edits (v3): v2's setting replacePointMarkers removed (v3 draws no '.' markers).

% AJB edits (v2): period used for MO-OOA and SR-OOA in panel A (points, fit lines and r).
% true = April-June 2023 only, as in Kavin's fits; false = the whole record
% (February 2023 - February 2024). See fix 4 in the header.
% AJB edits (v3): false (v2: true), so that every fit uses all months
ooaAprJunOnly = false;

% AJB edits (v3): how every species is drawn (by plotFitAndBinMeans, end of this file):
% its fit line plus its mean in each UCVF bin.
style.binEdges      = linspace(0, 1, 15);  % 15 bin edges = 14 UCVF bins, each 1/14 = 0.071 wide
style.lineWidth_pt  = 1;                   % width of the fit lines (MATLAB default: 0.5 pt)
style.markerSize_pt = 4;                   % size of the bin-mean markers

% AJB edits (v3): colours, RGB from 0 to 255, one row per species in legend order.
% Panel A: the colour-blind-friendly entries 61-73 of mySavedColors.m
% (Matlab-Research/startup), copied here so that the script runs without that file.
% Of all sets of six from those 13 colours, this one is the easiest to tell apart under
% simulated red-green colour blindness.
colorsA = [ 55  57 169          % AOS      - entry 64, eggplant
            99   0  25          % SO4      - entry 63, dark red
           238 172  47          % MO-OOA   - entry 65, orange
           145 127 197          % SR-OOA   - entry 71, purple
            75 114   0          % f44      - entry 67, broccoli green
             0 219 233] ./ 255; % FThickly - entry 70, electric blue
% Panels B and C: five more colour-blind-safe colours, each with at least 3.5:1
% contrast on white, from Okabe and Ito (2008), https://jfly.uni-koeln.de/color/, and
% Tol (2021), Colour Schemes, SRON/EPS/TN/09-002, https://personal.sron.nl/~pault/.
% Within each panel they stay distinct under simulated red-green colour blindness, and
% they differ from the panel A colours.
colorsB = [  0 114 178          % HOA            - Okabe-Ito blue (#0072B2)
           213  94   0          % rBC            - Okabe-Ito vermillion (#D55E00)
             0 153 136] ./ 255; % ultrafine salt - Tol vibrant teal (#009988)
colorsC = [136  34  85          % OH             - Tol muted wine (#882255)
           238  51 119] ./ 255; % July dD        - Tol vibrant magenta (#EE3377)
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

% AJB edits (v3): every species is drawn by plotFitAndBinMeans (end of this file): its
% fit line, which uses every sample, and its mean in each UCVF bin. The arguments are:
% UCVF, values, fit line x, fit line y, multiplier, offset, colour, style. This replaces
% Kavin's plot() and scatter() calls; the multipliers and offsets are his unless a
% comment says otherwise.

% ---------- Panel A ----------
nexttile(1,[4,3])
hold on
plotFitAndBinMeans(elavarasi.UCVF, elavarasi.AOS, xFit7, yFit7, 1, 0, colorsA(1,:), style)        % AOS
plotFitAndBinMeans(KAVIN_TABLE.UCVF, KAVIN_TABLE.SO4, xFit5, yFit5, 2, 3, colorsA(2,:), style)     % SO4

% AJB edits (v3): MO-OOA offset changed from Kavin's 2.5 to 0.
%   Why: v3 uses all months for MO-OOA (ooaAprJunOnly = false). Over all months the
%   MO-OOA values are on average about 3 times larger than in April-June 2023 (mean 0.77
%   vs. 0.27, in the units of KAVIN_TABLE). With Kavin's offset of 2.5, the MO-OOA line
%   therefore sat on top of the AOS line and ran off the top of the panel.
%   With offset 0, the line sits between SO4 and SR-OOA, where Kavin's April-June line was.
%   The offset only moves the line and markers up or down. It does not change the slope
%   or r, and the y axis has no numbers ("Arbitrary Units").
plotFitAndBinMeans(KAVIN_TABLE.UCVF(ooaRows), KAVIN_TABLE.("MO-OOA")(ooaRows), xFit3, yFit3, 5, 0, colorsA(3,:), style)   % MO-OOA

plotFitAndBinMeans(KAVIN_TABLE.UCVF(ooaRows), KAVIN_TABLE.("SR-OOA")(ooaRows), xFit6, yFit6, 7, 1.8, colorsA(4,:), style) % SR-OOA
plotFitAndBinMeans(KAVIN_TABLE.UCVF, KAVIN_TABLE.f44, xFit8, yFit8, 25, -3, colorsA(5,:), style)                          % f44
plotFitAndBinMeans(elavarasi.UCVF, elavarasi.AOSFthickly_90nm, xFit11, yFit11, 15, -4, colorsA(6,:), style)              % FThickly
% AJB note: no3index is not used
no3index = KAVIN_TABLE.NO3.*95 -29 < 6
hold off
% AJB edits (v3): top of the panel raised from 7 to 7.2, so that the highest AOS bin
% mean (at 7.0) is not cut in half
ylim([0 7.2])
ylabel("Arbitrary Units")
yticklabels("")
% AJB edits (v2): ooaLabel added to the MO-OOA and SR-OOA entries (fix 4). In v3 it is
% empty, because v3 uses all months.
lgnd = legend("AOS, "+ string7,"SO4, "+string5,"MO-OOA"+ooaLabel+", "+string3,"SR-OOA"+ooaLabel+", " + string6,"f44, "+string8,"FThickly, "+string11)
lgnd.FontSize = 5
% AJB edits: legend in the bottom-right corner, where it is in Kavin's figure
lgnd.Location = 'southeast';
text(0.03, 0.97*7,'A')
xticklabels([""])
xticks([0:0.1:1])

% ---------- Panel B ----------
% AJB edits (v2): HOA and rBC are not multiplied by 2 (fix 3)
nexttile(17,[2,3])
hold on
plotFitAndBinMeans(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.HOAMass_ugM_3_, xFit1, yFit1, 1, 0, colorsB(1,:), style)   % HOA
plotFitAndBinMeans(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.rBCMass_ugM_3_, xFit2, yFit2, 1, 0, colorsB(2,:), style)   % rBC
plotFitAndBinMeans(ultrafinetbl.avgUCVF, ultrafinetbl.avgUfineSalt, xFit14, yFit14, 1, 0, colorsB(3,:), style)        % ultrafine salt
hold off
ylabel("Arbitrary Units")
yticklabels("")
% AJB edits (v2): names of entries 1 and 2 swapped back, and entry 3 prints string14,
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

% ---------- Panel C ----------
% OH multiplied by 20 and dD shifted by +85, as in Kavin's figure
nexttile(25,[2,3])
hold on
plotFitAndBinMeans(OHvsUCVF.UCVF24, OHvsUCVF.OH_burst_25000_vv_dilution_ratio_uM, xFit10, yFit10, 20, 0, colorsC(1,:), style)   % OH
plotFitAndBinMeans(julyUCVF, julydD, xFit13, yFit13, 1, 85, colorsC(2,:), style)                                                 % July dD
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

% AJB edits (v3): v2's '.'-marker workaround removed (v3 draws no '.' markers).
%% 4. Save the figure
% AJB edits: this whole section is new. exportgraphics crops the figure to its
% contents, which removes the empty fourth tile column.

if saveFigure
    exportgraphics(fig, fullfile(outDir, [figName '.png']), 'Resolution', 300)
    exportgraphics(fig, fullfile(outDir, [figName '.pdf']), 'ContentType', 'vector')
    fprintf('Saved %s.png and %s.pdf in %s\n', figName, figName, outDir)
end

%% Local functions
% AJB edits (v3): new.

function plotFitAndBinMeans(ucvf, values, xFit, yFit, multiplier, offset, color, style)
% Draws one species in the current axes:
%   1. its least-squares fit line (xFit, yFit from section 2), which uses every sample;
%   2. its mean in each UCVF bin (bins: style.binEdges), as filled circles at the bin
%      centres.
% Both are multiplied by MULTIPLIER and shifted by OFFSET, so they line up.
% The fit line carries the circle marker too, with MarkerIndices = [] so that no
% circles are drawn on the line itself; the legend then shows line and circle together.
% The bin means are kept out of the legend.

% 1. Fit line
plot(xFit, yFit.*multiplier + offset, '-o', 'Color', color, ...
    'LineWidth', style.lineWidth_pt, 'MarkerFaceColor', color, ...
    'MarkerSize', style.markerSize_pt, 'MarkerIndices', [])

% 2. Bin means. Only rows where both UCVF and the value exist are used, as in the fit.
ok = ~isnan(ucvf) & ~isnan(values);
ucvfOk = ucvf(ok);
valuesOk = values(ok);
binNumber = discretize(ucvfOk, style.binEdges);   % which bin each sample falls in
nBins = numel(style.binEdges) - 1;
binMean = NaN(nBins, 1);                          % stays NaN for empty bins
for b = 1:nBins
    inBin = (binNumber == b);
    if any(inBin)
        binMean(b) = mean(valuesOk(inBin));
    end
end
binCentre = (style.binEdges(1:end-1) + style.binEdges(2:end))' / 2;
plot(binCentre, binMean.*multiplier + offset, 'o', 'Color', color, ...
    'MarkerFaceColor', color, 'MarkerSize', style.markerSize_pt, ...
    'HandleVisibility', 'off')
end
