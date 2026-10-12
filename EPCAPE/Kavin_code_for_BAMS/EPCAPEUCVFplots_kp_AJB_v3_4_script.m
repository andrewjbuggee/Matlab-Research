%% EPCAPE: UCVF correlation figure (Figure 8), v3-4 (FThickly from individual samples)
% AJB edits: plain-script copy of EPCAPEUCVFplots_kp_AJB_v3_4.mlx, made on 2026-10-08
% because the Live Editor hung at the export step (exportgraphics). Same code as the
% .mlx except that figureVisible = 'on' here, so the figure is shown in a window.
% Run this file (Run button, or run('EPCAPEUCVFplots_kp_AJB_v3_4_script.m')).
% It saves the same figure files as the .mlx: EPCAPEUCVFplots_kp_AJB_v3_4.png and .pdf.
%
% AJB edits (Andrew Buggee, 2026-10-12): version 3-4 is version 3-3
% (EPCAPEUCVFplots_kp_AJB_v3_3.mlx) with new data from Elavarasi Ravichandran (file
% "AOS rBC coating thickness, rBC Mass, and Fthickly Number.xlsx"): FThickly is now
% binned from the individual samples, with the bin median, an IQR band and a fit to the
% bin medians, like the other species; and the AOS series (rBC mass per particle) gets a
% band of plus/minus one standard deviation. See "Changes in version 3-4" below.
% (File names use v3_4 because MATLAB cannot run a file whose name contains a hyphen.)
% Version 3-3 adds the Pier dD and the spread in each bin; version 3-2 computes the fits
% and r from the bins; version 3-1 scales up the ultrafine salt; version 3 redraws
% version 2 with bins and colour-blind-friendly colours.
%
% How to find the changes: lines changed for version 3-4 are marked "AJB edits (v3-4):",
% for version 3-3 "AJB edits (v3-3):", for version 3-2 "AJB edits (v3-2):", for version 3-1 "AJB edits (v3-1):", and for
% version 3 "AJB edits (v3):". "AJB edits (v2):" marks the version 2 bug fixes, and
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
%% How to run this file
% # Folders: in "User settings" below, switch to your four folder lines (dataDir,
%   outDir, pierDir and aosDir) as described there. The six original input files must be
%   in dataDir. Kavin's original code read two of them from ".../EPCAPE/UCVF plot/" and
%   the other four from MATLAB's current folder, so copy those four into that folder if
%   needed. Elba Granados Escartin's file with both dD sites (dD_UCVF_soledad_pier.csv)
%   must be in pierDir, and Elavarasi Ravichandran's file "AOS rBC coating thickness, rBC
%   Mass, and Fthickly Number.xlsx" in aosDir.
% # Press Run. The figure is saved in outDir as EPCAPEUCVFplots_kp_AJB_v3_4.png and .pdf.
% # In the .mlx file the figure is hidden while it is made (figureVisible = 'off'); open
%   the saved PNG to see it. EPCAPEUCVFplots_kp_AJB_v3_4_script.m is the same code as a
%   plain script and shows the figure in a window.
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
%% Changes in version 3-1
% # Panel B: the ultrafine salt is multiplied by saltMultiplier (3; User settings).
%   Plotted at its real values it looked almost flat, although its correlation with UCVF
%   (r = -0.68) is the strongest in the panel. With the factor 3, its fit line drops by
%   about a fifth of the panel height across its UCVF range, instead of 6 %. It stays
%   below the rBC line except at the far right, where all three lines approach zero.
%   Like Kavin's panel A multipliers, this changes only how the salt is drawn: its r and
%   the legend are unchanged, and the y axis has no numbers ("Arbitrary Units").
%
%% Changes in version 3-2
% # Fits to bin means. For every variable except AOS and FThickly, the samples are first
%   averaged in the 14 UCVF bins (the same bins as the markers), and the fit line and r
%   are computed from those bin means, in the same style as Kavin's fits (new section 2c;
%   setting fitToBinMeans). The legends call this r "r_bin". AOS and FThickly keep the
%   fits on their own data, because the data owners already averaged them.
% # MO-OOA and SR-OOA months are set separately: MO-OOA uses April-June 2023
%   (moAprJunOnly = true, as in Kavin's file), SR-OOA uses all months
%   (srAprJunOnly = false). v2 to v3-1 used one setting for both.
% # MO-OOA offset back to Kavin's 2.5, because MO-OOA uses April-June again (it was 0
%   in v3 and v3-1, where MO-OOA used all months).
% # Section 2c prints a table of r from every sample and r from the bin means.
% # AOS and FThickly are drawn with their own points (18 and 21), not averaged again in
%   the UCVF bins (v3 and v3-1 averaged them for the markers). All other species are
%   shown by their bin means, as before.
% # Legend text 8 pt instead of Kavin's 5 pt (style.legendFontSize_pt), and the panel A
%   legend moved to the top-left corner. The top of panel A is raised to 8.4 (v3-1: 7.2)
%   so that the legend sits above the AOS data, and the letter "A" is placed in axes
%   units so that it stays in the corner.
% How to read r_bin: averaging removes the scatter within each bin, so r from bin means
% is higher than r from the samples. r_bin says how close the bin means lie to a
% straight line, not how well UCVF predicts a single sample, and r_bin^2 is not the
% fraction of the sample variance explained. Each bin counts equally in the fit,
% whether it holds 1 sample or 400.
%
%% Changes in version 3-3
% # Panel C: Pier dD added, next to the Mt. Soledad dD that was already in the figure.
%   The Pier UCVF and dD come from Elba Granados Escartin's file with both sites,
%   dD_UCVF_soledad_pier.csv (columns time, UCVF, soledad_dD, pier_dD; every 10 min
%   from 14 July 2023), so that the Pier dD and its UCVF come from one consistent file.
%   The Soledad dD stays as in Kavin's figure (dD_UCVFreal.xlsx). Both sites use July
%   2023, which for both files means 14-31 July (both records start on 14 July). The
%   Pier dD is binned and fitted like the other variables, and drawn with the same shift
%   (+85) as the Soledad dD, so the difference between the two sites can be read
%   directly.
% # Panel C legend names the site of each dD record ("07/2023"), and sits in the
%   bottom-left corner. OH is shifted up by 10 (Kavin: no shift), so that it does not
%   overlap the Pier dD; the y range is -24 to 23 (Kavin: -20 to 20).
% # Each bin is shown by its median (style.binStatistic = "median"; v3-2: the mean), and
%   the fit lines and r_bin are computed from the bin medians, so the lines follow the
%   markers. The median always lies inside the interquartile range and is not pulled up
%   by a few very high values.
% # Spread of the data in each bin, drawn behind every bin marker (settings style.spread
%   and style.spreadDisplay). Default: the middle half of the samples in the bin (25th to
%   75th percentile, the interquartile range), shown as a light shaded band. A label in
%   panel A, "IQR (25th - 75th pct.)", with a short pointer line (pointing up and
%   slightly left) from the top of the uppermost IQR band (SO4) at UCVF = 0.7 explains
%   it. Thin light bars, the standard deviation and the standard error are other
%   options.
%   The spread is multiplied by the same factor as the data. Bins with a single sample
%   show no spread. AOS and FThickly show none, because only their already-averaged
%   values are available.
% # Panel C colours: Pier dD in Tol bright green (#228833), which stays distinct from
%   the OH and Soledad colours, also for red-green colour blindness.
%
%% Changes in version 3-4
% New data from Elavarasi Ravichandran (Google Sheet shared 12 October 2026, file "AOS
% rBC coating thickness, rBC Mass, and Fthickly Number.xlsx"; CNW inlet at Mt. Soledad,
% April to July 2023, only times when the CVI was off). Its four tabs:
%   "Coating thickness": rBC coating thickness (nm), individual samples. Not in this
%       figure (Andrew's choice).
%   "rBC mass": rBC mass per particle (fg), already averaged in UCVF bins (bin mean and
%       standard deviation), April to July. This is the series that Kavin's figure calls
%       "AOS" (his file holds the same quantity for April).
%   "Fthickly" and "Fthickly < 90nm": FThickly, individual samples (the two tabs hold
%       identical numbers).
% Kavin's file (Petters and Liggio, "- April") holds the April data only: binning the
% April samples of the new FThickly tab reproduces every one of its values.
% # Period: April 2023 only, as in Kavin's figure (setting aosAprilOnly; false = April to
%   July). The legend says "(Apr 2023)".
% # FThickly is binned from the individual April samples (3236 samples) into the same 14
%   UCVF bins as the other species: bin medians as markers, IQR band, fit line and r_bin
%   from the bin medians (v3-3: Kavin's 21 pre-averaged points and a fit to them).
%   FThickly is drawn x8 - 2.4 (Kavin: x15 - 4), so that its IQR band sits at the bottom
%   of panel A instead of covering most of it; panel A now runs from -1.3 to 8.4 (v3-3:
%   0 to 8.4).
% # AOS (rBC mass per particle): only bin means and standard deviations exist, so no
%   median or IQR can be computed. It keeps Kavin's April bin means and their fit, and
%   gets a band of plus/minus one standard deviation (Kavin's column "Mass uncertainty",
%   which equals the "STD" column of the new file where the two hold the same data). The
%   AJB edits (v3-4): the series is named "rBC mass per particle (AOS)" in the legend
%   (Kavin: "AOS"), and a note at the bottom right of panel A says that its band is the
%   mean plus/minus one standard deviation, because all other bands are IQRs.
% # Andrew's edits to the v3-3 script are included: panel labels at the top left in
%   pnl_fnt = 15 pt (here the plain capital letters A, B, C, not bold, as Lynn asked;
%   v3-3: "(a)", "(b)", "(c)" in bold), the IQR label in 12 pt, the note on the bands in panel A (reworded
%   for v3-4), and "\deltaD, Mt. Soledad" / "\deltaD, Pier" in the panel C legend.
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
% To switch, remove the "% " at the start of Kavin's lines and add "% " at the start of
% Andrew's lines.
%   dataDir: folder with the six original input files listed below
%   outDir:  folder where the figure (PNG and PDF) is saved
%   pierDir: folder with Elba's file with both dD sites, dD_UCVF_soledad_pier.csv
%            (AJB edits (v3-3): new)
%   aosDir:  folder with Elavarasi Ravichandran's file "AOS rBC coating thickness, rBC
%            Mass, and Fthickly Number.xlsx" (AJB edits (v3-4): new)
% --- Andrew ---
dataDir = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper';
outDir  = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS';
pierDir = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed/UCVF - Elba Granados Escartin';
aosDir  = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed';
% --- Kavin (the folder used in Kavin's original EPCAPEUCVFplots_kp.mlx) ---
% dataDir = '/Users/2fast2kavin/Desktop/Atmospheric Aerosol Group/EPCAPE/UCVF plot';
% outDir  = dataDir;
% pierDir = dataDir;      % put dD_UCVF_soledad_pier.csv in the same folder
% aosDir  = dataDir;      % put "AOS rBC coating thickness, rBC Mass, and Fthickly Number.xlsx" in the same folder

% AJB edits: input files (all in dataDir)
ultrafineFile = fullfile(dataDir, 'ucvf vs ultrafine salt.xlsx');      % daily mean UCVF and ultrafine salt (panel B)
lanlFile      = fullfile(dataDir, 'LANL_UCVF_Correlations.xlsx');      % UCVF, HOA and rBC mass in ug m^-3 (panel B)
ohFile        = fullfile(dataDir, 'EPCAPE_overview_UCLAgraph.csv');    % OH burst (uM), RH and UCVF24 per sample (panel C)
kavinFile     = fullfile(dataDir, 'kavin_UCVF.mat');                   % KAVIN_TABLE: SR-OOA, MO-OOA, f44, NO3, SO4, UCVF (panel A)
dDFile        = fullfile(dataDir, 'dD_UCVFreal.xlsx');                  % time, dD and UCVF (panel C)
aqueousFile   = fullfile(dataDir, 'Cloud aqueous reaction data (Petters and Liggio) - April.csv');  % UCVF, AOS, AOS F_thickly 90 nm (panel A)
% AJB edits (v3-3): Elba Granados Escartin's dD file with both sites (time, UCVF,
% soledad_dD, pier_dD; every 10 min from 14 July 2023). Only its UCVF and pier_dD are
% used; the Soledad dD in the figure is the one in dDFile, as in Kavin's figure.
pierFile      = fullfile(pierDir, 'dD_UCVF_soledad_pier.csv');          % Pier dD and UCVF (panel C)
% AJB edits (v3-4): Elavarasi Ravichandran's AOS file (tabs "Coating thickness", "rBC
% mass", "Fthickly", "Fthickly < 90nm"; CNW inlet, Mt. Soledad, April-July 2023, CVI off)
aosFile       = fullfile(aosDir, 'AOS rBC coating thickness, rBC Mass, and Fthickly Number.xlsx');   % FThickly samples, rBC mass bins (panel A)

% AJB edits: figure output
saveFigure = true;          % save <figName>.png (300 dpi) and <figName>.pdf in outDir (set above)
figName    = 'EPCAPEUCVFplots_kp_AJB_v3_4';  % AJB edits (v3-4): v3_3 -> v3_4

% AJB edits: figure size in points, [width height]. Kavin's code does not set a size.
% This one reproduces the axes of his figure as measured on the image: 471 pt wide,
% panel A 272 pt tall, panels B and C 132 pt tall (with the default 10 pt tick labels).
% AJB edits: 'off' = the figure is made hidden; 'on' = it is shown in a window. The
% saved PNG and PDF are the same either way. Keep 'off' in the .mlx file: with a visible
% figure, the Live Editor (MATLAB R2024b) hangs when it saves the figure. The plain-script
% copy (..._script.m) uses 'on'.
figureVisible = 'on';   % AJB edits: 'on' in this plain-script copy

figSize_pt = [811 675];

% AJB edits (v3): v2's setting replacePointMarkers removed (v3 draws no '.' markers).

% AJB edits (v3-2): months used for MO-OOA and for SR-OOA in panel A (points, fit
% lines and r), now set separately (v2 to v3-1 used one setting, ooaAprJunOnly, for
% both). true = April-June 2023 only, as in Kavin's file; false = the whole record
% (February 2023 - February 2024).
moAprJunOnly = true;    % MO-OOA: April-June 2023, as in Kavin's file
srAprJunOnly = false;   % SR-OOA: all months

% AJB edits (v3-2): true = fit lines and r from the bin means (all variables except AOS
% and FThickly; section 2c). false = from every sample, as in v3-1.
fitToBinMeans = true;

% AJB edits (v3): how every species is drawn (by plotFitAndBinMeans, end of this file):
% its fit line plus its mean in each UCVF bin.
style.binEdges      = linspace(0, 1, 15);  % 15 bin edges = 14 UCVF bins, each 1/14 = 0.071 wide
style.lineWidth_pt  = 1;                   % width of the fit lines (MATLAB default: 0.5 pt)
style.markerSize_pt = 4;                   % size of the bin-mean markers
% AJB edits (v3-2): size of the legend text (Kavin: 5 pt, too small to read)
style.legendFontSize_pt = 8;
% AJB edits (v3-3): spread of the data in each UCVF bin, drawn behind the bin means.
% What is shown (style.spread):
%   "iqr"  = middle half of the samples: 25th to 75th percentile (interquartile range).
%            Recommended: many of these variables are skewed (a few very high values),
%            and the IQR follows that shape and never goes below the lowest sample.
%   "std"  = bin mean plus and minus one standard deviation (for skewed data this can
%            reach below zero, which concentrations cannot do)
%   "sem"  = bin mean plus and minus the standard error of the mean (std divided by the
%            square root of the number of samples): how well the bin mean is known, not
%            how much the samples spread; much narrower
%   "none" = no spread shown
% How it is shown (style.spreadDisplay):
%   "bars" = a thin, light vertical bar at each bin, without end caps (recommended: the
%            bars stay behind the markers and lines and do not hide the other species)
%   "band" = a light, see-through band joining the bins. Clean only when the species are
%            well apart (as in panel C); in panels A and B the bands overlap and the
%            figure gets muddy.
style.spread        = "iqr";
style.spreadDisplay = "band";  % AJB edits (v3-3): IQR shaded bands, chosen by Andrew
style.bandAlpha     = 0.2;     % band opacity: 0 = invisible, 1 = solid
style.barLineWidth_pt = 0.5;   % width of the bars (if spreadDisplay = "bars")
style.barLightening = 0.5;     % bar colour mixed with white: 0 = species colour, 1 = white
% AJB edits (v3-3): value that represents each bin, for the markers and for the fit lines
% and r_bin of section 2c: "median" (middle sample of the bin; it always lies inside the
% interquartile range and is not pulled up by a few very high values) or "mean".
style.binStatistic = "median";
% AJB edits (v3-3): label for the spread in panel A, with a short pointer line from the top
% edge of the uppermost band (SO4) at this UCVF, pointing up and slightly to the left.
spreadLabelUCVF = 0.7;

% AJB edits (v3-4): period of the AOS data (rBC mass per particle, "AOS"; FThickly).
% true = April 2023 only, as in Kavin's figure (Andrew's choice); false = April to July
% 2023, all the data in Elavarasi's file.
aosAprilOnly = true;

% AJB edits (v3-1): factor that the ultrafine salt (panel B) is multiplied by before it is
% drawn, so that its slope can be seen. 1 = real values (as in v3). With 5 its line sits
% just above the rBC line; with 10 just above the HOA line.
saltMultiplier = 3;

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
           238  51 119          % Soledad dD     - Tol vibrant magenta (#EE3377)
            34 136  51] ./ 255; % Pier dD        - Tol bright green (#228833), AJB edits (v3-3)
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
% AJB edits (v3-2): rows used for MO-OOA (moRows) and for SR-OOA (srRows) for points,
% fits and r, and the matching legend suffixes (moLabel, srLabel). Kavin selected
% April-June 2023 (dmoso 4-6) inside both fits. (v2 to v3-1: one set, ooaRows/ooaLabel.)
aprJun  = dmoso <= 6 & dmoso >= 4;          % April-June 2023
allRows = true(height(KAVIN_TABLE), 1);     % every month
if moAprJunOnly
    moRows = aprJun;   moLabel = " (Apr-Jun 2023)";
else
    moRows = allRows;  moLabel = "";
end
if srAprJunOnly
    srRows = aprJun;   srLabel = " (Apr-Jun 2023)";
else
    srRows = allRows;  srLabel = "";
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

% AJB edits (v3-3): Pier dD and UCVF from Elba's file with both sites, July 2023.
% The file starts with an unnamed index column, then time, UCVF, soledad_dD, pier_dD.
% ReadVariableNames makes sure the first line is read as the column names.
pierTable = readtable(pierFile, 'ReadVariableNames', true, 'Delimiter', ',');
% Times are text such as '2023-07-14 22:20:00+00:00' (UTC). Keep the first 19 characters,
% as for the Soledad times above. (If a MATLAB release reads them as datetimes directly,
% use them as they are, without the time zone.)
if isdatetime(pierTable.time)
    pierTime = pierTable.time;
    pierTime.TimeZone = '';
else
    pierTime = datetime(extractBetween(string(pierTable.time), 1, 19), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
end
% Keep July 2023, as for the Soledad dD (both records start on 14 July 2023).
inJuly = month(pierTime) == 7 & year(pierTime) == 2023;
pierdD   = pierTable.pier_dD(inJuly);      % Pier dD, July 2023
pierUCVF = pierTable.UCVF(inJuly);         % UCVF at the same times, from the same file

% AJB edits: "Cloud aqueous reaction data (Petters and Liggio) - April.csv" (current folder) -> aqueousFile
% AJB note: Petters and Liggio cloud aqueous reaction data, April (UCVF, AOS, AOSFthickly_90nm, ...)
elavarasi = readtable(aqueousFile)

% AJB edits (v3-4): AOS data for panel A, from Elavarasi Ravichandran's file.
% FThickly, individual samples (tab "Fthickly < 90nm"; the tab "Fthickly" holds the same
% numbers). Columns: Datetime, UCVF (-), Fthickly (-). Missing values are the text "NaN".
fthickTable = readtable(aosFile, 'Sheet', 'Fthickly < 90nm', 'Range', 'A:C', ...
    'ReadVariableNames', true, 'TreatAsMissing', 'NaN', 'VariableNamingRule', 'preserve');
if aosAprilOnly
    keepFthick = month(fthickTable{:,1}) == 4 & year(fthickTable{:,1}) == 2023;   % April 2023
    aosLabel = " (Apr 2023)";
    aosText  = "Apr 2023";      % AJB edits (v3-4): the same period, for the AOS legend entry
else
    keepFthick = true(height(fthickTable), 1);                                     % April-July 2023
    aosLabel = " (Apr-Jul 2023)";
    aosText  = "Apr-Jul 2023";
end
fthickUCVF = fthickTable{keepFthick, 2};   % UCVF of each FThickly sample
fthick     = fthickTable{keepFthick, 3};   % FThickly of each sample
% rBC mass per particle (fg), the series called "AOS". Only bin means and standard
% deviations exist (no individual samples).
if aosAprilOnly
    % Kavin's April file: columns AOS (Mass_fg) and Var3 (Mass uncertainty = standard
    % deviation; it equals the "STD" of Elavarasi's new file where both hold the same data)
    rbcMassUCVF = elavarasi.UCVF;
    rbcMass     = elavarasi.AOS;
    rbcMassStd  = elavarasi.Var3;
else
    % Elavarasi's new file, tab "rBC mass" (April-July): UCVF (mean of each bin), Mass (fg), STD
    massTable = readtable(aosFile, 'Sheet', 'rBC mass', 'Range', 'A:C', ...
        'ReadVariableNames', true, 'TreatAsMissing', 'NaN', 'VariableNamingRule', 'preserve');
    massTable = massTable(~isnan(massTable{:,1}), :);    % drop empty rows below the table
    rbcMassUCVF = massTable{:,1};
    rbcMass     = massTable{:,2};
    rbcMassStd  = massTable{:,3};
end
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
% AJB edits (v3-2): rows moRows (v2 to v3-1: ooaRows; Kavin: dmoso <= 6 & dmoso >= 4) in
% the four lines below, so that the period follows moAprJunOnly (Kavin's rows when true).
model = fitlm(KAVIN_TABLE.UCVF(moRows),KAVIN_TABLE.("MO-OOA")(moRows));
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit3 = linspace(min(KAVIN_TABLE.UCVF(moRows)), max(KAVIN_TABLE.UCVF(moRows)), 100);
yFit3 = feval(model, xFit3);

[r,p]=corrcoef(KAVIN_TABLE.UCVF(moRows), KAVIN_TABLE.("MO-OOA")(moRows),'rows','pairwise');
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
% AJB edits (v3-2): rows srRows (v2 to v3-1: ooaRows; Kavin: dmoso <= 6 & dmoso >= 4) in
% the four lines below, so that the period follows srAprJunOnly (Kavin's rows when true).
model = fitlm(KAVIN_TABLE.UCVF(srRows),KAVIN_TABLE.("SR-OOA")(srRows));
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit6 = linspace(min(KAVIN_TABLE.UCVF(srRows)), max(KAVIN_TABLE.UCVF(srRows)), 100);
yFit6 = feval(model, xFit6);

[r,p]=corrcoef(KAVIN_TABLE.UCVF(srRows), KAVIN_TABLE.("SR-OOA")(srRows),'rows','pairwise');
r=r(2,1);
p=p(2,1);
str1="r = ";
str2=string(round(r,2));
string6 =append(str1,str2);

% AJB note: panel A, AOS (fit 7): rBC mass per particle, bin means
% AJB edits (v3-4): rbcMassUCVF and rbcMass (section 1) instead of elavarasi.UCVF and
% elavarasi.AOS; with aosAprilOnly = true they are the same numbers.
model = fitlm(rbcMassUCVF,rbcMass);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit7 = linspace(min(rbcMassUCVF), max(rbcMassUCVF), 100);
yFit7 = feval(model, xFit7);

[r,p]=corrcoef(rbcMassUCVF, rbcMass,'rows','pairwise');
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

% AJB note: panel A, FThickly (fit 11)
% AJB edits (v3-4): fitted to the individual FThickly samples (fthickUCVF, fthick from
% section 1) instead of Kavin's pre-averaged points (elavarasi.AOSFthickly_90nm). With
% fitToBinMeans = true it is replaced by the fit to the bin medians in section 2c.
model = fitlm(fthickUCVF,fthick);
b = model.Coefficients.Estimate(1);
m = model.Coefficients.Estimate(2);
xFit11 = linspace(min(fthickUCVF), max(fthickUCVF), 100);
yFit11 = feval(model, xFit11);

[r,p]=corrcoef(fthickUCVF,fthick,'rows','pairwise');
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

% AJB edits (v3-3): panel C, Pier dD (fit 15), same style as Kavin's fits. With
% fitToBinMeans = true it is replaced by the fit to the bin means in section 2c.
model = fitlm(pierUCVF,pierdD);
xFit15 = linspace(min(pierUCVF), max(pierUCVF), 100);
yFit15 = feval(model, xFit15);
[r,p]=corrcoef(pierUCVF,pierdD,'rows','pairwise');
r=r(2,1);
string15=append("r = ",string(round(r,2)));

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
%% 2c. Fits to the bin means
% AJB edits (v3-2): this whole section is new. For every variable except AOS (and, until
% v3-3, FThickly), the samples are first summarised in the UCVF bins (style.binEdges, the same
% bins as the markers). AJB edits (v3-3): each bin is summarised by its median or its
% mean (style.binStatistic; median by default), the same value as the marker. The straight-line fit (fitlm) and r (corrcoef) are then computed
% from the bin means, one point per bin that holds samples, placed at the bin centre.
% The results replace the fits of section 2 for these variables (xFitN, yFitN and
% stringN are overwritten), so the figure code below does not change. AOS (fit 7) and
% FThickly (fit 11) keep the fits on their own data, because the data owners already
% averaged them. fitBinMeans is defined at the end of this file.
if fitToBinMeans
    % AJB edits (v3-3): Pier dD added (10th row); "July dD" renamed "Soledad dD"
    % AJB edits (v3-4): FThickly added (11th row)
    variable = ["SO4"; "MO-OOA" + moLabel; "SR-OOA" + srLabel; "f44"; "HOA"; "rBC"; "ultrafine salt"; "OH"; "Soledad dD (07/2023)"; "Pier dD (07/2023)"; "FThickly" + aosLabel];
    rSamples = NaN(11, 1);    % r from every sample (as in v3-1)
    rBinMeans = NaN(11, 1);   % r from the bin values (medians or means; used in the figure)
    nBins = NaN(11, 1);       % number of bins that hold samples
    [xFit5,  yFit5,  string5,  rBinMeans(1), rSamples(1), nBins(1)] = fitBinMeans(KAVIN_TABLE.UCVF, KAVIN_TABLE.SO4, style.binEdges, style.binStatistic);
    [xFit3,  yFit3,  string3,  rBinMeans(2), rSamples(2), nBins(2)] = fitBinMeans(KAVIN_TABLE.UCVF(moRows), KAVIN_TABLE.("MO-OOA")(moRows), style.binEdges, style.binStatistic);
    [xFit6,  yFit6,  string6,  rBinMeans(3), rSamples(3), nBins(3)] = fitBinMeans(KAVIN_TABLE.UCVF(srRows), KAVIN_TABLE.("SR-OOA")(srRows), style.binEdges, style.binStatistic);
    [xFit8,  yFit8,  string8,  rBinMeans(4), rSamples(4), nBins(4)] = fitBinMeans(KAVIN_TABLE.UCVF, KAVIN_TABLE.f44, style.binEdges, style.binStatistic);
    [xFit1,  yFit1,  string1,  rBinMeans(5), rSamples(5), nBins(5)] = fitBinMeans(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.HOAMass_ugM_3_, style.binEdges, style.binStatistic);
    [xFit2,  yFit2,  string2,  rBinMeans(6), rSamples(6), nBins(6)] = fitBinMeans(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.rBCMass_ugM_3_, style.binEdges, style.binStatistic);
    [xFit14, yFit14, string14, rBinMeans(7), rSamples(7), nBins(7)] = fitBinMeans(ultrafinetbl.avgUCVF, ultrafinetbl.avgUfineSalt, style.binEdges, style.binStatistic);
    [xFit10, yFit10, string10, rBinMeans(8), rSamples(8), nBins(8)] = fitBinMeans(OHvsUCVF.UCVF24, OHvsUCVF.OH_burst_25000_vv_dilution_ratio_uM, style.binEdges, style.binStatistic);
    [xFit13, yFit13, string13, rBinMeans(9), rSamples(9), nBins(9)] = fitBinMeans(julyUCVF, julydD, style.binEdges, style.binStatistic);
    [xFit15, yFit15, string15, rBinMeans(10), rSamples(10), nBins(10)] = fitBinMeans(pierUCVF, pierdD, style.binEdges, style.binStatistic);   % AJB edits (v3-3)
    [xFit11, yFit11, string11, rBinMeans(11), rSamples(11), nBins(11)] = fitBinMeans(fthickUCVF, fthick, style.binEdges, style.binStatistic);   % AJB edits (v3-4): FThickly
    % r from every sample and from the bin values, side by side (rounded to 2 decimals)
    % AJB edits (v3-3): column r_bins (v3-2: r_binMeans), as the bins can be medians
    rComparison = table(variable, round(rSamples, 2), round(rBinMeans, 2), nBins, ...
        'VariableNames', {'variable', 'r_samples', 'r_bins', 'binsWithData'})
end
%% 3. Figure

% AJB edits (v3-3, by Andrew): font size of the panel labels (A, B, C)
pnl_fnt = 15;

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

% AJB note: 8 x 4 grid of tiles. The panels span columns 1-3 only (column 4 stays empty
% and exportgraphics crops it away): A = rows 1-4, B = rows 5-6, C = rows 7-8.
tiledlayout(8,4,'TileSpacing','compact')

% AJB edits (v3): every species is drawn by plotFitAndBinMeans (end of this file): its
% fit line, which uses every sample, and its mean in each UCVF bin. The arguments are:
% UCVF, values, fit line x, fit line y, multiplier, offset, colour, style. This replaces
% Kavin's plot() and scatter() calls; the multipliers and offsets are his unless a
% comment says otherwise.

% ---------- Panel A ----------
% AJB edits (v3-2): the axes handle axA is kept, to place the legend (end of section 3)
axA = nexttile(1,[4,3]);
hold on
% AJB edits (v3-2): AOS is drawn with its own points (plotFitAndOwnPoints, end of this
% file), not averaged again in the UCVF bins: its owners already averaged it.
% AJB edits (v3-4): AOS from rbcMassUCVF and rbcMass (section 1), and a band of plus/minus
% one standard deviation (rbcMassStd), drawn like the other bands (style.spreadDisplay).
% It is a different measure of spread from the IQR of the other species, because only
% bin means and standard deviations exist for AOS; the legend says so.
plotFitAndOwnPoints(rbcMassUCVF, rbcMass, xFit7, yFit7, 1, 0, colorsA(1,:), style)       % AOS
if style.spread ~= "none"
    if style.spreadDisplay == "band"
        drawSpreadBand(rbcMassUCVF, rbcMass - rbcMassStd, rbcMass + rbcMassStd, colorsA(1,:), style.bandAlpha)
    else
        drawSpreadBars(rbcMassUCVF, rbcMass - rbcMassStd, rbcMass + rbcMassStd, colorsA(1,:), style)
    end
end
% AJB edits (v3-3): Kavin's SO4 factors (x2, +3) given names, so that the spread label
% below can find the top edge of the SO4 band
so4Multiplier = 2;
so4Offset = 3;
plotFitAndBinMeans(KAVIN_TABLE.UCVF, KAVIN_TABLE.SO4, xFit5, yFit5, so4Multiplier, so4Offset, colorsA(2,:), style)     % SO4

% AJB edits (v3-2): MO-OOA offset depends on the months used.
%   April-June only (moAprJunOnly = true, the default): Kavin's offset, 2.5.
%   All months (moAprJunOnly = false): offset 0, as in v3 and v3-1. Over all months the
%   MO-OOA values are on average about 3 times larger than in April-June 2023 (mean 0.77
%   vs. 0.27, in the units of KAVIN_TABLE), and with offset 2.5 the MO-OOA line sat on
%   top of the AOS line and ran off the top of the panel. Offset 0 puts it between SO4
%   and SR-OOA, where Kavin's April-June line was.
%   The offset only moves the line and markers up or down. It does not change the slope
%   or r, and the y axis has no numbers ("Arbitrary Units").
if moAprJunOnly
    moOffset = 2.5;
else
    moOffset = 0;
end
% AJB edits (v3-2): MO-OOA uses moRows and moOffset, SR-OOA uses srRows (v3-1: ooaRows)
plotFitAndBinMeans(KAVIN_TABLE.UCVF(moRows), KAVIN_TABLE.("MO-OOA")(moRows), xFit3, yFit3, 5, moOffset, colorsA(3,:), style)   % MO-OOA

plotFitAndBinMeans(KAVIN_TABLE.UCVF(srRows), KAVIN_TABLE.("SR-OOA")(srRows), xFit6, yFit6, 7, 1.8, colorsA(4,:), style)       % SR-OOA
plotFitAndBinMeans(KAVIN_TABLE.UCVF, KAVIN_TABLE.f44, xFit8, yFit8, 25, -3, colorsA(5,:), style)                          % f44
% AJB edits (v3-4): FThickly binned from its individual samples like the other species
% (v3-2 and v3-3: Kavin's pre-averaged points with plotFitAndOwnPoints).
% FThickly factors: Kavin used x15 - 4 for his pre-averaged April points. With the
% individual samples, the IQR of FThickly is 0.1 to 0.3 wide, so x15 made its band 1.5 to
% 4.5 units tall: it covered most of panel A and ran below the bottom. x8 - 2.4 (Andrew's
% choice from three options) keeps the FThickly band (about -1.1 to 2.0) at the bottom of
% the panel, mostly below the f44 band. Like Kavin's factors, this only changes where
% FThickly is drawn, not its slope sign or r.
fthickMultiplier = 8;
fthickOffset = -2.4;
plotFitAndBinMeans(fthickUCVF, fthick, xFit11, yFit11, fthickMultiplier, fthickOffset, colorsA(6,:), style)   % FThickly
% AJB edits (v3-3): label that explains the spread, with a short pointer line from the top
% edge of the uppermost IQR band in panel A (SO4; the AOS band above it is mean +/- std,
% see the legend) at UCVF = spreadLabelUCVF, pointing up and
% slightly to the left (north-north-west). The top edge is the 75th percentile of SO4
% in the bins on either side of that UCVF (linear interpolation, as the band is drawn),
% scaled like the SO4 data. The label text follows style.spread.
if style.spread ~= "none"
    [binCentreSO4, ~, ~, ~, ~, q75SO4] = binMeansOf(KAVIN_TABLE.UCVF, KAVIN_TABLE.SO4, style.binEdges);
    hasQ75 = ~isnan(q75SO4);
    pointerStart = [spreadLabelUCVF, interp1(binCentreSO4(hasQ75), q75SO4(hasQ75), spreadLabelUCVF)*so4Multiplier + so4Offset];
    pointerEnd = pointerStart + [-0.03, 1.2];     % up and slightly left, in data units
    plot([pointerStart(1) pointerEnd(1)], [pointerStart(2) pointerEnd(2)], '-', 'Color', [0.2 0.2 0.2], ...
        'LineWidth', 0.5, 'HandleVisibility', 'off')
    switch style.spread
        case "iqr", labelText = "IQR (25th - 75th pct.)";
        case "std", labelText = "\pm1 standard deviation";
        case "sem", labelText = "\pm1 standard error";
    end
    % label centred just above the end of the pointer line
    text(pointerEnd(1), pointerEnd(2), labelText, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', [0.2 0.2 0.2])   % AJB edits (v3-3, by Andrew): 12 pt
end
% AJB note: no3index is not used
no3index = KAVIN_TABLE.NO3.*95 -29 < 6
hold off
% AJB edits (v3-2): top of the panel raised to 8.4 (Kavin: 7; v3 and v3-1: 7.2), so
% that the legend in the top-left corner sits above the AOS line and points
% AJB edits (v3-4): bottom lowered to -1.3 (v3-3: 0) to make room for the FThickly band
ylim([-1.3 8.4])
ylabel("Arbitrary Units")
yticklabels("")
% AJB edits (v3-2): month labels moLabel and srLabel (v2 to v3-1: ooaLabel for both).
% With fitToBinMeans, the r of every entry except AOS and FThickly is r_bin (section 2c).
% AJB edits (v3-4): the AOS series is named "rBC mass per particle (AOS)" (Kavin: "AOS"),
% and the AOS and FThickly entries give the period (aosText, aosLabel). What the AOS band
% shows (mean +/- 1 std) is in the note at the bottom of the panel.
lgnd = legend("rBC mass per particle (AOS), "+aosText+", "+ string7,"SO4, "+string5,"MO-OOA"+moLabel+", "+string3,"SR-OOA"+srLabel+", " + string6,"f44, "+string8,"FThickly"+aosLabel+", "+string11)
% AJB edits (v3-2): larger legend text (Kavin: lgnd.FontSize = 5)
lgnd.FontSize = style.legendFontSize_pt;
% AJB edits (v3-2): legend in the top-left corner, where there is no data (v1 to v3-1:
% bottom right, where it is in Kavin's figure). It is moved a little to the right of
% the panel letter "A" at the end of section 3.
lgnd.Location = 'northwest';
% AJB edits (v3-2): panel letter placed in axes units, so that it stays in the top-left
% corner after the y range change (Kavin: text(0.03, 0.97*7,'A'))
% AJB edits (v3-3, by Andrew): panel label at the top left, pnl_fnt points (v3-2: 'A').
% AJB edits (v3-4): Lynn's request: a plain capital letter, not bold (v3-3: "(a)" in bold)
text(0.015, 0.95, "A", 'FontSize', pnl_fnt,...
    'Units', 'normalized')   % same place as the labels (b) and (c)
% AJB edits (v3-3, by Andrew): note on the bands. AJB edits (v3-4): reworded, because
% FThickly now has an IQR band (v3-3 text: "* No IQR bands for AOS and FThickly data.
% PIs provided binned data, rather than raw.")
% Right-aligned at the right edge of the panel, so that the text cannot run past it.
text(0.99, 0.05, "* AOS band: mean \pm 1 std (PIs provided binned data, rather than raw). All other bands: IQR.", 'FontSize', 7,...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'Color', [0.5020, 0.5020, 0.5020])
xticklabels([""])
xticks([0:0.1:1])

% ---------- Panel B ----------
% AJB edits (v2): HOA and rBC are not multiplied by 2 (fix 3)
nexttile(17,[2,3])
hold on
plotFitAndBinMeans(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.HOAMass_ugM_3_, xFit1, yFit1, 1, 0, colorsB(1,:), style)   % HOA
plotFitAndBinMeans(LANLdata_HOA_rBC.UCVF, LANLdata_HOA_rBC.rBCMass_ugM_3_, xFit2, yFit2, 1, 0, colorsB(2,:), style)   % rBC
% AJB edits (v3-1): ultrafine salt multiplied by saltMultiplier (v3: 1), so that its slope
% can be seen. Its r and legend are unchanged.
plotFitAndBinMeans(ultrafinetbl.avgUCVF, ultrafinetbl.avgUfineSalt, xFit14, yFit14, saltMultiplier, 0, colorsB(3,:), style)   % ultrafine salt
hold off
ylabel("Arbitrary Units")
yticklabels("")
% AJB edits (v2): names of entries 1 and 2 swapped back, and entry 3 prints string14,
% the ultrafine salt r, instead of string3, the MO-OOA r (fixes 1 and 2). Kavin's line:
% legend("rBC, "+ string1,"HOA, "+ string2,"NaCl2, "+ string3, FontSize=5)
% AJB edits (v3-2): larger legend text (Kavin: FontSize=5)
legend("HOA, "+ string1,"rBC, "+ string2,"NaCl2, "+ string14, FontSize=style.legendFontSize_pt)
% AJB edits: legend in the top-right corner, where it is in Kavin's figure
axB = gca;
axB.Legend.Location = 'northeast';
% AJB edits (v3-3, by Andrew): panel label at the top left, pnl_fnt points (Kavin: 'B').
% AJB edits (v3-4): Lynn's request: a plain capital letter, not bold (v3-3: "(b)" in bold)
text(0.015, 0.95, "B", 'FontSize', pnl_fnt,...
    'Units', 'normalized')   % same place as the labels (a) and (c)
xticklabels("")
ylim([0 0.6])
xticks([0:0.1:1])

% ---------- Panel C ----------
% OH multiplied by 20 and dD shifted by +85, as in Kavin's figure.
% AJB edits (v3-3): OH also shifted up by 10 (Kavin: no shift), so that it does not sit
% on top of the Pier dD; the Pier dD is added with the same +85 shift as the Soledad dD,
% so the two sites can be compared directly.
nexttile(25,[2,3])
hold on
plotFitAndBinMeans(OHvsUCVF.UCVF24, OHvsUCVF.OH_burst_25000_vv_dilution_ratio_uM, xFit10, yFit10, 20, 10, colorsC(1,:), style)  % OH
plotFitAndBinMeans(julyUCVF, julydD, xFit13, yFit13, 1, 85, colorsC(2,:), style)                                                 % Soledad dD
plotFitAndBinMeans(pierUCVF, pierdD, xFit15, yFit15, 1, 85, colorsC(3,:), style)                                                 % Pier dD (AJB edits (v3-3))
hold off
ylabel("Arbitrary Units")
yticklabels("")
% AJB edits (v3-3): the legend names the site of each dD record (Kavin: "7/2023 dD").
% "\delta" is shown as the Greek letter delta. Larger legend text (v3-2; Kavin: FontSize=5).
legend("OH, "+ string10, "\deltaD, Mt. Soledad (07/2023), "+string13, "\deltaD, Pier (07/2023), "+string15, ...   % AJB edits (v3-3, by Andrew): comma after \deltaD
    FontSize=style.legendFontSize_pt)
% AJB edits (v3-3): panel letter in axes units (Kavin: text(0.03, 0.95*40-20,'C')), and a
% y range with room for OH's new shift, the spread bands and the legend in the bottom-left
% corner (Kavin: ylim([-20 20]))
% AJB edits (v3-3, by Andrew): panel label at the top left, pnl_fnt points (Kavin: 'C').
% AJB edits (v3-4): Lynn's request: a plain capital letter, not bold (v3-3: "(c)" in bold)
text(0.015, 0.95, "C", 'FontSize', pnl_fnt,...
    'Units', 'normalized')   % same place as the labels (a) and (b)
ylim([-24 23])
xticks([0:0.1:1])
xticklabels([0:0.1:1])
% AJB edits: x-axis label (not in Kavin's figure)
xlabel("UCVF")

% AJB edits (v3-3): panel C legend in the bottom-left corner (v1 to v3-2: a position
% measured on Kavin's figure, a little higher). The dD records start at UCVF = 0.29, so
% this corner is empty below the dD spread bands.
axC = gca;
axC.Legend.Location = 'southwest';
drawnow

% AJB edits (v3-2): panel A legend: keep it in the top-left corner, but start it 6 % of
% the axes width from the left axis, so that it does not cover the panel letter "A".
% (The axes and its legend are both children of the tiled layout, so their Position
% values share one coordinate frame; drawnow above gives the final axes positions.)
lgnd.Position(1) = axA.Position(1) + 0.06*axA.Position(3);

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
%   1. its fit line (xFit, yFit from section 2, or from section 2c in v3-2);
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

% 2. Bin means (AJB edits (v3-2): computed by binMeansOf, shared with fitBinMeans)
[binCentre, binMean, binStd, binCount, binQ25, binQ75, binMedian] = binMeansOf(ucvf, values, style.binEdges);
% AJB edits (v3-3): the value drawn for each bin: median or mean (style.binStatistic)
if style.binStatistic == "median"
    binValue = binMedian;
else
    binValue = binMean;
end

% AJB edits (v3-3): 3. Spread of the samples in each bin (style.spread), as a band or as
% bars (style.spreadDisplay). lower and upper are the limits in data units; they are
% then multiplied and shifted like the data. Bins with fewer than 2 samples have NaN
% limits and show no spread.
switch style.spread
    case "iqr"                                   % 25th to 75th percentile
        lower = binQ25;
        upper = binQ75;
    case "std"                                   % mean +/- standard deviation
        lower = binMean - binStd;
        upper = binMean + binStd;
    case "sem"                                   % mean +/- standard error of the mean
        lower = binMean - binStd./sqrt(binCount);
        upper = binMean + binStd./sqrt(binCount);
    otherwise                                    % "none"
        lower = NaN(size(binMean));
        upper = NaN(size(binMean));
end
lowerPlot = min(lower.*multiplier, upper.*multiplier) + offset;   % min/max keep the order
upperPlot = max(lower.*multiplier, upper.*multiplier) + offset;   % if a multiplier is negative
if style.spreadDisplay == "band"
    drawSpreadBand(binCentre, lowerPlot, upperPlot, color, style.bandAlpha)
else
    drawSpreadBars(binCentre, lowerPlot, upperPlot, color, style)
end

% Bin markers (median or mean, see above), kept out of the legend
plot(binCentre, binValue.*multiplier + offset, 'o', 'Color', color, ...
    'MarkerFaceColor', color, 'MarkerSize', style.markerSize_pt, ...
    'HandleVisibility', 'off')
end


function drawSpreadBand(x, lowerY, upperY, color, alphaValue)
% AJB edits (v3-3): new. Light, see-through band between lowerY and upperY across the
% bins (x). A band is drawn for each run of neighbouring bins that all have a spread;
% bins without one (empty, or a single sample) interrupt the band. The band is kept out
% of the legend and moved behind the lines and markers (uistack).
hasSpread = ~isnan(lowerY) & ~isnan(upperY);
runStart = find(hasSpread & ~[false; hasSpread(1:end-1)]);   % first bin of each run
runEnd   = find(hasSpread & ~[hasSpread(2:end); false]);     % last bin of each run
for k = 1:numel(runStart)
    idx = runStart(k):runEnd(k);
    band = fill([x(idx); flipud(x(idx))], [lowerY(idx); flipud(upperY(idx))], color, ...
        'FaceAlpha', alphaValue, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(band, 'bottom')
end
end


function drawSpreadBars(x, lowerY, upperY, color, style)
% AJB edits (v3-3): new. Thin vertical bar from lowerY to upperY at each bin (x), in a
% lighter version of the species colour, without end caps. All bars of one species are
% one line object: NaN between the bars breaks the line.
lightColor = color + style.barLightening .* (1 - color);    % mix the colour with white
n = numel(x);
xs = reshape([x(:)'; x(:)'; NaN(1, n)], [], 1);
ys = reshape([lowerY(:)'; upperY(:)'; NaN(1, n)], [], 1);
plot(xs, ys, '-', 'Color', lightColor, 'LineWidth', style.barLineWidth_pt, 'HandleVisibility', 'off')
end


function plotFitAndOwnPoints(ucvf, values, xFit, yFit, multiplier, offset, color, style)
% AJB edits (v3-2): new. Like plotFitAndBinMeans, but draws the species' own data points
% instead of bin means. Used for AOS and FThickly, whose owners already averaged them.
% Both line and points are multiplied by MULTIPLIER and shifted by OFFSET.

% 1. Fit line (with the circle marker for the legend; MarkerIndices = [] draws no
%    circles on the line itself)
plot(xFit, yFit.*multiplier + offset, '-o', 'Color', color, ...
    'LineWidth', style.lineWidth_pt, 'MarkerFaceColor', color, ...
    'MarkerSize', style.markerSize_pt, 'MarkerIndices', [])

% 2. Own data points, same marker as the bin means elsewhere, kept out of the legend
plot(ucvf, values.*multiplier + offset, 'o', 'Color', color, ...
    'MarkerFaceColor', color, 'MarkerSize', style.markerSize_pt, ...
    'LineStyle', 'none', 'HandleVisibility', 'off')
end


function [binCentre, binMean, binStd, binCount, binQ25, binQ75, binMedian] = binMeansOf(ucvf, values, edges)
% AJB edits (v3-2): new (moved out of plotFitAndBinMeans so that the markers and the
% fits of section 2c use exactly the same bin means).
% Mean of VALUES in each UCVF bin; EDGES are the bin edges (15 edges = 14 bins).
% Only rows where both UCVF and the value exist are used, as in the linear fits.
% binCentre: centre of each bin; binMean: mean in each bin (NaN for an empty bin).
% AJB edits (v3-3): also returns binStd, the standard deviation of the samples in each
% bin (NaN for bins with fewer than 2 samples, where it is not defined), and binCount,
% the number of samples in each bin; binQ25 and binQ75, the 25th and 75th percentiles of
% the samples in each bin (NaN for bins with fewer than 2 samples); binMedian, the median
% of the samples in each bin (NaN for an empty bin).
ok = ~isnan(ucvf) & ~isnan(values);
ucvfOk = ucvf(ok);
valuesOk = values(ok);
binNumber = discretize(ucvfOk, edges);        % which bin each sample falls in
nBinsTotal = numel(edges) - 1;
binMean = NaN(nBinsTotal, 1);
binStd = NaN(nBinsTotal, 1);                  % AJB edits (v3-3)
binCount = zeros(nBinsTotal, 1);              % AJB edits (v3-3)
binQ25 = NaN(nBinsTotal, 1);                  % AJB edits (v3-3)
binQ75 = NaN(nBinsTotal, 1);                  % AJB edits (v3-3)
binMedian = NaN(nBinsTotal, 1);               % AJB edits (v3-3)
for b = 1:nBinsTotal
    inBin = (binNumber == b);
    binCount(b) = sum(inBin);
    if binCount(b) >= 1
        binMean(b) = mean(valuesOk(inBin));
        binMedian(b) = median(valuesOk(inBin));   % AJB edits (v3-3)
    end
    if binCount(b) >= 2                       % AJB edits (v3-3): std needs 2 or more samples
        binStd(b) = std(valuesOk(inBin));     % sample standard deviation (divides by n - 1)
        q = prctile(valuesOk(inBin), [25 75]);  % AJB edits (v3-3): quartiles of the bin
        binQ25(b) = q(1);
        binQ75(b) = q(2);
    end
end
binCentre = (edges(1:end-1) + edges(2:end))' / 2;
end


function [xFit, yFit, rText, rBins, rSamples, nBins] = fitBinMeans(ucvf, values, edges, statistic)
% AJB edits (v3-2): new. Averages VALUES in the UCVF bins (EDGES), then fits a straight
% line to the bin means, in the same style as Kavin's fits in section 2:
%   xFit, yFit : the fit line on 100 points between the first and last bin centre
%   rText      : legend text, e.g. "r_{bin} = 0.95" (shown as r with subscript "bin")
%   rBins      : r between bin centres and bin means
%   rSamples   : r from every sample (for comparison)
%   nBins      : number of bins that hold samples (points in the fit)
% AJB edits (v3-3): STATISTIC ("median" or "mean") picks the value that represents each
% bin; the default in the User settings is the median, the same value as the markers.
[binCentre, binMean, ~, ~, ~, ~, binMedian] = binMeansOf(ucvf, values, edges);
if statistic == "median"
    binValue = binMedian;
else
    binValue = binMean;
end
hasData = ~isnan(binValue);                    % bins that hold at least one sample
nBins = sum(hasData);
model = fitlm(binCentre(hasData), binValue(hasData));
xFit = linspace(min(binCentre(hasData)), max(binCentre(hasData)), 100);
yFit = feval(model, xFit);
r = corrcoef(binCentre(hasData), binValue(hasData));
rBins = r(2,1);
rText = "r_{bin} = " + string(round(rBins, 2));
r = corrcoef(ucvf, values, 'rows', 'pairwise');
rSamples = r(2,1);
end
