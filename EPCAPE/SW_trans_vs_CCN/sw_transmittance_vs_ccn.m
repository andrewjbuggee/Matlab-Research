%% EPCAPE: shortwave transmittance at the Scripps Pier vs. Pier CCN, grouped by cloud supersaturation
%
% DESCRIPTION:
%   Scatter plot of the surface SW transmittance (y) against the Pier CCN
%   concentration (x), averaged over consecutive averageWindow_min windows
%   (e.g. 10 min or 1 h), with the points colored by the
%   cloud's effective supersaturation, optionally in nSupersaturationBins
%   groups. The question behind it: at a given aerosol loading, does the
%   cloud let less sunlight through when more of the aerosol activates
%   (higher supersaturation -> more droplets -> brighter cloud, Twomey
%   1977)?
%
%   WHAT IS PLOTTED, AND FROM WHERE
%   y  T_SW = <SW down> / <clear-sky SW down> over the window, from ARM
%      RADFLUX at the Pier (Long & Ackerman 2000), using only the minutes
%      that pass every filter below (EPCAPE/SW_transmittance).
%   x  CCN at the CCN counter's ccnSetpoint_pct supersaturation, from the
%      Pier DMT CCN-200, corrected by Dedrick et al. (2026) for an inlet leak
%      (x 1.57; supplement Text S4): the CCN data Dedrick et al. used. With
%      ccnAtRetrievedS = true, x is instead the CCN spectrum evaluated at the
%      hour's retrieved supersaturation, i.e. the number of particles that
%      should have activated (for s_Nd this equals the FM-120 droplet
%      number by construction, N_d = CCN(s_Nd)).
%      WHY THE PIER: s_kappa, s_CCN, s_Nd,FM and s_Nd,VAP all use Scripps Pier
%      aerosol data (Dedrick et al. 2026, Sec. 2d; supplement Table S7). The
%      Mt. Soledad CCN-200 samples behind the GCVI during cloud events, i.e.
%      it measures droplet RESIDUALS (particles that already activated), was
%      used only for the extra s_Dact method (Text S10), starts at 0.2 % and
%      cycles in 2 h, and is not among the posted products on disk.
%   color  the effective cloud supersaturation s (%), either Dedrick's
%      posted hourly retrievals ("posted"; method and period chosen below)
%      or my own re-derivation of the cloud-scene method s_Nd,FM
%      ("AJB"; EPCAPE/supersaturation, not yet independently verified).
%      x and color are different supersaturations: ccnSetpoint_pct is the
%      chamber setting of the CCN counter (a fixed probe of the aerosol: how
%      many particles WOULD activate at that s), while the color is the s the
%      cloud actually reached, set by updraft and the droplets' competition
%      for vapor (Twomey 1959).
%
%   CASES AND FILTERS (set in Section 1)
%   caseSet "dedrick24" keeps the in-cloud windows of the 24 Mt. Soledad
%   cloud events of Dedrick et al. (2026, supplement Table S5). Then every
%   minute must pass the filters switched on below: daytime, LWP between
%   the bounds, no rain, single layer, liquid only (height proxy), and a
%   deck that persists at the Pier and at Mt. Soledad. What each filter
%   tests, and why, is in functions/analysis_filters_1min.m. The day is
%   then cut into averaging windows of averageWindow_min minutes; a window
%   enters the plot if >= minPassFraction of its minutes pass, and it has a
%   CCN value and a supersaturation. T_SW and LWP are averaged over the
%   passing minutes only. The CCN and s are hourly products, so for
%   windows shorter than an hour every window inside one hour carries that
%   hour's CCN and s (functions/hourly_values_on_windows.m).
%
%   WHY THE NUMBER OF POINTS IS SMALL (see the data-loss script in
%   EPCAPE/Dedrick_2026_coupled_cases)
%   * Mt. Soledad is in cloud mostly at night and in the early morning, so
%     the daytime cut removes most of each event. We need daytime data if
%     we wish to investigate SW transmittance with CCN.
%   * The Pier CCN counter came online on 23 May 2023: events 1-7 (late
%     April to 18 May) have no Pier CCN, and no s_CCN or s_Nd.
%   * The posted "SoledadCloudEvent" supersaturations are not the 24
%     events: 278 h from 12 Mar 2023 to 21 Jan 2024, of which ~141 h fall
%     inside the 24 windows, and none in events 23 and 24 (October).
%
%   CAUTION WHEN READING THE PLOT:
%   T_SW depends far more on liquid water path than on droplet number, so a
%   CCN or supersaturation signal is easily hidden by LWP variability; the
%   second panel (showLwpPanel) shows T_SW against LWP with the same
%   colors. With a few dozen hours from a handful of events, the points are
%   not independent (and short windows of one hour share one CCN and one
%   s), so only strong relationships can be detected.
%
% OUTPUTS (in ./output):
%   sw_transmittance_vs_ccn_<tag>.png/.pdf  the figure
%   sw_transmittance_vs_ccn_<tag>.csv       every plotted window
%   (<tag> = case set, s source and method, averaging window)
%   epcape_1min_inputs.mat                  cache of the 1-min inputs
%                                           (rebuildCache = true redoes it)
%
% DATA (Python project data/processed, and Dedrick's posted products):
%   ARM EPCAPE M1: RADFLUX (radflux1long c2), MWRRET v2 (mwrret2turn c1),
%     ceilometer (ceil b1), KAZR-ARSCL (arsclkazr1kollias c1), laser
%     disdrometer (ld b1); ARM S2 (Mt. Soledad): laser disdrometer. 
%   Mt. Soledad visibility (UCSD, msd_visibility_*.csv) and FM-120 fog
%     monitor 5-min files (Dalhousie: L. Robinson and R. Chang).
%   Berta, V., S. Han, A. Osawa, J. Dedrick, N. Maneenoi, and L. M. Russell,
%     2025: Aerosol quantification products at Mt. Soledad and Scripps Pier
%     during EPCAPE 2023-24. UC San Diego Library Digital Collections,
%     doi:10.6075/J0WM1DSC (corrected Pier CCN, supersaturations, low-cloud
%     periods).
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: Quantitative estimates of upwind
%     and on-site particle scavenging and retrieved cloud supersaturations
%     from urban coastal stratocumulus during EPCAPE. J. Atmos. Sci., 83,
%     1061-1082, doi:10.1175/JAS-D-26-0004.1, and supplement
%     doi:10.1175/JAS-D-26-0004.s1.
%   Long, C. N., and T. P. Ackerman, 2000: J. Geophys. Res., 105,
%     15609-15626, doi:10.1029/2000JD900077.
%   Twomey, S., 1959: The nuclei of natural cloud formation. Part II: The
%     supersaturation in natural clouds and the variation of cloud droplet
%     concentration. Geofis. Pura Appl., 43, 243-249.
%   Twomey, S., 1977: The influence of pollution on the shortwave albedo of
%     clouds. J. Atmos. Sci., 34, 1149-1152.
%   Rogers, R. R., and M. K. Yau, 1989: A Short Course in Cloud Physics,
%     3rd ed., Pergamon (CCN activation, pp. 95-96, 107-110).
%
% By Andrew John Buggee
%%

clear variables


%% 1. Settings: everything you may want to change is here

% --- which cloud cases --------------------------------------------------------
caseSet            = "dedrick24";   % "dedrick24" | "soledadInCloud" | "lowCloudPeriods" | "all"
cloudIntervalsOnly = false;          % dedrick24 only: keep just the 5-min intervals of each window that
                                     % meet the cloud condition (vis < 1 km, LWC > 0.01 g m-3)

% --- filters (true = apply) ---------------------------------------------------
applyDaytime     = true;  
szaMax_deg       = 80;      % sun higher than 10 deg

applyLWP         = true;   
lwpMin_gm2       = 20;      % LWP > 20 g m-2 ...
lwpMax_gm2       = 300;     % ... and < 300 g m-2 (MWRRET v2 at the Pier) - value used by Dedrick et al. (2026)

applyNoPrecip    = true;   
rainWindow_min   = 15;      % no rain within +/- rainWindow_min minutes ...
rainSites        = ["Pier" "Soledad"];   % ... at these disdrometers
   
requireArscl     = false;   % true: a minute without ARSCL data fails

% cloud filters
applySingleLayer       = true;
applyLiquid            = true;   
maxCloudTop_m          = 3000;    % liquid proxy: top and base below 3 km
maxCloudBase_m         = 3000;

% cloud persistance: tested minute by minute. A minute passes if, over the
% persistence window [minute - persistHalfWindow_min, minute + persistHalfWindow_min],
% cloud is detected at every site in persistenceSites for >= persistMinFraction
% of the minutes that have cloud-detection data, AND those minutes with data
% make up >= persistMinCoverage of the window (so a minute is not called
% "persistent" on the strength of 2 minutes of data in a 21-min window)
applyPersistence       = true;
persistenceSites       = ["Pier" "Soledad"];
persistHalfWindow_min  = 10;     % persistence window = minute +/- persistHalfWindow_min
persistMinFraction     = 0.9;    % cloud for >= 90 % of the window minutes with data ...
persistMinCoverage     = 0.5;    % ... and data for >= 50 % of the window


% --- averaging window: one plotted point per window ---------------------------------
% After every filter has been applied minute by minute, T_SW (ratio of the
% means), LWP, CCN and s are averaged over the PASSING minutes of
% consecutive windows [t, t + averageWindow_min), anchored at midnight UTC
% (10 -> :00, :10, :20, ...). This is independent of the persistence window:
% with averageWindow_min = 10, an hour in which every minute passes gives 6
% points. averageWindow_min = 60 gives the clock hours.
% NOTE: the Pier CCN and the supersaturations are hourly (the CCN counter
% measures each setpoint ~once per 1-h cycle), so every window inside one
% hour carries that hour's CCN and s: only T_SW and LWP vary within an hour.
% Use windows that divide the hour (10, 15, 20, 30, 60) or multiples of it.
averageWindow_min = 10;    % minutes
minPassFraction   = 0.5;   % a window is plotted if >= this fraction of its minutes pass every filter
                           % (unlike persistMinCoverage, which is about data inside the
                           % persistence test, this counts minutes passing ALL filters)


% --- x axis: CCN -----------------------------------------------------------------
ccnSetpoint_pct = 0.2;     % 0.1 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 (CCN-counter supersaturation, %)
ccnAtRetrievedS = false;   % true: x = CCN spectrum at the hour's retrieved s (activated CCN)
ccnLogAxis      = false;


% --- color: cloud supersaturation ------------------------------------------------
supersaturationSource = "posted";   % "posted" (Dedrick et al. 2026) | "AJB" (my s_Nd,FM; unverified)
supersaturationMethod = "sNd";      % posted: "sNd" | "sCCN" | "skappa" | "sARG00"
postedPeriods = ["SoledadCloudEvent" "LowCloudPeriod"];   % posted columns used, in order of preference
binBySupersaturation  = true;       % false: continuous color scale with a color bar
nSupersaturationBins  = 5;          % 5, 6, 7, ...
binEdgeRule           = "quantile"; % "quantile" (equal numbers of windows) | "equalWidth"

% --- figure and run control ---------------------------------------------------------
showLwpPanel = true;      % second panel: T_SW vs LWP, same colors
rebuildCache = false;     % true: re-read every instrument (~1-2 min) instead of the cache

% --- folders --------------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
epcapeDir = fileparts(scriptDir);
outDir    = fullfile(scriptDir, 'output');
processed = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed";
productsFolder = fullfile(processed, "Aerosol Quantification Products");

files.radflux      = fullfile(processed, "radflux_M1_20230215_20240214.nc");
files.mwrret2      = fullfile(processed, "lwp_mwrret2_M1_20230215_20240214.nc");
files.ceil         = fullfile(processed, "cbh_ceil_M1_20230215_20240214.nc");
files.arscl        = fullfile(processed, "cloud_arscl_M1_20230215_20240214.nc");
files.rain_pier    = fullfile(processed, "rain_ld_M1_20230215_20240214.nc");
files.rain_soledad = fullfile(processed, "rain_ld_S2_20230215_20240214.nc");
files.visibility   = [fullfile(processed, "msd_visibility_20230224-20230719.csv"), ...
                      fullfile(processed, "msd_visibility_20230719-20240214.csv")];
files.fm120_folder = fullfile(processed, "FM120", "Daily Averaged FM-120 FIles");

addpath(fullfile(scriptDir, 'functions'));
addpath(fullfile(epcapeDir, 'SW_transmittance'));
addpath(fullfile(epcapeDir, 'Dedrick_2026_coupled_cases'));
addpath(fullfile(epcapeDir, 'supersaturation', 'functions'));   % CCN spectrum + s retrieval + new_figure
if ~isfolder(outDir), mkdir(outDir); end

% collect the filter settings for analysis_filters_1min.m
F = struct('applyDaytime', applyDaytime, 'szaMax_deg', szaMax_deg, ...
           'applyLwp', applyLWP, 'lwpMin_gm2', lwpMin_gm2, 'lwpMax_gm2', lwpMax_gm2, ...
           'applyNoPrecip', applyNoPrecip, 'rainWindow_min', rainWindow_min, 'rainSites', rainSites, ...
           'applySingleLayer', applySingleLayer, 'requireArscl', requireArscl, ...
           'applyLiquid', applyLiquid, 'maxCloudTop_m', maxCloudTop_m, 'maxCloudBase_m', maxCloudBase_m, ...
           'applyPersistence', applyPersistence, 'persistenceSites', persistenceSites, ...
           'persistHalfWindow_min', persistHalfWindow_min, 'persistMinFraction', persistMinFraction, ...
           'persistMinCoverage', persistMinCoverage);

% fewest passing minutes for a window to be plotted
minPassMinutes = max(1, ceil(minPassFraction * averageWindow_min));

if caseSet == "lowCloudPeriods" && applyPersistence && any(lower(persistenceSites) == "soledad")
    warning(['Low-cloud periods are times when the cloud does NOT envelop Mt. Soledad, but the ' ...
             'Soledad persistence test needs Mt. Soledad in cloud: almost nothing will pass. ' ...
             'Use persistenceSites = "Pier" for this case set.']);
end



%% 2. The 1-min inputs (from the cache if it exists)

cacheFile = fullfile(outDir, 'epcape_1min_inputs.mat');
if isfile(cacheFile) && ~rebuildCache
    fprintf('Loading the 1-min inputs from %s\n', cacheFile);
    load(cacheFile, 'D', 'cloud5', 'fm');
else
    fprintf('Building the 1-min inputs (every instrument, ~1-2 min) ...\n');
    [D, cloud5, fm] = build_epcape_1min_inputs(files);
    save(cacheFile, 'D', 'cloud5', 'fm', 'files', '-v7.3');
end



%% 3. Hourly CCN spectra and supersaturations

[H, s_set_pct] = load_hourly_ccn_and_supersaturation(productsFolder, fm);



%% 4. Cases and filters, minute by minute

[in_case, event_id_1min] = case_set_mask_1min(D, caseSet, cloud5, H, cloudIntervalsOnly);
M = analysis_filters_1min(D, F);
pass = in_case & M.all;



%% 5. Window averages over the passing minutes, joined to CCN and supersaturation

% T_SW, LWP and SZA of each averaging window, from its passing minutes only
Hr = average_passing_minutes(D, pass, averageWindow_min, minPassMinutes, szaMax_deg);

% Dedrick event number of each window (the event of its passing minutes)
ev = event_id_1min;  ev(~pass) = NaN;
[~, ev_win] = window_mean(D.time, ev, averageWindow_min);   % all passing minutes of a window share one event
Hr.event_id = round(ev_win);

% the clock hour that holds each window's start: the CCN and s of a window
% shorter than an hour come from this hour
Hr.hour_utc = dateshift(Hr.time, 'start', 'hour');

% --- the supersaturation used for color (hourly, on the hourly grid H) ---
switch supersaturationSource
    case "posted"
        s_pct = nan(height(H), 1);
        for period = postedPeriods                   % first period with a value wins
            col = supersaturationMethod + "_" + period;
            fill_in = isnan(s_pct) & isfinite(H.(col));
            s_pct(fill_in) = H.(col)(fill_in);
        end
        sLabel = sprintf('Dedrick et al. (2026) posted s_{%s}', erase(supersaturationMethod, "s"));
    case "AJB"
        s_pct = H.sNdFM_AJB;
        sLabel = 's_{Nd,FM} (my re-derivation, unverified)';
    otherwise
        error('supersaturationSource must be "posted" or "AJB"');
end

% --- the CCN for the x axis (hourly, on the hourly grid H) ---
if ccnAtRetrievedS
    ccn = nan(height(H), 1);
    for i = find(isfinite(s_pct) & all(isfinite(H.ccn_cm3), 2))'
        [s_ext, c_ext] = extend_ccn_spectrum(H.ccn_cm3(i, :));   % spectrum extended to 0 % and 1.2 %
        ccn(i) = interp1(s_ext, c_ext, s_pct(i));
    end
    xLabel = 'Pier CCN at the retrieved supersaturation (cm^{-3})';
else
    col = find(abs(s_set_pct - ccnSetpoint_pct) < 1e-9);
    assert(~isempty(col), 'ccnSetpoint_pct must be one of %s', mat2str(s_set_pct));
    ccn = H.ccn_cm3(:, col);
    xLabel = sprintf('Pier CCN at %.1f %% supersaturation (cm^{-3})', ccnSetpoint_pct);
end

% --- hourly s and CCN onto the averaging windows (mean over each window's
% passing minutes; a window inside one hour simply gets that hour's value) ---
[Hr.s_pct, t_chk]   = hourly_values_on_windows(D.time, pass, H.time, s_pct, averageWindow_min);
[Hr.ccn_cm3, t_chk2] = hourly_values_on_windows(D.time, pass, H.time, ccn, averageWindow_min);
assert(isequal(t_chk, Hr.time) && isequal(t_chk2, Hr.time), 'window grids differ');



%% 6. How many windows survive each step

% windows with >= minPassMinutes minutes passing the case and the filters
% applied so far (a running AND, in the order of M.applied); "hours" is the
% same count expressed as hours of data (windows x averageWindow_min / 60)
count_windows = @(mask) nnz(window_count(D.time, mask, averageWindow_min) >= minPassMinutes);
step_name = "in case set (" + caseSet + ")";
step_windows = count_windows(in_case);
running = in_case;
for name = M.applied
    running = running & M.(name);
    step_name(end + 1, 1) = "+ " + name + ": " + M.description.(name); %#ok<SAGROW>
    step_windows(end + 1, 1) = count_windows(running); %#ok<SAGROW>
end
has_T   = isfinite(Hr.T_sw);
has_ccn = has_T & isfinite(Hr.ccn_cm3);
has_s   = has_ccn & isfinite(Hr.s_pct);
step_name    = [step_name; "+ valid SW and clear-sky SW"; "+ Pier CCN"; "+ supersaturation (" + supersaturationSource + ")"];
step_windows = [step_windows; nnz(has_T); nnz(has_ccn); nnz(has_s)];
step_hours   = step_windows * averageWindow_min / 60;
steps = table(step_name, step_windows, step_hours, ...
              'VariableNames', {'step', sprintf('windows_%gmin', averageWindow_min), 'hours'});
disp(steps)

P = Hr(has_s, :);          % the plotted windows
n_hours_P = numel(unique(P.hour_utc));   % distinct clock hours = distinct CCN and s values (windows <= 60 min)
n_events = numel(unique(P.event_id(isfinite(P.event_id))));
fprintf('Plotted: %d %g-min averages from %d distinct hours', height(P), averageWindow_min, n_hours_P);
if caseSet == "dedrick24"
    fprintf(', %d of the 24 events (%s)', n_events, mat2str(unique(P.event_id(isfinite(P.event_id)))'));
end
fprintf('\n');



%% 7. Supersaturation bins

if height(P) == 0
    error('No windows pass every step; relax a filter (see the table above).');
end
if binBySupersaturation
    switch binEdgeRule
        case "quantile"
            edges = prctile(P.s_pct, linspace(0, 100, nSupersaturationBins + 1));
        case "equalWidth"
            edges = linspace(min(P.s_pct), max(P.s_pct), nSupersaturationBins + 1);
        otherwise
            error('binEdgeRule must be "quantile" or "equalWidth"');
    end
    edges = unique(edges);                  % ties (repeated values) can merge quantile edges
    if numel(edges) - 1 < nSupersaturationBins
        warning('Only %d distinct bins could be formed from %d windows.', numel(edges) - 1, height(P));
    end
    P.s_bin = discretize(P.s_pct, edges);   % the last bin includes its right edge
    n_bins = numel(edges) - 1;
    binColors = ordinal_blue_ramp(n_bins);
end



%% 8. Figure

c_ink  = [0.043 0.043 0.043];    % #0b0b0b, primary text
c_ink2 = [0.322 0.318 0.306];    % #52514e, secondary text

n_panels = 1 + showLwpPanel;
fig = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 6.5 * n_panels 6]);
tl = tiledlayout(fig, 1, n_panels, 'TileSpacing', 'compact', 'Padding', 'compact');

x_vars = {'ccn_cm3', 'lwp_gm2'};
x_labels = {xLabel, 'Pier LWP over the same minutes (g m^{-2}, MWRRET v2)'};
panel_titles = {'(a) SW transmittance vs. CCN', '(b) SW transmittance vs. LWP'};

r_panel = nan(n_panels, 1);   % Pearson r of each panel (linear values, also with ccnLogAxis)
n_panel = zeros(n_panels, 1); % number of points behind it
for ip = 1:n_panels
    ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
    x = P.(x_vars{ip});

    % Pearson correlation coefficient between x and T_SW over every point of
    % the panel ('complete' skips a point without LWP when the LWP filter is off)
    R = corrcoef(x, P.T_sw, 'Rows', 'complete');
    r_panel(ip) = R(1, 2);
    n_panel(ip) = nnz(isfinite(x) & isfinite(P.T_sw));
    if binBySupersaturation
        h = gobjects(n_bins, 1);
        for b = 1:n_bins
            in = P.s_bin == b;
            h(b) = scatter(ax, x(in), P.T_sw(in), 45, binColors(b, :), 'filled', ...
                           'MarkerEdgeColor', 'w', 'LineWidth', 0.75, ...
                           'DisplayName', sprintf('%.3f-%.3f %%  (n = %d)', edges(b), edges(b + 1), nnz(in)));
        end
        if ip == 1
            % the legend sits in its own tile right of the panels, clear of the data
            lg = legend(ax, h, 'Box', 'off', 'FontSize', 9);
            lg.Layout.Tile = 'east';
            title(lg, sLabel, 'FontWeight', 'normal', 'FontSize', 9);
        end
    else
        scatter(ax, x, P.T_sw, 45, P.s_pct, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.75);
        colormap(ax, ordinal_blue_ramp(256));
        clim(ax, [min(P.s_pct) max(P.s_pct)]);              % same color scale in both panels
        if ip == n_panels
            cb = colorbar(ax);
            cb.Layout.Tile = 'east';
            cb.Label.String = [sLabel ' (%)'];
        end
    end
    if ip == 1 && ccnLogAxis
        set(ax, 'XScale', 'log');
    end
    xlim(ax, 'padded');                                     % keep edge points off the axis lines
    xlabel(ax, x_labels{ip});
    ylabel(ax, 'SW transmittance, T_{SW} (Pier, RADFLUX)');
    ylim(ax, [0 ceil(10 * (max(P.T_sw) + 0.05)) / 10]);   % from 0 to just above the largest value
    % r and n go in the panel subtitle, where they can never cover a point
    title(ax, panel_titles{ip}, sprintf('Pearson r = %.2f  (n = %d)', r_panel(ip), n_panel(ip)), ...
          'FontWeight', 'normal');
    ax.Subtitle.Color = c_ink2;
    ax.XColor = c_ink;  ax.YColor = c_ink;
end

% what was plotted, in words, under the panels
case_words = struct('dedrick24', 'the 24 Mt. Soledad cloud events of Dedrick et al. (2026), in-cloud windows', ...
                    'soledadInCloud', 'every hour with Mt. Soledad in cloud (vis < 1 km, LWC > 0.01 g m^{-3})', ...
                    'lowCloudPeriods', 'Dedrick et al. (2026) low-cloud periods at the Pier', ...
                    'all', 'all hours');
filter_words = [string(cellfun(@(nm) M.description.(nm), cellstr(M.applied), 'UniformOutput', false)), ...
                sprintf("%g-min averages of the passing minutes (>= %d of %g min pass)", ...
                        averageWindow_min, minPassMinutes, averageWindow_min)];
if averageWindow_min < 60
    filter_words(end + 1) = "CCN and s are hourly: the windows of one hour share them";
end
title(tl, sprintf('EPCAPE, %d %g-min averages from %d hours: %s', height(P), averageWindow_min, ...
                  n_hours_P, case_words.(caseSet)), 'FontSize', 12, 'Color', c_ink);
subtitle(tl, wrap_words("Filters: " + strjoin(filter_words, '; '), 70 + 60 * showLwpPanel), ...
         'FontSize', 8, 'Color', c_ink2, 'Interpreter', 'none');

tag = sprintf('%s_%s_%s', caseSet, supersaturationSource, supersaturationMethod);
if supersaturationSource == "AJB", tag = sprintf('%s_AJB_sNdFM', caseSet); end
tag = sprintf('%s_%gmin', tag, averageWindow_min);   % one set of outputs per averaging window
exportgraphics(fig, fullfile(outDir, "sw_transmittance_vs_ccn_" + tag + ".png"), 'Resolution', 200);
exportgraphics(fig, fullfile(outDir, "sw_transmittance_vs_ccn_" + tag + ".pdf"), 'ContentType', 'vector');
close(fig);



%% 9. The plotted windows, and simple statistics

out = timetable2table(P);
out.Properties.VariableNames{1} = 'window_start_utc';
writetable(out, fullfile(outDir, "sw_transmittance_vs_ccn_" + tag + ".csv"));

fprintf('Pearson r(T_SW, CCN) over the %d plotted windows: %.2f\n', n_panel(1), r_panel(1));
if showLwpPanel
    fprintf('Pearson r(T_SW, LWP) over the same windows:      %.2f\n', r_panel(2));
end
fprintf(['These windows are not independent: several per event, and for windows < 60 min the windows of one\n' ...
         'hour share one CCN and one s value (%d distinct hours). The usual significance test overstates the evidence.\n'], ...
        n_hours_P);
fprintf('Saved the figure and table to %s (tag %s)\n', outDir, tag);




%% Local functions

function n = window_count(t, mask, window_min)
% number of true minutes in each window [t, t + window_min) of the window_mean.m grid
m = double(mask);  m(~mask) = NaN;           % window_mean counts the finite samples
[~, ~, n] = window_mean(t, m, window_min);
end


function lines = wrap_words(text, max_chars)
% break a long string into lines of at most ~max_chars characters, at spaces
words = split(string(text), " ");
lines = strings(0, 1);
current = "";
for w = words'
    if strlength(current) + strlength(w) + 1 > max_chars && strlength(current) > 0
        lines(end + 1, 1) = current; %#ok<AGROW>
        current = w;
    else
        current = strtrim(current + " " + w);
    end
end
lines(end + 1, 1) = current;
end


function C = ordinal_blue_ramp(n)
% n colors, light to dark, from the reference sequential blue ramp
% (dataviz reference palette, steps 250 -> 700; step 250 is the lightest
% that keeps >= 2:1 contrast with a white background for ordinal marks)
anchors_hex = ["86b6ef" "6da7ec" "5598e7" "3987e5" "2a78d6" "256abf" "1c5cab" "184f95" "104281" "0d366b"];
A = zeros(numel(anchors_hex), 3);
for k = 1:numel(anchors_hex)
    A(k, :) = hex2dec([extractBetween(anchors_hex(k), 1, 2); extractBetween(anchors_hex(k), 3, 4); ...
                       extractBetween(anchors_hex(k), 5, 6)])' / 255;
end
if n == 1
    C = A(end, :);
    return
end
x = linspace(1, numel(anchors_hex), n);
C = interp1(1:numel(anchors_hex), A, x);
end
