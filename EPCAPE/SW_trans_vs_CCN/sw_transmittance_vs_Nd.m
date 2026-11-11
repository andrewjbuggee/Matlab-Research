%% EPCAPE: shortwave transmittance at the Scripps Pier vs. NDROP droplet number, colored by cloud supersaturation
%
% DESCRIPTION:
%   Companion to sw_transmittance_vs_ccn.m, with the CCN on the x axis
%   replaced by the cloud droplet number concentration N_d retrieved at the
%   Pier by the ARM NDROP value-added product. Scatter plot of the window-
%   averaged surface SW transmittance (y) against N_d (x), with the points
%   colored by the cloud's effective supersaturation (Dedrick's posted
%   hourly retrievals), optionally in nSupersaturationBins groups. Panel (b)
%   repeats T_SW against LWP, which explains most of the variance in T_SW.
%   The physical question: at a given LWP, does a cloud with more droplets
%   let less sunlight through (Twomey 1977), and does that relate to the
%   supersaturation the cloud reached?
%
%   WHAT IS PLOTTED, AND FROM WHERE
%   y  T_SW = <SW down> / <clear-sky SW down> over the window, from ARM
%      RADFLUX at the Pier (Long & Ackerman 2000), using only the minutes
%      that pass every filter AND have an NDROP value (EPCAPE/SW_transmittance).
%   x  N_d (cm-3), ARM NDROP (epcndropmfrsrM1.c1; Zhang et al. 2023): the
%      cloud-layer droplet number from the MFRSR optical depth and the MWR
%      LWP with an adiabatic cloud model, every 20 s, averaged to 1 min
%      (functions/read_ndrop_1min.m) and then over the same passing minutes
%      as T_SW. NDROP only exists in daylight (cos SZA >= 0.2) for overcast
%      (cloud fraction >= 0.9), single-layer, warm cloud, so the NDROP
%      requirement is itself a strict daytime + overcast screen.
%   color  Dedrick's posted hourly effective supersaturation (method and
%      period chosen below). Windows without a posted value are still
%      plotted, as hollow grey circles, because N_d does not need it.
%
%   WHY THERE IS MORE DATA THAN IN THE CCN PLOT
%   The Pier CCN counter came online on 23 May 2023, so events 1-7 had no
%   CCN; NDROP covers the whole campaign. N_d also has 20-s resolution, so
%   with 10-min windows every point has its own N_d (the CCN was hourly).
%
%   CAUTION WHEN READING THE PLOT
%   * N_d and T_SW are NOT independent measurements: NDROP derives N_d
%     from the cloud optical depth tau that the MFRSR infers from the SW
%     transmission, and for an adiabatic cloud N_d ~ tau^3 LWP^(-5/2)
%     (adiabatic relations reviewed by Grosvenor et al. 2018, Sec. 2). At fixed LWP, a cloud that
%     transmits less (larger tau) is assigned more droplets by
%     construction. A negative T_SW-N_d relation is therefore expected from
%     the retrieval alone, and must be read at fixed LWP.
%   * With postedPeriods including "LowCloudPeriod", the color can be
%     s_Nd,VAP, which Dedrick derived from this same NDROP N_d: then color
%     and x share an input. The default uses only "SoledadCloudEvent"
%     (s_Nd,FM, from the FM-120 droplet number at Mt. Soledad).
%   * Dedrick et al. (2026, Text S5) also removed NDROP values where the
%     1-min N_d exceeded the CN10 (CPC) concentration (5 % of the time).
%     CN10 is not on disk, so that screen is not applied here.
%   * The points are not independent (many per event).
%
% OUTPUTS (in ./output):
%   sw_transmittance_vs_Nd_<tag>.png/.pdf  the figure
%   sw_transmittance_vs_Nd_<tag>.csv       every plotted window
%   epcape_1min_inputs.mat                 shared cache of the 1-min inputs
%   ndrop_1min_<qc>.mat                    cache of the 1-min NDROP N_d
%
% DATA (Python project data/processed and data/arm, Dedrick's posted products):
%   ARM EPCAPE M1: NDROP (epcndropmfrsrM1.c1), RADFLUX, MWRRET v2,
%     ceilometer, KAZR-ARSCL, laser disdrometer; ARM S2: laser disdrometer.
%   Mt. Soledad visibility and FM-120 (case set and persistence test).
%   Berta, V., S. Han, A. Osawa, J. Dedrick, N. Maneenoi, and L. M. Russell,
%     2025: Aerosol quantification products at Mt. Soledad and Scripps Pier
%     during EPCAPE 2023-24. UC San Diego Library Digital Collections,
%     doi:10.6075/J0WM1DSC (supersaturations).
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: Quantitative estimates of upwind
%     and on-site particle scavenging and retrieved cloud supersaturations
%     from urban coastal stratocumulus during EPCAPE. J. Atmos. Sci., 83,
%     1061-1082, doi:10.1175/JAS-D-26-0004.1, and supplement
%     doi:10.1175/JAS-D-26-0004.s1.
%   Zhang, D., and coauthors, 2023: Evaluation of four ground-based
%     retrievals of cloud droplet number concentration in marine
%     stratocumulus with aircraft in situ measurements. Atmos. Meas. Tech.,
%     16, 5827-5846, doi:10.5194/amt-16-5827-2023.
%   Grosvenor, D. P., and coauthors, 2018: Remote sensing of droplet number
%     concentration in warm clouds: A review of the current state of
%     knowledge and perspectives. Rev. Geophys., 56, 409-453,
%     doi:10.1029/2017RG000593.
%   Long, C. N., and T. P. Ackerman, 2000: J. Geophys. Res., 105,
%     15609-15626, doi:10.1029/2000JD900077.
%   Twomey, S., 1977: The influence of pollution on the shortwave albedo of
%     clouds. J. Atmos. Sci., 34, 1149-1152.
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
applyDaytime     = true;    % apply a SZA filter on the data
szaMax_deg       = 80;      % sun higher than 10 deg from the horizon

applyLWP         = true;    % apply a LWP filter on the data
lwpMin_gm2       = 20;      % LWP > 20 g m-2 ...
lwpMax_gm2       = 300;     % ... and < 300 g m-2 (MWRRET v2 at the Pier) - value used by Dedrick et al. (2026)

applyNoPrecip    = true;    % filter the data for cases without precipitation
rainWindow_min   = 15;      % no rain within +/- rainWindow_min minutes ...
rainSites        = ["Pier" "Soledad"];   % ... at these disdrometers

requireArscl     = false;   % true: a minute without ARSCL data is ignored

% cloud filters
applySingleLayer       = true;    % filter for cloud scenes with only 1 cloud layer
applyLiquid            = true;    % filter for cloud scenes with liquid cloud only
maxCloudTop_m          = 3000;    % liquid proxy: top below 3 km
maxCloudBase_m         = 3000;    % base below 3 km

% cloud persistence (minute by minute; see sw_transmittance_vs_ccn.m)
applyPersistence       = true;
persistenceSites       = ["Pier" "Soledad"];
persistHalfWindow_min  = 10;     % persistence window = minute +/- persistHalfWindow_min
persistMinFraction     = 0.9;    % cloud for >= 90 % of the window minutes with data ...
persistMinCoverage     = 0.5;    % ... and data for >= 50 % of the window

% --- NDROP droplet number (x axis) ------------------------------------------------
ndropQc   = "strict";   % "strict": qc_drop_number_conc == 0 | "bad": reject only ARM "Bad" bits (+ bit 9)
ndLogAxis = false;      % log x axis

% --- averaging window: one plotted point per window ---------------------------------
% T_SW, LWP and N_d are averaged over the minutes of each window
% [t, t + averageWindow_min) that pass every filter AND have an NDROP value.
% N_d is 20-s data, so (unlike CCN) each window has its own N_d; the
% posted supersaturation is hourly, so windows of one hour share their s.
averageWindow_min = 10;    % minutes
minPassFraction   = 0.5;   % a window is plotted if >= this fraction of its minutes pass (incl. having N_d)

% --- color: cloud supersaturation ------------------------------------------------
supersaturationMethod = "sNd";      % posted: "sNd" | "sCCN" | "skappa" | "sARG00"
postedPeriods = "SoledadCloudEvent";   % add "LowCloudPeriod" to fill gaps (s_Nd,VAP uses NDROP N_d itself)
binBySupersaturation  = true;       % false: continuous color scale with a color bar
nSupersaturationBins  = 5;          % 5, 6, 7, ...
binEdgeRule           = "quantile"; % "quantile" (equal numbers of windows) | "equalWidth"

% --- figure and run control ---------------------------------------------------------
showLwpPanel = true;      % second panel: T_SW vs LWP, same colors
rebuildCache = false;     % true: re-read every instrument (incl. NDROP) instead of the caches

% --- folders --------------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
epcapeDir = fileparts(scriptDir);
outDir    = fullfile(scriptDir, 'output');
dataDir   = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data";
processed = fullfile(dataDir, "processed");
productsFolder = fullfile(processed, "Aerosol Quantification Products");
ndropFolder    = fullfile(dataDir, "arm", "epcndropmfrsrM1.c1");

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
addpath(fullfile(epcapeDir, 'supersaturation', 'functions'));   % posted products + new_figure
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



%% 2. The 1-min inputs (shared cache with sw_transmittance_vs_ccn.m) and NDROP

cacheFile = fullfile(outDir, 'epcape_1min_inputs.mat');
if isfile(cacheFile) && ~rebuildCache
    fprintf('Loading the 1-min inputs from %s\n', cacheFile);
    load(cacheFile, 'D', 'cloud5', 'fm');
else
    fprintf('Building the 1-min inputs (every instrument, ~1-2 min) ...\n');
    [D, cloud5, fm] = build_epcape_1min_inputs(files);
    save(cacheFile, 'D', 'cloud5', 'fm', 'files', '-v7.3');
end

ndCache = fullfile(outDir, "ndrop_1min_" + ndropQc + ".mat");
if isfile(ndCache) && ~rebuildCache
    load(ndCache, 'nd_cm3', 'n_nd', 't_nd');
    assert(isequal(t_nd, D.time), 'NDROP cache is on a different 1-min grid; set rebuildCache = true');
else
    fprintf('Reading NDROP (%s QC) ...\n', ndropQc);
    [nd_cm3, n_nd] = read_ndrop_1min(ndropFolder, D.time, ndropQc);
    t_nd = D.time;
    save(ndCache, 'nd_cm3', 'n_nd', 't_nd');
end
D.nd_cm3 = nd_cm3;



%% 3. Hourly posted supersaturations (the CCN columns of H are not used here)

H = load_hourly_ccn_and_supersaturation(productsFolder, fm);

s_pct = nan(height(H), 1);
for period = postedPeriods                   % first period with a value wins
    col = supersaturationMethod + "_" + period;
    fill_in = isnan(s_pct) & isfinite(H.(col));
    s_pct(fill_in) = H.(col)(fill_in);
end
sLabel = sprintf('Dedrick et al. (2026) posted s_{%s}', erase(supersaturationMethod, "s"));



%% 4. Cases and filters, minute by minute

[in_case, event_id_1min] = case_set_mask_1min(D, caseSet, cloud5, H, cloudIntervalsOnly);
M = analysis_filters_1min(D, F);
has_nd = isfinite(D.nd_cm3);
pass = in_case & M.all & has_nd;             % the minutes behind every plotted value



%% 5. Window averages over the passing minutes

% T_SW, LWP and SZA of each averaging window, from its passing minutes only
Hr = average_passing_minutes(D, pass, averageWindow_min, minPassMinutes, szaMax_deg);

% N_d over the same passing minutes
nd = D.nd_cm3;  nd(~pass) = NaN;
[t_chk, Hr.nd_cm3] = window_mean(D.time, nd, averageWindow_min);
assert(isequal(t_chk, Hr.time), 'window grids differ');

% Dedrick event number of each window
ev = event_id_1min;  ev(~pass) = NaN;
[~, ev_win] = window_mean(D.time, ev, averageWindow_min);
Hr.event_id = round(ev_win);
Hr.hour_utc = dateshift(Hr.time, 'start', 'hour');

% hourly s onto the windows (a window inside one hour gets that hour's value)
[Hr.s_pct, t_chk] = hourly_values_on_windows(D.time, pass, H.time, s_pct, averageWindow_min);
assert(isequal(t_chk, Hr.time), 'window grids differ');



%% 6. How many windows survive each step

count_windows = @(mask) nnz(window_count(D.time, mask, averageWindow_min) >= minPassMinutes);
step_name = "in case set (" + caseSet + ")";
step_windows = count_windows(in_case);
running = in_case;
for name = M.applied
    running = running & M.(name);
    step_name(end + 1, 1) = "+ " + name + ": " + M.description.(name); %#ok<SAGROW>
    step_windows(end + 1, 1) = count_windows(running); %#ok<SAGROW>
end
running = running & has_nd;
has_T = isfinite(Hr.T_sw) & isfinite(Hr.nd_cm3);
has_s = has_T & isfinite(Hr.s_pct);
step_name    = [step_name; "+ NDROP N_d (" + ndropQc + " QC)"; "+ valid SW and clear-sky SW = PLOTTED"; ...
                "  of which with a posted s (colored)"];
step_windows = [step_windows; count_windows(running); nnz(has_T); nnz(has_s)];
step_hours   = step_windows * averageWindow_min / 60;
steps = table(step_name, step_windows, step_hours, ...
              'VariableNames', {'step', sprintf('windows_%gmin', averageWindow_min), 'hours'});
disp(steps)

P = Hr(has_T, :);                            % the plotted windows
if height(P) == 0
    error('No windows pass every step; relax a filter (see the table above).');
end
ev_list = unique(P.event_id(isfinite(P.event_id)))';
fprintf('Plotted: %d %g-min averages (%d with a posted s)', height(P), averageWindow_min, nnz(isfinite(P.s_pct)));
if caseSet == "dedrick24"
    fprintf(', %d of the 24 events (%s)', numel(ev_list), mat2str(ev_list));
end
fprintf('\n');



%% 7. Supersaturation bins (over the windows that have a posted s)

with_s = isfinite(P.s_pct);
if binBySupersaturation && any(with_s)
    switch binEdgeRule
        case "quantile"
            edges = prctile(P.s_pct(with_s), linspace(0, 100, nSupersaturationBins + 1));
        case "equalWidth"
            edges = linspace(min(P.s_pct), max(P.s_pct), nSupersaturationBins + 1);
        otherwise
            error('binEdgeRule must be "quantile" or "equalWidth"');
    end
    edges = unique(edges);                   % ties can merge quantile edges
    if numel(edges) - 1 < nSupersaturationBins
        warning('Only %d distinct bins could be formed.', numel(edges) - 1);
    end
    P.s_bin = discretize(P.s_pct, edges);    % NaN for windows without s
    n_bins = numel(edges) - 1;
    binColors = ordinal_blue_ramp(n_bins);
end



%% 8. Figure

c_ink   = [0.043 0.043 0.043];    % #0b0b0b, primary text
c_ink2  = [0.322 0.318 0.306];    % #52514e, secondary text
c_nos   = [0.62 0.62 0.60];       % grey ring: no posted s

n_panels = 1 + showLwpPanel;
fig = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 6.5 * n_panels 6]);
tl = tiledlayout(fig, 1, n_panels, 'TileSpacing', 'compact', 'Padding', 'compact');

x_vars = {'nd_cm3', 'lwp_gm2'};
x_labels = {'Pier N_d over the same minutes (cm^{-3}, ARM NDROP)', ...
            'Pier LWP over the same minutes (g m^{-2}, MWRRET v2)'};
panel_titles = {'(a) SW transmittance vs. N_d', '(b) SW transmittance vs. LWP'};

r_panel = nan(n_panels, 1);
n_panel = zeros(n_panels, 1);
for ip = 1:n_panels
    ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
    x = P.(x_vars{ip});

    % Pearson r over every plotted window (colored or not)
    R = corrcoef(x, P.T_sw, 'Rows', 'complete');
    r_panel(ip) = R(1, 2);
    n_panel(ip) = nnz(isfinite(x) & isfinite(P.T_sw));

    % windows without a posted s: hollow grey, drawn first (underneath)
    h0 = scatter(ax, x(~with_s), P.T_sw(~with_s), 30, c_nos, 'LineWidth', 0.75, ...
                 'DisplayName', sprintf('no posted s  (n = %d)', nnz(~with_s)));
    if binBySupersaturation && any(with_s)
        h = gobjects(n_bins, 1);
        for b = 1:n_bins
            in = P.s_bin == b;
            h(b) = scatter(ax, x(in), P.T_sw(in), 45, binColors(b, :), 'filled', ...
                           'MarkerEdgeColor', 'w', 'LineWidth', 0.75, ...
                           'DisplayName', sprintf('%.3f-%.3f %%  (n = %d)', edges(b), edges(b + 1), nnz(in)));
        end
        if ip == 1
            lg = legend(ax, [h; h0], 'Box', 'off', 'FontSize', 9);
            lg.Layout.Tile = 'east';
            title(lg, sLabel, 'FontWeight', 'normal', 'FontSize', 9);
        end
    elseif any(with_s)
        scatter(ax, x(with_s), P.T_sw(with_s), 45, P.s_pct(with_s), 'filled', ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.75);
        colormap(ax, ordinal_blue_ramp(256));
        clim(ax, [min(P.s_pct) max(P.s_pct)]);
        if ip == n_panels
            cb = colorbar(ax);
            cb.Layout.Tile = 'east';
            cb.Label.String = [sLabel ' (%)'];
        end
    end
    if ip == 1 && ndLogAxis
        set(ax, 'XScale', 'log');
    end
    xlim(ax, 'padded');
    xlabel(ax, x_labels{ip});
    ylabel(ax, 'SW transmittance, T_{SW} (Pier, RADFLUX)');
    ylim(ax, [0 ceil(10 * (max(P.T_sw) + 0.05)) / 10]);
    title(ax, panel_titles{ip}, sprintf('Pearson r = %.2f  (n = %d)', r_panel(ip), n_panel(ip)), ...
          'FontWeight', 'normal');
    ax.Subtitle.Color = c_ink2;
    ax.XColor = c_ink;  ax.YColor = c_ink;
end

case_words = struct('dedrick24', 'the 24 Mt. Soledad cloud events of Dedrick et al. (2026), in-cloud windows', ...
                    'soledadInCloud', 'every hour with Mt. Soledad in cloud (vis < 1 km, LWC > 0.01 g m^{-3})', ...
                    'lowCloudPeriods', 'Dedrick et al. (2026) low-cloud periods at the Pier', ...
                    'all', 'all hours');
filter_words = [string(cellfun(@(nm) M.description.(nm), cellstr(M.applied), 'UniformOutput', false)), ...
                "NDROP N_d available (" + ndropQc + " QC)", ...
                sprintf("%g-min averages of the passing minutes (>= %d of %g min pass)", ...
                        averageWindow_min, minPassMinutes, averageWindow_min)];
title(tl, sprintf('EPCAPE, %d %g-min averages from %d events: %s', height(P), averageWindow_min, ...
                  numel(ev_list), case_words.(caseSet)), 'FontSize', 12, 'Color', c_ink);
subtitle(tl, wrap_words("Filters: " + strjoin(filter_words, '; '), 70 + 60 * showLwpPanel), ...
         'FontSize', 8, 'Color', c_ink2, 'Interpreter', 'none');

tag = sprintf('%s_posted_%s_%gmin', caseSet, supersaturationMethod, averageWindow_min);
exportgraphics(fig, fullfile(outDir, "sw_transmittance_vs_Nd_" + tag + ".png"), 'Resolution', 200);
exportgraphics(fig, fullfile(outDir, "sw_transmittance_vs_Nd_" + tag + ".pdf"), 'ContentType', 'vector');
close(fig);



%% 9. The plotted windows, and simple statistics

out = timetable2table(P);
out.Properties.VariableNames{1} = 'window_start_utc';
writetable(out, fullfile(outDir, "sw_transmittance_vs_Nd_" + tag + ".csv"));

fprintf('Pearson r(T_SW, N_d) over the %d plotted windows: %.2f\n', n_panel(1), r_panel(1));
if showLwpPanel
    fprintf('Pearson r(T_SW, LWP) over the same windows:      %.2f\n', r_panel(2));
end
fprintf(['N_d is retrieved from the MFRSR optical depth (i.e. from the transmitted SW), so T_SW and N_d\n' ...
         'are not independent; and the windows are not independent (several per event).\n']);
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
% n colors, light to dark, from the reference sequential blue ramp (steps 250 -> 700)
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
