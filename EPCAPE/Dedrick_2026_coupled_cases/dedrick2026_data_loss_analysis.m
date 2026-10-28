%% Where the data go: Dedrick et al. (2026)'s event selection, my cloud-scene filters, and what relaxing them gains
%
% DESCRIPTION:
%   The SW transmittance vs. CCN plot (EPCAPE/SW_trans_vs_CCN) starts from
%   the 24 Mt. Soledad cloud events of Dedrick et al. (2026) and then applies
%   cloud-scene filters (daytime, 20 < LWP < 250 g m-2, no rain, single
%   layer, liquid, persistent deck at the Pier and Mt. Soledad). This script
%   counts the data lost at every step, so the criteria worth relaxing can be
%   chosen with numbers in hand. It answers four questions:
%
%   A. How were the 24 events selected, and how many more would a looser
%      selection give? The GCVI-based selection of the paper (supplement
%      Table S4: 154 -> 95 -> 44 -> 24 events) is rebuilt from the posted
%      1-s GCVI file, the Mt. Soledad visibility and the FM-120, then re-run
%      with looser thresholds.
%   B. Event by event, how many hours of the 24 survive each step, and which
%      events end up in the plot?
%   C. Are my filters redundant with the event selection? For the minutes
%      inside the 24 in-cloud windows: how much each filter removes on its
%      own, how much it removes that NO other filter would ("unique"; ~0
%      means redundant), and how many plotted hours would come back if it
%      were dropped.
%   D. Which case sets give the most usable hours: the 24 events, all hours
%      with Mt. Soledad in cloud, Dedrick's low-cloud periods, or all hours?
%
%   WHAT DEDRICK'S EVENT CRITERIA DO AND DO NOT GUARANTEE (supplement Text
%   S7, Table S4)
%   The events require the GCVI on (visibility < 5 km, no precipitation
%   detected at Mt. Soledad) for > 1 h, cloud (visibility < 1 km, LWC > 0.01
%   g m-3) for > 2 h, and SEMS + FM-120 data. Nothing in that list tests
%   daytime, LWP, single-layer cloud, cloud phase, rain at the Pier, or the
%   cloud over the Pier. The single-layer, < 3 km, coupled criteria belong
%   to Dedrick's LOW-CLOUD PERIODS (cloud over the Pier, NOT enveloping Mt.
%   Soledad), and the 20-300 g m-2 LWP, daytime, overcast and single-layer
%   screens are those of the ARM NDROP droplet-number product (Text S5),
%   which Dedrick used only for the low-cloud periods. So my filters are not
%   implied by the 24 events; Part C measures how much they overlap anyway.
%
% OUTPUTS (in ./output):
%   dedrick2026_data_loss_tables.xlsx  one sheet per table printed below
%   dedrick2026_data_loss.png/.pdf     the summary figure
%
% DATA AND REFERENCES: see sw_transmittance_vs_ccn.m (same inputs and
%   cache), and
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement doi:10.1175/JAS-D-26-0004.s1
%     (Text S5, S7; Tables S4-S6).
%   Berta, V., and coauthors, 2025: Aerosol quantification products at Mt.
%     Soledad and Scripps Pier during EPCAPE 2023-24. UC San Diego Library
%     Digital Collections, doi:10.6075/J0WM1DSC (GCVI 1-s file, CCN,
%     supersaturations, low-cloud periods).
%
% By Andrew John Buggee
%%

clear variables


%% 1. Settings (the filter block is the same as in sw_transmittance_vs_ccn.m)

% --- my cloud-scene filters ----------------------------------------------------
applyDaytime     = true;   szaMax_deg     = 80;
applyLwp         = true;   lwpMin_gm2     = 20;     lwpMax_gm2 = 250;
applyNoPrecip    = true;   rainWindow_min = 30;     rainSites  = ["Pier" "Soledad"];
applySingleLayer = true;   requireArscl   = false;
applyLiquid      = true;   maxCloudTop_m  = 3000;   maxCloudBase_m = 3000;
applyPersistence = true;   persistenceSites      = ["Pier" "Soledad"];
                           persistHalfWindow_min = 30;
                           persistMinFraction    = 0.9;
                           persistMinCoverage    = 0.5;
minPassMinutesPerHour = 30;

% --- what a usable hour also needs ------------------------------------------------
ccnSetpoint_pct = 0.2;                                   % CCN on the x axis
postedMethod    = "sNd";                                 % Dedrick's posted retrieval used for colour
postedPeriods   = ["SoledadCloudEvent" "LowCloudPeriod"];

% --- Dedrick's event selection: the paper's thresholds, then looser ones ---------------
% columns: MinGcviDuration_h, MergeGap_h, MinCloudHours, MinDataCoverage
cascadeScenarios = table( ...
    ["paper criteria"; "cloud > 1 h"; "cloud > 0.5 h"; "GCVI > 0.5 h, cloud > 1 h"; "merge within 3 h"; "no data-coverage test"], ...
    [1; 1; 1; 0.5; 1; 1], [6; 6; 6; 6; 3; 6], [2; 1; 0.5; 1; 2; 2], [0.5; 0.5; 0.5; 0.5; 0.5; 0], ...
    'VariableNames', {'scenario', 'MinGcviDuration_h', 'MergeGap_h', 'MinCloudHours', 'MinDataCoverage'});

% --- folders ------------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
epcapeDir = fileparts(scriptDir);
outDir    = fullfile(scriptDir, 'output');
processed = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed";
productsFolder = fullfile(processed, "Aerosol Quantification Products");
cacheFile = fullfile(epcapeDir, 'SW_trans_vs_CCN', 'output', 'epcape_1min_inputs.mat');   % shared with the plot script

addpath(scriptDir);
addpath(fullfile(epcapeDir, 'SW_trans_vs_CCN', 'functions'));
addpath(fullfile(epcapeDir, 'SW_transmittance'));
addpath(fullfile(epcapeDir, 'supersaturation', 'functions'));
if ~isfolder(outDir), mkdir(outDir); end
xlsxFile = fullfile(outDir, 'dedrick2026_data_loss_tables.xlsx');
if isfile(xlsxFile), delete(xlsxFile); end

F = struct('applyDaytime', applyDaytime, 'szaMax_deg', szaMax_deg, ...
           'applyLwp', applyLwp, 'lwpMin_gm2', lwpMin_gm2, 'lwpMax_gm2', lwpMax_gm2, ...
           'applyNoPrecip', applyNoPrecip, 'rainWindow_min', rainWindow_min, 'rainSites', rainSites, ...
           'applySingleLayer', applySingleLayer, 'requireArscl', requireArscl, ...
           'applyLiquid', applyLiquid, 'maxCloudTop_m', maxCloudTop_m, 'maxCloudBase_m', maxCloudBase_m, ...
           'applyPersistence', applyPersistence, 'persistenceSites', persistenceSites, ...
           'persistHalfWindow_min', persistHalfWindow_min, 'persistMinFraction', persistMinFraction, ...
           'persistMinCoverage', persistMinCoverage);



%% 2. Inputs: the 1-min grid (cached), hourly CCN and supersaturations

if isfile(cacheFile)
    fprintf('Loading the 1-min inputs from %s\n', cacheFile);
    load(cacheFile, 'D', 'cloud5', 'fm');
else
    error(['Run EPCAPE/SW_trans_vs_CCN/sw_transmittance_vs_ccn.m once first: it builds the 1-min ' ...
           'cache this script reads (%s).'], cacheFile);
end
[H, s_set_pct] = load_hourly_ccn_and_supersaturation(productsFolder, fm);

% what a usable hour needs besides the filters: CCN at the setpoint and a supersaturation
ccn_h = H.ccn_cm3(:, abs(s_set_pct - ccnSetpoint_pct) < 1e-9);
s_posted = nan(height(H), 1);
for period = postedPeriods
    col = postedMethod + "_" + period;
    fill_in = isnan(s_posted) & isfinite(H.(col));
    s_posted(fill_in) = H.(col)(fill_in);
end
s_ajb = H.sNdFM_AJB;
need = struct('ccn', isfinite(ccn_h), 'posted', isfinite(s_posted), 'ajb', isfinite(s_ajb));

M = analysis_filters_1min(D, F);
filters = M.applied;
fprintf('Filters applied: %s\n', strjoin(filters, ', '));



%% Part A. Rebuild Dedrick's event selection (supplement Table S4), then loosen it

[segs, ~] = gcvi_sampling_segments(fullfile(productsFolder, "bb2743661p_36_1_CVI_EF.csv"));
[ev_paper, cascade_summary] = find_soledad_cloud_events(segs, cloud5);
fprintf('\nA1. Event selection rebuilt with the paper''s thresholds (cloud = 5-min intervals with vis < 1 km and LWC > 0.01):\n');
disp(cascade_summary)
paper = dedrick2026_soledad_cloud_events();
sel = ev_paper(ev_paper.selected, :);
found = intersect(paper.event_id, sel.dedrick_event_id);
fprintf(['The rebuilt events include %d of the paper''s 24 (missing: %s). The paper''s final step also needs SEMS\n' ...
         'data, which are not available here, so the rebuilt step 4 is an upper bound.\n'], numel(found), ...
        mat2str(setdiff(paper.event_id, found)'));
writetable(cascade_summary, xlsxFile, 'Sheet', 'A1 cascade');
writetable(ev_paper, xlsxFile, 'Sheet', 'A1 rebuilt events');

% looser thresholds
A2 = cascadeScenarios;
A2.events_step3 = zeros(height(A2), 1);   A2.cloud_hours_step3 = zeros(height(A2), 1);
A2.events_step4 = zeros(height(A2), 1);   A2.cloud_hours_step4 = zeros(height(A2), 1);
A2.paper_events_included = zeros(height(A2), 1);
for i = 1:height(A2)
    [e, s] = find_soledad_cloud_events(segs, cloud5, 'MinGcviDuration_h', A2.MinGcviDuration_h(i), ...
        'MergeGap_h', A2.MergeGap_h(i), 'MinCloudHours', A2.MinCloudHours(i), 'MinDataCoverage', A2.MinDataCoverage(i));
    A2.events_step3(i) = s.n_events(3);   A2.cloud_hours_step3(i) = s.cloud_hours(3);
    A2.events_step4(i) = s.n_events(4);   A2.cloud_hours_step4(i) = s.cloud_hours(4);
    A2.paper_events_included(i) = numel(intersect(paper.event_id, e.dedrick_event_id(e.selected)));
end
fprintf('\nA2. The event selection with looser thresholds (cloud hours = 5-min cloud intervals):\n');
disp(A2)
writetable(A2, xlsxFile, 'Sheet', 'A2 looser selection');



%% Part B. The 24 events, step by step

[in24, eid] = case_set_mask_1min(D, "dedrick24", cloud5, H);
in24_cloud  = case_set_mask_1min(D, "dedrick24", cloud5, H, true);
B = paper(:, {'event_id', 'in_start_utc', 'in_end_utc', 'in_window_hours'});
B.cloud_interval_h  = zeros(24, 1);   % hours of the window meeting the 5-min cloud condition
B.daytime_h         = zeros(24, 1);   % hours with >= minPass daytime minutes in the window
B.all_filters_h     = zeros(24, 1);   % ... passing every filter
B.with_T_h          = zeros(24, 1);   % ... and a valid hourly SW transmittance
B.with_ccn_h        = zeros(24, 1);   % ... and Pier CCN
B.with_posted_s_h   = zeros(24, 1);   % ... and Dedrick's posted s (plotted with "posted")
B.with_ajb_s_h      = zeros(24, 1);   % ... and my s_Nd,FM (plotted with "AJB")
B.posted_s_any_time_h = zeros(24, 1); % posted s hours inside the window, any time of day, no filters
for k = 1:24
    ink = eid == k;
    B.cloud_interval_h(k) = nnz(ink & in24_cloud) / 60;
    B.daytime_h(k) = n_hours(D.time, ink & M.daytime, minPassMinutesPerHour);
    r = usable_hours(D, ink & M.all, minPassMinutesPerHour, szaMax_deg, need);
    B.all_filters_h(k)   = r.n_filters;
    B.with_T_h(k)        = r.n_T;
    B.with_ccn_h(k)      = r.n_ccn;
    B.with_posted_s_h(k) = r.n_posted;
    B.with_ajb_s_h(k)    = r.n_ajb;
    in_hour = n_hours_mask(D.time, ink, 30);      % hours at least half inside the window
    B.posted_s_any_time_h(k) = nnz(in_hour & need.posted);
end
fprintf('\nB. The 24 events (hours; an hour counts if >= %d of its minutes pass):\n', minPassMinutesPerHour);
disp(B)
totals = varfun(@sum, B(:, 4:end));
disp(totals)
writetable(B, xlsxFile, 'Sheet', 'B 24 events');



%% Part C. Are the filters redundant with the 24 events?

base = usable_hours(D, in24 & M.all, minPassMinutesPerHour, szaMax_deg, need);
n_in = nnz(in24);
C = table(filters', 'VariableNames', {'filter'});
C.description         = strings(numel(filters), 1);
C.removed_alone_pct   = zeros(numel(filters), 1);   % % of the in-window minutes failing this filter
C.removed_unique_pct  = zeros(numel(filters), 1);   % ... failing ONLY this filter
C.plotted_h_if_dropped = zeros(numel(filters), 1);  % plotted hours (posted s) without this filter
C.gain_h              = zeros(numel(filters), 1);
for i = 1:numel(filters)
    f = filters(i);
    others = in24;
    for g = filters(filters ~= f)
        others = others & M.(g);
    end
    C.description(i)        = M.description.(f);
    C.removed_alone_pct(i)  = 100 * nnz(in24 & ~M.(f)) / n_in;
    C.removed_unique_pct(i) = 100 * nnz(others & ~M.(f)) / n_in;
    r = usable_hours(D, others, minPassMinutesPerHour, szaMax_deg, need);
    C.plotted_h_if_dropped(i) = r.n_posted;
    C.gain_h(i) = r.n_posted - base.n_posted;
end

% the two halves of the persistence test, separately
for site = ["Pier" "Soledad"]
    F1 = F;  F1.persistenceSites = site;
    M1 = analysis_filters_1min(D, F1);
    C = [C; {"persistence " + site + " only", M1.description.persistence, ...
             100 * nnz(in24 & ~M1.persistence) / n_in, NaN, NaN, NaN}]; %#ok<AGROW>
end
fprintf(['\nC. Redundancy of each filter inside the 24 in-cloud windows (%d minutes = %.0f h).\n' ...
         '   Plotted hours with every filter: %d (posted s), %d (my s_Nd,FM).\n'], n_in, n_in / 60, base.n_posted, base.n_ajb);
disp(C)
writetable(C, xlsxFile, 'Sheet', 'C redundancy');



%% Part D. Which case sets give the most usable hours?

% rebuilt events with looser thresholds (cloud > 1 h, no data-coverage test), as windows
[e_loose, ~] = find_soledad_cloud_events(segs, cloud5, 'MinCloudHours', 1, 'MinDataCoverage', 0);
e_loose = e_loose(e_loose.selected, :);
in_loose = in_windows(D.time, e_loose.start_utc, e_loose.end_utc);

in_soledad = case_set_mask_1min(D, "soledadInCloud", cloud5, H);
in_low     = case_set_mask_1min(D, "lowCloudPeriods", cloud5, H);
in_all     = true(height(D), 1);

% the Pier-only persistence test (a deck whose base is above Mt. Soledad cannot be seen there)
F_pier = F;  F_pier.persistenceSites = "Pier";
M_pier = analysis_filters_1min(D, F_pier);
F_none = F;  F_none.applyPersistence = false;
M_none = analysis_filters_1min(D, F_none);

scenario = ["24 events, all filters"; "24 events, cloud intervals only"; "24 events, Pier-only persistence"; ...
            "24 events, no persistence test"; sprintf("rebuilt events, cloud > 1 h, no coverage test (%d)", height(e_loose)); ...
            "every Mt. Soledad in-cloud minute"; "low-cloud periods, Pier-only persistence"; "all hours, Pier-only persistence"];
masks = {in24 & M.all, in24_cloud & M.all, in24 & M_pier.all, in24 & M_none.all, in_loose & M.all, ...
         in_soledad & M.all, in_low & M_pier.all, in_all & M_pier.all};
Dtab = table(scenario);
Dtab.hours_case = [n_hours(D.time, in24, minPassMinutesPerHour); n_hours(D.time, in24_cloud, minPassMinutesPerHour); ...
                   n_hours(D.time, in24, minPassMinutesPerHour); n_hours(D.time, in24, minPassMinutesPerHour); ...
                   n_hours(D.time, in_loose, minPassMinutesPerHour); n_hours(D.time, in_soledad, minPassMinutesPerHour); ...
                   n_hours(D.time, in_low, minPassMinutesPerHour); n_hours(D.time, in_all, minPassMinutesPerHour)];
for name = ["after_filters" "with_T" "with_ccn" "with_posted_s" "with_ajb_s" "with_either_s" "days_either"]
    Dtab.(name) = zeros(height(Dtab), 1);
end
for i = 1:numel(masks)
    r = usable_hours(D, masks{i}, minPassMinutesPerHour, szaMax_deg, need);
    Dtab.after_filters(i) = r.n_filters;
    Dtab.with_T(i)        = r.n_T;
    Dtab.with_ccn(i)      = r.n_ccn;
    Dtab.with_posted_s(i) = r.n_posted;
    Dtab.with_ajb_s(i)    = r.n_ajb;
    Dtab.with_either_s(i) = r.n_either;
    Dtab.days_either(i)   = r.n_days_either;
end
fprintf('\nD. Usable hours by case set (an hour counts if >= %d of its minutes pass):\n', minPassMinutesPerHour);
disp(Dtab)
writetable(Dtab, xlsxFile, 'Sheet', 'D case sets');



%% Figure: where the 24 events' hours go, and what dropping each filter would give back

c_bar  = [0.165 0.471 0.839];    % #2a78d6, categorical slot 1
c_bar2 = [0.922 0.408 0.204];    % #eb6834, categorical slot 2
c_ink  = [0.043 0.043 0.043];
c_ink2 = [0.322 0.318 0.306];

% (a) the running total for the 24 events
running = in24;
labels = "24 in-cloud windows";
vals = n_hours(D.time, running, minPassMinutesPerHour);
for f = filters
    running = running & M.(f);
    labels(end + 1) = "+ " + f; %#ok<SAGROW>
    vals(end + 1) = n_hours(D.time, running, minPassMinutesPerHour); %#ok<SAGROW>
end
labels = [labels "+ valid T_{SW}" "+ Pier CCN" "+ posted s"];
vals   = [vals base.n_T base.n_ccn base.n_posted];

fig = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 12 5.5]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'off');
barh(ax, 1:numel(vals), vals, 0.7, 'FaceColor', c_bar, 'EdgeColor', 'none');
for i = 1:numel(vals)
    text(ax, vals(i) + 3, i, sprintf('%d h', vals(i)), 'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', c_ink);
end
set(ax, 'YDir', 'reverse', 'YTick', 1:numel(vals), 'YTickLabel', labels, 'XGrid', 'on', ...
    'TickLabelInterpreter', 'tex', 'XColor', c_ink2, 'YColor', c_ink);
xlim(ax, [0 1.15 * max(vals)]);
xlabel(ax, sprintf('Hours with >= %d passing minutes', minPassMinutesPerHour));
title(ax, '(a) The 24 events: hours left after each step', 'FontWeight', 'normal');

% (b) plotted hours regained by dropping one filter: 24 events vs. every Soledad in-cloud minute
gain24 = C.gain_h(1:numel(filters));
gainS = zeros(numel(filters), 1);
baseS = usable_hours(D, in_soledad & M.all, minPassMinutesPerHour, szaMax_deg, need);
for i = 1:numel(filters)
    others = in_soledad;
    for g = filters(filters ~= filters(i))
        others = others & M.(g);
    end
    r = usable_hours(D, others, minPassMinutesPerHour, szaMax_deg, need);
    gainS(i) = r.n_either - baseS.n_either;
end
ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'off');
bh = barh(ax, 1:numel(filters), [gain24 gainS], 'grouped', 'EdgeColor', 'none');
bh(1).FaceColor = c_bar;  bh(2).FaceColor = c_bar2;
% value at the end of every bar, so a zero is visible as "0"
for j = 1:2
    for i = 1:numel(filters)
        text(ax, bh(j).YEndPoints(i) + 0.15, bh(j).XEndPoints(i), sprintf('%d', bh(j).YData(i)), ...
             'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', c_ink2);
    end
end
row_labels = "drop " + filters;
row_labels(filters == "daytime") = "drop daytime (T_{SW} needs sun anyway)";
set(ax, 'YDir', 'reverse', 'YTick', 1:numel(filters), 'YTickLabel', row_labels, 'XGrid', 'on', ...
    'XColor', c_ink2, 'YColor', c_ink);
xlim(ax, [0 1.15 * max([gain24; gainS; 1])]);
xlabel(ax, 'Plotted hours gained');
legend(ax, {sprintf('24 events, posted s (now %d h)', base.n_posted), ...
            sprintf('every Soledad in-cloud minute, posted or my s (now %d h)', baseS.n_either)}, ...
       'Location', 'southoutside', 'Box', 'off');
title(ax, '(b) Hours regained by dropping one filter', 'FontWeight', 'normal');

title(tl, 'EPCAPE: data loss along Dedrick et al. (2026) events and the cloud-scene filters', 'FontSize', 12, 'Color', c_ink);
exportgraphics(fig, fullfile(outDir, 'dedrick2026_data_loss.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(outDir, 'dedrick2026_data_loss.pdf'), 'ContentType', 'vector');
close(fig);
fprintf('\nSaved the tables (%s) and the figure to %s\n', xlsxFile, outDir);




%% Local functions

function in_hour = n_hours_mask(t, mask, min_minutes)
% true for each hour [h, h + 1) with at least min_minutes true minutes
[~, frac] = window_mean(t, double(mask), 60);
in_hour = round(60 * frac) >= min_minutes;
end


function n = n_hours(t, mask, min_minutes)
% number of hours with at least min_minutes true minutes
n = nnz(n_hours_mask(t, mask, min_minutes));
end


function r = usable_hours(D, pass, min_minutes, szaMax_deg, need)
% Hours with >= min_minutes passing minutes, then with a valid hourly SW
% transmittance over those minutes, Pier CCN at the setpoint, and a
% supersaturation (posted, mine, or either)
Hr = hourly_from_passing_minutes(D, pass, min_minutes, szaMax_deg);
ok_filters = n_hours_mask(D.time, pass, min_minutes);
ok_T = isfinite(Hr.T_sw);
r.n_filters = nnz(ok_filters);
r.n_T       = nnz(ok_T);
r.n_ccn     = nnz(ok_T & need.ccn);
r.n_posted  = nnz(ok_T & need.ccn & need.posted);
r.n_ajb     = nnz(ok_T & need.ccn & need.ajb);
either = ok_T & need.ccn & (need.posted | need.ajb);
r.n_either  = nnz(either);
r.n_days_either = numel(unique(dateshift(Hr.time(either), 'start', 'day')));
end


function in = in_windows(t, starts, ends)
% true for each minute [t, t + 1 min) whose middle lies in any window [start, end)
mid = t + seconds(30);
in = false(size(t));
for k = 1:numel(starts)
    in = in | (mid >= starts(k) & mid < ends(k));
end
end
