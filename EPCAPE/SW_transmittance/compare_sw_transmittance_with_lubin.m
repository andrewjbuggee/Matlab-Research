%% Compare my SW transmittance (RADFLUX) with Dan Lubin's hourly SW transmittance at the Scripps Pier
%
% DESCRIPTION:
%   Two hourly SW transmittances exist for the EPCAPE Pier:
%     MINE  = <SW down> / <clear-sky SW down> from ARM RADFLUX, both 1-min,
%             QC = 0, sun above 90 - szaMax_deg degrees, averaged over each
%             UTC hour (sw_transmittance.m with 60-min windows; the hourly
%             counterpart of the Python quantities.sw_transmittance). The
%             clear-sky value is RADFLUX's empirical fit to the site's own
%             clear periods (Long & Ackerman 2000).
%     DAN'S = his <MON>_SW_Boxplot.csv tables: measured SW down over a
%             clear-sky ANALYTICAL model (Atwater & Ball 1981, per Kavin's
%             table note), one value per PST hour (read_lubin_sw_transmittance.m).
%   Both use the same pyranometer data (the SKYRAD at the Pier) as far as is
%   known, so the differences should come mainly from (1) the clear-sky
%   model, (2) the time base and averaging, (3) the sun-elevation cut and
%   (4) QC. This script measures them.
%
%   WHAT IT DOES
%     1. reads both; mine is computed for two sun cuts: szaMax_deg (my
%        method's default, 80 deg) and 90 deg (every daylight minute, which
%        is closer to Dan's tables: they have values at very low sun)
%     2. TIME ALIGNMENT: correlates Dan's series against mine for whole-hour
%        shifts of -2 ... +2 h, and for my hours offset by 30 min. The best
%        alignment should be shift 0 with whole hours if Dan's column h means
%        [h, h + 1) PST.
%     3. FEBRUARY: Dan's FEB file pools both Februaries. The correlation with
%        mine is reported for the two possible row-to-date readings.
%     4. AGREEMENT on the matched hours: r, OLS fit, mean difference, median
%        ratio Dan / mine, by month and by solar zenith angle.
%     5. CLEAR-SKY HOURS: hours whose daylight minutes RADFLUX flagged as
%        clear by its SW method (clearsky_status = 1) for >= 90 % of the
%        minutes. There mine is ~1 by construction, so Dan's value there is
%        directly SW_clear,RADFLUX / SW_clear,Atwater-Ball: the ratio of the
%        two clear-sky models.
%     6. a 4-panel figure and a CSV of every matched hour.
%
% DATA:
%   * RADFLUX (radflux1long c2, EPCAPE M1), combined file in the Python
%     project's data/processed folder.
%   * Dan Lubin's tables in VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper
%     (8 months: APR, MAY, JUL-DEC). The other four months (JAN, FEB, MAR,
%     JUN) are only in Kavin_code_for_BAMS/Annual Echo Plot, which is read
%     second (the shared month, APR, is byte-identical in both).
%
% REFERENCES:
%   Long, C. N., and T. P. Ackerman, 2000: Identification of clear skies
%     from broadband pyranometer measurements and calculation of downwelling
%     shortwave cloud effects. J. Geophys. Res., 105, 15609-15626,
%     doi:10.1029/2000JD900077.
%   Atwater, M. A., and J. T. Ball, 1981: the clear-sky model named in
%     Kavin's note for Dan's tables (not checked here).
%   Python originals: comparisons/seasonal_averages/quantities.py
%     (sw_transmittance) and kavin.py (read_lubin_sw, lubin_sw_series,
%     hourly_ratio, LUBIN_UTC_OFFSET_H = 8).
%
% By Andrew John Buggee
%%

clear variables


%% 1. Settings

% --- folders and files -------------------------------------------------------
scriptDir   = fileparts(mfilename('fullpath'));
outDir      = fullfile(scriptDir, 'output');
processed   = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed";
radfluxFile = fullfile(processed, "radflux_M1_20230215_20240214.nc");
lubinFolders = ["/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper", ...
                "/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/EPCAPE/Kavin_code_for_BAMS/Annual Echo Plot"];

% new_figure.m (a figure with MATLAB's factory defaults, so startup.m's dark
% theme doesn't apply) lives with the supersaturation functions
addpath(fullfile(scriptDir, '..', 'supersaturation', 'functions'));

% --- my transmittance --------------------------------------------------------
szaMax_deg      = 80;     % my method's sun cut (Python sw_transmittance default)
minCountPerHour = 30;     % an hour needs >= 30 valid 1-min samples (Python hourly_ratio default)

% --- Dan's tables --------------------------------------------------------------
febMapping  = "chronological";   % "chronological" (Fig. 4 matching) or "dayOfMonth" (Python assumption)
utcOffset_h = 8;                 % Dan's hours are PST = UTC - 8

% --- clear-sky hours -----------------------------------------------------------
clearFractionMin = 0.9;   % an hour is "clear" if >= 90 % of its daylight minutes have clearsky_status = 1

if ~isfolder(outDir), mkdir(outDir); end
assert(isfile(radfluxFile), 'RADFLUX file not found: %s', radfluxFile);



%% 2. My hourly transmittance from RADFLUX

rad = read_radflux_sw(radfluxFile);
fprintf('RADFLUX: %d 1-min samples, %s to %s; SW down kept after QC = 0: %d of %d finite\n', ...
        numel(rad.time), string(rad.time(1)), string(rad.time(end)), rad.qc_info.n_kept, rad.qc_info.n_finite);

% 60-min windows anchored at midnight UTC and labelled by their start =
% UTC hour bins [h, h + 1), on minute START times
[T_mine80, t_hour, info80] = sw_transmittance(rad, 'SzaMax_deg', szaMax_deg, 'WindowMin', 60, 'MinCount', minCountPerHour);
[T_mine90, ~,      ~]      = sw_transmittance(rad, 'SzaMax_deg', 90,         'WindowMin', 60, 'MinCount', minCountPerHour);

% the same with hours offset by 30 min ([h - 30 min, h + 30 min)), for the
% alignment test: shifting every sample 30 min LATER puts the original
% interval [h - 30 min, h + 30 min) into the window labelled h
rad_off = rad;  rad_off.time = rad.time + minutes(30);
[T_mine90_off, t_hour_off] = sw_transmittance(rad_off, 'SzaMax_deg', 90, 'WindowMin', 60, 'MinCount', minCountPerHour);

% mean cos(SZA) of each hour's daylight minutes, for the SZA dependence
mu = rad.cos_sza;  mu(mu <= 0) = NaN;
[~, mu_hour] = window_mean(rad.time, mu, 60);
sza_hour_deg = acosd(mu_hour);

% fraction of each hour's daylight minutes that RADFLUX's SW method called clear
clr_status = double(ncread(radfluxFile, 'clearsky_status'));
is_clear = double(clr_status(:) == 1);  is_clear(~(rad.cos_sza > 0)) = NaN;
[~, clear_frac_hour] = window_mean(rad.time, is_clear, 60);

MINE = timetable(t_hour, T_mine80, T_mine90, sza_hour_deg, clear_frac_hour, info80.n_per_window, ...
                 'VariableNames', {'T_mine80', 'T_mine90', 'sza_deg', 'clear_frac', 'n_minutes80'});
fprintf('My hourly transmittance: %d hours (SZA < %g deg), %d hours (SZA < 90 deg)\n', ...
        nnz(isfinite(T_mine80)), szaMax_deg, nnz(isfinite(T_mine90)));



%% 3. Dan Lubin's hourly transmittance

[DAN, sources] = read_lubin_sw_transmittance(lubinFolders, 'FebMapping', febMapping, 'UtcOffset_h', utcOffset_h);
month_names = ["JAN" "FEB" "MAR" "APR" "MAY" "JUN" "JUL" "AUG" "SEP" "OCT" "NOV" "DEC"];
for m = 1:12
    if sources.file(m) == ""
        fprintf('  %s: no file\n', month_names(m));
    else
        [~, folderName] = fileparts(fileparts(sources.file(m)));
        fprintf('  %s: %4d hours, from %s\n', month_names(m), nnz(DAN.month == m), folderName);
    end
end



%% 4. Time alignment: which shift of Dan's series matches mine best?
% Dan's value labelled UTC hour h is compared with mine at h + shift. If
% the column label is the hour's start and the PST -> UTC offset is right,
% shift 0 on whole hours is best.

shifts_h = -2:2;
r_whole = nan(size(shifts_h));   n_whole = r_whole;
r_half  = nan(size(shifts_h));   n_half  = r_half;
for k = 1:numel(shifts_h)
    [r_whole(k), n_whole(k)] = corr_at(DAN.time + hours(shifts_h(k)), DAN.T_lubin, t_hour, T_mine90);
    [r_half(k),  n_half(k)]  = corr_at(DAN.time + hours(shifts_h(k)), DAN.T_lubin, t_hour_off, T_mine90_off);
end
alignment = table(shifts_h', r_whole', n_whole', r_half', n_half', 'VariableNames', ...
    {'shift_h', 'r_mine_on_[h,h+1)', 'n', 'r_mine_on_[h-30min,h+30min)', 'n_half'});
disp('Time alignment (Dan vs mine, SZA < 90 deg):'); disp(alignment)



%% 5. February: which row-to-date reading fits?

feb_check = table();
for mapping = ["chronological" "dayOfMonth"]
    D = read_lubin_sw_transmittance(lubinFolders, 'FebMapping', mapping, 'UtcOffset_h', utcOffset_h);
    D = D(D.month == 2, :);
    [r, n] = corr_at(D.time, D.T_lubin, t_hour, T_mine90);
    feb_check = [feb_check; table(mapping, r, n)]; %#ok<AGROW>
end
disp('FEB file, correlation with mine for the two row readings:'); disp(feb_check)



%% 6. Agreement on the matched hours

% join Dan's hours to mine (exact UTC hour starts)
[tf, loc] = ismember(DAN.time, MINE.Properties.RowTimes);
M = DAN(tf, :);
M.T_mine80   = MINE.T_mine80(loc(tf));
M.T_mine90   = MINE.T_mine90(loc(tf));
M.sza_deg    = MINE.sza_deg(loc(tf));
M.clear_frac = MINE.clear_frac(loc(tf));
M.ratio_dan_over_mine = M.T_lubin ./ M.T_mine80;

both = isfinite(M.T_lubin) & isfinite(M.T_mine80);
fprintf('\nMatched hours (Dan and mine with SZA < %g deg): %d\n', szaMax_deg, nnz(both));
overall = agreement_stats(M.T_mine80(both), M.T_lubin(both));
disp(struct2table(overall))

% by month (calendar month of Dan's file)
byMonth = table();
for m = 1:12
    k = both & M.month == m;
    if nnz(k) < 5, continue, end
    s = agreement_stats(M.T_mine80(k), M.T_lubin(k));
    byMonth = [byMonth; [table(month_names(m), 'VariableNames', {'month'}), struct2table(s)]]; %#ok<AGROW>
end
disp('By month:'); disp(byMonth)

% by solar zenith angle
sza_edges = [0 30 40 50 60 70 80];
bySza = table();
for k = 1:numel(sza_edges) - 1
    in = both & M.sza_deg >= sza_edges(k) & M.sza_deg < sza_edges(k + 1);
    if nnz(in) < 5, continue, end
    s = agreement_stats(M.T_mine80(in), M.T_lubin(in));
    bySza = [bySza; [table(string(sprintf('%d-%d', sza_edges(k), sza_edges(k + 1))), 'VariableNames', {'sza_deg'}), struct2table(s)]]; %#ok<AGROW>
end
disp('By hourly-mean solar zenith angle:'); disp(bySza)

% clear-sky hours: the ratio of the two clear-sky models
clear_h = both & M.clear_frac >= clearFractionMin;
fprintf('Clear-sky hours (>= %.0f %% of daylight minutes flagged clear by RADFLUX): %d\n', 100 * clearFractionMin, nnz(clear_h));
fprintf('  mine:  median %.3f (IQR %.3f-%.3f)  [~1 by construction]\n', prctile(M.T_mine80(clear_h), [50 25 75]));
fprintf('  Dan''s: median %.3f (IQR %.3f-%.3f)  = SW_clear,RADFLUX / SW_clear,Atwater-Ball\n', prctile(M.T_lubin(clear_h), [50 25 75]));



%% 7. Figure: four panels

c_mine = [0.165 0.471 0.839];    % #2a78d6, categorical slot 1 (blue)
c_dan  = [0.922 0.408 0.204];    % #eb6834, categorical slot 2 (orange)
c_ink  = [0.043 0.043 0.043];    % #0b0b0b, primary text
c_ink2 = [0.322 0.318 0.306];    % #52514e, secondary text and reference lines

fig = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 11 8.5]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Dan vs mine, every matched hour
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
scatter(ax, M.T_mine80(both), M.T_lubin(both), 12, c_dan, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
plot(ax, [0 1.3], [0 1.3], '-', 'Color', c_ink2, 'LineWidth', 1);
plot(ax, [0 1.3], overall.intercept + overall.slope * [0 1.3], '--', 'Color', c_ink, 'LineWidth', 1.5);
axis(ax, [0 1.3 0 1.3]);  axis(ax, 'square');
xlabel(ax, 'My T_{SW} (RADFLUX SW / clear-sky SW)');  ylabel(ax, 'Dan Lubin''s T_{SW}');
title(ax, sprintf('(a) Every matched hour, n = %d', overall.n), 'FontWeight', 'normal');
legend(ax, {'hours', '1:1', sprintf('OLS: y = %.2f x + %.2f, r = %.3f', overall.slope, overall.intercept, overall.r)}, ...
       'Location', 'northwest', 'Box', 'off');

% (b) ratio vs solar zenith angle
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
scatter(ax, M.sza_deg(both), M.ratio_dan_over_mine(both), 12, c_dan, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
scatter(ax, M.sza_deg(clear_h), M.ratio_dan_over_mine(clear_h), 16, c_mine, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
yline(ax, 1, '-', 'Color', c_ink2, 'LineWidth', 1);
ylim(ax, [0.5 1.3]);  xlim(ax, [0 80]);
xlabel(ax, 'Hourly-mean solar zenith angle (deg)');  ylabel(ax, 'Dan / mine');
title(ax, '(b) Ratio vs. sun elevation', 'FontWeight', 'normal');
legend(ax, {'all matched hours', sprintf('clear-sky hours (n = %d)', nnz(clear_h))}, 'Location', 'southwest', 'Box', 'off');

% (c) monthly median ratio and IQR
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
mm = find(ismember(month_names, byMonth.month));
q = nan(numel(mm), 3);
for k = 1:numel(mm)
    q(k, :) = prctile(M.ratio_dan_over_mine(both & M.month == mm(k)), [25 50 75]);
end
errorbar(ax, mm, q(:, 2), q(:, 2) - q(:, 1), q(:, 3) - q(:, 2), 'o', 'Color', c_dan, ...
         'MarkerFaceColor', c_dan, 'MarkerSize', 6, 'LineWidth', 1.5, 'CapSize', 0);
yline(ax, 1, '-', 'Color', c_ink2, 'LineWidth', 1);
xlim(ax, [0.5 12.5]);  ylim(ax, [0.5 1.3]);
xticks(ax, 1:12);  xticklabels(ax, month_names);
xlabel(ax, 'Calendar month of Dan''s file');  ylabel(ax, 'Dan / mine');
title(ax, '(c) Monthly median ratio and IQR', 'FontWeight', 'normal');

% (d) an example week, both series
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
t0 = datetime(2023, 7, 1);  t1 = t0 + days(7);
% put the week on a full hourly grid so that night-time gaps (NaN) break the lines
W = retime(M(M.time >= t0 & M.time < t1, {'T_mine80', 'T_lubin'}), (t0:hours(1):t1 - hours(1))', 'fillwithmissing');
plot(ax, W.Properties.RowTimes, W.T_mine80, 'o-', 'Color', c_mine, 'MarkerFaceColor', c_mine, 'MarkerSize', 4, 'LineWidth', 1.5);
plot(ax, W.Properties.RowTimes, W.T_lubin,  's-', 'Color', c_dan,  'MarkerFaceColor', c_dan,  'MarkerSize', 4, 'LineWidth', 1.5);
ylim(ax, [0 1.3]);  ylabel(ax, 'T_{SW}');  xlabel(ax, 'UTC');
title(ax, '(d) Example week, 1-7 July 2023', 'FontWeight', 'normal');
legend(ax, {'mine (RADFLUX)', 'Dan Lubin'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

title(tl, 'EPCAPE Scripps Pier hourly SW transmittance: RADFLUX clear-sky ratio vs. Dan Lubin''s tables', 'FontSize', 13);

exportgraphics(fig, fullfile(outDir, 'compare_sw_transmittance_with_lubin.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(outDir, 'compare_sw_transmittance_with_lubin.pdf'), 'ContentType', 'vector');
close(fig);



%% 8. Save every matched hour

out = timetable2table(M(both, {'T_lubin', 'T_mine80', 'T_mine90', 'ratio_dan_over_mine', 'sza_deg', 'clear_frac', 'month', 'hour_pst'}));
out.Properties.VariableNames{1} = 'hour_start_utc';
writetable(out, fullfile(outDir, 'compare_sw_transmittance_with_lubin_matched_hours.csv'));
fprintf('\nSaved the figure and %d matched hours to %s\n', height(out), outDir);




%% Local functions

function [r, n] = corr_at(t_a, a, t_b, b)
% Pearson r between series a and b at the time stamps they share exactly
[tf, loc] = ismember(t_a, t_b);
x = a(tf);  y = b(loc(tf));
ok = isfinite(x) & isfinite(y);
n = nnz(ok);
if n < 3
    r = NaN;
    return
end
R = corrcoef(x(ok), y(ok));
r = R(1, 2);
end


function s = agreement_stats(x_mine, y_dan)
% n, Pearson r, OLS slope and intercept of Dan on mine, mean of each, mean
% difference (Dan - mine), RMS difference, median ratio Dan / mine
s.n = numel(x_mine);
R = corrcoef(x_mine, y_dan);
s.r = R(1, 2);
c = polyfit(x_mine, y_dan, 1);
s.slope = c(1);
s.intercept = c(2);
s.mean_mine = mean(x_mine);
s.mean_dan = mean(y_dan);
s.mean_diff = mean(y_dan - x_mine);
s.rms_diff = sqrt(mean((y_dan - x_mine).^2));
s.median_ratio = median(y_dan ./ x_mine);
end
