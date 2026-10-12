%% plot_missing_EF.m
% Shows that, in some periods, Abbey Williams' 15-min file (EPCAPE_15mindata.nc) flags the
% GCVI as ON while the 1-s GCVI record (msd_cvi_ef_calc.csv) is empty. In those periods her
% enhancement-factor (EF) code uses EF = 1, so the AMS "residual" concentrations were never
% divided by the enhancement factor.
%
% ONE FIGURE PER TIME WINDOW, 4 panels on a shared time axis (UTC):
%   a) CVIon_class from the 15-min file (left axis: off / changing / on) and the 1-s EF from
%      msd_cvi_ef_calc.csv (right axis). Where the flag is "on" but there are no EF dots, the
%      GCVI logged nothing.
%   b) Abbey's EF code (from her email of 8 Oct 2026), re-run here: the EF used for each 15-min
%      interval. It is the mean of the 1-s EF, or EF = 1 when the interval has no 1-s record
%      at all (her comment: "set EF to 1 (no enhancement) for times when CVI is off").
%   c) The AMS mass of the four species in the 15-min file, i.e. the values Kavin averaged.
%      Dashed lines show the red-shaded intervals divided by a typical EF (the median of the
%      intervals that do have an EF), i.e. roughly what the residuals should have been.
%   d) Proof that the EF was or was not applied: the 15-min file's AMS_mass divided by the
%      CE-corrected AMS from the UC San Diego Library file, averaged over the same 15 min.
%      About 1/EF (~0.2) where the EF was applied; exactly 1 where it was not (and, as it should
%      be, 1 whenever the CVI was off and the AMS sampled ambient air).
%   Shading in every panel: blue = CVI on and a 1-s EF record exists;
%                           red  = CVI flagged on, but no 1-s EF record (so EF = 1 was used).
%
% WHY THESE WINDOWS (see "cases" in Section 1)
%   * 18-20 March 2023 (overview) and 20 March 2023 (zoom). The 1-s EF record stops on 18 Mar
%     at 14:26 and only resumes on 20 Mar at 16:34. On 18 March the 15-min file flags the CVI
%     OFF from 14:30 on, so EF = 1 is right there. The problem appears on 20 March,
%     08:00-14:30, when the flag turns ON inside the same gap (12 haze intervals with AMS data).
%   * 30 April 2023. The 1-s EF record stops at 00:31 and resumes at 16:32. The 15-min file
%     flags the CVI ON from 08:30 to 16:30 (27 cloud + 4 haze intervals with AMS data). The AMS
%     drops from the ambient level (~0.7 ug m-3 sulfate) to ~0.03 when it switches to the CVI,
%     then climbs back to ambient-like levels: undivided CVI concentrations.
%   The same check over the whole campaign finds 77 CVI-on intervals with no EF record
%   (20 Mar, 30 Apr, 6 Nov, 27 Dec 2023); it is printed in the Command Window.
%
% DATA FILES (Section 1)
%   EPCAPE_15mindata.nc        Abbey's 15-min merged file (from Kavin): CVIon_class,
%                              ClrCldHaz_class, AMS_mass (5 species, ug m-3), 15-min intervals
%                              labelled by their START time.
%   msd_cvi_ef_calc.csv        Michael's 1-s GCVI file, attached to Abbey's email. It is
%                              byte-for-byte identical (same MD5, e1344be9...) to the UC San
%                              Diego Library copy, bb2743661p_36_1_CVI_EF.csv.
%   AMS_CE_corrected_msptof.nc the CE-corrected AMS from the UC San Diego Library ("Aerosol
%                              Quantification Products"), ~5.5-min samples, NOT divided by EF.
%
% Uses read_netcdf_time.m and read_15min_file.m from the readers/ folder of this project.
% Written for MATLAB R2024b and later (needs xregion, R2023a+); base MATLAB only.
% Andrew Buggee, 2026-10-12.

clear; clc;

%% 1. Settings: EDIT THESE PATHS ------------------------------------------------------------
kavin_folder = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper";
file_15min   = fullfile(kavin_folder, "EPCAPE_15mindata.nc");
file_cvi_ef  = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed/derived/seasonal_averages_check/msd_cvi_ef_calc.csv";
file_ams_ce  = "/Users/andrewbuggee/Documents/VS_CODE/Python-Research/EPCAPE/data/processed/Aerosol Quantification Products/AMS_CE_corrected_msptof.nc";
output_folder = fullfile(fileparts(mfilename('fullpath')), 'output', 'figures');

% Time windows to plot (UTC). Add or edit rows: {short tag, title text, start, end}
cases = {
    "18to20Mar2023", "18-20 March 2023", datetime(2023, 3, 18, 8, 0, 0), datetime(2023, 3, 20, 20, 0, 0)
    "20Mar2023",     "20 March 2023",    datetime(2023, 3, 20, 5, 0, 0), datetime(2023, 3, 20, 18, 0, 0)
    "30Apr2023",     "30 April 2023",    datetime(2023, 4, 29, 22, 0, 0), datetime(2023, 4, 30, 20, 0, 0)
};

% Colors. Species = the first four categorical slots of the reference palette (fixed order);
% states use neutral gray, the series blue, and the "critical" status red.
ink        = [0.043 0.043 0.043];          % #0b0b0b  text and single-series marks
ink_2      = [0.322 0.318 0.306];          % #52514e  secondary text
gray_mark  = [0.62 0.62 0.60];             % 1-s EF dots behind the 15-min EF
species_color = [0.165 0.471 0.839;        % #2a78d6  sulfate
                 0.922 0.408 0.204;        % #eb6834  organics
                 0.106 0.686 0.478;        % #1baf7a  nitrate
                 0.929 0.631 0.000];       % #eda100  ammonium
blue_state = [0.165 0.471 0.839];          % CVI on, EF record exists
red_state  = [0.816 0.231 0.231];          % #d03b3b  CVI flagged on, no EF record
off_state  = [0.55 0.55 0.53];             % CVI off (ambient)

% Make the readers in this project visible (read_netcdf_time, read_15min_file)
code_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(code_folder, 'readers'), fullfile(code_folder, 'helpers'));
if ~isfolder(output_folder), mkdir(output_folder); end

%% 2. Abbey's 15-min file -----------------------------------------------------------------------
% F.time (interval start), F.class (ClrCldHaz_class), F.ams_ugm3 (intervals x 5 species),
% F.species = nitrate, sulfate, organics, ammonium, chloride (checked by the reader)
F = read_15min_file(file_15min);
cvi_on = double(ncread(file_15min, 'CVIon_class'));   % 1 = on, 0 = off, 0.5 = changing state
t15 = F.time;                                          % 15-min interval STARTS
n15 = numel(t15);
plot_species = ["sulfate", "organics", "nitrate", "ammonium"];
[~, col] = ismember(plot_species, F.species);         % columns of F.ams_ugm3 for these species
ams15 = F.ams_ugm3(:, col);                            % (intervals x 4), ug m-3

%% 3. The 1-s EF record, and Abbey's EF code re-run ---------------------------------------------
% Read only the time and EF columns of Michael's 1-s file (2.5 million rows; ~20 s).
fprintf('Reading the 1-s GCVI file (about 20 s) ...\n');
opts = detectImportOptions(file_cvi_ef, 'CommentStyle', '#', 'Delimiter', ',');
opts.SelectedVariableNames = {'msd__cvi_v1__datetime', 'msd__cvi_v1__EF_cnttemp_round'};
opts = setvartype(opts, 'msd__cvi_v1__datetime', 'datetime');
opts = setvaropts(opts, 'msd__cvi_v1__datetime', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
opts = setvartype(opts, 'msd__cvi_v1__EF_cnttemp_round', 'double');
G = readtable(file_cvi_ef, opts);
t_1s  = G.msd__cvi_v1__datetime;
ef_1s = G.msd__cvi_v1__EF_cnttemp_round;

% --- Abbey's code, step by step (her lines quoted as comments) ---
%   cvi_EF(isinf(cvi_EF)) = NaN;  cvi_EF(cvi_EF>50) = NaN;  cvi_EF(cvi_EF<0) = NaN;  cvi_EF(cvi_EF==0) = NaN;
ef_clean = ef_1s;
ef_clean(isinf(ef_clean) | ef_clean > 50 | ef_clean < 0 | ef_clean == 0) = NaN;
%   cviEF_mean = mean(cvi_EF,'omitnan');  cviEF_SD = std(cvi_EF,'omitnan');
%   cviEF_outlier_class = 1 where cvi_EF > mean + 3 SD or < mean - 3 SD
ef_mean = mean(ef_clean, 'omitnan');
ef_sd   = std(ef_clean, 'omitnan');
outlier = ef_clean > ef_mean + 3 * ef_sd | ef_clean < ef_mean - 3 * ef_sd;
%   for i = 1:length(EPCAPE_15min)-1
%       ind = find(cvi_time >= EPCAPE_15min(i) & cvi_time < EPCAPE_15min(i+1));
%       if ~isempty(ind): cviEF_15min(i) = mean(cvi_EF(ind),'omitnan'); flag if any outlier
%       else:             cviEF_15min(i) = 1;   % "set EF to 1 (no enhancement) for times when CVI is off"
% Her loop runs find() 35,041 times over 2.5 million rows. The same result, much faster:
% give every 1-s sample the number of the 15-min interval it falls in, then sum per interval.
k = floor((t_1s - t15(1)) / minutes(15)) + 1;          % interval number of each 1-s sample
use = k >= 1 & k <= n15 - 1;                           % her loop stops at the second-to-last interval
n_records = accumarray(k(use), 1, [n15 1]);            % 1-s rows in each interval (any value)
good = use & ~isnan(ef_clean);
ef_sum = accumarray(k(good), ef_clean(good), [n15 1]);
n_good = accumarray(k(good), 1, [n15 1]);
ef15 = ef_sum ./ n_good;                               % mean(...,'omitnan'); NaN if every value is NaN
ef15(n_records == 0) = 1;                              % her rule: no record -> EF = 1 ("CVI off")
ef15(n15) = NaN;                                       % never set by her loop
ef15_flag = accumarray(k(use), double(outlier(use)), [n15 1]) > 0;   % only feeds cviEF_15min_noflags

% The two states that matter, per 15-min interval
on_with_ef = cvi_on == 1 & n_records > 0;
on_no_ef   = cvi_on == 1 & n_records == 0;
ef_typical = median(ef15(on_with_ef), 'omitnan');      % ~4.96, used for the dashed lines in panel c

%% 4. CE-corrected AMS (Library file) averaged over the same 15-min intervals --------------------
t_ams = read_netcdf_time(file_ams_ce, 'time');
k_ams = floor((t_ams - t15(1)) / minutes(15)) + 1;
ce15 = NaN(n15, numel(plot_species));
for j = 1:numel(plot_species)
    v = double(ncread(file_ams_ce, plot_species(j)));
    ok = k_ams >= 1 & k_ams <= n15 & isfinite(v);
    ce15(:, j) = accumarray(k_ams(ok), v(ok), [n15 1]) ./ accumarray(k_ams(ok), 1, [n15 1]);
end
% Panel d: one ratio per interval (the EF divides every species alike, so take the median)
ratio = median(ams15 ./ ce15, 2, 'omitnan');

%% 5. Campaign-wide summary in the Command Window ---------------------------------------------------
fprintf('\nEF statistics used by Abbey''s outlier test: mean %.2f, SD %.2f\n', ef_mean, ef_sd);
fprintf(['15-min intervals with an outlier-flagged EF: %d; with the CVI on and AMS data: %d ' ...
         '(flags only enter her cviEF_15min_noflags)\n'], sum(ef15_flag), ...
        sum(ef15_flag & cvi_on == 1 & any(isfinite(ams15), 2)));
fprintf('Median EF of the CVI-on intervals that have a 1-s record: %.2f\n', ef_typical);
fprintf('\nCVI flagged ON (CVIon_class = 1) but NO 1-s EF record, so EF = 1 was used: %d intervals (%.2f h)\n', ...
    sum(on_no_ef), sum(on_no_ef) / 4);
days_no_ef = unique(dateshift(t15(on_no_ef), 'start', 'day'));
for d = days_no_ef(:)'
    in_day = on_no_ef & dateshift(t15, 'start', 'day') == d;
    has_ams = in_day & any(isfinite(ams15), 2);
    fprintf('   %s: %2d intervals (cloud %2d, haze %2d), %2d with AMS data; their AMS_mass / CE-corrected AMS = %s\n', ...
        string(d, 'dd MMM yyyy'), sum(in_day), sum(in_day & F.class == 1), sum(in_day & F.class == 2), sum(has_ams), ...
        ratio_text(ratio(has_ams)));
end
fprintf('   (CVI-on intervals WITH an EF record: AMS_mass / CE-corrected AMS has median %.3f = 1/%.2f)\n', ...
    median(ratio(on_with_ef), 'omitnan'), 1 / median(ratio(on_with_ef), 'omitnan'));

%% 6. One figure per time window ---------------------------------------------------------------------
for c = 1:size(cases, 1)
    [tag, title_text, t0, t1] = cases{c, :};
    win15 = t15 >= t0 - minutes(15) & t15 <= t1;        % intervals overlapping the window
    win1s = t_1s >= t0 & t_1s <= t1;
    tw = t15(win15);
    t_steps = [tw; tw(end) + minutes(15)];              % stairs need the end of the last interval

    fig = figure('Color', 'w', 'Position', [80 60 1100 1150]);
    if ~isMATLABReleaseOlderThan("R2025a")
        theme(fig, "light");                            % R2025a+ may default to the dark theme
    end
    tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: the CVI is flagged on, but the GCVI enhancement factor is missing', title_text), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', ink);
    axs = gobjects(4, 1);

    % --- a) CVI flag (left) and 1-s EF (right) ---
    axs(1) = nexttile(tl); ax = axs(1); hold(ax, 'on');
    shade_states(ax, t15, on_with_ef & win15, on_no_ef & win15, blue_state, red_state, true);
    yyaxis(ax, 'left');
    stairs(ax, t_steps, [cvi_on(win15); cvi_on(find(win15, 1, 'last'))], '-', 'Color', ink, 'LineWidth', 2, ...
        'DisplayName', 'CVIon\_class (15-min file)');
    ylim(ax, [-0.15 1.9]); yticks(ax, [0 0.5 1]); yticklabels(ax, {'off', 'changing', 'on'});  % headroom for the legend
    ylabel(ax, 'CVIon_class (15-min file)', 'Interpreter', 'none', 'Color', ink);
    ax.YAxis(1).Color = ink;
    yyaxis(ax, 'right');
    plot(ax, t_1s(win1s), ef_1s(win1s), '.', 'Color', blue_state, 'MarkerSize', 8, ...
        'DisplayName', '1-s EF (msd\_cvi\_ef\_calc.csv)');
    ylim(ax, [0 12]); ylabel(ax, 'EF, 1-s record', 'Color', blue_state);
    ax.YAxis(2).Color = blue_state;
    mark_record_gaps(ax, t_1s, t0, t1, t15(on_no_ef), ink_2);
    legend(ax, 'Location', 'north', 'Box', 'off', 'FontSize', 9, 'TextColor', ink, 'NumColumns', 2, 'AutoUpdate', 'off');
    title(ax, 'a) EPCAPE-15mindata.nc flags the CVI; the 1-s GCVI file logs the enhancement factor', ...
        'FontWeight', 'normal', 'Color', ink);

    % --- b) Abbey's code: EF used for each 15-min interval ---
    axs(2) = nexttile(tl); ax = axs(2); hold(ax, 'on');
    shade_states(ax, t15, on_with_ef & win15, on_no_ef & win15, blue_state, red_state, false);
    plot(ax, t_1s(win1s), ef_1s(win1s), '.', 'Color', gray_mark, 'MarkerSize', 8, 'DisplayName', '1-s EF');
    stairs(ax, t_steps, [ef15(win15); ef15(find(win15, 1, 'last'))], '-', 'Color', blue_state, 'LineWidth', 2, ...
        'DisplayName', 'EF used by Abbey''s code (15-min)');
    yline(ax, 1, ':', 'EF = 1: no 1-s record ("CVI off" in her code)', 'Color', red_state, 'LineWidth', 1.5, ...
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize',...
        9, 'HandleVisibility', 'off');
    ylim(ax, [0 9]); ylabel(ax, 'Enhancement factor', 'Color', ink);
    legend(ax, 'Location', 'north', 'Box', 'off', 'FontSize', 9, 'TextColor', ink, 'NumColumns', 2, 'AutoUpdate', 'off');
    title(ax, 'b) Abbey''s EF code re-run: the 15-min mean of the 1-s EF, or EF = 1 when there is no record', ...
        'FontWeight', 'normal', 'Color', ink);

    % --- c) AMS mass in the 15-min file (what Kavin averaged), and what it would be divided by EF ---
    axs(3) = nexttile(tl); ax = axs(3); hold(ax, 'on');
    shade_states(ax, t15, on_with_ef & win15, on_no_ef & win15, blue_state, red_state, false);
    for j = 1:numel(plot_species)
        y = ams15(win15, j);
        stairs(ax, t_steps, [y; y(end)], '-', 'Color', species_color(j, :), 'LineWidth', 2, ...
            'DisplayName', char(plot_species(j)));
        y_div = NaN(size(y));
        red_here = on_no_ef(win15);
        y_div(red_here) = y(red_here) / ef_typical;     % what the residual would be with a typical EF
        stairs(ax, t_steps, [y_div; y_div(end)], '--', 'Color', species_color(j, :), 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
    plot(ax, NaT, NaN, '--', 'Color', ink_2, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('red intervals divided by a typical EF (%.2f)', ef_typical));
    ylabel(ax, 'AMS mass (\mug m^{-3})', 'Color', ink);
    ylim(ax, [0 1.3 * max(ams15(win15, :), [], 'all')]);  % headroom for the legend
    legend(ax, 'Location', 'north', 'Box', 'off', 'FontSize', 9, 'NumColumns', 5, 'TextColor', ink, 'AutoUpdate', 'off');
    title(ax, 'c) AMS mass in the 15-min file (the values Kavin averaged); behind the CVI it should be the residual ÷ EF', ...
        'FontWeight', 'normal', 'Color', ink);

    % --- d) Was the EF applied? 15-min AMS_mass / CE-corrected AMS ---
    axs(4) = nexttile(tl); ax = axs(4); hold(ax, 'on');
    shade_states(ax, t15, on_with_ef & win15, on_no_ef & win15, blue_state, red_state, false);
    yline(ax, 1, '-', 'not divided by EF', 'Color', ink_2, 'LabelHorizontalAlignment', 'right', 'FontSize', 9, ...
        'HandleVisibility', 'off');
    yline(ax, 1 / ef_typical, '--', sprintf('divided by EF (1/%.2f)', ef_typical), 'Color', ink_2, ...
        'LabelHorizontalAlignment', 'right', 'FontSize', 9, 'HandleVisibility', 'off');
    tc = tw + minutes(7.5);                             % plot each interval at its middle
    r = ratio(win15);
    s_off = cvi_on(win15) == 0; s_ef = on_with_ef(win15); s_noef = on_no_ef(win15);
    plot(ax, tc(s_off),  r(s_off),  'o', 'MarkerSize', 6, 'MarkerFaceColor', off_state, 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'CVI off (ambient air)');
    plot(ax, tc(s_ef),   r(s_ef),   'o', 'MarkerSize', 8, 'MarkerFaceColor', blue_state, 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'CVI on, EF record');
    plot(ax, tc(s_noef), r(s_noef), 'o', 'MarkerSize', 8, 'MarkerFaceColor', red_state, 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'CVI on, no EF record');
    ylim(ax, [0 1.75]); ylabel(ax, {'15-min AMS ÷', 'CE-corrected AMS'}, 'Color', ink);
    legend(ax, 'Location', 'north', 'Box', 'off', 'FontSize', 9, 'TextColor', ink, 'NumColumns', 3, 'AutoUpdate', 'off');
    title(ax, 'd) Was the EF applied? Ratio of the 15-min file''s AMS to the Library''s CE-corrected AMS (median of 4 species)', ...
        'FontWeight', 'normal', 'Color', ink);
    xlabel(ax, 'Time (UTC)', 'Color', ink);

    % --- shared time axis ---
    linkaxes(axs, 'x');
    for a = axs'
        xlim(a, [t0 t1]);
        a.XAxis.TickLabelFormat = 'dd-MMM HH:mm';
        a.Color = 'w';                                  % white plot area (batch mode/dark theme default to dark)
        grid(a, 'on'); a.GridColor = [0.15 0.15 0.15]; a.GridAlpha = 0.12; a.Box = 'off';
        a.FontSize = 10; a.XColor = ink_2;
        a.TitleHorizontalAlignment = 'left';
        a.Toolbar.Visible = 'off';                      % keep the axes toolbar out of the exported image
    end
    for a = axs(2:4)'
        a.YColor = ink_2;                               % panel a keeps its two coloured y-axes
    end
    for a = axs(1:3)'
        a.XTickLabel = {};
    end

    png = fullfile(output_folder, "missing_EF_" + tag + ".png");
    exportgraphics(fig, png, 'Resolution', 200);
    exportgraphics(fig, replace(png, ".png", ".pdf"), 'ContentType', 'image', 'Resolution', 300);
    fprintf('Saved %s (and .pdf)\n', png);
end

%% Local functions -------------------------------------------------------------------------------
function shade_states(ax, t15, on_with_ef, on_no_ef, blue_state, red_state, with_legend)
% SHADE_STATES  Shade the 15-min intervals with the CVI on: blue where an EF record exists,
% red where there is none. Consecutive intervals are merged into one band.
for s = 1:2
    if s == 1, mask = on_with_ef; color = blue_state; alpha = 0.10; name = 'CVI on, EF record';
    else,      mask = on_no_ef;   color = red_state;  alpha = 0.18; name = 'CVI flagged on, no EF record (EF = 1 used)';
    end
    [starts, ends] = block_edges(t15, mask);
    if isempty(starts), continue; end
    h = xregion(ax, starts, ends, 'FaceColor', color, 'FaceAlpha', alpha, 'EdgeColor', 'none');  % one band per block
    set(h, 'HandleVisibility', 'off');
    if with_legend
        h(1).HandleVisibility = 'on';                   % one legend entry per state
        h(1).DisplayName = name;
    end
end
end

function [starts, ends] = block_edges(t15, mask)
% BLOCK_EDGES  Start and end times of each run of consecutive true 15-min intervals.
mask = mask(:);
d = diff([false; mask; false]);
starts = t15(d(1:end-1) == 1);
ends   = t15(find(d(2:end) == -1)) + minutes(15);
end

function mark_record_gaps(ax, t_1s, t0, t1, t_on_no_ef, color)
% MARK_RECORD_GAPS  Dashed vertical lines where the 1-s EF record stops and resumes, for the gaps
% (> 1 h) that contain at least one interval with the CVI flagged on (t_on_no_ef = their starts).
gap = find(diff(t_1s) > hours(1));
for g = gap(:)'
    stop = t_1s(g); resume = t_1s(g + 1);
    if ~any(t_on_no_ef >= stop & t_on_no_ef < resume)
        continue                                         % an ordinary CVI-off gap: not marked
    end
    if stop >= t0 && stop <= t1
        xline(ax, stop, '--', "1-s record stops " + string(stop, 'dd MMM HH:mm'), 'Color', color, ...
            'LabelVerticalAlignment', 'top', 'FontSize', 9, 'HandleVisibility', 'off');
    end
    if resume >= t0 && resume <= t1
        xline(ax, resume, '--', "resumes " + string(resume, 'dd MMM HH:mm'), 'Color', color, ...
            'LabelVerticalAlignment', 'top', 'FontSize', 9, 'HandleVisibility', 'off');
    end
    if stop < t0 && resume > t0                         % the gap began before the window
        text(ax, 0.005, 0.72, "no 1-s EF record since " + string(stop, 'dd MMM HH:mm'), 'Units', 'normalized', ...
            'Color', color, 'FontSize', 9, 'HorizontalAlignment', 'left');
    end
    if resume > t1 && stop < t1                         % the gap ends after the window
        text(ax, 0.995, 0.72, "no 1-s EF record until " + string(resume, 'dd MMM HH:mm'), 'Units', 'normalized', ...
            'Color', color, 'FontSize', 9, 'HorizontalAlignment', 'right');
    end
end
end

function txt = ratio_text(r)
% RATIO_TEXT  "min-max" of the ratios, or "no AMS data".
r = r(isfinite(r));
if isempty(r)
    txt = "no AMS data";
else
    txt = sprintf('%.3f-%.3f', min(r), max(r));
end
end
