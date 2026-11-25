function fig = plot_fig10_comparison(SS, MINE, fig_path, layout)
% PLOT_FIG10_COMPARISON  Fig. 10 of Dedrick et al. (2026) with my re-derived retrievals overlaid.
%
%   fig = plot_fig10_comparison(SS, MINE, fig_path)
%   fig = plot_fig10_comparison(SS, MINE, fig_path, layout)
%
%   Jeramy Dedrick's posted hourly retrievals are drawn as in plot_fig10_recreation.m
%   (same bins, greys, fonts, axis limits), and my re-derivations are drawn on top in
%   blue, so each panel compares the two. Code and settings are copied from
%   plot_fig10_recreation.m wherever they apply; if the style of that figure changes,
%   change it here too.
%
%   layout  '5panel' (default) or '9panel'
%
%   '5panel': the paper's layout.
%     (a)-(d): histograms of the 1-h retrievals (s_CCN, s_kappa, s_Nd, s_ARG00), each
%       with both periods. Posted: solid stairs, light grey = Pier low-cloud periods,
%       dark grey = Mt. Soledad cloud events (the paper's style). Mine: broken stairs
%       (line style mine_style), light blue = low-cloud periods, dark blue = Soledad
%       cloud events. Where two histograms coincide, the broken blue line runs on top
%       of the grey one, so both stay visible. The mean +/- std of each series is
%       printed in its own colour, posted first.
%     (e): the box-plot panel (below), in column 3 under the legend text.
%   '9panel': the same histograms split by period, so each panel compares only two
%     series (posted and mine, one period). A 2 x 4 grid of histograms: top row
%     (a)-(d) = low-cloud periods, bottom row (e)-(h) = Soledad cloud events, methods
%     in the paper's order (s_CCN, s_kappa, s_Nd, s_ARG00). The s_Nd panels are
%     labelled s_{Nd,VAP} (low-cloud) and s_{Nd,FM} (Soledad): the two versions the
%     paper combines in its panel (c). (i) is the box-plot panel, identical to (e) of
%     '5panel', in a 5th column under the legend text.
%   Box-plot panel: four boxes per retrieval: posted low-cloud, mine low-cloud, posted
%     Soledad, mine Soledad (top to bottom). It spans more rows than in
%     plot_fig10_recreation.m to make room for the 16 boxes. Outliers are not drawn
%     (whiskers still end at the furthest point within 1.5 IQR, as in MATLAB's boxchart).
%
%   Inputs
%     SS        timetable of Jeramy's posted retrievals with CORRECTED labels (variables
%               sCCN_*, skappa_*, sNd_*, sARG00_* for *_LowCloudPeriod and *_SoledadCloudEvent),
%               supersaturation in %
%     MINE      timetable with the same variables and rows as SS, holding my
%               re-derivation (%) at Jeramy's posted hours and NaN elsewhere. A variable
%               that is all NaN (not re-derived) gets no blue line or blue mean +/- std,
%               and "not re-derived" in its box-plot slot (and, in '9panel', in its
%               histogram panel).
%     fig_path  where to save the PNG ('' = do not save)
%     layout    '5panel' or '9panel' (optional; default '5panel', so older calls are unchanged)
%
%   Colors of my retrievals: blue was chosen because it contrasts with both greys
%   and with the red setpoint lines. Worst-case separation under simulated colour
%   blindness (Machado et al. 2009 matrices, OKLab distance x 100, my check) is ~12
%   between each blue and its grey, which is low for tritan vision; the broken line
%   style is the second cue that keeps the two apart.
%
%   As in plot_fig10_recreation.m, the figure is drawn INVISIBLE and saved with
%   exportgraphics (exporting a visible figure has hung MATLAB R2024b). The caller
%   closes it and shows the saved PNG.

if nargin < 4 || isempty(layout)
    layout = '5panel';
end
layout = validatestring(layout, {'5panel', '9panel'});

% ---- settings shared with plot_fig10_recreation.m (keep in step) -------------------
light = [0.65 0.65 0.65];          % the paper's light grey (posted, low-cloud periods)
dark  = [0.24 0.24 0.24];          % the paper's dark grey (posted, Soledad cloud events)
periods = {'LowCloudPeriod', 'SoledadCloudEvent'};
period_names = {'Low Cloud Periods', 'Soledad Cloud Events'};
colors  = {light, dark};
methods = {'sCCN', 'skappa', 'sNd', 'sARG00'};
titles  = {'s_{CCN}', 's_\kappa', 's_{Nd,VAP} & s_{Nd,FM}', 's_{ARG00}'};   % '5panel': both periods
titles_by_period = {'s_{CCN}', 's_\kappa', 's_{Nd,VAP}', 's_{ARG00}';    % '9panel': low-cloud row
                    's_{CCN}', 's_\kappa', 's_{Nd,FM}',  's_{ARG00}'};   % '9panel': Soledad row
ccn_based = [true false true false];
edges = 0:0.05:1.2;                % 0.05 % bins (see plot_fig10_recreation.m)

tick_fs  = 20;   % tick-label font size (points), all panels
label_fs = 20;   % axis-label font size (points), all panels
fnt_sz_mean_std = 18;
fnt_sz_panel_labels = 20;
fnt_sz_panelE_labels = 12;
ln_wdth_jeramy = 6;

% ---- settings for my retrievals (new in this figure) ---------------------------------
mine_light = [0.00 0.59 1.00];     % #0096FF, my low-cloud periods
mine_dark  = [0.00 0.16 0.78];     % #0028C8, my Soledad cloud events
mine_colors = {mine_light, mine_dark};
mine_style = ':';                  % dotted, so a blue line on top of an identical grey one shows both
mine_width = 2.2;
switch mine_style                  % the style's name, for the legend text
    case '--', mine_style_name = 'dashed';
    case ':',  mine_style_name = 'dotted';
    case '-.', mine_style_name = 'dash-dot';
    otherwise, mine_style_name = 'solid';
end
ink = [0.11 0.11 0.11];

% ---- layout --------------------------------------------------------------------------
% Both layouts use a 6-row grid; each histogram spans 3 rows. In the last column the
% legend text takes rows 1-2 and the box-plot panel takes rows 3-6 (4/6 of the height).
% Tile numbers run along rows: tile = n_col * (row - 1) + column.
% P lists the histogram panels: its tile, method k, the period(s) js it shows, its
% label, its y-label (empty = none), whether it gets an x-label, and the height y0 of
% its first mean +/- std line (normalized units).
P = struct('tile', {}, 'k', {}, 'js', {}, 'label', {}, 'ylab', {}, 'xlab', {}, 'y0', {});
switch layout
    case '5panel'
        % 6 x 3 grid: (a) (b) in rows 1-3, (c) (d) in rows 4-6, columns 1-2
        n_col = 3; fig_pos = [100 100 1250 1000];
        hist_tiles = [1 2 10 11];
        for k = 1:4
            if any(k == [1 3]), ylab = {'number of hours'}; else, ylab = {}; end
            if k == 3, y0 = 0.80; else, y0 = 0.90; end   % panel (c) starts lower, below its long title
            P(end + 1) = struct('tile', hist_tiles(k), 'k', k, 'js', [1 2], ...
                'label', sprintf('(%c) %s', 'a' + k - 1, titles{k}), ...
                'ylab', {ylab}, 'xlab', k >= 3, 'y0', y0);   %
        end
        box_letter = '(e)';
    case '9panel'
        % 6 x 5 grid: low-cloud row (a)-(d) in rows 1-3 and Soledad row (e)-(h) in rows
        % 4-6, columns 1-4; legend and box plots (i) in column 5. The figure is wider so
        % each histogram is about as wide as in '5panel' (~400 px).
        n_col = 5; fig_pos = [100 100 2100 1000];
        for j = 1:2
            for k = 1:4
                if k == 1
                    % first column: the period name (bold) above 'number of hours' names the row
                    ylab = {['\bf' period_names{j}], '\rmnumber of hours'};
                else
                    ylab = {};
                end
                P(end + 1) = struct('tile', n_col * 3 * (j - 1) + k, 'k', k, 'js', j, ...
                    'label', sprintf('(%c) %s', 'a' + 4 * (j - 1) + k - 1, titles_by_period{j, k}), ...
                    'ylab', {ylab}, 'xlab', j == 2, 'y0', 0.90);   %
            end
        end
        box_letter = '(i)';
end
% Larger text and tighter spacing in '9panel': its figure is 2100 px wide, so it is
% shown scaled down more than '5panel', which makes the same point sizes look smaller.
% (factor chosen by eye; 1 = the shared sizes above)
switch layout
    case '5panel', font_scale = 1.0; tile_spacing = 'compact';
    case '9panel', font_scale = 1.2; tile_spacing = 'tight';
end
tick_fs  = round(font_scale * tick_fs);
label_fs = round(font_scale * label_fs);
fnt_sz_mean_std      = round(font_scale * fnt_sz_mean_std);
fnt_sz_panel_labels  = round(font_scale * fnt_sz_panel_labels);
fnt_sz_panelE_labels = round(font_scale * fnt_sz_panelE_labels);
leg_fs  = 18;                      % legend text: not scaled, larger sizes run its two columns together
note_fs = round(font_scale * 12);  % source note under the legend
title_fs = round(font_scale * 20); % figure title

legend_tile = n_col;               % row 1, last column
box_tile = 3 * n_col;              % row 3, last column
span_hist = [3 1];                 % rows x columns spanned by each histogram

fig = new_figure('Visible', 'off', 'Position', fig_pos);   % invisible: see the note above
tl = tiledlayout(fig, 6, n_col, 'TileSpacing', tile_spacing, 'Padding', 'compact');

% ---- histograms -----------------------------------------------------------------------
for p = 1:numel(P)
    k = P(p).k;
    ax = nexttile(tl, P(p).tile, span_hist); hold(ax, 'on'); box(ax, 'off');
    % Tick-label size lives on each axes' rulers; setting it also resets the axis
    % label, so the rulers go first and the labels below override them.
    ax.XAxis.FontSize = tick_fs;
    ax.YAxis.FontSize = tick_fs;
    max_count = 0;                 % tallest bin drawn, to check the fixed y-limits below
    n_js = numel(P(p).js);         % 2 in '5panel', 1 in '9panel'

    for jj = 1:n_js
        j = P(p).js(jj);
        col = [methods{k} '_' periods{j}];

        % posted (Jeramy): as in plot_fig10_recreation.m
        v = SS.(col);
        v = min(v(isfinite(v)), 1.2);   % values above 1.2 % go into the last bin
        h = histogram(ax, v, 'BinEdges', edges, 'DisplayStyle', 'stairs', 'EdgeColor', colors{j}, ...
                      'LineWidth', ln_wdth_jeramy);
        max_count = max([max_count, h.Values]);
        % mean +/- std lines: posted series first (lines 1..n_js), then mine
        % (lines n_js+1..2 n_js). Lines stay at fixed positions even when a series is
        % missing, so a colour is always on the same line.
        text(ax, 0.75, P(p).y0 - 0.08 * (jj - 1), sprintf('%.2f±%.2f', mean(v), std(v)), 'Units', 'normalized', ...
             'HorizontalAlignment', 'right', 'Color', colors{j}, 'FontSize', fnt_sz_mean_std);

        % mine: same bins, broken blue line on top
        w = MINE.(col);
        w = min(w(isfinite(w)), 1.2);
        y_mine = P(p).y0 - 0.08 * (n_js + jj - 1);
        if ~isempty(w)
            h = histogram(ax, w, 'BinEdges', edges, 'DisplayStyle', 'stairs', 'EdgeColor', mine_colors{j}, ...
                          'LineWidth', mine_width, 'LineStyle', mine_style);
            max_count = max([max_count, h.Values]);
            text(ax, 0.75, y_mine, sprintf('%.2f±%.2f', mean(w), std(w)), 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', 'Color', mine_colors{j}, 'FontSize', fnt_sz_mean_std);
        elseif n_js == 1
            % '9panel' shows one period per panel, so say why there is no blue line
            text(ax, 0.75, y_mine, 'not re-derived', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'Color', mine_colors{j}, 'FontSize', fnt_sz_mean_std - 4, 'FontAngle', 'italic');
        end
    end

    if ccn_based(k)
        xline(ax, [0.1 1.0], ':', 'Color', [0.84 0.15 0.16], 'LineWidth', 2.5);
    end
    text(ax, 0.03, 0.96, P(p).label, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'FontSize', fnt_sz_panel_labels, 'Color', ink);

    % axis limits as in plot_fig10_recreation.m, the same in both layouts so the two
    % figures compare directly. They clip the posted dark-grey s_ARG00 peak (151 h in
    % the first bin) at 130 h, as the recreated figure does.
    xlim(ax, [0 1.2]);
    if k == 1 || k == 2
        ylim(ax, [0, 70]);
    else
        ylim(ax, [0, 130]);
    end

    % The S_ARG00 retrieval is much higher than what is plotted in Jeramy's
    % paper
    if p == 8
        ylim([0. 180])
    end

    if max_count > ax.YLim(2)
        fprintf('plot_fig10_comparison: panel %s, tallest bin (%d h) is above the y-limit (%d h)\n', ...
                P(p).label(1:3), max_count, ax.YLim(2));
    end

    if ~isempty(P(p).ylab)
        ylabel(ax, P(p).ylab, 'FontSize', label_fs);
    end
    if P(p).xlab
        xlabel(ax, 'supersaturation (%)', 'FontSize', label_fs);
    end
    xticks(ax, 0:0.2:1.2);
    grid(ax, 'on');
end

% ---- legend text (last column, rows 1-2): posted in the left column, mine in the right ----
axl = nexttile(tl, legend_tile, [2 1]); axis(axl, 'off');
text(axl, 0.02, 0.95, 'Jeramy, posted (solid)', 'Color', ink, 'FontSize', leg_fs, 'FontWeight', 'bold', 'Units', 'normalized');
text(axl, 0.02, 0.80, period_names{1}, 'Color', light, 'FontSize', leg_fs, 'FontWeight', 'bold', 'Units', 'normalized');
text(axl, 0.02, 0.67, period_names{2}, 'Color', dark, 'FontSize', leg_fs, 'FontWeight', 'bold', 'Units', 'normalized');
text(axl, 0.65, 0.95, ['mine (' mine_style_name ')'], 'Color', ink, 'FontSize', leg_fs, 'FontWeight', 'bold', 'Units', 'normalized');
text(axl, 0.65, 0.80, period_names{1}, 'Color', mine_light, 'FontSize', leg_fs, 'FontWeight', 'bold', 'Units', 'normalized');
text(axl, 0.65, 0.67, period_names{2}, 'Color', mine_dark, 'FontSize', leg_fs, 'FontWeight', 'bold', 'Units', 'normalized');

% Text describing data and a few caveats
% text(axl, 0.02, 0.30, {'Posted: UC San Diego Library, doi:10.6075/J0WM1DSC,', ...
%      'with the skappa / sCCN column labels corrected.', ...
%      'Mine: this script (Parts 5-8), at Jeramy''s posted hours.', ...
%      'Low-cloud s_{Nd} is not re-derived.'}, ...
%      'Color', [0.32 0.32 0.31], 'FontSize', note_fs, 'Units', 'normalized');



% ---- box plots (last column, rows 3-6): the same in both layouts --------------------------
% Four boxes per retrieval, one unit apart, with a gap of one unit between retrievals:
%   s_CCN 19-16, s_kappa 14-11, s_Nd 9-6, s_ARG00 4-1 (y positions, top to bottom).
% Within a retrieval: posted low-cloud, mine low-cloud, posted Soledad, mine Soledad, so
% each posted box sits directly above my box for the same period.
axe = nexttile(tl, box_tile, [4 1]); hold(axe, 'on');
for i = 1:4
    top = 20 - 5 * (i - 1) - 1;    % y of the first box of retrieval i
    for j = 1:2
        col = [methods{i} '_' periods{j}];
        pos_posted = top - 2 * (j - 1);   % slots 1 and 3
        pos_mine   = pos_posted - 1;      % slots 2 and 4

        v = SS.(col); v = v(isfinite(v));
        boxchart(axe, pos_posted * ones(size(v)), v, 'Orientation', 'horizontal', 'BoxFaceColor', colors{j}, ...
                 'BoxWidth', 0.8, 'WhiskerLineColor', [0 0 0], 'MarkerStyle', 'none');   % no outliers

        w = MINE.(col); w = w(isfinite(w));
        if ~isempty(w)
            boxchart(axe, pos_mine * ones(size(w)), w, 'Orientation', 'horizontal', 'BoxFaceColor', mine_colors{j}, ...
                     'BoxWidth', 0.8, 'WhiskerLineColor', [0 0 0], 'MarkerStyle', 'none');
        else
            text(axe, 0.02, pos_mine, 'not re-derived', 'Color', mine_colors{j}, 'FontSize', fnt_sz_panelE_labels, ...
                 'FontAngle', 'italic', 'VerticalAlignment', 'middle');
        end
    end
end
yticks(axe, [2.5 7.5 12.5 17.5]);   % centre of each group of four
yticklabels(axe, {'s_{ARG00}', 's_{Nd}', 's_\kappa', 's_{CCN}'});
axe.YAxisLocation = 'right';
axe.XAxis.FontSize = tick_fs;
axe.YAxis.FontSize = label_fs;      % the method names act as labels here, so label size

% group separators, as in the paper: Hoppel minimum (s_CCN, s_kappa), cloud scene (s_Nd), parcel (s_ARG00)
groups = {20.0, 'Hoppel Minimum'; 10.0, 'Cloud-Scene'; 5.0, 'Parcel'};
for g = 1:3
    yline(axe, groups{g, 1}, ':', 'Color', [0.54 0.54 0.52]);
    text(axe, 1.18, groups{g, 1} - 0.1, groups{g, 2}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'FontSize', fnt_sz_panelE_labels, 'Color', [0.32 0.32 0.31]);
end
xlim(axe, [0 1.2]);
ylim(axe, [0 20.8]);

xlabel(axe, 'supersaturation (%)', 'FontSize', label_fs);   % after the ruler, so it is not reset
text(axe, 0.02, 0.995, box_letter, 'Units', 'normalized', 'VerticalAlignment', 'top', 'Color', ink, ...
     'FontSize', fnt_sz_panel_labels);
grid(axe, 'on');

title(tl, 'Recreation of Fig. 10 from Dedrick et al. (2026): Jeramy''s posted 1-h retrievals (grey) and mine (blue)',...
    'FontSize', title_fs);
if ~isempty(fig_path)
    exportgraphics(fig, fig_path, 'Resolution', 200);
end

end
