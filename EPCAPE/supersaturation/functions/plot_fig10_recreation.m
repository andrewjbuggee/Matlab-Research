function fig = plot_fig10_recreation(SS, fig_path)
% PLOT_FIG10_RECREATION  Recreate Fig. 10 of Dedrick et al. (2026) from the posted hourly retrievals.
%
%   fig = plot_fig10_recreation(SS, fig_path)
%
%   Panels (a)-(d): histograms of the 1-h supersaturation retrievals (s_CCN,
%   s_kappa, s_Nd, s_ARG00); light grey = Scripps Pier low-cloud periods, dark grey =
%   Mt. Soledad cloud events; dotted red lines = the CCN counter's lowest and highest
%   setpoints (0.1 and 1.0 %), drawn on the two CCN-based panels as in the paper.
%   Panel (e): box plots of all eight.
%
%   Inputs
%     SS        timetable of the retrievals with CORRECTED labels (variables
%               sCCN_*, skappa_*, sNd_*, sARG00_* for *_LowCloudPeriod and *_SoledadCloudEvent)
%     fig_path  where to save the PNG ('' = do not save)
%
%   Bins are 0.05 % wide: they reproduce the paper's light-grey peak heights in
%   (a)-(c) (47, 62, 78 h vs ~47, 62, 80 h read off the paper). The paper's dark-grey
%   peaks in (c)-(d) suggest MATLAB's automatic, per-series bin widths, not reproduced.
%
%   The figure is drawn INVISIBLE ('Visible', 'off') and saved with exportgraphics.
%   Exporting a VISIBLE figure has hung MATLAB R2024b for > 10 min (in -batch runs and
%   in the Live Editor; EPCAPE BAMS Fig. 7 and this script, Oct 2026), while the same
%   figure drawn invisibly saves in ~1 s. The caller closes it (close(fig)) and, to see
%   it, displays the saved PNG. It uses factory graphics defaults (new_figure), so
%   startup.m's dark theme does not apply.

light = [0.65 0.65 0.65];          % the paper's light grey
dark  = [0.24 0.24 0.24];          % the paper's dark grey
periods = {'LowCloudPeriod', 'SoledadCloudEvent'};
colors  = {light, dark};
methods = {'sCCN', 'skappa', 'sNd', 'sARG00'};
titles  = {'(a) s_{CCN}', '(b) s_\kappa', '(c) s_{Nd,VAP} & s_{Nd,FM}', '(d) s_{ARG00}'};
ccn_based = [true false true false];
edges = 0:0.05:1.2;

fig = new_figure('Visible', 'off', 'Position', [100 100 1250 760]);   % invisible: see the note above
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tiles = [1 2 4 5];
for k = 1:4
    ax = nexttile(tl, tiles(k)); hold(ax, 'on'); box(ax, 'off');
    for j = 1:2
        v = SS.([methods{k} '_' periods{j}]);
        v = min(v(isfinite(v)), 1.2);   % values above 1.2 % go into the last bin
        histogram(ax, v, 'BinEdges', edges, 'DisplayStyle', 'stairs', 'EdgeColor', colors{j}, 'LineWidth', 2.2);
        text(ax, 0.97, 0.90 - 0.08 * (j - 1), sprintf('%.2f±%.2f', mean(v), std(v)), 'Units', 'normalized', ...
             'HorizontalAlignment', 'right', 'Color', colors{j}, 'FontSize', 11);
    end
    if ccn_based(k)
        xline(ax, [0.1 1.0], ':', 'Color', [0.84 0.15 0.16], 'LineWidth', 1.2);
    end
    text(ax, 0.03, 0.96, titles{k}, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', [0.11 0.11 0.11]);
    xlim(ax, [0 1.2]);
    if any(k == [1 3]), ylabel(ax, 'number of hours'); end
    if k >= 3, xlabel(ax, 'supersaturation (%)'); end
    grid(ax, 'on');
end

% legend text, as in the paper
axl = nexttile(tl, 3); axis(axl, 'off');
text(axl, 0.05, 0.80, 'Low Cloud Periods', 'Color', light, 'FontSize', 15, 'FontWeight', 'bold', 'Units', 'normalized');
text(axl, 0.05, 0.66, 'Soledad Cloud Events', 'Color', dark, 'FontSize', 15, 'FontWeight', 'bold', 'Units', 'normalized');
text(axl, 0.05, 0.38, {'Recreated from the posted hourly retrievals', '(UC San Diego Library, doi:10.6075/J0WM1DSC),', ...
     'with the skappa / sCCN column labels corrected.'}, 'Color', [0.32 0.32 0.31], 'FontSize', 9, 'Units', 'normalized');

% (e) box plots: boxes at y = 12, 11 (s_CCN), 9, 8 (s_kappa), 6, 5 (s_Nd), 3, 2 (s_ARG00); light above dark
axe = nexttile(tl, 6); hold(axe, 'on');
for i = 1:4
    for j = 1:2
        v = SS.([methods{i} '_' periods{j}]);
        v = v(isfinite(v));
        pos = 12 - 3 * (i - 1) - (j - 1);
        boxchart(axe, pos * ones(size(v)), v, 'Orientation', 'horizontal', 'BoxFaceColor', colors{j}, ...
                 'MarkerColor', colors{j}, 'BoxWidth', 0.8, 'WhiskerLineColor', [0 0 0]);
    end
end
yticks(axe, [2.5 5.5 8.5 11.5]);
yticklabels(axe, {'s_{ARG00}', 's_{Nd}', 's_\kappa', 's_{CCN}'});
axe.YAxisLocation = 'right';
groups = {13.0, 'Hoppel Minimum'; 7.0, 'Cloud-Scene'; 4.0, 'Parcel'};
for g = 1:3
    yline(axe, groups{g, 1}, ':', 'Color', [0.54 0.54 0.52]);
    text(axe, 1.18, groups{g, 1} - 0.1, groups{g, 2}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'FontSize', 9, 'Color', [0.32 0.32 0.31]);
end
xlim(axe, [0 1.2]); ylim(axe, [1 13.6]);
xlabel(axe, 'supersaturation (%)');
text(axe, 0.02, 0.98, '(e)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'Color', [0.11 0.11 0.11]);
grid(axe, 'on');

title(tl, 'Recreation of Dedrick et al. (2026) Fig. 10: 1-h retrievals of cloud supersaturation');
if ~isempty(fig_path)
    exportgraphics(fig, fig_path, 'Resolution', 200);
end

end
