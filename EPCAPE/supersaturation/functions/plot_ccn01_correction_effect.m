function fig = plot_ccn01_correction_effect(D, fig_path, method_text)
% PLOT_CCN01_CORRECTION_EFFECT  Box plots of how the CCN-based retrievals change with Sourita's 0.1 % correction.
%
%   fig = plot_ccn01_correction_effect(D, fig_path, method_text)
%
%   method_text  how s_new was made, printed on the second title line, e.g.
%                '0.1 % setpoint scaled by gauss / measured'. Called once per
%                version of the correction (Part 11 draws two).
%
%   Left: absolute difference |s_new - s_old| per hour, in units of supersaturation
%   (%), for the four CCN-based retrievals. Right: the relative change
%   100 (s_new - s_old) / s_old, in % of s_old. s_old uses the posted Pier CCN;
%   s_new uses the posted CCN with the 0.1 % setpoint replaced by Sourita Saha's
%   Gaussian-corrected value, in one of two ways (Part 11 of the live script;
%   method_text says which). s_kappa and s_ARG00 are not
%   shown: they do not use the CCN counter, so their difference is zero by
%   construction.
%
%   Boxes follow MATLAB's boxchart: median line, 25th-75th percentile box,
%   whiskers to the furthest point within 1.5 IQR, outliers as dots. Colours follow Fig. 10 of Dedrick et al. (2026):
%   light grey = Pier low-cloud periods, dark grey = Mt. Soledad cloud events.
%   The median of each box is printed under its tick label.
%
%   Inputs
%     D         struct array, one element per box, with fields
%               label   (text under the box), period ('LowCloudPeriod' or 'SoledadCloudEvent'),
%               s_old, s_new (column vectors, %, same hours, finite only)
%     fig_path  where to save the PNG ('' = do not save)
%
%   Drawn INVISIBLE and saved with exportgraphics, as in plot_fig10_recreation.m
%   (exporting a visible figure has hung MATLAB R2024b). The caller closes it and
%   shows the saved PNG.

light = [0.65 0.65 0.65];   % Fig. 10 light grey (Pier low-cloud periods)
dark  = [0.24 0.24 0.24];   % Fig. 10 dark grey (Mt. Soledad cloud events)
ink   = [0.11 0.11 0.11];

fig = new_figure('Visible', 'off', 'Position', [100 100 1150 470]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

panels = {'abs', 'rel'};
for p = 1:2
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'off');
    pos = [1 2 4 5];   % two methods, two periods each, with a gap between methods
    for k = 1:numel(D)
        if strcmp(panels{p}, 'abs')
            y = abs(D(k).s_new - D(k).s_old);                         % supersaturation units (%)
        else
            y = 100 * (D(k).s_new - D(k).s_old) ./ D(k).s_old;        % % of s_old
        end
        if strcmp(D(k).period, 'LowCloudPeriod'), c = light; else, c = dark; end
        boxchart(ax, pos(k) * ones(size(y)), y, 'BoxFaceColor', c, 'MarkerColor', c, ...
                 'BoxWidth', 0.7, 'WhiskerLineColor', ink, 'MarkerStyle', '.');
        D(k).med(p) = median(y);   % kept for the label below
    end
    % tick labels: the retrieval on the first line, its median on the second
    % (inside the axes, median text collided with the panel title)
    if strcmp(panels{p}, 'abs'), fmt = 'median %.4f'; else, fmt = 'median %+.2f %%'; end
    tick_text = cell(1, numel(D));
    for k = 1:numel(D)
        tick_text{k} = [D(k).label '\newline' sprintf(fmt, D(k).med(p))];
    end
    xticks(ax, pos); xticklabels(ax, tick_text); xlim(ax, [0.2 5.8]);
    ax.TickLabelInterpreter = 'tex';
    if strcmp(panels{p}, 'abs')
        ylabel(ax, '|s_{new} - s_{old}|  (supersaturation, %)');
        title(ax, 'Absolute change', 'FontWeight', 'normal');
    else
        yline(ax, 0, '-', 'Color', [0.54 0.54 0.52]);
        ylabel(ax, '(s_{new} - s_{old}) / s_{old}  (%)');
        title(ax, 'Relative change', 'FontWeight', 'normal');
    end
end

% legend for the two greys (invisible patches)
ax1 = nexttile(tl, 1);
h1 = patch(ax1, NaN, NaN, light, 'EdgeColor', 'none', 'DisplayName', 'Pier low-cloud periods');
h2 = patch(ax1, NaN, NaN, dark, 'EdgeColor', 'none', 'DisplayName', 'Mt. Soledad cloud events');
legend(ax1, [h1 h2], 'Location', 'northeast');

title(tl, {'Effect of the Gaussian-corrected CCN at the 0.1 % setpoint on the CCN-based supersaturation retrievals', ...
           ['s_{old}: posted Pier CCN;  s_{new}: ' method_text]}, ...
      'FontSize', 11);

if ~isempty(fig_path)
    exportgraphics(fig, fig_path, 'Resolution', 200);
end

end
