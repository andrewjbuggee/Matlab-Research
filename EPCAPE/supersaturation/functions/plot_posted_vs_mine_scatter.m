function fig = plot_posted_vs_mine_scatter(SS, X, fig_path, x_text, notes)
% PLOT_POSTED_VS_MINE_SCATTER  Jeramy's posted retrievals against mine, hour by hour, in a 2 x 4 grid of scatter plots.
%
%   fig = plot_posted_vs_mine_scatter(SS, X, fig_path, x_text, notes)
%
%   One panel per method and period, laid out as the '9panel' version of
%   plot_fig10_comparison.m: top row (a)-(d) = Pier low-cloud periods, bottom row
%   (e)-(h) = Mt. Soledad cloud events; columns s_CCN, s_kappa, s_Nd, s_ARG00.
%   Each marker is ONE hour: Jeramy's posted value of that method and period (y)
%   against my value of the same method, period and hour (x). Rows of SS and X are
%   the same hours, so the pairing is by row; hours where either is missing are
%   left out.
%
%   Each panel shows the 1:1 line (dashed) and the ordinary least-squares fit of
%   y on x (solid), with the fit's slope, intercept, R^2 and the number of hours.
%   R^2 is the coefficient of determination of that fit, which for a straight line
%   fitted by least squares equals the squared Pearson correlation (a standard
%   identity).
%
%   Inputs
%     SS        timetable of the posted retrievals with CORRECTED labels (Part 2),
%               variables sCCN_*, skappa_*, sNd_*, sARG00_* for *_LowCloudPeriod and
%               *_SoledadCloudEvent, %
%     X         timetable with the same rows and variables: my retrievals, %
%               (NaN where there is none)
%     fig_path  where to save the PNG ('' = do not save)
%     x_text    what x is, printed on the second title line
%     notes     struct of short notes printed in a panel, field = variable name
%               (e.g. notes.sNd_LowCloudPeriod = '...'); [] for none
%
%   Drawn INVISIBLE and saved with exportgraphics, like the other figures of the
%   live script (exporting a visible figure has hung MATLAB R2024b); the caller
%   closes it and shows the PNG (show_saved_png).

if nargin < 5 || isempty(notes), notes = struct(); end

periods = {'LowCloudPeriod', 'SoledadCloudEvent'};
period_names = {'Low Cloud Periods', 'Soledad Cloud Events'};
methods = {'sCCN', 'skappa', 'sNd', 'sARG00'};
titles = {'s_{CCN}', 's_\kappa', 's_{Nd,VAP}', 's_{ARG00}';
          's_{CCN}', 's_\kappa', 's_{Nd,FM}',  's_{ARG00}'};
marker_colors = {[0.00 0.59 1.00], [0.00 0.16 0.78]};   % the blues of plot_fig10_comparison.m (low cloud, Soledad)
ink = [0.11 0.11 0.11];
fit_color = [0.84 0.37 0.00];   % Okabe & Ito (2008) vermillion #D55E00, for the regression line
fs = 16;                        % font size (points)

fig = new_figure('Visible', 'off', 'Position', [100 100 1900 1000]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:2
    for k = 1:4
        ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        ax.FontSize = fs;
        col = [methods{k} '_' periods{j}];
        y = SS.(col); x = X.(col);
        ok = isfinite(x) & isfinite(y);
        x = x(ok); y = y(ok);
        label = sprintf('(%c) %s', 'a' + 4 * (j - 1) + k - 1, titles{j, k});

        if numel(x) >= 3
            lim = [0, 1.05 * max([x; y])];
            h_pts = scatter(ax, x, y, 22, marker_colors{j}, 'filled', 'MarkerFaceAlpha', 0.6, ...
                            'DisplayName', 'one hour');
            h_11 = plot(ax, lim, lim, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.5, 'DisplayName', '1:1');
            c = polyfit(x, y, 1);                          % least-squares fit y = c(1) x + c(2)
            y_fit = polyval(c, x);
            R2 = 1 - sum((y - y_fit) .^ 2) / sum((y - mean(y)) .^ 2);   % coefficient of determination
            h_fit = plot(ax, lim, polyval(c, lim), '-', 'Color', fit_color, 'LineWidth', 2, ...
                         'DisplayName', 'least-squares fit');
            text(ax, 0.97, 0.04, {sprintf('y = %.2f x %+.3f', c(1), c(2)), sprintf('R^2 = %.2f', R2), ...
                 sprintf('n = %d h', numel(x))}, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'bottom', 'FontSize', fs - 1, 'Color', ink);
            xlim(ax, lim); ylim(ax, lim);
            axis(ax, 'square');
            if j == 1 && k == 1
                lg = legend(ax, [h_pts h_11 h_fit], 'NumColumns', 3, 'FontSize', fs);
                lg.Layout.Tile = 'north';
            end
        else
            axis(ax, 'square'); xlim(ax, [0 1]); ylim(ax, [0 1]);
            text(ax, 0.5, 0.5, 'not re-derived', 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                 'FontAngle', 'italic', 'FontSize', fs, 'Color', marker_colors{j});
        end
        text(ax, 0.03, 0.97, label, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fs + 2, 'Color', ink);
        if isfield(notes, col)
            text(ax, 0.03, 0.85, notes.(col), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'FontSize', fs - 4, 'FontAngle', 'italic', 'Color', [0.32 0.32 0.31]);
        end
        if k == 1
            ylabel(ax, {['\bf' period_names{j}], '\rmposted, Jeramy (%)'});
        end
        if j == 2
            xlabel(ax, 'mine (%)');
        end
    end
end
title(tl, {'Jeramy''s posted 1-h retrievals against mine, hour by hour', x_text}, 'FontSize', fs + 2);

if ~isempty(fig_path)
    exportgraphics(fig, fig_path, 'Resolution', 200);
end

end
