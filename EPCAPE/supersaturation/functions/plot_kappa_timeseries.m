function fig = plot_kappa_timeseries(K, kappa_const, fig_path)
% PLOT_KAPPA_TIMESERIES  Hourly, daily and weekly means of the posted EPCAPE kappa series.
%
%   fig = plot_kappa_timeseries(K, kappa_const, fig_path)
%
%   Three panels, one per averaging interval: (a) hourly, (b) daily, (c) weekly
%   means of each series in K (output of read_kappa_product). Each mean is the
%   plain mean of the samples in the interval [t, t + interval), ignoring NaN; an
%   interval without samples is NaN, so the lines break across data gaps. Weeks
%   are MATLAB's calendar weeks (retime 'weekly', starting on Sunday).
%   A dashed black line marks kappa_const, the constant kappa used so far for
%   s_kappa (0.37, the paper's campaign mean at the Pier, supplement Text S10).
%
%   Colours: the Okabe & Ito (2008) colour-blind-safe palette ("Color Universal
%   Design", https://jfly.uni-koeln.de/color/). The Pier series, the one that goes
%   with the Pier Hoppel minimum in s_kappa, is drawn last (on top) and thickest.
%
%   Inputs
%     K            timetable from read_kappa_product (kappa, unitless)
%     kappa_const  the constant kappa to mark (unitless)
%     fig_path     where to save the PNG ('' = do not save)
%
%   Drawn INVISIBLE and saved with exportgraphics, as the other figures of the live
%   script (exporting a visible figure has hung MATLAB R2024b). The caller closes
%   it and shows the saved PNG (show_saved_png).

% series to draw, in drawing order (the Pier series last, so it is on top)
vars   = {'kappa_sol_f44', 'kappa_sol_OC', 'kappa_sol_01', 'kappa_pier_01'};
labels = {'Mt. Soledad, \kappa_{org} from f_{44}', 'Mt. Soledad, \kappa_{org} from O/C', ...
          'Mt. Soledad, \kappa_{org} = 0.10', 'Scripps Pier, \kappa_{org} = 0.10'};
colors = {[0.000 0.620 0.451], ...   % Okabe-Ito bluish green #009E73
          [0.800 0.475 0.655], ...   % Okabe-Ito reddish purple #CC79A7
          [0.902 0.624 0.000], ...   % Okabe-Ito orange #E69F00
          [0.000 0.447 0.698]};      % Okabe-Ito blue #0072B2
widths = {[0.5 1.2 1.8], [0.5 1.2 1.8], [0.5 1.2 1.8], [0.9 2.0 2.6]};   % line width per panel (hourly, daily, weekly)

steps  = {'hourly', 'daily', 'weekly'};
titles = {'(a) hourly means', '(b) daily means', '(c) weekly means'};
fs = 14;                           % font size (points)

fig = new_figure('Visible', 'off', 'Position', [100 100 1400 1000]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:3
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    ax.FontSize = fs;
    A = retime(K(:, vars), steps{p}, @(x) mean(x, 'omitnan'));   % mean of each interval; NaN if empty
    h = gobjects(numel(vars), 1);
    for k = 1:numel(vars)
        y = A.(vars{k});
        if p == 3
            h(k) = plot(ax, A.time, y, '-o', 'Color', colors{k}, 'LineWidth', widths{k}(p), ...
                        'MarkerSize', 4, 'MarkerFaceColor', colors{k}, 'DisplayName', labels{k});
        else
            h(k) = plot(ax, A.time, y, '-', 'Color', colors{k}, 'LineWidth', widths{k}(p), 'DisplayName', labels{k});
        end
    end
    hc = yline(ax, kappa_const, '--', 'Color', [0 0 0], 'LineWidth', 1.8, ...
               'DisplayName', sprintf('\\kappa = %.2f (constant used so far)', kappa_const));
    ylabel(ax, '\kappa');
    ylim(ax, [0 1]);
    xlim(ax, [min(K.time) max(K.time)]);
    text(ax, 0.01, 0.95, titles{p}, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fs + 2);
    if p == 1
        % one legend above all panels, so the three x axes stay aligned
        lg = legend(ax, [h(end:-1:1); hc], 'NumColumns', 3, 'FontSize', fs - 1);
        lg.Layout.Tile = 'north';
    end
end
title(tl, 'Composition-based hygroscopicity \kappa at EPCAPE (bb2743661p\_37\_1\_kappa.nc)', 'FontSize', fs + 2);

if ~isempty(fig_path)
    exportgraphics(fig, fig_path, 'Resolution', 200);
end

end
