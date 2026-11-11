%% Timeline of the 24 Mt. Soledad cloud events, aligned on the day the cloud arrives (PDT, horizontal)
%
% DESCRIPTION:
%   Alternative to plot_dedrick2026_event_timing.m (same data and colors).
%   One horizontal row per event, time running left to right like a
%   timeline, with every event aligned on the local day its IN-CLOUD
%   window starts:
%
%   x  local time, Pacific Daylight Time (PDT = UTC - 7 h; all 24 events
%      fall inside daylight saving time), as a clock that repeats 00-24 h,
%      counted from local midnight of the day the in-cloud window starts
%      ("Cloud day"). The day before is "Day -1", the days after "Day +1",
%      "Day +2". Dark lines mark local midnight.
%   y  event number (event 1 at the top), with the date of the cloud day.
%
%   WHY ALIGN ON THE CLOUD DAY: the in-cloud windows then line up in one
%   band, so it is easy to read at which time of day the cloud arrives and
%   leaves, and how that sits relative to sunrise and sunset; the pre- and
%   post-cloud windows fall to the left and right of it. With one day ~1.1
%   inch wide, a window of a few hours is still several millimetres long.
%
%   THE THREE WINDOWS (supplement Table S5)
%     pre-cloud   grey, stippled     ambient sampling before the event
%     in-cloud    white (a cloud)    the GCVI event: Mt. Soledad in cloud
%     post-cloud  grey, diagonal     ambient sampling after the event
%   DAY AND NIGHT: light yellow where the sun is up at Mt. Soledad on that
%   event's dates, semi-transparent black where it is down (SZA > 90.833
%   deg, the NOAA sunrise/sunset convention), computed every 5 min.
%
% OUTPUTS (in ../output):
%   dedrick2026_event_timing_PDT_alignedCloudDay.png  (500 dpi)
%   dedrick2026_event_timing_PDT_alignedCloudDay.pdf  (vector)
%
% DATA:
%   Event windows: dedrick2026_soledad_cloud_events.m, transcribed from
%     Table S5 of the supplement (UTC).
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement Table S5,
%     doi:10.1175/JAS-D-26-0004.s1.
%   NOAA solar position equations: see solar_zenith_deg.m.
%
% DEPENDENCIES: solar_zenith_deg.m, draw_pattern_rect.m (this folder),
%   EPCAPE/Dedrick_2026_coupled_cases (event table),
%   EPCAPE/supersaturation/functions (new_figure.m)
%
% By Andrew John Buggee
%%

clear variables


%% 1. Settings

utcOffset_h  = -7;             % PDT = UTC - 7 h
lat_deg      = 32.8397;        % Mt. Soledad, La Jolla (approx.), degrees north
lon_deg      = -117.2519;      % degrees east (negative = west)
szaNight_deg = 90;         % sun below the horizon for SZA > this

% colors and transparency (same scheme as the other two versions)
c_cloudGrey = [0.55 0.55 0.55];
c_pattern   = [0.15 0.15 0.15];
c_inCloud   = [1 1 1];
c_edge      = [0.25 0.25 0.25];
c_day       = [1.00 0.96 0.78];
c_night     = [0 0 0];
a_night     = 0.22;
a_stage     = 0.65;
a_in        = 0.90;
c_midnight  = [0.15 0.15 0.15];
c_ink       = [0.043 0.043 0.043];

rowHeight     = 0.84;             % background height of each row (rows are 1 apart)
barHeight     = 0.56;             % height of the window bars
patSpacing_in = 0.035;

% figure layout (inches)
figSize_in = [7.5 6.6];
axPos_in   = [1.05 0.6 6.25 5.2];

% folders
scriptDir = fileparts(mfilename('fullpath'));
outDir    = fullfile(fileparts(scriptDir), 'output');
epcapeDir = fileparts(fileparts(scriptDir));
addpath(scriptDir);
addpath(fullfile(epcapeDir, 'Dedrick_2026_coupled_cases'));
addpath(fullfile(epcapeDir, 'supersaturation', 'functions'));
if ~isfolder(outDir), mkdir(outDir); end



%% 2. The 24 events in PDT, as hours since local midnight of the cloud day

E = dedrick2026_soledad_cloud_events();
n_ev = height(E);
to_pdt = @(t) t + hours(utcOffset_h);
win.pre  = [to_pdt(E.pre_start_utc)  to_pdt(E.pre_end_utc)];
win.in   = [to_pdt(E.in_start_utc)   to_pdt(E.in_end_utc)];
win.post = [to_pdt(E.post_start_utc) to_pdt(E.post_end_utc)];

cloudDay = dateshift(win.in(:, 1), 'start', 'day');       % local midnight of each event's cloud day
x_of = @(t, k) hours(t - cloudDay(k));

% axis limits: the nearest 6 h beyond the earliest and latest window
x_min = 6 * floor(min(hours(win.pre(:, 1) - cloudDay)) / 6);
x_max = 6 * ceil(max(hours(win.post(:, 2) - cloudDay)) / 6);



%% 3. Figure and axes

fig = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 figSize_in]);
ax = axes(fig, 'Units', 'inches', 'Position', axPos_in);
hold(ax, 'on');
xlim(ax, [x_min x_max]);
ylim(ax, [0.4, n_ev + 0.6]);
set(ax, 'YDir', 'reverse');
in_per_x = axPos_in(3) / diff(xlim(ax));
in_per_y = axPos_in(4) / diff(ylim(ax));

for k = 1:n_ev
    yb = k + [-1 1] * rowHeight / 2;
    yw = k + [-1 1] * barHeight / 2;

    % ---- day and night on this event's own dates ----
    t = (cloudDay(k) + hours(x_min):minutes(5):cloudDay(k) + hours(x_max))';
    night = solar_zenith_deg(t - hours(utcOffset_h), lat_deg, lon_deg) > szaNight_deg;
    edges = [1; find(diff(night)) + 1; numel(t)];
    for r = 1:numel(edges) - 1
        x0 = x_of(t(edges(r)), k);  x1 = x_of(t(edges(r + 1)), k);
        if night(edges(r))
            draw_pattern_rect(ax, x0, x1, yb(1), yb(2), c_night, a_night, "none", 1, 1, 1, c_pattern);
        else
            draw_pattern_rect(ax, x0, x1, yb(1), yb(2), c_day, 1, "none", 1, 1, 1, c_pattern);
        end
    end

    % ---- the three windows ----
    draw_pattern_rect(ax, x_of(win.pre(k, 1), k), x_of(win.pre(k, 2), k), yw(1), yw(2), ...
                      c_cloudGrey, a_stage, "stipple", in_per_x, in_per_y, patSpacing_in, c_pattern);
    draw_pattern_rect(ax, x_of(win.post(k, 1), k), x_of(win.post(k, 2), k), yw(1), yw(2), ...
                      c_cloudGrey, a_stage, "stripes", in_per_x, in_per_y, 1.5 * patSpacing_in, c_pattern);
    p = draw_pattern_rect(ax, x_of(win.in(k, 1), k), x_of(win.in(k, 2), k), yw(1), yw(2), ...
                          c_inCloud, a_in, "none", 1, 1, 1, c_pattern);
    p.EdgeColor = c_edge;  p.LineWidth = 0.5;
end

% midnight lines, and day names above the axes centred on the visible part of each day
for xd = 24 * ceil(x_min / 24):24:x_max
    if xd > x_min && xd < x_max
        plot(ax, [xd xd], ylim(ax), '-', 'Color', c_midnight, 'LineWidth', 1);
    end
end
for xd = 24 * floor(x_min / 24):24:x_max - 1
    n = xd / 24;
    if n == 0, name = 'Cloud day'; else, name = sprintf('Day %+d', n); end
    xc = (max(xd, x_min) + min(xd + 24, x_max)) / 2;
    text(ax, xc, 0.25, name, 'Interpreter', 'latex', 'FontSize', 10, 'Color', c_ink, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');
end



%% 4. Axes

xt = x_min:6:x_max;
set(ax, 'XTick', xt, 'XTickLabel', compose('%02d', mod(xt, 24)'), ...
        'YTick', 1:n_ev, 'YTickLabel', compose('%d \\ (%s)', (1:n_ev)', string(cloudDay, 'd MMM')), ...
        'TickDir', 'out', 'TickLength', [0.006 0.006], 'Layer', 'top', 'Box', 'off', 'FontSize', 9, ...
        'TickLabelInterpreter', 'latex', 'XColor', c_ink, 'YColor', c_ink, 'LineWidth', 0.75);
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = x_min:1:x_max;
xlabel(ax, 'Local time (PDT, UTC$-$7) [h]', 'Interpreter', 'latex', 'FontSize', 10);
ylabel(ax, 'Cloud event (cloud day, 2023)', 'Interpreter', 'latex', 'FontSize', 10);



%% 5. Legend

draw_legend(fig, [axPos_in(1), axPos_in(2) + axPos_in(4) + 0.35, axPos_in(3), 0.32], ...
            patSpacing_in, c_cloudGrey, c_pattern, c_inCloud, c_edge, c_day, c_night, ...
            a_stage, a_in, a_night, c_ink);

tag = 'dedrick2026_event_timing_PDT_alignedCloudDay';
exportgraphics(fig, fullfile(outDir, [tag '.png']), 'Resolution', 500);
exportgraphics(fig, fullfile(outDir, [tag '.pdf']), 'ContentType', 'vector');
close(fig);
fprintf('Saved %s.png/.pdf to %s\n', tag, outDir);




%% Local functions

function draw_legend(fig, pos_in, sp, c_grey, c_pat, c_in, c_edge, c_day, c_night, a_stage, a_in, a_night, c_ink)
% one row of swatches above the axes, in an invisible axes whose data units are inches
lax = axes(fig, 'Units', 'inches', 'Position', pos_in, 'Visible', 'off');
hold(lax, 'on');
xlim(lax, [0 pos_in(3)]);  ylim(lax, [0 pos_in(4)]);
items = {'Pre-cloud', 'In-cloud', 'Post-cloud', 'Day', 'Night'};
sw = 0.30;  sh = 0.16;  pitch = pos_in(3) / numel(items);
for i = 1:numel(items)
    x0 = (i - 1) * pitch + 0.1;  x1 = x0 + sw;  y0 = (pos_in(4) - sh) / 2;  y1 = y0 + sh;
    switch items{i}
        case 'Pre-cloud'
            draw_pattern_rect(lax, x0, x1, y0, y1, c_grey, a_stage, "stipple", 1, 1, sp, c_pat);
        case 'In-cloud'
            p = draw_pattern_rect(lax, x0, x1, y0, y1, c_in, a_in, "none", 1, 1, sp, c_pat);
            p.EdgeColor = c_edge;  p.LineWidth = 0.5;
        case 'Post-cloud'
            draw_pattern_rect(lax, x0, x1, y0, y1, c_grey, a_stage, "stripes", 1, 1, 1.5 * sp, c_pat);
        case 'Day'
            p = draw_pattern_rect(lax, x0, x1, y0, y1, c_day, 1, "none", 1, 1, sp, c_pat);
            p.EdgeColor = [0.7 0.7 0.7];  p.LineWidth = 0.5;
        case 'Night'
            draw_pattern_rect(lax, x0, x1, y0, y1, c_night, a_night, "none", 1, 1, sp, c_pat);
    end
    text(lax, x1 + 0.07, (y0 + y1) / 2, items{i}, 'Interpreter', 'latex', 'FontSize', 9, ...
         'Color', c_ink, 'VerticalAlignment', 'middle');
end
end
