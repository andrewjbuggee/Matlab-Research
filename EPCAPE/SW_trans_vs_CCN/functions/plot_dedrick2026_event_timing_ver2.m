%% Timeline (PDT) of the pre-, in- and post-cloud windows of the 24 Mt. Soledad cloud events, over day and night
%
% DESCRIPTION:
%   Version 2 of plot_dedrick2026_event_timing.m. Time runs along the x
%   axis, from the start of event 1's pre-cloud window to the end of event
%   24's post-cloud window, in Pacific Daylight Time (PDT = UTC - 7 h): all
%   24 events (26 Apr - 10 Oct 2023) fall inside daylight saving time
%   (12 Mar - 5 Nov 2023), so PDT is the local clock time. Each event has
%   its own row (event 1 at the top).
%
%   THE THREE WINDOWS (supplement Table S5)
%     pre-cloud   grey, stippled     ambient sampling before the event
%     in-cloud    white (a cloud)    the GCVI event: Mt. Soledad in cloud
%     post-cloud  grey, diagonal     ambient sampling after the event
%   Within one event the three windows follow one another in time and never
%   overlap: the gap between the end of pre-cloud and the start of in-cloud
%   is 0-7 h, and between in-cloud and post-cloud 0-6 h. (In version 1 the
%   windows only appeared to overlap because each was wrapped onto a 0-24 h
%   time-of-day axis, and every event lasts 26-48 h in total.) So all three
%   are drawn on the event's row, semi-transparent. Consecutive EVENTS can
%   overlap in time (e.g. event 5's post-cloud window runs into event 6's
%   pre-cloud window), which shows up as bars on neighbouring rows sharing
%   the same time.
%
%   DAY AND NIGHT BACKGROUND
%   Light yellow where the sun is up at Mt. Soledad, semi-transparent black
%   where it is down (solar zenith angle > 90.833 deg: the NOAA sunrise and
%   sunset convention, upper limb of the sun with standard refraction),
%   computed every 5 min over the whole time range.
%
%   TIME AXIS
%   timeAxis = "continuous": one unbroken axis from 26 Apr to 10 Oct 2023
%     (167 days). Each event is only ~1-2 days long, so the bars are narrow.
%   timeAxis = "gapsRemoved": the same continuous time inside each group of
%     events, but the days with no event between groups are cut out (marked
%     by a double slash on the axis), so each day is ~7x wider.
%
% OUTPUTS (in ../output):
%   dedrick2026_event_timeline_PDT_<timeAxis>.png  (500 dpi)
%   dedrick2026_event_timeline_PDT_<timeAxis>.pdf  (vector)
%
% DATA:
%   Event windows: dedrick2026_soledad_cloud_events.m, transcribed from
%     Table S5 of the supplement (UTC).
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: Quantitative estimates of upwind
%     and on-site particle scavenging and retrieved cloud supersaturations
%     from urban coastal stratocumulus during EPCAPE. J. Atmos. Sci., 83,
%     1061-1082, doi:10.1175/JAS-D-26-0004.1; supplement Table S5,
%     doi:10.1175/JAS-D-26-0004.s1.
%   NOAA Global Monitoring Laboratory, "General Solar Position
%     Calculations" (Spencer 1971 Fourier series),
%     https://gml.noaa.gov/grad/solcalc/solareqns.PDF.
%   Spencer, J. W., 1971: Fourier series representation of the position of
%     the Sun. Search, 2, 172.
%
% DEPENDENCIES: EPCAPE/Dedrick_2026_coupled_cases (event table),
%   EPCAPE/supersaturation/functions (new_figure.m)
%
% By Andrew John Buggee
%%

clear variables


%% 1. Settings

timeAxis     = "continuous";   % "continuous" | "gapsRemoved"
gapPad_h     = 6;              % gapsRemoved: hours kept before/after each group of events

utcOffset_h  = -7;             % PDT = UTC - 7 h
lat_deg      = 32.8397;        % Mt. Soledad, La Jolla (approx.), degrees north
lon_deg      = -117.2519;      % degrees east (negative = west)
szaNight_deg = 90.833;         % sun below the horizon for SZA > this

% colors and transparency
c_cloudGrey = [0.55 0.55 0.55];   % pre- and post-cloud fill
c_pattern   = [0.15 0.15 0.15];   % stipple dots and diagonal stripes
c_inCloud   = [1 1 1];            % in-cloud: white
c_edge      = [0.25 0.25 0.25];   % thin outline so the white bars show on the yellow
c_day       = [1.00 0.96 0.78];   % light yellow
c_night     = [0 0 0];            % black, drawn semi-transparent
a_night     = 0.22;               % night opacity
a_stage     = 0.65;               % opacity of the pre-/post-cloud fill
a_in        = 0.90;               % opacity of the in-cloud fill
c_ink       = [0.043 0.043 0.043];

barHeight   = 0.62;               % bar height (rows are 1 apart)
patSpacing_in = 0.035;            % spacing of stipple dots and stripes (inches on the page)

% figure layout (inches)
figSize_in = [7.5 6.2];
axPos_in   = [0.85 0.75 6.45 4.75];    % [left bottom width height] of the main axes

% folders
scriptDir = fileparts(mfilename('fullpath'));            % .../SW_trans_vs_CCN/functions
outDir    = fullfile(fileparts(scriptDir), 'output');    % .../SW_trans_vs_CCN/output
epcapeDir = fileparts(fileparts(scriptDir));
addpath(fullfile(epcapeDir, 'Dedrick_2026_coupled_cases'));
addpath(fullfile(epcapeDir, 'supersaturation', 'functions'));   % new_figure.m
if ~isfolder(outDir), mkdir(outDir); end



%% 2. The 24 events in PDT, and the time -> x mapping

E = dedrick2026_soledad_cloud_events();
n_ev = height(E);
to_pdt = @(t) t + hours(utcOffset_h);
win.pre  = [to_pdt(E.pre_start_utc)  to_pdt(E.pre_end_utc)];
win.in   = [to_pdt(E.in_start_utc)   to_pdt(E.in_end_utc)];
win.post = [to_pdt(E.post_start_utc) to_pdt(E.post_end_utc)];

t_first = min(win.pre(:, 1));
t_last  = max(win.post(:, 2));

% x is in hours. "continuous": hours since t_first. "gapsRemoved": the
% groups of events (all windows, padded by gapPad_h and merged where they
% touch) are placed side by side with the empty days between them removed
switch timeAxis
    case "continuous"
        seg = [t_first t_last];
    case "gapsRemoved"
        seg = merge_intervals([win.pre(:, 1) - hours(gapPad_h), win.post(:, 2) + hours(gapPad_h)]);
    otherwise
        error('timeAxis must be "continuous" or "gapsRemoved"');
end
seg_len_h = hours(seg(:, 2) - seg(:, 1));
seg_x0    = [0; cumsum(seg_len_h(1:end - 1))];     % x where each segment starts
x_max     = sum(seg_len_h);
t2x = @(t) time_to_x(t, seg, seg_x0);



%% 3. Figure and axes (fixed size in inches, so the patterns can be drawn at 45 deg)

fig = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 figSize_in]);
fig.PaperPositionMode = 'auto';
ax = axes(fig, 'Units', 'inches', 'Position', axPos_in);
hold(ax, 'on');
xlim(ax, [0 x_max]);
ylim(ax, [0.4, n_ev + 0.6]);
set(ax, 'YDir', 'reverse');

% inches per data unit, for the patterns
in_per_x = axPos_in(3) / x_max;
in_per_y = axPos_in(4) / diff(ylim(ax));



%% 4. Day and night background, continuous over each segment

for s = 1:size(seg, 1)
    t = (seg(s, 1):minutes(5):seg(s, 2))';
    sza = solar_zenith_deg(t - hours(utcOffset_h), lat_deg, lon_deg);   % SZA needs UTC
    night = sza > szaNight_deg;
    k = [1; find(diff(night)) + 1; numel(t) + 1];                      % runs of day or night
    for r = 1:numel(k) - 1
        x0 = t2x(t(k(r)));
        x1 = t2x(min(t(min(k(r + 1), numel(t))), seg(s, 2)));
        if r == numel(k) - 1, x1 = t2x(seg(s, 2)); end
        if night(k(r))
            fill_rect(ax, x0, x1, 0.4, n_ev + 0.6, c_night, a_night);
        else
            fill_rect(ax, x0, x1, 0.4, n_ev + 0.6, c_day, 1);
        end
    end
end

% faint row separators
for k = 1:n_ev - 1
    plot(ax, [0 x_max], [k k] + 0.5, '-', 'Color', [1 1 1 0.6], 'LineWidth', 0.5);
end



%% 5. The three windows of every event

for k = 1:n_ev
    y0 = k - barHeight / 2;  y1 = k + barHeight / 2;

    % pre-cloud: grey, stippled
    x0 = t2x(win.pre(k, 1));  x1 = t2x(win.pre(k, 2));
    fill_rect(ax, x0, x1, y0, y1, c_cloudGrey, a_stage);
    stipple(ax, x0, x1, y0, y1, in_per_x, in_per_y, patSpacing_in, c_pattern);

    % post-cloud: grey, diagonal stripes
    x0 = t2x(win.post(k, 1));  x1 = t2x(win.post(k, 2));
    fill_rect(ax, x0, x1, y0, y1, c_cloudGrey, a_stage);
    stripes(ax, x0, x1, y0, y1, in_per_x, in_per_y, 1.5 * patSpacing_in, c_pattern);

    % in-cloud: white, outlined
    x0 = t2x(win.in(k, 1));  x1 = t2x(win.in(k, 2));
    p = fill_rect(ax, x0, x1, y0, y1, c_inCloud, a_in);
    p.EdgeColor = c_edge;  p.LineWidth = 0.5;
end



%% 6. Axes, ticks and the axis breaks

set(ax, 'YTick', 1:n_ev, 'YTickLabel', compose('%d', (1:n_ev)'), 'TickDir', 'out', ...
        'TickLength', [0.006 0.006], 'Layer', 'top', 'Box', 'off', 'FontSize', 9, ...
        'TickLabelInterpreter', 'latex', 'XColor', c_ink, 'YColor', c_ink, 'LineWidth', 0.75);
ylabel(ax, 'Mt.\ Soledad cloud event (Table S5)', 'Interpreter', 'latex', 'FontSize', 10);

switch timeAxis
    case "continuous"
        % a tick at the first of every month, minor ticks every 7 days
        months = dateshift(t_first, 'start', 'month', 'next'):calmonths(1):t_last;
        set(ax, 'XTick', t2x(months), 'XTickLabel', cellstr(string(months, 'd MMM')));
        ax.XAxis.MinorTick = 'on';
        ax.XAxis.MinorTickValues = t2x(dateshift(t_first, 'start', 'day', 'next'):days(7):t_last);
        xlabel(ax, 'Date in 2023 (PDT, UTC$-$7)', 'Interpreter', 'latex', 'FontSize', 10);
    case "gapsRemoved"
        % a tick at every local midnight inside a segment, labelled with the
        % date at the first midnight of each segment; // marks each cut
        ticks = [];  labels = {};
        for s = 1:size(seg, 1)
            mids = dateshift(seg(s, 1), 'start', 'day', 'next'):days(1):seg(s, 2);
            ticks = [ticks, t2x(mids)]; %#ok<AGROW>
            lab = repmat({''}, 1, numel(mids));
            if ~isempty(mids), lab{1} = char(string(mids(1), 'd MMM')); end
            labels = [labels, lab]; %#ok<AGROW>
        end
        set(ax, 'XTick', ticks, 'XTickLabel', labels, 'XTickLabelRotation', 90);
        for s = 2:size(seg, 1)
            xb = seg_x0(s);
            plot(ax, [xb xb], [0.4, n_ev + 0.6], '-', 'Color', 'w', 'LineWidth', 2.5);   % visible cut
            text(ax, xb, n_ev + 0.6, '//', 'Interpreter', 'latex', 'FontSize', 9, 'Color', c_ink, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Clipping', 'off');
        end
        xlabel(ax, 'Date in 2023 (PDT, UTC$-$7); days without events removed', ...
               'Interpreter', 'latex', 'FontSize', 10);
end



%% 7. Legend (drawn by hand, so the stipple and stripes appear in it)

lgPos_in = [axPos_in(1), axPos_in(2) + axPos_in(4) + 0.12, axPos_in(3), 0.32];
lax = axes(fig, 'Units', 'inches', 'Position', lgPos_in, 'Visible', 'off');
hold(lax, 'on');
xlim(lax, [0 lgPos_in(3)]);  ylim(lax, [0 lgPos_in(4)]);    % data units = inches
items = {'Pre-cloud', 'In-cloud', 'Post-cloud', 'Day', 'Night'};
sw = 0.30;  sh = 0.16;  pitch = lgPos_in(3) / numel(items);
for i = 1:numel(items)
    x0 = (i - 1) * pitch + 0.15;  x1 = x0 + sw;  y0 = (lgPos_in(4) - sh) / 2;  y1 = y0 + sh;
    switch items{i}
        case 'Pre-cloud'
            fill_rect(lax, x0, x1, y0, y1, c_cloudGrey, a_stage);
            stipple(lax, x0, x1, y0, y1, 1, 1, patSpacing_in, c_pattern);
        case 'In-cloud'
            p = fill_rect(lax, x0, x1, y0, y1, c_inCloud, a_in);
            p.EdgeColor = c_edge;  p.LineWidth = 0.5;
        case 'Post-cloud'
            fill_rect(lax, x0, x1, y0, y1, c_cloudGrey, a_stage);
            stripes(lax, x0, x1, y0, y1, 1, 1, 1.5 * patSpacing_in, c_pattern);
        case 'Day'
            p = fill_rect(lax, x0, x1, y0, y1, c_day, 1);
            p.EdgeColor = [0.7 0.7 0.7];  p.LineWidth = 0.5;
        case 'Night'
            fill_rect(lax, x0, x1, y0, y1, c_night, a_night);
    end
    text(lax, x1 + 0.07, (y0 + y1) / 2, items{i}, 'Interpreter', 'latex', 'FontSize', 9, ...
         'Color', c_ink, 'VerticalAlignment', 'middle');
end

tag = "dedrick2026_event_timeline_PDT_" + timeAxis;
exportgraphics(fig, fullfile(outDir, tag + ".png"), 'Resolution', 500);
exportgraphics(fig, fullfile(outDir, tag + ".pdf"), 'ContentType', 'vector');
close(fig);
fprintf('Saved %s.png/.pdf to %s\n', tag, outDir);




%% Local functions

function x = time_to_x(t, seg, seg_x0)
% hours along the (possibly cut) time axis; t must lie inside a segment
x = nan(size(t));
for s = 1:size(seg, 1)
    in = t >= seg(s, 1) & t <= seg(s, 2);
    x(in) = seg_x0(s) + hours(t(in) - seg(s, 1));
end
end


function m = merge_intervals(iv)
% merge overlapping [start end] rows (datetime), sorted by start
iv = sortrows(iv, 1);
m = iv(1, :);
for k = 2:size(iv, 1)
    if iv(k, 1) <= m(end, 2)
        m(end, 2) = max(m(end, 2), iv(k, 2));
    else
        m(end + 1, :) = iv(k, :); %#ok<AGROW>
    end
end
end


function p = fill_rect(ax, x0, x1, y0, y1, color, alpha)
% filled rectangle without an edge (patch, so it can be transparent)
p = patch(ax, [x0 x1 x1 x0], [y0 y0 y1 y1], color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end


function stipple(ax, x0, x1, y0, y1, in_per_x, in_per_y, spacing_in, color)
% dots on a staggered grid with spacing_in inches between them on the page
X0 = x0 * in_per_x;  X1 = x1 * in_per_x;  Y0 = y0 * in_per_y;  Y1 = y1 * in_per_y;   % to inches
xs = []; ys = [];
rows = (Y0 + spacing_in / 2):spacing_in:Y1;
for r = 1:numel(rows)
    off = mod(r, 2) * spacing_in / 2;                        % stagger every other row
    cols = (X0 + spacing_in / 4 + off):spacing_in:X1;
    xs = [xs, cols]; %#ok<AGROW>
    ys = [ys, repmat(rows(r), 1, numel(cols))]; %#ok<AGROW>
end
if isempty(xs), return, end
plot(ax, xs / in_per_x, ys / in_per_y, '.', 'Color', color, 'MarkerSize', 2.5);
end


function stripes(ax, x0, x1, y0, y1, in_per_x, in_per_y, spacing_in, color)
% 45-degree stripes (on the page) clipped to the rectangle; all in one line object
X0 = x0 * in_per_x;  X1 = x1 * in_per_x;  Y0 = y0 * in_per_y;  Y1 = y1 * in_per_y;   % to inches
xs = []; ys = [];
for c = (X0 - (Y1 - Y0)):spacing_in:X1           % lines X - (Y - Y0) = c
    xa = max(X0, c);           xb = min(X1, c + (Y1 - Y0));
    if xb <= xa, continue, end
    xs = [xs, xa, xb, NaN]; %#ok<AGROW>
    ys = [ys, Y0 + (xa - c), Y0 + (xb - c), NaN]; %#ok<AGROW>
end
if isempty(xs), return, end
plot(ax, xs / in_per_x, ys / in_per_y, '-', 'Color', color, 'LineWidth', 0.4);
end


function sza_deg = solar_zenith_deg(t_utc, lat_deg, lon_deg)
% Solar zenith angle (deg) from the NOAA "General Solar Position
% Calculations" (Spencer 1971 Fourier series); ~0.1-0.2 deg accuracy
doy  = day(t_utc, 'dayofyear');
hr   = hour(t_utc) + minute(t_utc) / 60 + second(t_utc) / 3600;   % UTC hour of day
n_yr = 365 + (mod(year(t_utc), 4) == 0);
g = 2 * pi ./ n_yr .* (doy - 1 + (hr - 12) / 24);                 % fractional year (rad)
eqtime = 229.18 * (0.000075 + 0.001868 * cos(g) - 0.032077 * sin(g) ...
                   - 0.014615 * cos(2 * g) - 0.040849 * sin(2 * g));          % equation of time (min)
decl = 0.006918 - 0.399912 * cos(g) + 0.070257 * sin(g) - 0.006758 * cos(2 * g) ...
       + 0.000907 * sin(2 * g) - 0.002697 * cos(3 * g) + 0.00148 * sin(3 * g);  % declination (rad)
tst = hr * 60 + eqtime + 4 * lon_deg;            % true solar time (min), UTC input
ha_deg = tst / 4 - 180;                          % hour angle (deg)
cos_z = sind(lat_deg) .* sin(decl) + cosd(lat_deg) .* cos(decl) .* cosd(ha_deg);
sza_deg = acosd(min(max(cos_z, -1), 1));
end
