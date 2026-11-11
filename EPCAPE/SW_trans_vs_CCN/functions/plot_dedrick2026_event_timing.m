%% Time of day (PDT) of the pre-, in- and post-cloud windows of the 24 Mt. Soledad cloud events, over several days
%
% DESCRIPTION:
%   For each of the 24 Mt. Soledad cloud events of Dedrick et al. (2026,
%   supplement Table S5), shows when the pre-cloud, in-cloud and post-cloud
%   windows occurred, on top of the night and day of the event's own dates.
%
%   x  event number, 1-24 (Table S5 order)
%   y  local time, Pacific Daylight Time (PDT = UTC - 7 h; all 24 events,
%      26 Apr - 10 Oct 2023, fall inside daylight saving time, 12 Mar -
%      5 Nov 2023). The axis is a clock that runs 0-24 h and then repeats:
%      "Day 1" is the local calendar day on which the event's pre-cloud
%      window starts, "Day 2" the next day, and so on. Every event lasts
%      26-48 h from the start of pre-cloud to the end of post-cloud, so the
%      events reach into Day 2 or Day 3. A dark line marks each local
%      midnight, and the day number is written beside each 24 h block.
%
%   THE THREE WINDOWS (supplement Table S5), in one column per event
%     pre-cloud   grey, stippled     ambient sampling before the event
%     in-cloud    white (a cloud)    the GCVI event: Mt. Soledad in cloud
%     post-cloud  grey, diagonal     ambient sampling after the event
%   Within one event the windows follow one another and never overlap (gaps
%   of 0-7 h between them), so they share a column; they are drawn
%   semi-transparent so the day/night shading shows through.
%
%   DAY AND NIGHT BACKGROUND
%   Light yellow where the sun is up at Mt. Soledad on that event's dates,
%   semi-transparent black where it is down (solar zenith angle > 90.833
%   deg, the NOAA sunrise/sunset convention), computed every 5 min.
%
% OUTPUTS (in ../output):
%   dedrick2026_event_timing_PDT.png  (500 dpi)
%   dedrick2026_event_timing_PDT.pdf  (vector)
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
szaNight_deg = 90.833;         % sun below the horizon for SZA > this

% colors and transparency (same scheme as plot_dedrick2026_event_timing_ver2.m)
c_cloudGrey = [0.55 0.55 0.55];   % pre- and post-cloud fill
c_pattern   = [0.15 0.15 0.15];   % stipple dots and diagonal stripes
c_inCloud   = [1 1 1];            % in-cloud: white
c_edge      = [0.25 0.25 0.25];   % thin outline so the white bars show on the yellow
c_day       = [1.00 0.96 0.78];   % light yellow
c_night     = [0 0 0];            % black, drawn semi-transparent
a_night     = 0.22;               % night opacity
a_stage     = 0.65;               % opacity of the pre-/post-cloud fill
a_in        = 0.90;               % opacity of the in-cloud fill
c_midnight  = [0.15 0.15 0.15];   % line at each local midnight
c_ink       = [0.043 0.043 0.043];

eventWidth    = 0.84;             % width of each event's background (events are 1 apart)
barWidth      = 0.56;             % width of the window bars
patSpacing_in = 0.035;            % spacing of stipple dots and stripes (inches on the page)

% figure layout (inches)
figSize_in = [7.5 7.6];
axPos_in   = [0.75 0.65 6.2 6.3];      % [left bottom width height] of the main axes

% folders
scriptDir = fileparts(mfilename('fullpath'));            % .../SW_trans_vs_CCN/functions
outDir    = fullfile(fileparts(scriptDir), 'output');    % .../SW_trans_vs_CCN/output
epcapeDir = fileparts(fileparts(scriptDir));
addpath(scriptDir);
addpath(fullfile(epcapeDir, 'Dedrick_2026_coupled_cases'));
addpath(fullfile(epcapeDir, 'supersaturation', 'functions'));   % new_figure.m
if ~isfolder(outDir), mkdir(outDir); end



%% 2. The 24 events in PDT, as hours since local midnight of Day 1

E = dedrick2026_soledad_cloud_events();
n_ev = height(E);
to_pdt = @(t) t + hours(utcOffset_h);
win.pre  = [to_pdt(E.pre_start_utc)  to_pdt(E.pre_end_utc)];
win.in   = [to_pdt(E.in_start_utc)   to_pdt(E.in_end_utc)];
win.post = [to_pdt(E.post_start_utc) to_pdt(E.post_end_utc)];

day1 = dateshift(win.pre(:, 1), 'start', 'day');          % local midnight starting Day 1 of each event
y_of = @(t, k) hours(t - day1(k));                         % hours since then

% enough whole days for the longest event
n_days = ceil(max(hours(win.post(:, 2) - day1)) / 24);
y_max  = 24 * n_days;



%% 3. Figure and axes (fixed size in inches, so the patterns can be drawn at 45 deg)

fig = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 figSize_in]);
ax = axes(fig, 'Units', 'inches', 'Position', axPos_in);
hold(ax, 'on');
xlim(ax, [0.4, n_ev + 0.6]);
ylim(ax, [0 y_max]);
in_per_x = axPos_in(3) / diff(xlim(ax));
in_per_y = axPos_in(4) / diff(ylim(ax));

for k = 1:n_ev
    xc = E.event_id(k);
    xb = xc + [-1 1] * eventWidth / 2;                    % background edges
    xw = xc + [-1 1] * barWidth / 2;                      % bar edges

    % ---- day and night on this event's own dates ----
    t = (day1(k):minutes(5):day1(k) + hours(y_max))';
    night = solar_zenith_deg(t - hours(utcOffset_h), lat_deg, lon_deg) > szaNight_deg;   % SZA needs UTC
    edges = [1; find(diff(night)) + 1; numel(t)];         % runs of day or night
    for r = 1:numel(edges) - 1
        y0 = y_of(t(edges(r)), k);  y1 = y_of(t(edges(r + 1)), k);
        if night(edges(r))
            draw_pattern_rect(ax, xb(1), xb(2), y0, y1, c_night, a_night, "none", 1, 1, 1, c_pattern);
        else
            draw_pattern_rect(ax, xb(1), xb(2), y0, y1, c_day, 1, "none", 1, 1, 1, c_pattern);
        end
    end

    % ---- the three windows ----
    draw_pattern_rect(ax, xw(1), xw(2), y_of(win.pre(k, 1), k), y_of(win.pre(k, 2), k), ...
                      c_cloudGrey, a_stage, "stipple", in_per_x, in_per_y, patSpacing_in, c_pattern);
    draw_pattern_rect(ax, xw(1), xw(2), y_of(win.post(k, 1), k), y_of(win.post(k, 2), k), ...
                      c_cloudGrey, a_stage, "stripes", in_per_x, in_per_y, 1.5 * patSpacing_in, c_pattern);
    p = draw_pattern_rect(ax, xw(1), xw(2), y_of(win.in(k, 1), k), y_of(win.in(k, 2), k), ...
                          c_inCloud, a_in, "none", 1, 1, 1, c_pattern);
    p.EdgeColor = c_edge;  p.LineWidth = 0.5;
end

% a line at every local midnight, and the day number beside each 24 h block
for d = 1:n_days - 1
    plot(ax, xlim(ax), 24 * d * [1 1], '-', 'Color', c_midnight, 'LineWidth', 1);
end
for d = 1:n_days
    text(ax, n_ev + 0.75, 24 * (d - 0.5), sprintf('Day %d', d), 'Interpreter', 'latex', ...
         'FontSize', 10, 'Color', c_ink, 'Rotation', 90, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', 'Clipping', 'off');
end



%% 4. Axes: the clock repeats 0-24 in every day

yt = 0:6:y_max;
ytl = compose('%02d:00', mod(yt, 24)');
set(ax, 'XTick', 1:n_ev, 'YTick', yt, 'YTickLabel', ytl, 'TickDir', 'out', ...
        'TickLength', [0.006 0.006], 'Layer', 'top', 'Box', 'off', 'FontSize', 9, ...
        'TickLabelInterpreter', 'latex', 'XColor', c_ink, 'YColor', c_ink, 'LineWidth', 0.75);
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0:1:y_max;
xlabel(ax, 'Mt.\ Soledad cloud event (Dedrick et al.\ 2026, Table S5)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel(ax, 'Local time (PDT, UTC$-$7); Day 1 = day pre-cloud starts', 'Interpreter', 'latex', 'FontSize', 10);



%% 5. Legend (drawn by hand, so the stipple and stripes appear in it)

draw_legend(fig, [axPos_in(1), axPos_in(2) + axPos_in(4) + 0.15, axPos_in(3), 0.32], ...
            patSpacing_in, c_cloudGrey, c_pattern, c_inCloud, c_edge, c_day, c_night, ...
            a_stage, a_in, a_night, c_ink);

exportgraphics(fig, fullfile(outDir, 'dedrick2026_event_timing_PDT.png'), 'Resolution', 500);
exportgraphics(fig, fullfile(outDir, 'dedrick2026_event_timing_PDT.pdf'), 'ContentType', 'vector');
close(fig);
fprintf('Saved dedrick2026_event_timing_PDT.png/.pdf to %s\n', outDir);




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
