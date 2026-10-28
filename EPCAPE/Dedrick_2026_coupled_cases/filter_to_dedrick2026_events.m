%% Keep only the data that fall inside the 24 Mt. Soledad cloud events of Dedrick et al. (2026)
%
% DESCRIPTION:
%   Takes ANY time-stamped data (a datetime vector, a timetable, or a table
%   with a datetime column) and finds the samples that fall inside the
%   pre-cloud, in-cloud or post-cloud windows of the 24 Mt. Soledad cloud
%   events (Dedrick et al. 2026, supplement Table S5; the windows come from
%   dedrick2026_soledad_cloud_events.m). It returns a logical mask, the
%   event number of each kept sample, and the filtered data.
%
%   SAMPLES THAT ARE AVERAGES OVER AN INTERVAL
%   Many EPCAPE products are averages: the posted supersaturations and the
%   low-cloud flag are hourly, the Pier CCN is 30-min, the FM-120 is 5-min.
%   Such a sample covers an interval, not an instant, and may straddle the
%   start or end of a window. Set 'SampleLength' to the averaging interval
%   and 'TimeLabel' to where the time stamp sits in that interval; a sample
%   is then kept when at least 'MinOverlap' of its interval lies inside a
%   window. With the default SampleLength = 0 a sample is kept when its
%   time stamp lies inside a window, start included and end excluded.
%
%   KEEPING ONLY THE CLOUD INTERVALS OF AN EVENT
%   The in-cloud window is the whole GCVI event, haze included: the paper's
%   "214 cloud hours" are the total length of the 24 windows (see
%   dedrick2026_soledad_cloud_events.m). To keep only the 5-min intervals
%   inside a window that meet the cloud condition (visibility < 1 km and
%   FM-120 LWC > 0.01 g m-3; supplement Table S4), pass the classification
%   from soledad_cloud_condition_5min.m as 'SoledadCloud5min'; the in-cloud
%   windows are then narrowed to those intervals (137 h in total).
%
%   HOW THE OVERLAP IS COMPUTED
%   For each event the qualifying time is a set of disjoint intervals
%   [a_k, b_k). The cumulative qualifying time C(t), i.e. the total time in
%   those intervals before t, is piecewise linear: slope 1 inside an
%   interval and flat between intervals. So the time any sample interval
%   [s, e) shares with the set is simply C(e) - C(s), found by linear
%   interpolation through the interval end points. This is exact, needs no
%   loop over samples, and works for millions of samples.
%
% INPUTS:
%   data - one of
%            * datetime vector (UTC, no time zone, as all EPCAPE readers here)
%            * timetable (its row times are used)
%            * table with a datetime variable (name it with 'TimeVariable')
%          Datetimes that carry a time zone are converted to UTC first.
%
%   Optional name-value pairs:
%   'Period'           "inCloud" (default) | "preCloud" | "postCloud" |
%                      "allWindows" (pre, in and post together)
%   'Events'           event numbers to keep (default 1:24)
%   'SampleLength'     duration of each sample's averaging interval
%                      (default seconds(0) = instantaneous samples).
%                      Example: hours(1) for the hourly supersaturations.
%   'TimeLabel'        "start" (default) | "center" | "end": where the time
%                      stamp sits in its interval. EPCAPE files here are
%                      labelled by the interval START, except ARM RADFLUX,
%                      whose 1-min values are labelled by the END.
%   'MinOverlap'       fraction (0-1] of a sample's interval that must lie
%                      inside a window for the sample to be kept (default 0.5)
%   'SoledadCloud5min' timetable from soledad_cloud_condition_5min.m. If
%                      given, the in-cloud windows are reduced to the 5-min
%                      intervals with cloud == true. Ignored for pre- and
%                      post-cloud windows (those are out-of-cloud periods).
%   'TimeVariable'     name of the datetime variable when data is a table
%
% OUTPUTS:
%   keep      - logical column, true for samples inside a window
%   event_id  - event number (1-24) of each sample, NaN where keep is false.
%               In "allWindows" mode a sample can sit in two events (event
%               n's post-cloud window is often event n+1's pre-cloud window,
%               e.g. events 1 and 2); it gets the event with the larger
%               overlap, the lower number on a tie.
%   data_out  - the input with only the kept samples. For a timetable or
%               table, a variable dedrick2026_event_id is added.
%
% EXAMPLE:
%   SS = read_posted_product('.../bb2743661p_34_1_supersaturation.txt');   % hourly
%   [keep, id, SS_events] = filter_to_dedrick2026_events(SS, 'SampleLength', hours(1));
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement doi:10.1175/JAS-D-26-0004.s1
%     (Text S7, Tables S4 and S5).
%
% By Andrew John Buggee
%%

function [keep, event_id, data_out] = filter_to_dedrick2026_events(data, varargin)


% -------------------------------------------------------------------
% ------------------------ Parse the inputs -------------------------
% -------------------------------------------------------------------

p = inputParser;
p.FunctionName = 'filter_to_dedrick2026_events';
addParameter(p, 'Period', "inCloud", @(x) any(strcmpi(string(x), ["inCloud" "preCloud" "postCloud" "allWindows"])));
addParameter(p, 'Events', 1:24, @(x) isnumeric(x) && all(ismember(x, 1:24)));
addParameter(p, 'SampleLength', seconds(0), @(x) isduration(x) && isscalar(x) && x >= 0);
addParameter(p, 'TimeLabel', "start", @(x) any(strcmpi(string(x), ["start" "center" "end"])));
addParameter(p, 'MinOverlap', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'SoledadCloud5min', [], @(x) isempty(x) || istimetable(x));
addParameter(p, 'TimeVariable', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
opt = p.Results;

% The time stamps of the data, whatever form they come in
if isdatetime(data)
    t = data(:);
elseif istimetable(data)
    t = data.Properties.RowTimes;
elseif istable(data)
    if strlength(string(opt.TimeVariable)) == 0
        is_time = varfun(@isdatetime, data, 'OutputFormat', 'uniform');
        if nnz(is_time) ~= 1
            error('filter_to_dedrick2026_events: the table has %d datetime variables; name one with ''TimeVariable''.', nnz(is_time));
        end
        opt.TimeVariable = data.Properties.VariableNames{is_time};
    end
    t = data.(opt.TimeVariable);
else
    error('filter_to_dedrick2026_events: data must be a datetime vector, a timetable or a table.');
end

% Event times are UTC with no time zone. Convert zoned input to match.
if ~isempty(t.TimeZone)
    t.TimeZone = 'UTC';
    t.TimeZone = '';
end
t = t(:);



% -------------------------------------------------------------------
% ---------- The interval [s, e) that each sample covers ------------
% -------------------------------------------------------------------

% Work in seconds since 1970-01-01, which makes the overlap arithmetic
% simple and exact to well below a second.
t_s = posixtime(t);
L_s = seconds(opt.SampleLength);

if L_s == 0
    % An instantaneous sample is treated as a 1-ms interval starting at its
    % time stamp, so "overlap >= 50 %" means "start of window <= t < end".
    L_eff = 1e-3;
    s = t_s;
else
    L_eff = L_s;
    switch lower(string(opt.TimeLabel))
        case "start",  s = t_s;
        case "center", s = t_s - L_s/2;
        case "end",    s = t_s - L_s;
    end
end
e = s + L_eff;



% -------------------------------------------------------------------
% ---------------- The windows of each chosen event -----------------
% -------------------------------------------------------------------

events = dedrick2026_soledad_cloud_events();
events = events(ismember(events.event_id, opt.Events), :);

switch lower(string(opt.Period))
    case "precloud"
        win_start = events.pre_start_utc;   win_end = events.pre_end_utc;
    case "incloud"
        win_start = events.in_start_utc;    win_end = events.in_end_utc;
    case "postcloud"
        win_start = events.post_start_utc;  win_end = events.post_end_utc;
    case "allwindows"
        win_start = events.pre_start_utc;   win_end = events.post_end_utc;   % pre, in, post and the gaps between them
end

% Optional: reduce the in-cloud windows to the 5-min cloud intervals
use_cloud_intervals = ~isempty(opt.SoledadCloud5min) && strcmpi(string(opt.Period), "inCloud");
if use_cloud_intervals
    if ~ismember('cloud', opt.SoledadCloud5min.Properties.VariableNames)
        error('filter_to_dedrick2026_events: SoledadCloud5min must have a logical variable "cloud" (see soledad_cloud_condition_5min.m).');
    end
    cloud_start_s = posixtime(opt.SoledadCloud5min.Properties.RowTimes(opt.SoledadCloud5min.cloud));
    cloud_len_s   = 300;   % each row of the classification is a 5-min interval [t, t + 5 min)
end



% -------------------------------------------------------------------
% ----- Overlap of every sample with every event (one per loop) -----
% -------------------------------------------------------------------

n_samples = numel(t_s);
best_overlap = zeros(n_samples, 1);   % largest fraction of each sample inside any event
event_id = nan(n_samples, 1);

for k = 1:height(events)

    a = posixtime(win_start(k));   % window start (s)
    b = posixtime(win_end(k));     % window end (s)

    % The qualifying intervals of this event, as rows [start end]
    if use_cloud_intervals
        in_win = cloud_start_s >= a - cloud_len_s & cloud_start_s < b;
        iv = [cloud_start_s(in_win), cloud_start_s(in_win) + cloud_len_s];
        iv(:, 1) = max(iv(:, 1), a);   % clip the first and last 5-min interval to the window
        iv(:, 2) = min(iv(:, 2), b);
        iv = iv(iv(:, 2) > iv(:, 1), :);
        if isempty(iv)
            continue
        end
        iv = merge_intervals(iv);
    else
        iv = [a b];
    end

    % Only samples that could touch this event need the calculation
    near = e > iv(1, 1) & s < iv(end, 2);
    if ~any(near)
        continue
    end

    overlap_s = qualifying_time(iv, e(near)) - qualifying_time(iv, s(near));
    frac = overlap_s / L_eff;

    % Keep the event with the largest overlap (strictly larger, so ties go
    % to the earlier, lower-numbered event)
    idx = find(near);
    better = frac > best_overlap(idx);
    best_overlap(idx(better)) = frac(better);
    event_id(idx(better)) = events.event_id(k);
end

keep = best_overlap >= opt.MinOverlap - 1e-9;   % tolerance for float round-off at exactly 50 %
event_id(~keep) = NaN;



% -------------------------------------------------------------------
% -------------------- The filtered data, same type -----------------
% -------------------------------------------------------------------

if isdatetime(data)
    data_out = data(keep);
else
    data_out = data(keep, :);
    data_out.dedrick2026_event_id = event_id(keep);
end


end




% ===================================================================
% ========================= Local functions =========================
% ===================================================================

function C = qualifying_time(iv, x)
% Total time (s) inside the disjoint, sorted intervals iv = [a_k b_k]
% before each time x. C is piecewise linear in x with knots at every a_k
% and b_k: it rises with slope 1 inside an interval and stays flat between
% intervals, so linear interpolation through the knots is exact.
lengths = iv(:, 2) - iv(:, 1);
C_at_a = [0; cumsum(lengths(1:end-1))];   % C at each interval start
C_at_b = cumsum(lengths);                 % C at each interval end
knots_t = reshape(iv', [], 1);            % a1 b1 a2 b2 ...
knots_C = reshape([C_at_a C_at_b]', [], 1);
C = interp1(knots_t, knots_C, x, 'linear');
C(x <= iv(1, 1))   = 0;                   % before the first interval
C(x >= iv(end, 2)) = C_at_b(end);         % after the last interval
end


function iv = merge_intervals(iv)
% Merge touching or overlapping intervals so that every knot in
% qualifying_time is strictly increasing (interp1 needs distinct knots).
iv = sortrows(iv, 1);
out = iv(1, :);
for k = 2:size(iv, 1)
    if iv(k, 1) <= out(end, 2)
        out(end, 2) = max(out(end, 2), iv(k, 2));   % touches or overlaps the last one: extend it
    else
        out(end + 1, :) = iv(k, :); %#ok<AGROW>
    end
end
iv = out;
end
