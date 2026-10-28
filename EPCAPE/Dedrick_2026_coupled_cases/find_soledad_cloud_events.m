%% Rebuild the Mt. Soledad cloud-event selection of Dedrick et al. (2026) with adjustable criteria
%
% DESCRIPTION:
%   Re-runs, step by step, the selection that turned the GCVI record into
%   the 24 Mt. Soledad cloud events (supplement Table S4), with every
%   threshold an input, so the effect of relaxing a criterion can be seen:
%
%     step 1  GCVI sampling continuously for > MinGcviDuration_h (1 h)
%                                             paper: 154 periods
%     step 2  periods within MergeGap_h (6 h) of each other merged into one
%             event                           paper: 95 events, 342 cloud h, 204 haze h
%     step 3  cloud duration > MinCloudHours (2 h), where cloud = the 5-min
%             intervals of soledad_cloud_condition_5min.m (visibility < 1 km
%             and FM-120 LWC > 0.01 g m-3)
%                                             paper: 44 events, 272 cloud h, 94 haze h
%     step 4  data for >= MinDataCoverage (50 %) of the event
%                                             paper: 24 events, 214 cloud h, 42 haze h
%
%   WHAT CANNOT BE REPRODUCED HERE
%   * Step 4 needs SEMS (residual size distribution) AND FM-120 coverage.
%     No SEMS data are in the local copy, so only the FM-120 coverage is
%     tested; the step-4 count here is therefore an UPPER bound.
%   * Table S4 does not say whether "cloud duration > 2 h" is the total of
%     the cloud intervals or one continuous stretch. The total is used; it
%     gives 46 events and 318 cloud hours at step 3 (paper: 44 and 272),
%     whereas the event WINDOW length > 2 h would give 76 events.
%   * The posted GCVI file has no sampling between 00:32 and 16:32 UTC on
%     30 Apr 2023, so the paper's event 3 cannot be found from it.
%   * The algorithm is Williams et al. (2026, submitted); details beyond the
%     table (e.g. how a segment's start is defined, how gaps in the
%     visibility record are treated) may differ.
%   The output marks which rebuilt event matches which of the paper's 24
%   (largest overlap of the in-cloud windows), so the reconstruction can be
%   judged against Table S5.
%
% INPUTS:
%   segs   - table from gcvi_sampling_segments.m (start_utc, end_utc, duration_h)
%   cloud5 - timetable from soledad_cloud_condition_5min.m
%
%   Optional name-value pairs (defaults = the paper's criteria):
%   'MinGcviDuration_h'  1
%   'MergeGap_h'         6
%   'MinCloudHours'      2
%   'MinDataCoverage'    0.5   fraction of the event's 5-min intervals
%                              with an FM-120 value
%
% OUTPUTS:
%   events  - table, one row per merged event (step 2), with
%             start_utc, end_utc, window_h, n_gcvi_periods, cloud_hours,
%             haze_hours, fm120_coverage, pass_cloud (step 3),
%             pass_coverage (step 4), selected (passes 3 and 4),
%             dedrick_event_id (the matching Table S5 event, NaN if none),
%             start_diff_min, end_diff_min (vs. Table S5)
%   summary - table: step, criterion, n_events, cloud_hours, haze_hours, and
%             the paper's numbers for comparison
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement Text S7, Tables S4 and S5.
%
% By Andrew John Buggee
%%

function [events, summary] = find_soledad_cloud_events(segs, cloud5, varargin)


% -------------------------------------------------------------------
% ------------------------ Parse the inputs -------------------------
% -------------------------------------------------------------------

p = inputParser;
addParameter(p, 'MinGcviDuration_h', 1);
addParameter(p, 'MergeGap_h', 6);
addParameter(p, 'MinCloudHours', 2);
addParameter(p, 'MinDataCoverage', 0.5);
parse(p, varargin{:});
opt = p.Results;



% -------------------------------------------------------------------
% -------- Step 1: GCVI periods longer than MinGcviDuration_h -------
% -------------------------------------------------------------------

long = segs(segs.duration_h > opt.MinGcviDuration_h, :);
n_step1 = height(long);



% -------------------------------------------------------------------
% ------------ Step 2: merge periods within MergeGap_h --------------
% -------------------------------------------------------------------

% a period joins the event before it if it starts within MergeGap_h of
% that event's end
ev_start = long.start_utc(1);  ev_end = long.end_utc(1);  n_periods = 1;
starts = NaT(0, 1);  ends = NaT(0, 1);  counts = zeros(0, 1);
for k = 2:height(long)
    if long.start_utc(k) - ev_end <= hours(opt.MergeGap_h)
        ev_end = max(ev_end, long.end_utc(k));
        n_periods = n_periods + 1;
    else
        starts(end + 1, 1) = ev_start;  ends(end + 1, 1) = ev_end;  counts(end + 1, 1) = n_periods; %#ok<AGROW>
        ev_start = long.start_utc(k);  ev_end = long.end_utc(k);  n_periods = 1;
    end
end
starts(end + 1, 1) = ev_start;  ends(end + 1, 1) = ev_end;  counts(end + 1, 1) = n_periods;

events = table(starts, ends, hours(ends - starts), counts, ...
    'VariableNames', {'start_utc', 'end_utc', 'window_h', 'n_gcvi_periods'});



% -------------------------------------------------------------------
% ---- Cloud and haze hours, and FM-120 coverage, of each event -----
% -------------------------------------------------------------------

% a 5-min interval [t, t + 5 min) belongs to an event when it overlaps the
% event's window (Table S5 counts whole intervals, see
% dedrick2026_soledad_cloud_events.m)
t5 = cloud5.Properties.RowTimes;
n_ev = height(events);
events.cloud_hours    = zeros(n_ev, 1);
events.haze_hours     = zeros(n_ev, 1);
events.fm120_coverage = zeros(n_ev, 1);
for k = 1:n_ev
    in = t5 + minutes(5) > events.start_utc(k) & t5 < events.end_utc(k);
    events.cloud_hours(k)    = nnz(cloud5.cloud(in)) * 5 / 60;
    events.haze_hours(k)     = nnz(cloud5.haze(in)) * 5 / 60;
    events.fm120_coverage(k) = mean(cloud5.has_fm120(in));
end



% -------------------------------------------------------------------
% --------------------- Steps 3 and 4 ------------------------------
% -------------------------------------------------------------------

events.pass_cloud    = events.cloud_hours > opt.MinCloudHours;
events.pass_coverage = events.fm120_coverage >= opt.MinDataCoverage;
events.selected      = events.pass_cloud & events.pass_coverage;



% -------------------------------------------------------------------
% ------------- Match each rebuilt event to Table S5 ----------------
% -------------------------------------------------------------------

paper = dedrick2026_soledad_cloud_events();
events.dedrick_event_id = nan(n_ev, 1);
events.start_diff_min   = nan(n_ev, 1);
events.end_diff_min     = nan(n_ev, 1);
for k = 1:n_ev
    overlap_h = hours(min(events.end_utc(k), paper.in_end_utc) - max(events.start_utc(k), paper.in_start_utc));
    [best, j] = max(overlap_h);
    if best > 0
        events.dedrick_event_id(k) = paper.event_id(j);
        events.start_diff_min(k)   = minutes(events.start_utc(k) - paper.in_start_utc(j));
        events.end_diff_min(k)     = minutes(events.end_utc(k) - paper.in_end_utc(j));
    end
end



% -------------------------------------------------------------------
% ------------------------- Step summary ----------------------------
% -------------------------------------------------------------------

step = (1:4)';
criterion = [ ...
    sprintf("GCVI sampling > %g h", opt.MinGcviDuration_h); ...
    sprintf("periods within %g h merged", opt.MergeGap_h); ...
    sprintf("+ cloud > %g h (vis < 1 km, LWC > 0.01)", opt.MinCloudHours); ...
    sprintf("+ FM-120 data >= %g %% of event (SEMS not tested)", 100 * opt.MinDataCoverage)];
n_events = [n_step1; n_ev; nnz(events.pass_cloud); nnz(events.selected)];
cloud_hours = [NaN; sum(events.cloud_hours); sum(events.cloud_hours(events.pass_cloud)); sum(events.cloud_hours(events.selected))];
haze_hours  = [NaN; sum(events.haze_hours);  sum(events.haze_hours(events.pass_cloud));  sum(events.haze_hours(events.selected))];
paper_n_events    = [154; 95; 44; 24];
paper_cloud_hours = [NaN; 342; 272; 214];
paper_haze_hours  = [NaN; 204; 94; 42];
summary = table(step, criterion, n_events, cloud_hours, haze_hours, paper_n_events, paper_cloud_hours, paper_haze_hours);


end
