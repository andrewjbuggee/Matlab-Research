%% The 24 Mt. Soledad cloud events of Dedrick et al. (2026), with their pre-, in- and post-cloud windows
%
% DESCRIPTION:
%   Returns the 24 "Mt. Soledad Cloud Events" analysed by Dedrick et al.
%   (2026) as a table, one row per event. The times are transcribed
%   verbatim from Table S5 of the paper's supplement (p. 14 of the
%   supplement PDF), which gives them in UTC. They were copied from the
%   PDF's text layer with pdftotext (not retyped), and the 24 cloud-hour
%   values add up to 213.7 h, the "214 h" quoted in the paper (p. 1064).
%
%   WHAT AN EVENT IS (supplement Text S7 and Table S4):
%     At Mt. Soledad (250 m MSL) the aerosol inlets switched to the
%     ground-based counterflow virtual impactor (GCVI), which samples the
%     residual particles of cloud droplets larger than ~9 um, whenever the
%     visibility was below 5000 m and no precipitation was detected
%     (Dedrick et al. 2026, p. 1064). Starting from these GCVI periods:
%       1. GCVI sampling continuously for > 1 h ............. 154 events
%       2. events within 6 h of each other merged .............. 95 events
%       3. cloud duration > 2 h, with "cloud" = every 5-min
%          interval with visibility < 1 km (no 1-s value above
%          1 km in the interval) and FM-120 LWC > 0.01 g m-3 .... 44 events
%       4. SEMS and FM-120 data available for >= 50 % of the
%          event ................................................ 24 events
%     The 24 survivors total 214 cloud hours, 9 +/- 5 h per event. The
%     algorithm is credited to Williams et al. (2026, submitted).
%
%   THE THREE WINDOWS OF EACH EVENT (Table S5):
%     pre-cloud   ambient (out-of-cloud) sampling before the event
%     in-cloud    the GCVI event itself: Mt. Soledad is inside the cloud
%     post-cloud  ambient sampling after the event
%   NOTE ON "CLOUD HOURS" (checked 27 Oct 2026): the Table S5 caption says
%   the cloud and haze hours count the 5-min intervals that meet the cloud
%   or haze condition. But for all 24 events cloud_hours equals the length
%   of the in-cloud window to within +/- 5 min (event 6: window 13.33 h,
%   cloud 13.33 h; event 24: 24.31 h vs 24.25 h), and the 24 windows add up
%   to 214.2 h, the paper's "214 h". So cloud_hours is in effect the
%   duration of the GCVI event, haze included (cloud + haze hours exceed
%   the window for most events). Counting only the 5-min intervals that
%   meet the cloud condition (soledad_cloud_condition_5min.m) gives 137 h
%   for the 24 windows. Table S4's step totals (342, 272 h), on the
%   other hand, are consistent with counts of cloud intervals (see
%   find_soledad_cloud_events.m). Use filter_to_dedrick2026_events.m with
%   'SoledadCloud5min' to keep only the cloud intervals of each window.
%
% INPUTS:
%   none
%
% OUTPUTS:
%   events - table with one row per event (24 rows) and the columns
%       event_id          event number 1-24 as in Table S5
%       pre_start_utc     start of the pre-cloud window  (datetime, UTC)
%       pre_end_utc       end   of the pre-cloud window
%       in_start_utc      start of the in-cloud window (GCVI event start)
%       in_end_utc        end   of the in-cloud window
%       post_start_utc    start of the post-cloud window
%       post_end_utc      end   of the post-cloud window
%       cloud_hours       Table S5 "cloud hours" (= the window length to
%                         within 5 min; see the note above)
%       haze_hours        Table S5 hours meeting the haze condition
%       in_window_hours   length of the in-cloud window (h), computed here
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: Quantitative estimates of upwind
%     and on-site particle scavenging and retrieved cloud supersaturations
%     from urban coastal stratocumulus during EPCAPE. J. Atmos. Sci., 83,
%     1061-1082, doi:10.1175/JAS-D-26-0004.1. Supplement:
%     doi:10.1175/JAS-D-26-0004.s1 (Text S7, Tables S4 and S5).
%
% By Andrew John Buggee
%%

function events = dedrick2026_soledad_cloud_events()


% -------------------------------------------------------------------
% ------------- Table S5 of the supplement, verbatim (UTC) ----------
% -------------------------------------------------------------------

% Columns: event, pre-cloud start, pre-cloud end, in-cloud start,
% in-cloud end, post-cloud start, post-cloud end, cloud hours, haze hours
S5 = {
     1, '26-Apr-2023 13:05:00', '27-Apr-2023 00:45:00', '27-Apr-2023 01:03:49', '27-Apr-2023 16:53:22', '27-Apr-2023 17:30:00', '27-Apr-2023 23:30:00', 15.83, 1.50
     2, '27-Apr-2023 17:30:00', '27-Apr-2023 23:30:00', '28-Apr-2023 05:07:53', '28-Apr-2023 17:36:44', '28-Apr-2023 18:45:00', '29-Apr-2023 00:45:00', 12.50, 2.75
     3, '30-Apr-2023 00:45:00', '30-Apr-2023 08:05:00', '30-Apr-2023 08:51:11', '30-Apr-2023 16:34:40', '30-Apr-2023 22:50:00', '01-May-2023 04:30:00',  7.67, 0.83
     4, '13-May-2023 17:45:00', '14-May-2023 05:00:00', '14-May-2023 05:40:11', '14-May-2023 14:44:26', '14-May-2023 17:10:00', '15-May-2023 02:40:00',  9.00, 3.83
     5, '15-May-2023 19:40:00', '16-May-2023 06:40:00', '16-May-2023 07:39:45', '16-May-2023 16:42:56', '16-May-2023 17:00:00', '17-May-2023 04:40:00',  9.08, 1.92
     6, '16-May-2023 18:50:00', '17-May-2023 06:45:00', '17-May-2023 06:45:38', '17-May-2023 20:05:25', '17-May-2023 22:10:00', '18-May-2023 08:05:00', 13.33, 1.75
     7, '17-May-2023 22:20:00', '18-May-2023 08:40:00', '18-May-2023 10:18:17', '18-May-2023 15:19:09', '18-May-2023 16:25:00', '19-May-2023 03:15:00',  5.00, 2.00
     8, '08-Jun-2023 21:55:00', '09-Jun-2023 09:20:00', '09-Jun-2023 09:51:21', '09-Jun-2023 12:21:02', '09-Jun-2023 13:25:00', '10-Jun-2023 00:20:00',  2.50, 0.08
     9, '16-Jun-2023 23:50:00', '17-Jun-2023 11:45:00', '17-Jun-2023 11:47:46', '17-Jun-2023 17:19:32', '17-Jun-2023 20:00:00', '18-Jun-2023 01:50:00',  5.50, 1.25
    10, '17-Jun-2023 20:00:00', '18-Jun-2023 01:50:00', '18-Jun-2023 03:38:57', '18-Jun-2023 12:09:53', '18-Jun-2023 15:00:00', '19-Jun-2023 00:05:00',  8.50, 3.42
    11, '01-Jul-2023 22:45:00', '02-Jul-2023 03:35:00', '02-Jul-2023 10:44:01', '02-Jul-2023 17:04:55', '02-Jul-2023 21:40:00', '03-Jul-2023 03:15:00',  6.33, 1.58
    12, '02-Jul-2023 21:40:00', '03-Jul-2023 03:15:00', '03-Jul-2023 07:11:14', '03-Jul-2023 19:17:36', '03-Jul-2023 20:25:00', '04-Jul-2023 03:05:00', 12.08, 1.08
    13, '03-Jul-2023 20:25:00', '04-Jul-2023 03:05:00', '04-Jul-2023 03:52:33', '04-Jul-2023 18:09:29', '04-Jul-2023 18:25:00', '05-Jul-2023 06:05:00', 14.25, 1.00
    14, '11-Jul-2023 19:25:00', '12-Jul-2023 07:20:00', '12-Jul-2023 07:20:29', '12-Jul-2023 14:48:25', '12-Jul-2023 14:55:00', '13-Jul-2023 02:45:00',  7.42, 0.50
    15, '12-Jul-2023 23:25:00', '13-Jul-2023 11:20:00', '13-Jul-2023 11:20:48', '13-Jul-2023 14:33:33', '13-Jul-2023 15:20:00', '14-Jul-2023 02:30:00',  3.17, 0.33
    16, '13-Jul-2023 19:00:00', '14-Jul-2023 04:35:00', '14-Jul-2023 06:56:42', '14-Jul-2023 11:17:58', '14-Jul-2023 14:50:00', '14-Jul-2023 23:15:00',  4.33, 1.00
    17, '19-Jul-2023 21:15:00', '20-Jul-2023 03:50:00', '20-Jul-2023 08:26:14', '20-Jul-2023 14:04:33', '20-Jul-2023 15:15:00', '21-Jul-2023 02:00:00',  5.58, 2.58
    18, '18-Aug-2023 19:50:00', '19-Aug-2023 04:10:00', '19-Aug-2023 07:45:09', '19-Aug-2023 15:57:41', '19-Aug-2023 16:00:00', '20-Aug-2023 03:55:00',  8.17, 1.83
    19, '23-Aug-2023 14:10:00', '23-Aug-2023 19:40:00', '23-Aug-2023 23:40:25', '24-Aug-2023 14:44:57', '24-Aug-2023 14:45:00', '25-Aug-2023 02:40:00', 15.00, 3.92
    20, '25-Aug-2023 00:40:00', '25-Aug-2023 12:30:00', '25-Aug-2023 12:38:51', '25-Aug-2023 18:09:49', '25-Aug-2023 18:20:00', '26-Aug-2023 06:05:00',  5.50, 0.75
    21, '25-Aug-2023 19:30:00', '26-Aug-2023 06:20:00', '26-Aug-2023 07:26:24', '26-Aug-2023 20:25:08', '27-Aug-2023 00:40:00', '27-Aug-2023 08:25:00', 13.00, 2.58
    22, '03-Sep-2023 00:40:00', '03-Sep-2023 12:35:00', '03-Sep-2023 12:39:46', '03-Sep-2023 14:56:12', '03-Sep-2023 15:00:00', '04-Sep-2023 02:55:00',  2.33, 0.08
    23, '03-Oct-2023 19:55:00', '04-Oct-2023 07:50:00', '04-Oct-2023 07:52:47', '04-Oct-2023 11:17:04', '04-Oct-2023 11:20:00', '04-Oct-2023 23:15:00',  3.42, 0.25
    24, '08-Oct-2023 04:30:00', '08-Oct-2023 10:35:00', '08-Oct-2023 16:25:57', '09-Oct-2023 16:44:37', '09-Oct-2023 17:30:00', '10-Oct-2023 04:40:00', 24.25, 5.08
    };



% -------------------------------------------------------------------
% ------------------ Convert the text to a table --------------------
% -------------------------------------------------------------------

% All times are UTC (Table S5 caption). The datetimes carry no time zone,
% the same convention as every other EPCAPE reader in this repository.
fmt = 'dd-MMM-yyyy HH:mm:ss';
to_time = @(col) datetime(S5(:, col), 'InputFormat', fmt, 'Locale', 'en_US');

event_id        = cell2mat(S5(:, 1));
pre_start_utc   = to_time(2);
pre_end_utc     = to_time(3);
in_start_utc    = to_time(4);
in_end_utc      = to_time(5);
post_start_utc  = to_time(6);
post_end_utc    = to_time(7);
cloud_hours     = cell2mat(S5(:, 8));
haze_hours      = cell2mat(S5(:, 9));

% Length of the in-cloud (GCVI) window
in_window_hours = hours(in_end_utc - in_start_utc);

events = table(event_id, pre_start_utc, pre_end_utc, in_start_utc, in_end_utc, ...
               post_start_utc, post_end_utc, cloud_hours, haze_hours, in_window_hours);



% -------------------------------------------------------------------
% ----------------------- Consistency checks ------------------------
% -------------------------------------------------------------------

% The windows of each event must be ordered pre <= in <= post, and the
% cloud hours must add up to the 214 h stated in the paper (p. 1064).
assert(all(events.pre_start_utc < events.pre_end_utc & events.pre_end_utc <= events.in_start_utc ...
           & events.in_start_utc < events.in_end_utc & events.in_end_utc <= events.post_start_utc ...
           & events.post_start_utc < events.post_end_utc), ...
       'dedrick2026_soledad_cloud_events: a window in Table S5 is out of order');
assert(abs(sum(events.cloud_hours) - 214) < 0.5, ...
       'dedrick2026_soledad_cloud_events: cloud hours do not add up to 214 h');
% Table S5's cloud hours equal the window length to within one 5-min
% interval (see the note in the header); a larger difference would mean a
% transcription error in the times
assert(all(abs(events.cloud_hours - events.in_window_hours) <= 5/60 + 1e-6), ...
       'dedrick2026_soledad_cloud_events: cloud hours differ from the window length by > 5 min');


end
