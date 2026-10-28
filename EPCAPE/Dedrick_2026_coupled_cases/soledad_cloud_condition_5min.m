%% Classify every 5-min interval at Mt. Soledad as cloud, haze or out-of-cloud (Dedrick et al. 2026, Table S4)
%
% DESCRIPTION:
%   The cloud definition behind the Mt. Soledad cloud events (supplement
%   Table S4, "Cloud Conditions", "any 5 min point where"):
%     cloud         visibility < 1 km and FM-120 LWC > 0.01 g m-3, and the
%                   1-s visibility "could not exceed 1 km at any point
%                   during the 5 min average"
%     haze          1 km < visibility < 5 km and LWC < 0.01 g m-3
%     out-of-cloud  visibility > 5 km
%   For the year they report 305 cloud hours, 180 haze hours and 7005
%   out-of-cloud hours.
%
%   HOW EACH TEST IS EVALUATED HERE
%     * cloud: the MAXIMUM 1-s visibility in the interval < CloudVisMax_m
%       (1000 m), which is the "could not exceed 1 km at any point" rule,
%       AND the FM-120 5-min LWC > CloudLwcMin_gm3 (0.01 g m-3).
%     * haze and out-of-cloud: Table S4 does not say which 5-min statistic of
%       the visibility is meant; the MEDIAN is used (my choice). haze =
%       CloudVisMax_m < median < HazeVisMax_m and LWC <= CloudLwcMin_gm3
%       (or no FM-120 value); out-of-cloud = median > HazeVisMax_m.
%     * an interval needs >= MinVisSamples valid 1-s visibility samples (a
%       full interval holds ~300); fewer -> every class is false and
%       has_vis is false. Without an FM-120 value an interval cannot be
%       cloud (has_fm120 is false).
%
%   WHY THESE THRESHOLDS
%   Visibility < 1 km is the WMO distinction between fog and mist (Dedrick
%   et al. 2026, supplement Text S7, citing WMO 2008); at 250 m MSL, being
%   "in fog" means being inside the stratocumulus deck. LWC > 0.01 g m-3
%   requires real cloud droplets (2-50 um) in the FM-120, which separates
%   cloud from haze of swollen aerosol (low visibility, little liquid).
%
% INPUTS:
%   vis - struct from read_soledad_visibility.m (time, m), ~1-s, UTC
%   fm  - struct from read_fm120_daily_files.m (time = 5-min interval START,
%         lwc_gm3, nd_cm3), UTC
%
%   Optional name-value pairs:
%   'CloudVisMax_m'    1000   (m)
%   'HazeVisMax_m'     5000   (m)
%   'CloudLwcMin_gm3'  0.01   (g m-3)
%   'MinVisSamples'    125    (the Python visibility_on_intervals default)
%   'TimeRange'        [datetime(2023,2,15) datetime(2024,2,15)], [start end)
%
% OUTPUTS:
%   C - timetable, one row per 5-min interval [t, t + 5 min) (row time =
%       interval start, UTC), with variables
%       vis_max_m, vis_median_m, n_vis, lwc_gm3, nd_cm3,
%       has_vis, has_fm120, cloud, haze, out_of_cloud (logical)
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement Text S7 and Table S4 (the
%     algorithm is credited to Williams et al. 2026, submitted).
%   Python counterpart of the interval statistics: quantities.py
%     (visibility_on_intervals).
%
% By Andrew John Buggee
%%

function C = soledad_cloud_condition_5min(vis, fm, varargin)


% -------------------------------------------------------------------
% ------------------------ Parse the inputs -------------------------
% -------------------------------------------------------------------

p = inputParser;
addParameter(p, 'CloudVisMax_m', 1000);
addParameter(p, 'HazeVisMax_m', 5000);
addParameter(p, 'CloudLwcMin_gm3', 0.01);
addParameter(p, 'MinVisSamples', 125);
addParameter(p, 'TimeRange', [datetime(2023, 2, 15) datetime(2024, 2, 15)]);
parse(p, varargin{:});
opt = p.Results;

% the 5-min grid: every interval start in [TimeRange(1), TimeRange(2))
t0 = opt.TimeRange(1);
n_int = floor(minutes(opt.TimeRange(2) - t0) / 5);
t_grid = t0 + minutes(5 * (0:n_int - 1)');



% -------------------------------------------------------------------
% ------- Visibility: max, median and count in each interval --------
% -------------------------------------------------------------------

% interval index of every 1-s sample (interval [t, t + 5 min) holds t)
k = floor(seconds(vis.time - t0) / 300) + 1;
ok = isfinite(vis.m) & k >= 1 & k <= n_int;
k = k(ok);
v = vis.m(ok);

n_vis        = accumarray(k, 1, [n_int 1], @sum, 0);
vis_max_m    = accumarray(k, v, [n_int 1], @max, NaN);
vis_median_m = accumarray(k, v, [n_int 1], @median, NaN);

has_vis = n_vis >= opt.MinVisSamples;
vis_max_m(~has_vis)    = NaN;
vis_median_m(~has_vis) = NaN;



% -------------------------------------------------------------------
% ----------------- FM-120 LWC and N_d on the grid ------------------
% -------------------------------------------------------------------

% the FM-120 files are labelled by interval start, on the same 5-min grid
lwc_gm3 = nan(n_int, 1);
nd_cm3  = nan(n_int, 1);
kf = round(minutes(fm.time - t0) / 5) + 1;
inside = kf >= 1 & kf <= n_int;
lwc_gm3(kf(inside)) = fm.lwc_gm3(inside);
nd_cm3(kf(inside))  = fm.nd_cm3(inside);
has_fm120 = isfinite(lwc_gm3);



% -------------------------------------------------------------------
% --------------------------- Classify -------------------------------
% -------------------------------------------------------------------

cloud = has_vis & has_fm120 & vis_max_m < opt.CloudVisMax_m & lwc_gm3 > opt.CloudLwcMin_gm3;

low_lwc = ~has_fm120 | lwc_gm3 <= opt.CloudLwcMin_gm3;
haze = has_vis & ~cloud & vis_median_m > opt.CloudVisMax_m & vis_median_m < opt.HazeVisMax_m & low_lwc;

out_of_cloud = has_vis & vis_median_m > opt.HazeVisMax_m;

C = timetable(t_grid, vis_max_m, vis_median_m, n_vis, lwc_gm3, nd_cm3, has_vis, has_fm120, ...
              cloud, haze, out_of_cloud);
C.Properties.DimensionNames{1} = 'time';


end
