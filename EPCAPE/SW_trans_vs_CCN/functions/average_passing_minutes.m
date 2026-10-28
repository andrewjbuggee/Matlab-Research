%% Window-averaged SW transmittance, LWP and sun angle from the minutes that pass every filter
%
% DESCRIPTION:
%   Turns the 1-min grid into consecutive, non-overlapping averaging
%   windows [t, t + window_min) (UTC, labelled by the window start). The
%   windows are anchored at midnight (window_mean.m), so 10-min windows
%   start at :00, :10, :20, ... and 60-min windows are the clock hours.
%   In each window, over ONLY the minutes that pass every filter:
%     * T_sw: sw_transmittance.m (SKYRAD SW / RADFLUX clear-sky SW, ratio of
%       the window's means), so cloud-free or rainy minutes inside an
%       otherwise good window do not enter the ratio
%     * lwp_gm2: mean MWRRET v2 LWP of the same minutes
%     * sza_deg: solar zenith angle of the mean cos(SZA) of those minutes
%     * n_pass: how many passing minutes had valid SW and clear-sky SW
%   A window is kept (finite T_sw) only if at least min_pass_minutes of its
%   minutes pass every filter AND have valid SW and clear-sky values.
%
%   The persistence filter is evaluated minute by minute before this step
%   (analysis_filters_1min.m), so the averaging window and the persistence
%   window are independent: e.g. a +/- 30 min persistence test with 10-min
%   averages gives six points for an hour in which every minute passes.
%
% INPUTS:
%   D                - 1-min timetable (build_epcape_1min_inputs.m)
%   pass             - logical, one per minute: in the case set AND passing
%                      every filter
%   window_min       - averaging window length (minutes), e.g. 10 or 60
%   min_pass_minutes - fewest passing minutes for a window to count
%   szaMax_deg       - the daytime cut handed to sw_transmittance.m
%
% OUTPUTS:
%   W - timetable, one row per window (row time = window start, UTC):
%       T_sw, lwp_gm2, sza_deg, n_pass
%
% DEPENDENCIES: EPCAPE/SW_transmittance (sw_transmittance.m, window_mean.m)
%
% By Andrew John Buggee
%%

function W = average_passing_minutes(D, pass, window_min, min_pass_minutes, szaMax_deg)


% the 1-min grid in the form sw_transmittance.m expects (times = minute start)
rad.time         = D.time;
rad.time_end     = D.time + minutes(1);
rad.sw_Wm2       = D.sw_Wm2;
rad.sw_clear_Wm2 = D.sw_clear_Wm2;
rad.cos_sza      = D.cos_sza;

[T_sw, t_win, info] = sw_transmittance(rad, 'SzaMax_deg', szaMax_deg, 'WindowMin', window_min, ...
                                       'MinCount', min_pass_minutes, 'Keep', pass);

% LWP and sun angle over the same passing minutes
lwp = D.lwp_gm2;  lwp(~pass) = NaN;
mu  = D.cos_sza;  mu(~pass)  = NaN;
[t_l, lwp_w] = window_mean(D.time, lwp, window_min);
[~,   mu_w]  = window_mean(D.time, mu, window_min);

% window_mean starts at the first sample's day in both cases, so the windows line up
assert(isequal(t_l, t_win), 'average_passing_minutes: window grids differ');

W = timetable(t_win, T_sw, lwp_w, acosd(mu_w), info.n_per_window, ...
              'VariableNames', {'T_sw', 'lwp_gm2', 'sza_deg', 'n_pass'});
W.Properties.DimensionNames{1} = 'time';


end
