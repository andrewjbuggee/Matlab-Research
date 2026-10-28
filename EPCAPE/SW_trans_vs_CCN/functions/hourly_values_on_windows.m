%% Put hourly values (CCN, supersaturation) onto the averaging windows, over the passing minutes
%
% DESCRIPTION:
%   The Pier CCN spectra and the cloud supersaturations are hourly: the
%   CCN-200 steps through its six setpoints in a ~1-h cycle (~10 min each;
%   Dedrick et al. 2026, p. 1063), so one setpoint is measured about once
%   an hour, and Dedrick et al. (2026) averaged every supersaturation
%   retrieval to 1 h. Neither can be resolved more finely than an hour.
%
%   To pair them with the window-averaged T_SW, every minute takes the
%   value of the hour [h, h + 1) that holds it, and the window value is the
%   mean over the window's PASSING minutes, the same minutes behind T_SW:
%     * windows of <= 60 min that divide the hour (10, 15, 20, 30, 60):
%       every window lies inside one hour and gets that hour's value, so all
%       the windows of one hour share the same CCN and supersaturation
%     * longer windows, or windows that straddle an hour boundary: the mean
%       of the hourly values, weighted by the passing minutes in each hour
%   A window is NaN if any of its passing minutes falls in an hour without
%   a value, so an average never silently covers part of the window.
%
% INPUTS:
%   t_min      - the 1-min grid (datetime column, minute start, UTC)
%   pass       - logical, one per minute: in the case set AND passing every filter
%   t_hour     - hour starts of the hourly values (datetime column, UTC,
%                consecutive hours)
%   x_hour     - hourly values (column, same length as t_hour; NaN = none)
%   window_min - averaging window length (minutes)
%
% OUTPUTS:
%   x_win - one value per window, on the window_mean.m grid of t_min
%   t_win - the window starts
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1 (Pier CCN-200 cycle, p. 1063; hourly
%     supersaturations, Sec. 2d).
%
% DEPENDENCIES: EPCAPE/SW_transmittance/window_mean.m
%
% By Andrew John Buggee
%%

function [x_win, t_win] = hourly_values_on_windows(t_min, pass, t_hour, x_hour, window_min)


% the hour that holds each minute (hour 1 = t_hour(1))
k = floor(minutes(t_min - t_hour(1)) / 60) + 1;
x_min = nan(size(t_min));
ok = k >= 1 & k <= numel(t_hour);
x_min(ok) = x_hour(k(ok));
x_min(~pass) = NaN;                          % only the passing minutes count

% passing minutes per window, and those with an hourly value
p = double(pass);  p(~pass) = NaN;
[t_win, ~, n_pass] = window_mean(t_min, p, window_min);
[~, x_win, n_val]  = window_mean(t_min, x_min, window_min);

x_win(n_val < n_pass) = NaN;                 % part of the passing minutes has no value


end
