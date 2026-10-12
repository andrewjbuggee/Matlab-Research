function t_floor = floor_time(t, step_min)
% FLOOR_TIME  Round times DOWN to the start of their step_min-minute interval.
%
%   t_floor = floor_time(t, step_min)
%
%   Example: with step_min = 5,  12:03:41 -> 12:00:00  and  12:05:00 -> 12:05:00.
%   Intervals are counted from midnight (00:00, 00:05, 00:10, ...), which is
%   how the FM-120 5-min files and the 15-min file label their intervals
%   (each interval [t, t + step) is labelled by its START time t).
%
%   Inputs   t         datetime vector
%            step_min  interval length in minutes (e.g. 5, 15, 60)
%   Output   t_floor   datetime vector, same size as t
%
%   How: convert to seconds since 1970-01-01 ("posix time"), round down to a
%   multiple of the step, convert back.

step_s  = 60 * step_min;
sec     = posixtime(t);                      % seconds since 1970-01-01 00:00
t_floor = datetime(floor(sec / step_s) * step_s, 'ConvertFrom', 'posixtime');

end
