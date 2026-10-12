function is_near = near_rain(rain_times, t, window_min)
% NEAR_RAIN  True for each time that lies within ±window_min minutes of any rain.
%
%   is_near = near_rain(rain_times, t, window_min)
%
%   Why: a microwave radiometer's LWP is far too high while its radome is
%   wet, because liquid on the radome emits at the instrument's frequencies
%   (23.8 and 31.4 GHz). Every LWP sample within ±30 min of rain at the
%   Pier disdrometer is therefore removed.
%
%   Inputs
%     rain_times  datetime vector of the disdrometer samples with rain rate > 0
%     t           datetime vector of the LWP samples
%     window_min  half-width of the window in minutes (30)
%   Output
%     is_near     logical vector, true where at least one rain sample falls in
%                 [t - window_min, t + window_min]
%
%   How: with the rain times sorted, the number of rain samples inside the
%   window [lo, hi] is  count(<= hi) - count(< lo)  (count_sorted.m).

rain_s = sort(posixtime(rain_times(:)));   % rain sample times in seconds, sorted
t_s    = posixtime(t(:));                  % LWP sample times in seconds
lo     = t_s - 60 * window_min;            % start of each window
hi     = t_s + 60 * window_min;            % end of each window

n_rain_in_window = count_sorted(rain_s, hi, 'at_or_below') - count_sorted(rain_s, lo, 'below');
is_near = reshape(n_rain_in_window > 0, size(t));

end
