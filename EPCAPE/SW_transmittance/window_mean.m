%% Average a time series over consecutive, non-overlapping windows of a fixed length
%
% DESCRIPTION:
%   MATLAB version of window_mean (and grid_mean) in the Python module
%   comparisons/seasonal_averages/quantities.py.
%
%   * window_min > 0: the mean of the finite samples in each window, labelled
%     by the window's START time. Windows are anchored at midnight (00:00) of
%     the first sample's day, so 5-min windows start at :00, :05, :10, ...
%     Every window between the first and the last sample is returned; windows
%     with fewer than min_count finite samples are NaN (not 0). This is what
%     pandas' Series.resample(...).mean() does in the Python version.
%   * window_min == 0: no averaging. The finite samples are returned as they
%     are, at their own time stamps.
%
%   A sample at time t belongs to the window [w, w + window) that holds t.
%
% INPUTS:
%   t          - time stamps (datetime column, UTC), any order
%   x          - values (numeric column, same length as t; NaN = missing)
%   window_min - window length in minutes, any non-negative number
%                (fractional allowed, e.g. 7.5); 0 = no averaging
%   min_count  - (optional, default 1) fewest finite samples a window needs
%
% OUTPUTS:
%   t_win - window start times (datetime column); for window_min == 0, the
%           time stamps of the finite samples
%   x_win - window means (NaN where a window has < min_count finite samples)
%   n_win - number of finite samples in each window
%
% REFERENCES:
%   Python original: VS_CODE/Python-Research/EPCAPE/comparisons/
%     seasonal_averages/quantities.py (grid_mean, window_mean), and its
%     test test_window_mean_zero_is_every_sample_and_lengths_are_arbitrary.
%
% By Andrew John Buggee
%%

function [t_win, x_win, n_win] = window_mean(t, x, window_min, min_count)


% -------------------------------------------------------------------
% ------------------------- Check the inputs ------------------------
% -------------------------------------------------------------------

if nargin < 4 || isempty(min_count)
    min_count = 1;
end
if ~isscalar(window_min) || ~isfinite(window_min) || window_min < 0
    error('window_mean: window_min must be a finite number >= 0 (minutes), got %g', window_min);
end
t = t(:);
x = double(x(:));
if numel(t) ~= numel(x)
    error('window_mean: t and x must have the same number of elements');
end



% -------------------------------------------------------------------
% ------------------- window_min = 0: no averaging ------------------
% -------------------------------------------------------------------

if window_min == 0
    ok = isfinite(x);
    t_win = t(ok);
    x_win = x(ok);
    n_win = ones(nnz(ok), 1);
    return
end



% -------------------------------------------------------------------
% ---------- window_min > 0: mean of each fixed window --------------
% -------------------------------------------------------------------

% Seconds since midnight of the first sample's day. Every window index is
% then floor(seconds / window length), so the windows start at midnight.
origin = dateshift(min(t), 'start', 'day');
sec = seconds(t - origin);
w_s = 60 * window_min;
k = floor(sec / w_s + 1e-9) + 1;          % window index of each sample (1 = first window); the
                                          % 1e-9 guards against float round-off at a boundary
n_windows = max(k);

ok = isfinite(x);
sum_win = accumarray(k(ok), x(ok), [n_windows 1], @sum, 0);
n_win   = accumarray(k(ok), 1,     [n_windows 1], @sum, 0);

% Return the windows from the one holding the first sample to the one
% holding the last (pandas' resample spans the same range)
first = min(k);
sum_win = sum_win(first:end);
n_win   = n_win(first:end);

x_win = sum_win ./ n_win;
x_win(n_win < min_count) = NaN;          % too few samples: NaN, not 0
t_win = origin + seconds(((first:n_windows)' - 1) * w_s);


end
