%% Hourly ratio of two time series, <num> / <den>, over the samples where both are valid
%
% DESCRIPTION:
%   MATLAB version of hourly_ratio in the Python module
%   comparisons/seasonal_averages/kavin.py. It gives one SW transmittance per
%   hour, the resolution of Dan Lubin's hourly tables, from the 1-min RADFLUX
%   irradiances:
%       T(h) = mean(num over hour h) / mean(den over hour h)
%   where both means use the same samples: num and den finite and den > 0.
%   Hour bins are [h, h + 1 h) on the time stamps given, so pass the START of
%   each sample (read_radflux_sw.m: rad.time). Hours with fewer than
%   min_count usable samples are NaN; with 1-min data, min_count = 30 (the
%   Python default) asks for at least half the hour.
%
%   The hours returned are every hour from the one holding the first sample
%   to the one holding the last (pandas' resample("1h") spans the same).
%
% INPUTS:
%   t         - sample START times (datetime column, UTC)
%   num       - numerator values, e.g. measured SW down (W m-2)
%   den       - denominator values, e.g. clear-sky SW down (W m-2)
%   min_count - (optional, default 30) fewest usable samples per hour
%
% OUTPUTS:
%   ratio  - hourly ratio (NaN where fewer than min_count usable samples)
%   t_hour - hour start times (datetime column, UTC)
%   n      - usable samples in each hour
%
% REFERENCES:
%   Python original: VS_CODE/Python-Research/EPCAPE/comparisons/
%     seasonal_averages/kavin.py (hourly_ratio).
%
% By Andrew John Buggee
%%

function [ratio, t_hour, n] = hourly_ratio(t, num, den, min_count)


if nargin < 4 || isempty(min_count)
    min_count = 30;
end

% usable samples: both values finite and a positive denominator
ok = isfinite(num) & isfinite(den) & (den > 0);
num(~ok) = NaN;
den(~ok) = NaN;

% hourly means of the same samples (window_mean anchors the hours at
% midnight and labels each by its start)
[t_hour, num_avg, n] = window_mean(t, num, 60);
[~,      den_avg]    = window_mean(t, den, 60);

ratio = num_avg ./ den_avg;
ratio(n < min_count) = NaN;


end
