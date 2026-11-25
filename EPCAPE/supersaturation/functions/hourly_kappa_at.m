function [kappa_h, filled] = hourly_kappa_at(t_hour, kh_time, kh_kappa, max_fill_h)
% HOURLY_KAPPA_AT  Hourly-mean kappa at a list of hours, with short gaps filled from the nearest hour.
%
%   [kappa_h, filled] = hourly_kappa_at(t_hour, kh_time, kh_kappa, max_fill_h)
%
%   Used for s_kappa (Part 7 of the live script): the Hoppel-minimum diameter of
%   each hour is the mean of the hour's 5-min size-distribution fits over
%   [t, t + 1 h), so kappa is averaged over the same window. The hourly means
%   (kh_time, kh_kappa) come from the Pier series of bb2743661p_37_1_kappa.nc,
%   one sample about every 30 min (Part 3: retime(..., 'hourly', mean)).
%
%   An hour with no kappa sample takes the hourly mean of the nearest hour that
%   has one, if that hour is at most max_fill_h hours away; if the nearest hours
%   before and after are equally far, their mean is used. Otherwise kappa is NaN
%   and so is s_kappa. This fill rule is my choice (Dedrick et al. 2026 do not
%   state theirs). Against the posted s_kappa it reproduces Jeramy's values as well
%   as linear interpolation in time and better than interpolating the samples to
%   the mid-hour or taking the nearest sample (my test of 24 Nov 2026), and every
%   posted s_kappa hour without a sample has one in the hour before or after.
%
%   Inputs
%     t_hour      hours wanted (datetime column, start of each hour)
%     kh_time     start times of the hourly means (datetime column, regular 1-h grid)
%     kh_kappa    hourly-mean kappa (unitless), NaN for hours without a sample
%     max_fill_h  largest gap (hours) that may be filled; 0 = no filling
%   Outputs
%     kappa_h     kappa at t_hour (unitless; NaN where none within max_fill_h)
%     filled      true where kappa_h was taken from a neighbouring hour

assert(all(diff(kh_time) == hours(1)), 'kh_time must be a regular hourly grid');
kappa_h = NaN(numel(t_hour), 1);
filled = false(numel(t_hour), 1);

% index of each wanted hour on the hourly grid (NaN if outside it)
idx = round(hours(t_hour - kh_time(1))) + 1;
inside = idx >= 1 & idx <= numel(kh_time);

for i = find(inside)'
    k = idx(i);
    if isfinite(kh_kappa(k))
        kappa_h(i) = kh_kappa(k);                        % measured in this hour
        continue
    end
    for g = 1:max_fill_h                                 % look 1, 2, ... hours away
        before = NaN; after = NaN;
        if k - g >= 1,               before = kh_kappa(k - g); end
        if k + g <= numel(kh_kappa), after  = kh_kappa(k + g); end
        if isfinite(before) || isfinite(after)
            kappa_h(i) = mean([before after], 'omitnan');   % both sides equally near -> their mean
            filled(i) = true;
            break
        end
    end
end

end
