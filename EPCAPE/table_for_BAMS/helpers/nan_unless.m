function factor = nan_unless(sel)
% NAN_UNLESS  1 where sel is true, NaN where it is false.
%
%   factor = nan_unless(sel)
%
%   Multiplying data by this factor keeps the selected values and turns all
%   others into NaN, which season_stats.m then ignores. Example:
%       values .* nan_unless(class == 1)     % keep only the cloud intervals

factor = ones(size(sel));
factor(~sel) = NaN;

end
