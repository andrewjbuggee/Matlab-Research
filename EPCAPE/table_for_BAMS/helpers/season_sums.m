function stats = season_sums(values, masks)
% SEASON_SUMS  Sum of the values in each season (used for sampling HOURS).
%
%   stats = season_sums(values, masks)
%
%   Inputs
%     values  vector, e.g. the hours each interval represents (0 or 5/60 h)
%     masks   logical matrix (samples x 5) from sheet_season_masks.m or
%             calendar_quarter_masks.m
%   Output (struct, every field 1 x 5: EPCAPE, Spring, Summer, Fall, Winter)
%     stats.mean  the SUM over the season (called "mean" only so that every
%                 result has the same fields; format_value prints it alone)
%     stats.std   NaN (a total has no spread)
%     stats.n     number of non-zero values that were added up

values = double(values(:));
values(~isfinite(values)) = 0;     % a missing value adds nothing

stats.mean = zeros(1, 5);
stats.std  = nan(1, 5);
stats.n    = zeros(1, 5);
for k = 1:5
    stats.mean(k) = sum(values(masks(:, k)));
    stats.n(k)    = sum(values(masks(:, k)) ~= 0);
end

end
