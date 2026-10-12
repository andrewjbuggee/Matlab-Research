function stats = season_stats(values, masks)
% SEASON_STATS  Mean, standard deviation and number of samples in each season.
%
%   stats = season_stats(values, masks)
%
%   Inputs
%     values  vector of samples (NaN samples are ignored)
%     masks   logical matrix (samples x 5) from sheet_season_masks.m or
%             calendar_quarter_masks.m
%   Output (a struct; every field is 1 x 5, one value per season, in the
%   order EPCAPE, Spring, Summer, Fall, Winter)
%     stats.mean  arithmetic mean of the finite samples in the season
%     stats.std   SAMPLE standard deviation (normalised by N - 1, MATLAB's
%                 default for std); NaN when the season has fewer than 2 samples
%     stats.n     number of finite samples used
%
%   Every sample has the same weight. So the campaign value is the mean of
%   all samples, NOT the mean of the four seasonal means.

values = double(values(:));
stats.mean = nan(1, 5);
stats.std  = nan(1, 5);
stats.n    = zeros(1, 5);

for k = 1:5
    v = values(masks(:, k));     % samples in season k
    v = v(isfinite(v));          % drop NaN (missing, or removed by a filter)
    stats.n(k) = numel(v);
    if numel(v) >= 1
        stats.mean(k) = mean(v);
    end
    if numel(v) >= 2
        stats.std(k) = std(v);   % N - 1 normalisation
    end
end

end
