function n = count_sorted(sorted_values, x, rule)
% COUNT_SORTED  How many of sorted_values are below (or at or below) each x.
%
%   n = count_sorted(sorted_values, x, 'below')        n(i) = number of sorted_values <  x(i)
%   n = count_sorted(sorted_values, x, 'at_or_below')  n(i) = number of sorted_values <= x(i)
%
%   Why this is useful: for a SORTED list of sample times, the samples
%   inside a window [lo, hi] are numbers  count(<lo) + 1  to  count(<=hi).
%   With that, a mean over millions of windows can be computed without a
%   loop (see kavin_window_mean in calc_rows31_32_FM120_hours.m and
%   near_rain.m).
%
%   Inputs
%     sorted_values  numeric vector sorted in increasing order (e.g. times in
%                    seconds); repeated values are allowed
%     x              numeric vector of query points (any order)
%     rule           'below' or 'at_or_below'
%   Output
%     n              vector the size of x
%
%   How: the unique query points split the number line into bins;
%   discretize puts every sorted value into its bin, and a running total of
%   the bin counts gives the count below each query point.

sorted_values = sorted_values(:);
[xu, ~, back] = unique(x(:));                   % xu is sorted and has no repeats
edges = [-Inf; xu; Inf];

switch rule
    case 'below'
        % bins [edge_k, edge_k+1): a value equal to xu(k) lands in the bin AFTER xu(k)
        bin = discretize(sorted_values, edges, 'IncludedEdge', 'left');
    case 'at_or_below'
        % bins (edge_k, edge_k+1]: a value equal to xu(k) lands in the bin BEFORE xu(k)
        bin = discretize(sorted_values, edges, 'IncludedEdge', 'right');
    otherwise
        error('rule must be ''below'' or ''at_or_below''');
end
bin = bin(~isnan(bin));
counts = accumarray(bin, 1, [numel(xu) + 1, 1]);  % how many values in each bin
running = cumsum(counts);                          % running(k) = values in bins 1..k
n = running(back);                                 % count below / at-or-below each x
n = reshape(n, size(x));

end
