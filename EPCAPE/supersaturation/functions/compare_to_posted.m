function out = compare_to_posted(ours, posted, tolerance_pct, label)
% COMPARE_TO_POSTED  How closely a re-derived retrieval matches Jeramy's posted values (same hours).
%
%   out = compare_to_posted(ours, posted, tolerance_pct, label)
%
%   Inputs
%     ours, posted   column vectors, %, one element per posted hour (NaN where ours is missing)
%     tolerance_pct  "agrees" if |ours - posted| < tolerance_pct (percentage points)
%     label          text for the printed line
%   Output struct: n_posted, n_valid, r (Pearson correlation), median_diff_pct, frac_within
%   (also printed)

ok = isfinite(ours) & isfinite(posted);
d = ours(ok) - posted(ok);
R = corrcoef(ours(ok), posted(ok));
out.n_posted = numel(posted);
out.n_valid = nnz(ok);
out.r = R(1, 2);
out.median_diff_pct = median(d);
out.frac_within = mean(abs(d) < tolerance_pct);
fprintf('%s: re-derived %d of %d posted hours, r = %.3f, median difference %+.4f %%, within ±%g %%: %.0f %%\n', ...
        label, out.n_valid, out.n_posted, out.r, out.median_diff_pct, tolerance_pct, 100 * out.frac_within);

end
