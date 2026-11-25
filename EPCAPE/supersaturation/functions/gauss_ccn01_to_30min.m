function out = gauss_ccn01_to_30min(t_s, meas_cm3, gauss_cm3, r2, r2_min, t30)
% GAUSS_CCN01_TO_30MIN  Put Sourita Saha's 1-s Gaussian-corrected 0.1 % CCN onto the posted 30-min rows.
%
%   out = gauss_ccn01_to_30min(t_s, meas_cm3, gauss_cm3, r2, r2_min, t30)
%
%   Why stages, not clock bins. The Pier CCN-200 measures the 0.1 % setpoint in
%   one ~10-min stage per cycle, and in Sourita's 1-s file the cycle repeats every
%   ~73 min (e.g. 1 Jun 2023: stages start 00:12, 01:25, 02:38, ...), not every
%   60 min (my measurement from the file; the paper gives a 1-h cycle at the Pier,
%   Dedrick et al. 2026, p. 1063, and a 1-1.5 h duty cycle in general, p. 1062).
%   Each posted 30-min row holds one stage's average, but a stage often straddles
%   two clock bins. So the samples are first grouped into stages (a new stage
%   starts after a gap > 60 s), and each 30-min row gets the stage that has the
%   most samples in [t, t + 30 min).
%
%   Why a ratio. The posted 0.1 % values come from ARM's stage-averaged product
%   (epcccn2colaavgM1), whose averaging window is not the full 600 s that the 1-s
%   file holds, so 1.57 x the stage mean of the 1-s "measured" CCN does not equal
%   the posted value exactly (Part 11 of the live script shows the spread). The
%   stage ratio gauss / measured isolates Sourita's correction from that
%   averaging difference; the live script multiplies the posted 0.1 % value by it.
%   The 1.57 leak factor cancels in the ratio (Sourita: both the measured and the
%   Gaussian-corrected CCN need the 1.57; email of 20 Oct 2026). (My choice of
%   method, not from Sourita or Jeramy.)
%
%   Samples with fit R^2 < r2_min are dropped from both sums: Sourita reported
%   rejecting fits with R^2 < 0.8 (email of 20 Oct 2026), but the file still
%   holds 1,720 such 1-s samples.
%
%   Inputs
%     t_s        1-s sample times (datetime column, UTC)
%     meas_cm3   measured CCN at 0.1 %, cm-3 (no 1.57 factor)
%     gauss_cm3  Gaussian-corrected CCN at 0.1 %, cm-3 (no 1.57 factor)
%     r2         R^2 of each sample's Gaussian fit
%     r2_min     samples with r2 < r2_min are not used
%     t30        start times of the posted 30-min rows (datetime column, regular 30-min grid)
%   Output struct (each field numel(t30) x 1; NaN where no stage falls in the row)
%     ratio            stage sum(gauss) / sum(measured) = ratio of the stage means
%     meas_stage_cm3   stage mean of the measured CCN, cm-3 (no 1.57 factor)
%     gauss_stage_cm3  stage mean of the Gaussian-corrected CCN, cm-3 (no 1.57 factor)
%     n_samples        number of 1-s samples in that stage (after the R^2 cut)
%     n_r2_dropped     scalar: number of 1-s samples dropped by the R^2 cut

% 0. keep valid samples, in time order
keep = isfinite(meas_cm3) & isfinite(gauss_cm3) & r2 >= r2_min;
out.n_r2_dropped = nnz(r2 < r2_min);
[t_s, order] = sort(t_s(keep));
meas_cm3 = meas_cm3(keep);  meas_cm3 = meas_cm3(order);
gauss_cm3 = gauss_cm3(keep); gauss_cm3 = gauss_cm3(order);

% 1. stages: a new stage starts after a gap of more than 60 s between samples
stage = cumsum([true; diff(t_s) > seconds(60)]);   % stage number of each sample
n_stage = stage(end);
sum_meas  = accumarray(stage, meas_cm3, [n_stage 1]);    % cm-3 s
sum_gauss = accumarray(stage, gauss_cm3, [n_stage 1]);   % cm-3 s
n_in_stage = accumarray(stage, 1, [n_stage 1]);

% 2. which 30-min row each sample falls in: [t30(k), t30(k) + 30 min)
assert(all(diff(t30) == minutes(30)), 't30 must be a regular 30-min grid');
row = discretize(t_s, [t30; t30(end) + minutes(30)]);   % NaN outside the grid
ok = ~isnan(row);

% 3. for each row, the stage with the most samples in it
counts = sparse(row(ok), stage(ok), 1, numel(t30), n_stage);   % rows x stages
[n_max, best] = max(counts, [], 2);
has = full(n_max) > 0;
best = full(best(has));

out.ratio = NaN(numel(t30), 1);
out.meas_stage_cm3 = NaN(numel(t30), 1);
out.gauss_stage_cm3 = NaN(numel(t30), 1);
out.n_samples = zeros(numel(t30), 1);
out.ratio(has) = sum_gauss(best) ./ sum_meas(best);
out.meas_stage_cm3(has) = sum_meas(best) ./ n_in_stage(best);
out.gauss_stage_cm3(has) = sum_gauss(best) ./ n_in_stage(best);
out.n_samples(has) = n_in_stage(best);

end
