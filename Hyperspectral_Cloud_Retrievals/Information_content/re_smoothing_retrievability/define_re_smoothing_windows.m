function windows_m = define_re_smoothing_windows()
% Vertical moving-average window widths (full width, in METERS) applied to the
% in-situ effective-radius profile to test the retrievability of vertical
% droplet-size structure.
%
% These define how aggressively the high-frequency vertical variability in
% r_e(z) is removed. The "raw" (unsmoothed) profile is always computed in
% ADDITION to these windows, so the number of versions per profile is
% 1 + numel(windows_m).
%
% Context for choosing values: the in-situ profiles are sampled at dz ~ 5 m,
% and a typical marine boundary-layer cloud in these data is ~200-400 m deep
% (median ~40 levels). So, roughly:
%   25 m  -> removes only the finest-scale fluctuations
%   75 m  -> mild smoothing
%   125 m -> moderate smoothing
%   250 m -> heavy smoothing (approaches a near-homogeneous r_e profile)
%
% *** SINGLE SOURCE OF TRUTH ***
% Change the values here and they propagate to the preview script, the
% reflectance calculation, and the aggregate figure. Re-run preview_re_smoothing
% to eyeball the effect before submitting the supercomputer jobs.
%
% By Andrew John Buggee

windows_m = [25, 50, 75, 100, 150];   % meters - full vertical width of moving-average window

end
