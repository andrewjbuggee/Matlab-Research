function [t_hour, ccn_hour] = hourly_ccn_spectrum(t30, ccn30_cm3, max_fill_h)
% HOURLY_CCN_SPECTRUM  Hourly CCN spectra (all six setpoints) from the Pier 30-min CCN file.
%
%   [t_hour, ccn_hour] = hourly_ccn_spectrum(t30, ccn30_cm3, max_fill_h)
%
%   CREDIT: method of Jeramy L. Dedrick (Dedrick et al. 2026, J. Atmos. Sci., Sec. 2d, pp. 1065-1066). His code
%   for this step is not public; this is a reconstruction from the paper, checked against his posted
%   retrievals (UC San Diego Library, doi:10.6075/J0WM1DSC). It may differ from his code in details.
%
%   The Pier CCN-200 steps through its supersaturation setpoints (0.1, 0.2,
%   0.4, 0.6, 0.8, 1.0 %), about 10 min each, so a 30-min average holds only
%   1-3 of the 6 setpoints. An hourly spectrum therefore needs a rule for the
%   missing setpoints. Dedrick et al. (2026) do not state theirs; this one
%   reproduces their posted retrievals best (checked in the live script):
%     1. average each setpoint over the hour [t, t + 1 h), ignoring NaN;
%     2. where a setpoint has no value in an hour, interpolate it linearly in
%        time between the neighbouring hours that have one, filling at most
%        max_fill_h hours of each gap (the first ones; gaps at the ends are not
%        filled). This copies pandas' interpolate(method="time",
%        limit=max_fill_h, limit_area="inside") used in the Python notebook;
%     3. keep only hours with all six setpoints.
%
%   Physics and sources: each setpoint gives one point of the aerosol's activity
%   spectrum, the number of particles activated as a function of supersaturation,
%   measured by stepping a chamber through controlled supersaturations (Rogers &
%   Yau 1989, p. 95). The Pier DMT CCN-200 ran 0.1-1.0 % in a 1-h cycle (Dedrick et
%   al. 2026, p. 1063); the posted counts are already multiplied by 1.57 for an inlet
%   leak (Dedrick et al. 2026, p. 1063; supplement Text S4, pp. 4-5). Uncertainty:
%   22 +/- 5 % on average, up to 60 % at 0.1-0.2 % (supplement Text S3, p. 4).
%
%   Inputs
%     t30         30-min time stamps (datetime column)
%     ccn30_cm3   CCN concentrations, cm-3 (N x 6; columns = setpoints)
%     max_fill_h  hours of each gap that may be filled (Dedrick reproduction: 2)
%   Outputs
%     t_hour      start of each hour with a complete spectrum (datetime column)
%     ccn_hour    the hourly spectra, cm-3 (M x 6)

% 1. hourly mean of each setpoint (retime puts the samples in [t, t + 1 h) into hour t)
TT = array2timetable(ccn30_cm3, 'RowTimes', t30);
TH = retime(TT, 'hourly', @(x) mean(x, 'omitnan'));   % mean([]) = NaN for an hour without data
ccn = TH{:, :};

% 2. fill interior gaps of each setpoint, at most max_fill_h hours per gap
for k = 1:size(ccn, 2)
    ccn(:, k) = fill_interior_gaps(ccn(:, k), max_fill_h);
end

% 3. complete spectra only
complete = all(isfinite(ccn), 2);
t_hour = TH.Properties.RowTimes(complete);
ccn_hour = ccn(complete, :);

end


function y = fill_interior_gaps(x, limit)
% Linear interpolation across runs of NaN that have a value on both sides (equally spaced samples).
% Only the first `limit` NaNs of each run are filled (pandas: limit=limit, limit_area="inside").
y = x;
good = find(~isnan(x));
for k = 1:numel(good) - 1
    i0 = good(k);
    i1 = good(k + 1);
    for j = i0 + 1 : min(i1 - 1, i0 + limit)
        y(j) = x(i0) + (x(i1) - x(i0)) * (j - i0) / (i1 - i0);
    end
end
end
