function s_pct = retrieve_at_hours_pct(t, n_target_cm3, t_hour, ccn_hour_cm3)
% RETRIEVE_AT_HOURS_PCT  CCN-based supersaturation retrieval (%) at a list of hours, from one set of hourly spectra.
%
%   s_pct = retrieve_at_hours_pct(t, n_target_cm3, t_hour, ccn_hour_cm3)
%
%   CREDIT: method of Jeramy L. Dedrick (Dedrick et al. 2026, J. Atmos. Sci., Sec. 2d, pp. 1065-1066).
%   A loop over retrieve_supersaturation_pct, so the same hours can be retrieved
%   with two versions of the hourly CCN spectra (e.g. with and without Sourita
%   Saha's Gaussian correction of the 0.1 % setpoint) and differenced.
%
%   Inputs
%     t             hours (datetime column, start of each hour)
%     n_target_cm3  target number at those hours, cm-3: N_d for s_Nd, N(D > D_HM) for s_CCN
%     t_hour, ccn_hour_cm3   hourly spectra (output of hourly_ccn_spectrum)
%   Output
%     s_pct         retrieved supersaturation, % (NaN where the hour has no complete
%                   spectrum or the spectrum never reaches the target)

s_pct = NaN(numel(t), 1);
for i = 1:numel(t)
    spec_cm3 = ccn_spectrum_at_hour(t(i), t_hour, ccn_hour_cm3);   % 1 x 6, cm-3
    s_pct(i) = retrieve_supersaturation_pct(n_target_cm3(i), spec_cm3);
end

end
