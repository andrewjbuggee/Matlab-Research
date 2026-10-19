function mean_s_pct = mean_snd_with_scaled_ccn01(t, nd_cm3, factor, t_hour, ccn_hour_cm3)
% MEAN_SND_WITH_SCALED_CCN01  Mean cloud-scene supersaturation (%) if the CCN count at 0.1 % were scaled.
%
%   mean_s_pct = mean_snd_with_scaled_ccn01(t, nd_cm3, factor, t_hour, ccn_hour_cm3)
%
%   CREDIT: method of Jeramy L. Dedrick (Dedrick et al. 2026, J. Atmos. Sci., Sec. 2d, pp. 1065-1066). His code
%   for this step is not public; this is a reconstruction from the paper, checked against his posted
%   retrievals (UC San Diego Library, doi:10.6075/J0WM1DSC). It may differ from his code in details.
%
%   Re-retrieves s_Nd at each hour t(i) with the CCN count at the 0.1 % setpoint
%   multiplied by `factor` (the other setpoints unchanged), and returns the mean.
%   Used to show what a correction of the 0.1 % CCN value (Sourita Saha's) would
%   do. Below 0.1 % the retrieval is s = 0.1 % x N_d / CCN(0.1 %), so those
%   retrievals scale exactly as 1 / factor (Dedrick et al. 2026, pp. 1065-1066).
%
%   Inputs
%     t             hours (datetime column)
%     nd_cm3        droplet number at those hours, cm-3
%     factor        multiplier for the 0.1 % CCN count (1 = as posted)
%     t_hour, ccn_hour_cm3   hourly spectra (hourly_ccn_spectrum)
%   Output
%     mean_s_pct    mean retrieved supersaturation, % (hours that give no retrieval are skipped)

s_pct = NaN(numel(t), 1);
for i = 1:numel(t)
    spec = ccn_spectrum_at_hour(t(i), t_hour, ccn_hour_cm3);
    spec(1) = spec(1) * factor;
    s_pct(i) = retrieve_supersaturation_pct(nd_cm3(i), spec);
end
mean_s_pct = mean(s_pct, 'omitnan');

end
