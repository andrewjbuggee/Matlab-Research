function spec_cm3 = ccn_spectrum_at_hour(t, t_hour, ccn_hour_cm3)
% CCN_SPECTRUM_AT_HOUR  The hourly Pier CCN spectrum at one hour, or NaNs if that hour has none.
%
%   spec_cm3 = ccn_spectrum_at_hour(t, t_hour, ccn_hour_cm3)
%
%   Inputs
%     t             the hour wanted (datetime, start of the hour)
%     t_hour        start of each hour with a complete spectrum (output of hourly_ccn_spectrum)
%     ccn_hour_cm3  the hourly spectra, cm-3 (M x 6; setpoints 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 %)
%   Output
%     spec_cm3      1 x 6, cm-3 (NaN(1, 6) if hour t has no complete spectrum)

k = find(t_hour == t, 1);
if isempty(k)
    spec_cm3 = NaN(1, 6);
else
    spec_cm3 = ccn_hour_cm3(k, :);
end

end
