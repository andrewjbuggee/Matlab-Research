function s_pct = retrieve_supersaturation_pct(n_target_cm3, ccn6_cm3)
% RETRIEVE_SUPERSATURATION_PCT  Supersaturation (%) where the CCN spectrum reaches a target number.
%
%   s_pct = retrieve_supersaturation_pct(n_target_cm3, ccn6_cm3)
%
%   The core of the CCN-based retrievals of Dedrick et al. (2026, pp. 1065-1066).
%   Physical basis: in a rising parcel the supersaturation peaks (S_max) near
%   cloud base, every nucleus whose critical supersaturation is below S_max
%   activates, and the droplet concentration is then fixed (Rogers & Yau 1989,
%   pp. 107-110), so N_d = CCN(S_max). Observed droplet spectra confirm the close
%   link between the CCN spectrum and the droplets (Rogers & Yau 1989, p. 96).
%     * cloud-scene (s_Nd):  n_target = droplet number N_d. Because the droplets
%       are the particles that activated, N_d = CCN(S_max), so the S at which the
%       measured CCN spectrum equals N_d is the cloud's effective supersaturation.
%     * Hoppel minimum (s_CCN): n_target = number of particles larger than the
%       Hoppel minimum diameter (the particles that activated in recent clouds).
%
%   The spectrum is extended to (0 %, 0 cm-3) and to 1.2 % (extend_ccn_spectrum)
%   and treated as straight lines between its 8 points. The function walks the 7
%   segments from 0 % upward and returns the first crossing. NaN if an input is
%   missing or the spectrum never reaches n_target (that is, s > 1.2 %).
%
%   Below 0.1 % the spectrum is the line from (0, 0) to (0.1 %, CCN(0.1 %)), so
%   there s = 0.1 % * n_target / CCN(0.1 %): the answer depends on the CCN count
%   at 0.1 % alone.
%
%   Inputs
%     n_target_cm3  scalar, cm-3
%     ccn6_cm3      CCN at the six setpoints, cm-3
%   Output
%     s_pct         supersaturation, % (NaN if not retrievable)

s_pct = NaN;
if ~isfinite(n_target_cm3) || any(~isfinite(ccn6_cm3))
    return
end
[s_ext, c_ext] = extend_ccn_spectrum(ccn6_cm3);
for k = 1:numel(s_ext) - 1
    c0 = c_ext(k);
    c1 = c_ext(k + 1);
    if n_target_cm3 >= min(c0, c1) && n_target_cm3 <= max(c0, c1)   % the target lies on this segment
        if c1 == c0
            s_pct = s_ext(k);
        else
            s_pct = s_ext(k) + (n_target_cm3 - c0) * (s_ext(k + 1) - s_ext(k)) / (c1 - c0);
        end
        return
    end
end

end
