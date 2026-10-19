function [s_ext_pct, ccn_ext_cm3] = extend_ccn_spectrum(ccn6_cm3)
% EXTEND_CCN_SPECTRUM  The 8-point CCN spectrum of Dedrick et al. (2026, Sec. 2d).
%
%   [s_ext_pct, ccn_ext_cm3] = extend_ccn_spectrum(ccn6_cm3)
%
%   CREDIT: method of Jeramy L. Dedrick (Dedrick et al. 2026, J. Atmos. Sci., Sec. 2d, pp. 1065-1066). His code
%   for this step is not public; this is a reconstruction from the paper, checked against his posted
%   retrievals (UC San Diego Library, doi:10.6075/J0WM1DSC). It may differ from his code in details.
%
%   The Pier CCN counter measured at 0.1-1.0 % only. To retrieve
%   supersaturations outside that range, the paper extends the spectrum:
%     * below 0.1 %: a straight line to 0 cm-3 at 0 % (no particle activates at 0 %);
%     * above 1.0 %: the last segment (0.8 -> 1.0 %) continued linearly to 1.2 %
%       ("no limit on the CCN number concentration at 1.2 %").
%   Retrievals above 1.2 % (or below 0) are discarded (Dedrick et al. 2026,
%   pp. 1065-1066). The (0 %, 0) end point agrees with the power-law activity
%   spectrum N_c = C s^k, which is zero at s = 0 (Rogers & Yau 1989, p. 95,
%   eq. 6.20); the straight line between 0 and 0.1 % is the paper's choice, not a
%   physical law (a power law with k < 1 would be curved there).
%
%   Input
%     ccn6_cm3     CCN at the six setpoints 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 %, cm-3 (1 x 6 or 6 x 1)
%   Outputs
%     s_ext_pct    [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2] %
%     ccn_ext_cm3  the CCN at those points, cm-3 (1 x 8)

s_set_pct = [0.1 0.2 0.4 0.6 0.8 1.0];
c = ccn6_cm3(:)';
slope = (c(6) - c(5)) / (s_set_pct(6) - s_set_pct(5));   % cm-3 per % between 0.8 and 1.0 %
s_ext_pct   = [0, s_set_pct, 1.2];
ccn_ext_cm3 = [0, c, c(6) + slope * (1.2 - 1.0)];

end
