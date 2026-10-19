function s_c_pct = critical_supersaturation_pct(d_dry_m, kappa)
% CRITICAL_SUPERSATURATION_PCT  Critical supersaturation (%) of a dry particle, kappa-Kohler approximation.
%
%   s_c_pct = critical_supersaturation_pct(d_dry_m, kappa)
%
%   For kappa > 0.2 and s_c < 1 % the peak of the kappa-Kohler curve is, to a good
%   approximation (Petters & Kreidenweis 2007, p. 1964, eq. 10, written for kappa),
%
%       ln S_c = sqrt( 4 A^3 / (27 kappa Dd^3) )
%
%   so bigger (Dd) and more soluble (kappa) particles activate at LOWER
%   supersaturation: s_c ~ Dd^(-3/2) kappa^(-1/2).
%
%   This is the textbook critical saturation ratio S* = 1 + sqrt(4 a^3 / (27 b))
%   (Rogers & Yau 1989, p. 89, eq. 6.8) written in diameters: with D = 2r, a = A/2
%   and b = kappa r_d^3 (the dilute solute term of the kappa-Kohler equation),
%   4 a^3 / (27 b) = 4 A^3 / (27 kappa Dd^3) (algebra only). Larger nuclei have
%   lower critical supersaturations (Rogers & Yau 1989, p. 89, Table 6.2).
%
%   This follows Jeramy Dedrick's kappakohler.m, option 1, which returns
%   100 * ln(1 + sqrt(...)). For s_c << 1 that equals 100 * sqrt(...) and
%   100 * (exp(sqrt(...)) - 1) to better than 0.1 % of the value.
%
%   Inputs
%     d_dry_m  dry diameter, m (array)
%     kappa    hygroscopicity (scalar or same size as d_dry_m)
%   Output
%     s_c_pct  critical supersaturation, % (same size as d_dry_m)

c = kohler_constants();
s_c_pct = 100 .* log(1 + sqrt(4 .* c.A_m .^ 3 ./ (27 .* kappa .* d_dry_m .^ 3)));

end
