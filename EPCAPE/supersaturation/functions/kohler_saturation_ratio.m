function S = kohler_saturation_ratio(d_wet_m, d_dry_m, kappa)
% KOHLER_SATURATION_RATIO  Equilibrium saturation ratio over a solution droplet (kappa-Kohler theory).
%
%   S = kohler_saturation_ratio(d_wet_m, d_dry_m, kappa)
%
%   Petters & Kreidenweis (2007, Atmos. Chem. Phys. 7, 1961-1971, p. 1962, eq. 6):
%
%       S(D) = (D^3 - Dd^3) / (D^3 - Dd^3 (1 - kappa)) * exp(A / D)
%
%   * the first factor is the solute (Raoult) term: it is < 1 and LOWERS the
%     saturation ratio a droplet needs, because dissolved material lowers the
%     vapour pressure over the solution;
%   * the second factor is the curvature (Kelvin) term: it is > 1 and RAISES it,
%     because molecules escape a strongly curved surface more easily.
%   Their product has a maximum, the critical saturation ratio S_c: air with a
%   supersaturation above S_c - 1 lets the particle grow without limit
%   ("activate") into a cloud droplet.
%
%   Textbook background (Rogers & Yau 1989): the Kelvin equation e_s(r) =
%   e_s(inf) exp(2 sigma / (r rho_L R_v T)) (p. 84, eq. 6.1); Raoult's law for the
%   solute (pp. 86-87, eq. 6.4); their combination, the Kohler curve
%   S = 1 + a/r - b/r^3 (p. 88, eq. 6.6); and activation at the curve's peak, below
%   which the droplet is a stable "haze particle" (p. 89).
%
%   Inputs
%     d_wet_m  wet (droplet) diameter, m (scalar or array; must be > d_dry_m)
%     d_dry_m  dry particle diameter, m (scalar)
%     kappa    hygroscopicity parameter, dimensionless (scalar)
%   Output
%     S        saturation ratio e / e_s (dimensionless; supersaturation % = 100 (S - 1))

c = kohler_constants();
d3  = d_wet_m .^ 3;
dd3 = d_dry_m .^ 3;
S = (d3 - dd3) ./ (d3 - dd3 .* (1 - kappa)) .* exp(c.A_m ./ d_wet_m);

end
