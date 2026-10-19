function d_crit_m = critical_diameter_m(s_pct, kappa)
% CRITICAL_DIAMETER_M  Smallest dry diameter (m) that activates at supersaturation s (kappa-Kohler).
%
%   BASED ON JERAMY L. DEDRICK'S GITHUB CODE: option 2 of his kappakohler.m
%   (https://github.com/jdedrick95/kappakohler, 2024).
%   Please cite his code and Dedrick et al. (2026, J. Atmos. Sci., doi:10.1175/JAS-D-26-0004.1) if you use
%   it. His repository has no licence: ask him before sharing this function outside the group.
%
%   d_crit_m = critical_diameter_m(s_pct, kappa)
%
%   The inverse of critical_supersaturation_pct (kappakohler.m, option 2; Petters &
%   Kreidenweis 2007, p. 1964, eq. 10 solved for the dry diameter):
%
%       Dc = ( 4 A^3 / (27 kappa (ln S)^2) )^(1/3),   S = 1 + s/100
%
%   Inputs
%     s_pct  supersaturation, % (array)
%     kappa  hygroscopicity (scalar or same size)
%   Output
%     d_crit_m  critical dry diameter, m

c = kohler_constants();
d_crit_m = (4 .* c.A_m .^ 3 ./ (27 .* kappa .* log(1 + s_pct ./ 100) .^ 2)) .^ (1 / 3);

end
