function smax_pct = smax_arg2000_pct(modes, w_ms, T_K, p_Pa)
% SMAX_ARG2000_PCT  Maximum parcel supersaturation (%) from the Abdul-Razzak & Ghan (2000) parameterization.
%
%   BASED ON JERAMY L. DEDRICK'S GITHUB CODE: the structure and constants of his ParcelActivation_ARG2000.m
%   (https://github.com/jdedrick95/ParcelActivation_ARG2000, 2025),
%   with the corrections described in smax_posted_arg2000_pct.m.
%   Please cite his code and Dedrick et al. (2026, J. Atmos. Sci., doi:10.1175/JAS-D-26-0004.1) if you use
%   it. His repository has no licence: ask him before sharing this function outside the group.
%
%   smax_pct = smax_arg2000_pct(modes, w_ms, T_K, p_Pa)
%
%   A rising parcel cools at a rate set by the updraft w, so its supersaturation
%   rises; once droplets activate they take up vapour, and the supersaturation
%   peaks at S_max. Abdul-Razzak, Ghan & Rivera-Carpio (1998, J. Geophys. Res.
%   103, 6123-6131) and Abdul-Razzak & Ghan (2000, J. Geophys. Res. 105,
%   6837-6844) fit that balance into one formula (ARG2000 eq. 6):
%
%     S_max = [ sum_i (1/S_m,i^2) ( f_i (zeta/eta_i)^(3/2) + g_i (S_m,i^2/(eta_i + 3 zeta))^(3/4) ) ]^(-1/2)
%
%     S_m,i  critical supersaturation of mode i's median particle: (2/sqrt(kappa_i)) (A/(3 a_i))^(3/2),
%            with a_i the mode's median RADIUS and A = 2 M_w sigma / (R T rho_w) (radius form)
%     f_i = 0.5 exp(2.5 (ln sigma_i)^2),  g_i = 1 + 0.25 ln sigma_i      (mode-width terms)
%     eta_i = (alpha w / G)^(3/2) / (2 pi rho_w gamma N_i)               (updraft "supply" per particle)
%     zeta  = (2/3) A (alpha w / G)^(1/2)                                (Kelvin-effect term)
%     alpha = g M_w L / (c_p R T^2) - g M_a / (R T)                      (ARG98 eq. 11; L in J kg-1, c_p in J kg-1 K-1)
%     gamma = R T / (e_s M_w) + M_w L^2 / (c_p p M_a T)                  (ARG98 eq. 12)
%     G     = 1 / ( rho_w R T / (e_s D_v M_w) + L rho_w / (k_a T) (L M_w / (R T) - 1) )   (growth coefficient)
%
%   More particles -> smaller eta -> lower S_max (competition for vapour).
%   Faster updraft -> larger eta -> higher S_max.
%
%   The physics behind the terms (Rogers & Yau 1989, "A Short Course in Cloud Physics"):
%     * supersaturation budget of a rising parcel, dS/dt = Q1 dz/dt - Q2 dchi/dt:
%       cooling in ascent produces supersaturation, condensation removes it
%       (p. 106, eqs. 7.21-7.24). With R' = R/M_a and epsilon = M_w/M_a, their Q1 is
%       exactly alpha above, and their Q2 is gamma times the air density (algebra only).
%       Their Fig. 7.2 (p. 106) plots Q1 = 7.6e-4 to 5.0e-4 m-1 for -30 to 20 C.
%     * G = 1/(F_k + F_d), the condensational growth coefficient of r dr/dt =
%       (S - 1)/(F_k + F_d) (pp. 102-104, eqs. 7.17-7.18 and xi_1 on p. 104);
%       k_a and D_v are close to Table 7.1 (p. 103) at 10-20 C.
%     * A = 2 M_w sigma/(R T rho_w) is the radius-form Kelvin coefficient a of the
%       Kohler curve (p. 88, eq. 6.6); S_m,i is that curve's critical supersaturation
%       for the mode's median particle (p. 89, eq. 6.8, with b = kappa a_i^3).
%     * the supersaturation peaks near cloud base and fixes the droplet number
%       (pp. 107-110); Twomey's S_max = 3.6 [1.6e-3 U^(3/2) / C]^(1/(k+2)) is the
%       simpler power-law-spectrum version (pp. 95-96, eq. 6.22).
%   The equation forms (radius form A, the square root of the SUM over modes, L and
%   c_p per kg) were cross-checked against the independent pyrcel implementation
%   (github.com/darothen/pyrcel, pyrcel/activation/_arg2000.py, accessed 19 Oct 2026).
%
%   Inputs (SI units)
%     modes  K x 4 matrix, one row per mode: [N (m-3), D_g (m, median DIAMETER), sigma_g, kappa]
%     w_ms   updraft speed at cloud base, m s-1
%     T_K    cloud-base temperature, K
%     p_Pa   cloud-base pressure, Pa
%   Output
%     smax_pct  maximum supersaturation, %
%
%   Constants match Jeramy Dedrick's ParcelActivation_ARG2000.m where it states
%   them (surface tension 0.076 N m-1 for the solution, as in Ghan et al. 2011),
%   but L and c_p are per kg as the formulas require.

g_ms2    = 9.81;        % gravity, m s-2
M_w      = 0.018;       % molar mass of water, kg mol-1
M_a      = 0.02897;     % molar mass of dry air, kg mol-1
R        = 8.314;       % universal gas constant, J mol-1 K-1
rho_w    = 1000;        % density of water, kg m-3
L_v      = 2.5e6;       % latent heat of vaporization, J kg-1
c_p      = 1005.7;      % specific heat of air at constant pressure, J kg-1 K-1
D_v      = 2.55e-5;     % diffusivity of water vapour in air, m2 s-1
k_a      = 2.40e-2;     % thermal conductivity of air, J m-1 s-1 K-1
sigma_sa = 0.076;       % surface tension of the solution droplet, N m-1

e_s   = 611.094 * exp(17.625 * (T_K - 273.15) / (T_K - 273.15 + 243.04));   % Magnus formula, Pa
alpha = g_ms2 * M_w * L_v / (c_p * R * T_K^2) - g_ms2 * M_a / (R * T_K);      % m-1: cooling drives S up
gamma = R * T_K / (e_s * M_w) + M_w * L_v^2 / (c_p * p_Pa * M_a * T_K);       % m3 kg-1: condensation drives S down
G     = 1 / (rho_w * R * T_K / (e_s * D_v * M_w) + ...
             L_v * rho_w / (k_a * T_K) * (L_v * M_w / (R * T_K) - 1));        % m2 s-1
A     = 2 * M_w * sigma_sa / (R * T_K * rho_w);                               % Kelvin parameter, RADIUS form, m
zeta  = (2 / 3) * A * sqrt(alpha * w_ms / G);

total = 0;
for i = 1:size(modes, 1)
    N = modes(i, 1); a = modes(i, 2) / 2; sg = modes(i, 3); kap = modes(i, 4);
    Sm  = 2 / sqrt(kap) * (A / (3 * a))^1.5;                       % critical supersaturation of the median radius
    eta = (alpha * w_ms / G)^1.5 / (2 * pi * rho_w * gamma * N);
    f   = 0.5 * exp(2.5 * log(sg)^2);
    gg  = 1 + 0.25 * log(sg);
    total = total + (f * (zeta / eta)^1.5 + gg * (Sm^2 / (eta + 3 * zeta))^0.75) / Sm^2;
end
smax_pct = 100 / sqrt(total);   % square root of the SUM over modes (ARG2000 eq. 6)

end
