function smax_pct = smax_posted_arg2000_pct(modes, w_ms, T_K, p_Pa)
% SMAX_POSTED_ARG2000_PCT  Jeramy Dedrick's posted ParcelActivation_ARG2000.m, ported line by line.
%
%   THIS IS JERAMY L. DEDRICK'S GITHUB CODE, translated line by line: his ParcelActivation_ARG2000.m
%   (https://github.com/jdedrick95/ParcelActivation_ARG2000, 2025).
%   Please cite his code and Dedrick et al. (2026, J. Atmos. Sci., doi:10.1175/JAS-D-26-0004.1) if you use
%   it. His repository has no licence: ask him before sharing this function outside the group.
%
%   smax_pct = smax_posted_arg2000_pct(modes, w_ms, T_K, p_Pa)
%
%   Same inputs as smax_arg2000_pct (modes = [N m-3, D_g m, sigma_g, kappa],
%   rows = coarse, accumulation, Aitken, as in his function). This reproduces
%   the posted code (github.com/jdedrick95/ParcelActivation_ARG2000, dated
%   14 Jul 2025) INCLUDING the following, so the two can be compared:
%     * certain: the Aitken mode is given the accumulation-mode diameter
%       (d_aero_mode2 = mode2_param(2) inside the Aitken block);
%     * likely: S_m and zeta use the diameter-form Kelvin parameter (4 sigma)
%       with the mode RADIUS (ARG2000 uses the radius form, 2 sigma, with the radius);
%     * likely: alpha and gamma use molar latent heat and heat capacity
%       (J mol-1, J mol-1 K-1) in formulas written per kg;
%     * likely: the square roots are taken per mode and then summed
%       (ARG2000 eq. 6 takes the square root of the sum).
%   Status of the "likely" items after checking (19 Oct 2026): CONFIRMED.
%     * alpha and gamma: Rogers & Yau (1989, p. 106, eqs. 7.23-7.24) write them with L
%       in J kg-1 and c_p in J kg-1 K-1; the textbook form reproduces their Fig. 7.2
%       (Q1 = 7.6e-4 to 5.0e-4 m-1 for -30 to 20 C), while this code gives alpha =
%       2.8e-4 m-1 at 285 K, below the plotted range.
%     * the Kelvin coefficient: ARG2000's A is the radius form, Rogers & Yau's
%       a = 2 sigma/(rho_L R_v T) (p. 88), used with the mode radius; so is pyrcel's.
%     * the sum: pyrcel takes smax = 1/sqrt(sum over modes), as ARG2000 eq. 6.
%   It may not be the exact version used for the paper. Do not use it for new work.

lv  = 4.5e4;              % latent heat, J mol-1 (as posted)
Ma  = 28.97e-3;           % molar mass of air, kg mol-1
Mw  = 18.02e-3;           % molar mass of water, kg mol-1
R   = 8.3144621;          % J mol-1 K-1
cp  = 1005.7 * Ma;        % J mol-1 K-1 (as posted: cp = cp .* Ma)
Dv  = 2.55e-5;            % m2 s-1
g   = 9.81;
kv  = 2.40e-2;
rho_l = 1e3;
sigma_sa = 0.076;

es = 6.1094 * exp((17.625 * (T_K - 273.15)) / ((T_K - 273.15) + 243.04)) * 100;   % Pa
G  = ((rho_l * R * T_K) / (Mw * Dv * es)) + (((rho_l * lv) / (Mw * kv * T_K)) * ((lv / (R * T_K)) - 1));
A_kohler = (4 * sigma_sa * Mw) / (rho_l * R * T_K);        % diameter form...

d1 = modes(1, 2); d2 = modes(2, 2);
d3 = d2;                                                     % ...Aitken diameter copied from mode 2 (as posted)
r = [d1, d2, d3] / 2;                                        % radii
kap = modes(:, 4)';
Sm = sqrt((4 .* A_kohler^3) ./ (27 .* kap .* r.^3));        % ...with the radius

alpha = ((g * Mw * lv) / (cp * R * T_K^2)) - (g * Ma) / (R * T_K);
gamma = ((R * T_K) / (es * Mw)) + ((Mw * lv^2) / (cp * p_Pa * Ma * T_K));
f   = 0.5 .* exp(2.5 .* (log(modes(:, 3)').^2));
gi  = 1 + (0.25 .* log(modes(:, 3)'));
eta = ((alpha .* w_ms) ./ (1 ./ G)).^(3/2) ./ (2 .* pi .* rho_l .* gamma .* modes(:, 1)');
zeta = (2/3) .* sqrt((alpha .* w_ms) ./ (1 ./ G)) .* A_kohler;

smax_pct = 100 / sum(((1 ./ (Sm.^2)) .* ((f .* ((zeta ./ eta).^(3/2))) + ...
                     (gi .* ((Sm.^2) ./ (eta + 3 .* zeta)).^(3/4)))).^(1/2), 'omitnan');   % sum of roots (as posted)

end
