function [E_down, E_up, prof] = slab_emission(T_mean_K, eps_slab, lay, P)
% SLAB_EMISSION  Up- and downward longwave emission of one atmospheric slab,
%                for an isothermal slab or a slab with an internal lapse rate.
%
%   [E_down, E_up, prof] = slab_emission(T_mean_K, eps_slab, lay, P)
%
%   T_mean_K   the slab's prognostic (mass-mean) temperature
%   eps_slab   its broadband emissivity (gases + cloud, already combined)
%   lay        the layer record from P.atm.layers (dz_m, p_mid_Pa, profile, lapse)
%   prof       diagnostics: Gamma_K_m, T_base_K, T_top_K
%
% ISOTHERMAL ('atm_profile' = 'isothermal'): the slab emits eps*B(T_mean) both
% ways, B from planck_emit (full or linearised).  This is the grey-slab
% atmosphere of Pierrehumbert (2010, sec. 4.3) and of the first version of this
% model.
%
% LAPSE-RATE PROFILE ('atm_profile' = 'lapse'): the temperature inside the slab is
% T(z) = T_mean + Gamma (z_mid - z), with Gamma from lapse_rate ('dry', 'moist'
% or a number; negative = inversion).  The slab is a homogeneous grey absorber
% of optical depth tau = -ln(1 - eps) spread uniformly over its depth dz, so the
% emission reaching its BASE from height z is attenuated by exp(-kappa z),
% kappa = tau/dz (Schwarzschild equation for a non-scattering layer, Goody &
% Yung 1989, ch. 2; the diffusivity factor is already inside eps):
%
%     E_down = int_0^dz  B(T(z)) kappa exp(-kappa z)        dz
%     E_up   = int_0^dz  B(T(z)) kappa exp(-kappa (dz - z)) dz
%
% For Gamma = 0 both integrals give B(T_mean)(1 - exp(-tau)) = eps B(T_mean),
% the isothermal result (checked in test T6).  With Gamma > 0 the base is warmer
% than the mean, so the slab emits MORE downward and LESS upward than an
% isothermal slab at the same mean temperature; a surface-based inversion does
% the opposite.  This is the only place the lapse rate enters the radiation;
% the turbulent fluxes use the profile's temperature at the reference height
% (seb_column_run, step 3).
%
% The integrals are evaluated on 41 nodes with B(T(z)) taken LINEAR between the
% nodes and the exponential kernel integrated EXACTLY on each interval (the
% "linear-in-tau" scheme of Clough, Iacono & Moncet 1992, JGR 97, 15761-15785,
% eq. 13):  with a = kappa h and h the node spacing,
%     int_0^h [B_i + (B_i+1 - B_i) u/h] kappa e^{-kappa u} du
%         = B_i (1 - e^{-a}) + (B_i+1 - B_i) [1 - (1 + a) e^{-a}] / a.
% An isothermal slab therefore returns eps B(T_mean) to round-off for ANY eps,
% and the black limit eps -> 1 returns B at the emitting face (no NaN at
% tau = inf, which a quadrature of the kernel gives).  The remaining error is
% the curvature of B between nodes, < 1e-6 relative here.  (Replaced a 41-point
% trapezoid rule, 0.09 % low at eps = 0.98 and NaN at eps = 1: 26 Oct 2026 review.)

Tref = lay.T_bg_K;                               % linearisation reference for this layer
if strcmpi(P.opt.atm_profile, 'isothermal') || eps_slab <= 0 || ~isfinite(lay.dz_m) || lay.dz_m <= 0
    B = planck_emit(T_mean_K, P, Tref);
    E_down = eps_slab * B;  E_up = E_down;
    prof.Gamma_K_m = 0;  prof.T_base_K = T_mean_K;  prof.T_top_K = T_mean_K;
    return
end

kind = lay.lapse;  if ischar(kind) && strcmpi(kind, 'auto'), kind = P.opt.lapse_kind; end
Gamma = lapse_rate(kind, T_mean_K, lay.p_mid_Pa, P);         % K/m, >0 cooling upward
dz    = lay.dz_m;
z     = linspace(0, dz, 41)';                                % height above the slab base
T_z   = T_mean_K + Gamma * (dz/2 - z);                       % mass-mean = mid-layer value
B_z   = planck_emit(T_z, P, Tref);
tau   = min(-log(1 - min(eps_slab, 1)), 700);                % optical depth, capped so eps = 1 stays finite
kappa = tau / dz;
E_down = kernel_integral(B_z, kappa, dz);                    % reaches the base (z = 0)
E_up   = kernel_integral(flipud(B_z), kappa, dz);            % reaches the top (z = dz): same integral, profile reversed
prof.Gamma_K_m = Gamma;  prof.T_base_K = T_z(1);  prof.T_top_K = T_z(end);
end

% ------------------------------------------------------------------------
function E = kernel_integral(B, kappa, dz)
% int_0^dz B(z) kappa exp(-kappa z) dz, B given at equally spaced nodes from
% z = 0 and linear between them, the exponential integrated exactly per interval.
n  = numel(B) - 1;  h = dz / n;  a = kappa * h;
em = -expm1(-a);                                             % 1 - e^{-a}, accurate for small a
w  = exp(-kappa * h * (0:n-1)');                             % e^{-kappa z_i} at the lower node of each interval
dB = diff(B(:));
E  = sum(w .* (B(1:n) * em + dB * (em - a * exp(-a)) / a));  % (em - a e^{-a})/a = [1 - (1 + a) e^{-a}]/a
end
