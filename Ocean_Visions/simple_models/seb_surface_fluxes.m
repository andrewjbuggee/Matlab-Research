function F = seb_surface_fluxes(T_sk_K, T_la_K, P)
% SEB_SURFACE_FLUXES  Bulk sensible and latent heat fluxes between the skin and
%                     the sub-cloud air.  Both are POSITIVE DOWNWARD (into the
%                     surface), the opposite of [S25] eq. (1) and the same as
%                     ERA5 / [LSM].  Say which convention you are in.
%
%   F = seb_surface_fluxes(T_sk_K, T_la_K, P)
%
%   F.SH_down_W_m2   sensible heat flux into the surface
%   F.LH_down_W_m2   latent heat flux into the surface (deposition > 0,
%                    sublimation < 0)
%   F.lambda_SH      d(SH_down)/d(T_la - T_sk) = rho c_p C_H U    [W m-2 K-1]
%   F.lambda_LH      d(LH_down)/dT_sk, the latent damping           [W m-2 K-1]
%   F.rho_air_kg_m3, F.q_la_kg_kg, F.q_sk_kg_kg   diagnostics
%
% BULK AERODYNAMIC FORMULAE (Monin-Obukhov similarity collapsed to one
% coefficient; e.g. [A10] for winter sea ice):
%
%   SH_down = rho_a c_p C_H U (T_la - T_sk)
%   LH_down = rho_a L_s C_E U (q_la - q_sat,ice(T_sk)),   q_la = RH_ice q_sat,ice(T_la)
%
% The product rho*c_p*C_H*U is the coupling coefficient lambda_SH of the
% linearised skin balance.  The ERA5 regressions over Barrow pack ice give 9.9
% W m-2 K-1 [LSM]; the defaults (rho = 1.3, C_H = 1.5e-3, U = 5) give 9.8.
%
% LIMITATION. C_H is held constant.  Over winter sea ice the real coefficient
% falls sharply with increasing stability (Grachev et al. 2007, BLM 124,
% 315-333; [A10]), so a strongly stable surface layer -- exactly what a cloud
% glaciation produces -- would couple more weakly than this. That makes the
% modelled SH response an UPPER bound.

c = P.const;  a = P.atm;

% air density at the surface pressure and sub-cloud temperature (ideal gas)
rho_air = a.p_sfc_Pa / (c.R_dry_J_kg_K * T_la_K);

% humidities over ice [MK05]
[~, q_sat_la] = sat_vapor_ice(T_la_K, a.p_sfc_Pa, c.eps_molar);
[~, q_sat_sk] = sat_vapor_ice(T_sk_K, a.p_sfc_Pa, c.eps_molar);
q_la = a.RH_ice_la * q_sat_la;

lambda_SH = rho_air * c.cp_air_J_kg_K * a.C_H * a.U_m_s;      % W m-2 K-1
SH_down   = lambda_SH * (T_la_K - T_sk_K);
LH_down   = rho_air * c.L_sub_J_kg * a.C_E * a.U_m_s * (q_la - q_sat_sk);

% latent damping coefficient: d(LH_down)/dT_sk = -rho L C_E U dq_sat/dT, with
% dq_sat/dT from Clausius-Clapeyron, q L_s/(R_v T^2), R_v = R_d/eps
R_v = c.R_dry_J_kg_K / c.eps_molar;
dq_dT = q_sat_sk * c.L_sub_J_kg / (R_v * T_sk_K^2);
lambda_LH = rho_air * c.L_sub_J_kg * a.C_E * a.U_m_s * dq_dT;

F.SH_down_W_m2 = SH_down;   F.LH_down_W_m2 = LH_down;
F.lambda_SH = lambda_SH;    F.lambda_LH = lambda_LH;
F.rho_air_kg_m3 = rho_air;  F.q_la_kg_kg = q_la;  F.q_sk_kg_kg = q_sat_sk;
end
