function [S, F] = ice_slab_step(S, forcing, P, dt_s)
% ICE_SLAB_STEP  One time step of the SIMPLEST ice representation: a single slab
%                at one temperature, with no separate skin -- 'ice_model' = 'slab'.
%
%   [S, F] = ice_slab_step(S, forcing, P, dt_s)
%
% This is the ice of the simplest cartoon in the 26 Oct 2026 notes: one box with
% a temperature T_i and a thickness h.  Longwave, sensible and latent heat act on
% the slab directly; conduction connects the slab centre to the ice-ocean
% interface held at the freezing point; growth is the Stefan condition at that
% interface.  It is the ice-model analogue of the single-slab energy balance of
% Budyko (1969) / North et al. (1981): one heat capacity, one temperature.
%
%     rho c h  dT_i/dt = eps (LWD - B(T_i)) + F_extra + SH_down + LH_down + F_bot_up   (slab budget)
%     F_bot_up         = k (T_f - T_i) / (h/2)                                          (Fourier, half-slab)
%     rho L   dh/dt    = F_bot_up - F_ocean_up                                          (Stefan)
%
% with B the (full or linearised) Planck emission.  New ice forms at T_f and is
% mixed into the slab conserving enthalpy.  The step is linearly implicit in T_i
% (the Planck, turbulent and conductive slopes go into the denominator), so it
% is stable for any dt.
%
% WHAT IT CANNOT DO: there is no skin, so the emitting temperature is the slab
% mean -- warmer than a real surface in winter -- and no internal gradient, so
% the response to a forcing is spread over the whole slab's heat capacity at
% once (no delay, no 1/sqrt(t) admittance).  Compare with ice_column_step.  The
% hierarchy exists so this difference is a number, not an argument.
%
% Sign conventions as in ice_column_step: SH_down, LH_down, F_extra positive INTO
% the ice; F_bot_up, F_ocean_up positive UPWARD; dh/dt positive for growth.
% F_cond_top_up is reported as minus the net surface flux so that the same
% partition bookkeeping (dE_ice/dt = F_bot_up - F_cond_top_up) holds for both
% ice models.

ice = P.ice;  T_f = P.ocean.T_f_K;
T_i = S.T_ice_K(1);  h = S.h_m;
[k_l, rho_c, L_eff] = ice_properties(T_i, P);
C_slab = rho_c * h;                                   % J m-2 K-1

[E_i, dEdT] = planck_emit(T_i, P, ice.T_sk_init_K);
SF = seb_surface_fluxes(T_i, forcing.T_air_K, P);
K_bot = k_l / (h/2);                                  % conductance base -> slab centre
F_bot_up = K_bot * (T_f - T_i);
net = ice.emis_sfc*(forcing.LWD_sfc_W_m2 - E_i) + forcing.F_extra_down_W_m2 ...
      + SF.SH_down_W_m2 + SF.LH_down_W_m2 + F_bot_up;
% linearly implicit update: every term that falls with T_i goes in the denominator
damp  = ice.emis_sfc*dEdT + SF.lambda_SH + SF.lambda_LH + K_bot;
T_new = T_i + dt_s * net / (C_slab + dt_s*damp);
melt = 0;
if T_new > ice.T_melt_sfc_K                           % melt regime [S25]
    T_new = ice.T_melt_sfc_K;
    melt  = net - damp*(T_new - T_i) - C_slab*(T_new - T_i)/dt_s;   % surplus AFTER warming the slab to T_melt (heat-capacity term added, 26 Oct 2026 review)
end
% The fluxes the rest of the column sees are the ones the implicit update
% actually applied -- the old-state values corrected by their slopes times the
% temperature change -- so that the energy the slab gained is exactly the
% energy the atmosphere and the ocean lost (conservation to round-off).
dT = T_new - T_i;
E_new    = E_i + dEdT*dT;
SH_new   = SF.SH_down_W_m2 - SF.lambda_SH*dT;
LH_new   = SF.LH_down_W_m2 - SF.lambda_LH*dT;
F_bot_up = K_bot * (T_f - T_new);                     % basal flux at the new temperature
dhdt_bot = (F_bot_up - forcing.F_ocean_up_W_m2) / (ice.rho_kg_m3 * L_eff);
dhdt_top = -melt / (ice.rho_kg_m3 * L_eff);
if strcmpi(P.opt.ice_growth, 'fixed'), dhdt_bot = 0; dhdt_top = 0; end   % constant thickness, as ERA5
dh = (dhdt_bot + dhdt_top) * dt_s;
h_new = h + dh;
if h_new < ice.h_min_m
    error('ice_slab_step:vanished', 'ice thickness fell below h_min (%.3f m)', ice.h_min_m);
end
% new ice forms at T_f: mix it into the slab conserving enthalpy
if dhdt_bot > 0
    T_new = (T_new*h + T_f*dhdt_bot*dt_s) / (h + dhdt_bot*dt_s);
end

S.T_ice_K = T_new;
S.T_sk_K  = T_new;                                    % no skin: the slab IS the surface
S.h_m     = h_new;

F.LWU_sfc_W_m2       = ice.emis_sfc*E_new + (1 - ice.emis_sfc)*forcing.LWD_sfc_W_m2;
F.SH_down_W_m2       = SH_new;
F.LH_down_W_m2       = LH_new;
F.F_cond_top_up_W_m2 = -(ice.emis_sfc*(forcing.LWD_sfc_W_m2 - E_new) + forcing.F_extra_down_W_m2 ...
                         + SH_new + LH_new) + melt;
F.F_cond_bot_up_W_m2 = F_bot_up;
F.melt_W_m2          = melt;
F.dhdt_bot_m_s       = dhdt_bot;
F.dhdt_top_m_s       = dhdt_top;   F.E_runoff_J_m2 = rho_c * (T_new - T_f) * (max(-dhdt_top*dt_s, 0) + max(-dhdt_bot*dt_s, 0));   % sensible heat leaving with melted ice (top and base); energy diagnostic
F.lambda_SH_W_m2_K   = SF.lambda_SH;
F.lambda_LH_W_m2_K   = SF.lambda_LH;
F.lambda_LW_W_m2_K   = ice.emis_sfc*dEdT;
F.C_H = SF.C_H;  F.zeta = SF.zeta;  F.Ri_b = SF.Ri_b;
F.L_eff_J_kg         = L_eff;
F.newton_iter        = 0;
F.skin_residual_W_m2 = 0;
F.E_ice_sens_J_m2    = rho_c * h_new * (T_new - T_f);
end
