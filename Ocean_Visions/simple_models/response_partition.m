function D = response_partition(Rp, Rc)
% RESPONSE_PARTITION  Perturbed minus control, and the partition of the LWD forcing.
%
%   D = response_partition(R_perturbed, R_control)
%
% At every instant the skin balance (or the slab budget) gives
%     dLWD = dLWU - dSH_down - dLH_down - dF_cond_top_up
% so the four shares f_LWU + f_SH + f_LH + f_G = 1.  The ice energy budget
%     F_cond_top_up = F_ocean_up + rho L dh/dt - dE_ice/dt
% splits the extra conduction out of the ice into extra FREEZING (latent heat
% released at the base) and extra COOLING of the ice (sensible energy removed),
% plus any change in the ocean heat flux.  Cumulative energies are in MJ m-2;
% negative = energy removed from the surface.  dh_zapalac is the bound in which
% every W m-2 of forcing freezes water, dh = |E_forcing| / (rho L) [Z23].
rho_i = Rp.P.ice.rho_kg_m3;
dt_s = (Rp.t_days(2) - Rp.t_days(1)) * 86400;
D.t_h   = Rp.t_days * 24;
D.dLWD  = Rp.LWD_sfc_W_m2 - Rc.LWD_sfc_W_m2;
D.dLWU  = Rp.LWU_sfc_W_m2 - Rc.LWU_sfc_W_m2;
D.dSH   = Rp.SH_down_W_m2 - Rc.SH_down_W_m2;
D.dLH   = Rp.LH_down_W_m2 - Rc.LH_down_W_m2;
D.dG    = Rp.F_cond_top_up_W_m2 - Rc.F_cond_top_up_W_m2;
D.dT_sk = Rp.T_sk_K - Rc.T_sk_K;
D.dT_air = Rp.T_atm_K(1,:) - Rc.T_atm_K(1,:);
D.dh    = Rp.h_m - Rc.h_m;
D.f_LWU =  D.dLWU ./ D.dLWD;
D.f_SH  = -D.dSH  ./ D.dLWD;
D.f_LH  = -D.dLH  ./ D.dLWD;
D.f_G   = -D.dG   ./ D.dLWD;
cs = @(x) cumsum(x) * dt_s * 1e-6;
D.E_forcing = cs(D.dLWD);
D.E_LWU     = cs(D.dLWU);
D.E_SH      = -cs(D.dSH);
D.E_LH      = -cs(D.dLH);
D.E_G       = -cs(D.dG);
D.E_growth  = -rho_i * Rp.L_eff_J_kg * D.dh * 1e-6;
D.E_icecool = (Rp.E_ice_sens_J_m2 - Rc.E_ice_sens_J_m2) * 1e-6;
D.E_ocean   = -cs(Rp.F_ocean_up_W_m2 - Rc.F_ocean_up_W_m2);
D.dh_zapalac = -D.E_forcing * 1e6 / (rho_i * Rp.L_eff_J_kg);
D.alpha_air = D.dT_air(end) / D.dT_sk(end);
end
