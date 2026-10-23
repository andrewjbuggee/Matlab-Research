function [k_l, rho_c_l, L_eff] = ice_properties(T_K, P)
% ICE_PROPERTIES  Layer conductivity, volumetric heat capacity, effective latent heat.
%
%   [k_l, rho_c_l, L_eff] = ice_properties(T_K, P)
%
%   T_K      layer temperatures (any shape); outputs are column vectors
%   k_l      conductivity            [W m-1 K-1]
%   rho_c_l  volumetric heat capacity [J m-3 K-1]
%   L_eff    latent heat of fusion to use in the Stefan condition [J kg-1]
%
% DEFAULT (P.ice.use_brine = false): pure-ice constants, k = 2.03 W m-1 K-1
% (Untersteiner 1961; the IFS value [IFS eq. 8.149]) and rho*c = 917*2106 =
% 1.93e6 J m-3 K-1 [BL99 c0], L = 3.34e5 J/kg (fresh).
%
% BRINE OPTION (P.ice.use_brine = true), the corrections of Bitz & Lipscomb
% (1999) [BL99] for sea ice of bulk salinity S [ppt], with T in degC (< 0):
%   k_i   = k0 + beta*S/T                 eq. (6), after Untersteiner (1961):
%                                          brine pockets lower the conductivity,
%                                          most strongly near the melting point
%   c_i   = c0 + L0*mu*S/T^2              eq. (4): the heat that goes into
%                                          melting/freezing brine pockets as the
%                                          ice warms/cools shows up as a large
%                                          apparent heat capacity near 0 degC
%   L_eff = L0*(1 - T_m/T_f), T_m = -mu*S  the latent part of eq. (3) evaluated
%                                          at the ice base: ice that already
%                                          holds brine needs less heat removed
%                                          to form (~12 % less at S = 4 ppt)
% The 1/T singularity is guarded at -0.1 degC.
ice = P.ice;  c = P.const;
n = numel(T_K);
if ~ice.use_brine
    k_l     = ice.k_W_m_K      * ones(n, 1);
    rho_c_l = ice.rho_c_J_m3_K * ones(n, 1);
    L_eff   = c.L_fus_J_kg;
else
    T_degC  = min(T_K(:) - c.T0_K, -0.1);
    k_l     = ice.k_W_m_K + ice.beta_k_W_m_ppt * ice.S_ice_ppt ./ T_degC;
    c_l     = ice.c_J_kg_K + c.L_fus_J_kg * ice.mu_K_ppt * ice.S_ice_ppt ./ T_degC.^2;
    rho_c_l = ice.rho_kg_m3 * c_l;
    T_m     = -ice.mu_K_ppt * ice.S_ice_ppt;      % liquidus temperature [degC]
    T_fC    = P.ocean.T_f_K - c.T0_K;
    L_eff   = c.L_fus_J_kg * (1 - T_m / T_fC);
end
end
