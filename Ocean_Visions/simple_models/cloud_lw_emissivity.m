function [eps_c, k_abs_m2_g] = cloud_lw_emissivity(WP_g_m2, phase, P, r_e_um)
% CLOUD_LW_EMISSIVITY  Broadband longwave emissivity of a cloud layer.
%
%   [eps_c, k_abs] = cloud_lw_emissivity(WP_g_m2, phase, P)
%   [eps_c, k_abs] = cloud_lw_emissivity(WP_g_m2, phase, P, r_e_um)
%
%   WP_g_m2  condensed water path [g m-2] (LWP for liquid, IWP for ice)
%   phase    'liquid' | 'ice' | 'none'
%   P        parameter struct from seb_column_params
%   r_e_um   effective radius [um]; defaults to P.cloud.r_e_liq_um / r_e_ice_um
%
% THE PHYSICS.  A cloud is a grey absorber in the thermal infrared with
% negligible scattering [Z23 sec. "Radiative forcing"; St78], so its emissivity
% follows Beer's law in the water path:
%
%       eps_c = 1 - exp( -k_abs * WP )                       [St78 eq. 13]
%
% k_abs is a MASS absorption coefficient for hemispheric flux [m2 g-1].
%
% LIQUID.  Stephens (1978) fitted k_abs = 0.130 m2 g-1 for downward emission
% (0.158 upward) to detailed calculations for stratiform water clouds [St78].
% Zapalac (2023) uses 0.156 [Z23]. The downward value is used here because the
% surface forcing is what this model is about.  At LWP = 30 g m-2 it gives
% eps = 0.98, which is why ~30 g m-2 is where the surface stops responding to
% more liquid -- the "black" threshold that Shupe & Intrieri (2004) and [Z23]
% quote.
%
% WHY k_abs SCALES AS 1/r_e.  In the large-particle (geometric-optics) limit the
% absorption cross-section of a sphere is Q_abs*pi*r^2 with Q_abs -> 1, and its
% mass is (4/3)*pi*r^3*rho, so the cross-section per unit mass is
%
%       k_beam = 3*Q_abs / (4*rho*r_e)                      (per unit mass)
%
% and the flux (hemispheric) coefficient is larger by the diffusivity factor
% 1.66 (Elsasser 1942; Goody & Yung 1989 sec. 2.4): k_abs = 1.66*k_beam.
% CHECK: for liquid, r_e = 10 um, rho = 1000 kg/m3:  1.66 * 3/(4*1e6 g/m3*1e-5 m)
% = 0.125 m2 g-1, within 4 % of the Stephens value. So the r_e dependence used
% here is the geometric-optics scaling ANCHORED at the Stephens fit:
%
%       k_abs(liquid, r_e) = 0.130 * (10 um / r_e)
%
% This is an extrapolation away from r_e = 10 um; Q_abs is not exactly 1 for
% 10-um drops across the 8-12 um window, so the agreement above is partly
% fortuitous.  (Relevant to the month-to-month shift of the LWP at which DLR
% saturates seen in the Utqiagvik data -- a larger r_e moves the asymptote to
% larger LWP.)
%
% ICE.  No published fit is used; the coefficient is the geometric-optics
% estimate above with rho_ice = 917 kg/m3 and r_e(ice) = 30 um by default:
%
%       k_abs(ice) = 1.66 * 3*Q_abs / (4*rho_ice*r_e)  = 0.045 m2 g-1 at 30 um
%
% i.e. about one third of the liquid value per gram.  Treat it as an ESTIMATE
% (the choice of r_e for ice is the main uncertainty; for a proper
% parameterization see Ebert & Curry 1992, JGR 97, 3831-3836, or Fu et al.
% 1998). The glaciation experiments do not hinge on it, because the dominant
% effect of glaciation in the LES of [V22] is the LOSS OF CONDENSED MASS to
% precipitation (TWC 0.170 -> 0.009 g/kg), not the per-gram change.

if nargin < 4 || isempty(r_e_um)
    switch lower(phase)
        case 'liquid', r_e_um = P.cloud.r_e_liq_um;
        case 'ice',    r_e_um = P.cloud.r_e_ice_um;
        otherwise,     r_e_um = NaN;
    end
end

switch lower(phase)
    case 'liquid'
        % Stephens (1978) anchor at 10 um, geometric-optics 1/r_e scaling
        k_abs_m2_g = P.cloud.k_abs_liq_ref_m2_g * (10.0 / r_e_um);
    case 'ice'
        % geometric optics: 1.66 * 3 Q / (4 rho r_e); convert rho to g m-3 and r_e to m
        rho_g_m3   = P.cloud.rho_ice_kg_m3 * 1e3;
        k_abs_m2_g = P.cloud.diffusivity * 3 * P.cloud.Q_abs_ice ...
                     / (4 * rho_g_m3 * r_e_um * 1e-6);
    case 'none'
        k_abs_m2_g = 0;
    otherwise
        error('cloud_lw_emissivity: unknown phase "%s"', phase);
end

eps_c = 1 - exp(-k_abs_m2_g .* WP_g_m2);
end
