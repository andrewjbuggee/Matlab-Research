function Gamma_K_m = lapse_rate(kind, T_K, p_Pa, P)
% LAPSE_RATE  Temperature lapse rate -dT/dz inside an atmospheric layer.
%
%   Gamma = lapse_rate('dry',   T_K, p_Pa, P)   g / c_p                  = 9.76 K km-1
%   Gamma = lapse_rate('moist', T_K, p_Pa, P)   saturated (pseudo-)adiabatic
%   Gamma = lapse_rate(value,   T_K, p_Pa, P)   prescribed, K m-1
%
% SIGN: Gamma > 0 means temperature DECREASES with height (the usual
% convention, Wallace & Hobbs 2006, sec. 3.4).  A surface-based inversion, the
% normal state over winter sea ice, is Gamma < 0 and is given as a number.
%
% DRY ADIABATIC: Gamma_d = g / c_p, from the first law for an adiabatic parcel in
% hydrostatic balance (Wallace & Hobbs 2006, eq. 3.53).
%
% SATURATED ADIABATIC (Wallace & Hobbs 2006, sec. 3.5.3; AMS Glossary,
% "saturation-adiabatic lapse rate"):
%
%     Gamma_s = g * (1 + L r_s / (R_d T)) / (c_p + L^2 r_s eps / (R_d T^2))
%
% with r_s the saturation mixing ratio at (T, p) over LIQUID water (a mixed-phase
% cloud is held at water saturation while its liquid persists; the vapour
% pressure over supercooled water is Murphy & Koop 2005, eq. 10) and L the
% latent heat of vaporisation.  At 250 K and 950 hPa r_s = 0.7 g kg-1 and
% Gamma_s = 8.9 K km-1 -- almost the dry value, because so little vapour is
% available to condense at these temperatures.  The moist/dry distinction
% therefore matters much less here than in the tropics; the option is kept
% because the request asked for it and because it is what a cloud-topped,
% well-mixed boundary layer follows.

c = P.const;
if isnumeric(kind)
    Gamma_K_m = kind;
    return
end
switch lower(kind)
    case 'dry'
        Gamma_K_m = c.g_m_s2 / c.cp_air_J_kg_K;
    case 'moist'
        % saturation vapour pressure over supercooled liquid water [MK05 eq. 10]
        e_w = exp(54.842763 - 6763.22./T_K - 4.210*log(T_K) + 0.000367*T_K ...
              + tanh(0.0415*(T_K - 218.8)) .* (53.878 - 1331.22./T_K - 9.44523*log(T_K) + 0.014025*T_K));
        r_s = c.eps_molar * e_w ./ (p_Pa - e_w);          % saturation mixing ratio [kg/kg]
        L   = c.L_vap_J_kg;
        Rd  = c.R_dry_J_kg_K;
        Gamma_K_m = c.g_m_s2 * (1 + L*r_s./(Rd*T_K)) ./ (c.cp_air_J_kg_K + L^2*r_s*c.eps_molar./(Rd*T_K.^2));
    case 'isothermal'
        Gamma_K_m = 0;
    otherwise
        error('lapse_rate: unknown kind "%s"', kind);
end
end
