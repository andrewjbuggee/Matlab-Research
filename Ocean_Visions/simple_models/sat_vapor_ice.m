function [e_i_Pa, q_i_kg_kg] = sat_vapor_ice(T_K, p_Pa, eps_molar)
% SAT_VAPOR_ICE  Saturation vapour pressure and specific humidity over ICE.
%
%   [e_i_Pa, q_i_kg_kg] = sat_vapor_ice(T_K, p_Pa, eps_molar)
%
% Murphy & Koop (2005), QJRMS 131, 1539-1565, eq. (7), valid for T > 110 K:
%
%   ln e_i [Pa] = 9.550426 - 5723.265/T + 3.53068 ln T - 0.00728332 T
%
% This is the reference formulation for ice (it reproduces the ITS-90 triple
% point, 611.657 Pa at 273.16 K, to better than 0.1 %). The surface over sea
% ice sublimates/deposits, so the skin saturation is over ice, and the sub-cloud
% air over winter pack ice is itself near ice saturation [A02].
%
% Specific humidity from vapour pressure (exact for an ideal mixture):
%   q = eps e / (p - (1 - eps) e)

e_i_Pa    = exp(9.550426 - 5723.265./T_K + 3.53068*log(T_K) - 0.00728332*T_K);
q_i_kg_kg = eps_molar .* e_i_Pa ./ (p_Pa - (1 - eps_molar).*e_i_Pa);
end
