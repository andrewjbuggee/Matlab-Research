function [C_H, C_M, zeta, Ri_b] = most_exchange_coefficient(U_m_s, Tv_air_K, Tv_sfc_K, z_ref_m, z0m_m, z0h_m, P)
% MOST_EXCHANGE_COEFFICIENT  Surface-layer transfer coefficients from
%   Monin-Obukhov similarity, as implemented in the IFS (Cy50r1, Part IV, Ch. 3).
%
%   [C_H, C_M, zeta, Ri_b] = most_exchange_coefficient(U, Tv_air, Tv_sfc, z_ref, z0m, z0h, P)
%
% The sensible heat flux is then H = rho c_p C_H |U| (theta_air - T_sk)  [IFS eq. 8.8].
%
% THE PHYSICS.  In the surface layer the turbulent fluxes are constant with
% height and the dimensionless wind and temperature gradients are universal
% functions of zeta = z/L, L the Obukhov length (Monin & Obukhov 1954).
% Integrating those gradient functions from the roughness lengths to the
% reference height gives the bulk coefficients [IFS eqs. 3.16-3.17]:
%
%   C_M = kappa^2 / [ ln((z+z0m)/z0m) - Psi_M(zeta) + Psi_M(z0m zeta/(z+z0m)) ]^2
%   C_H = kappa^2 / { [ ... same ... ] [ ln((z+z0m)/z0h) - Psi_H(zeta) + Psi_H(z0h zeta/(z+z0m)) ] }
%
% with the stability functions of IFS Section 3.2.2:
%   unstable (zeta < 0): Dyer & Hicks gradient functions, Paulson (1970) integrals
%       Psi_M = pi/2 - 2 atan(x) + ln[(1+x)^2 (1+x^2)/8],  Psi_H = 2 ln[(1+x^2)/2],
%       x = (1 - 16 zeta)^(1/4)                                             [eq. 3.22]
%   stable (zeta > 0): Holtslag & de Bruin (1988) profile functions
%       Psi_M = -b (zeta - c/d) exp(-d zeta) - a zeta - b c/d
%       Psi_H = -b (zeta - c/d) exp(-d zeta) - (1 + 2a zeta/3)^1.5 - b c/d + 1
%       a = 1, b = 2/3, c = 5, d = 0.35                                     [eq. 3.24]
%
% zeta itself depends on the fluxes, so it is found from the bulk Richardson
% number [IFS eq. 3.26]  Ri_b = g z (Tv_air - Tv_sfc) / (Tv |U|^2)  by solving
% the implicit relation of IFS eq. (3.25) -- here by bisection on zeta, which
% is monotone in Ri_b.  As in the IFS, zeta is capped at 5 (Section 3.2, "the
% ratio z/L is limited to 5").
%
% WHAT IT DOES TO THE COUPLING.  Neutral: C_H = 1.6e-3 for z = 10 m, z0m = 1.3e-3,
% z0h = 1.3e-4 (the IFS ice-cap roughness lengths, Table 8.3, index 12).  A 4-K
% inversion over 10 m at U = 5 m/s gives Ri_b = 0.06 and C_H 40 % lower; the
% flux itself SATURATES near 32 W m-2 for 8-16 K inversions and falls slightly
% (-0.2 W m-2 K-1) before the zeta = 5 cap makes it linear again -- so SH is not
% monotone in T_sk; the skin solve still has a unique root because the Planck
% and conduction slopes (> 3 and > 20 W m-2 K-1) dominate.  Reference height:
% this is the coefficient between the surface and z_ref; feeding it a 400-m slab
% mean overstates the flux unless the slab carries a lapse-rate profile ('lapse').
%
% Sources: IFS Documentation Cy50r1 Part IV, Sections 3.2.1-3.2.3 (eqs. 3.15-3.27);
% Dyer (1974); Paulson (1970); Holtslag & de Bruin (1988), J. Appl. Meteor. 27,
% 689-704; Monin & Obukhov (1954).

kap = 0.4;                                     % von Karman constant
g   = P.const.g_m_s2;
U   = max(U_m_s, 0.1);                         % avoid the singular calm limit
z   = z_ref_m;
Tv  = 0.5*(Tv_air_K + Tv_sfc_K);
Ri_b = g * z * (Tv_air_K - Tv_sfc_K) / (Tv * U^2);           % [eq. 3.26]

lnm = log((z + z0m_m)/z0m_m);
lnh = log((z + z0m_m)/z0h_m);
fM = @(zt) lnm - psiM(zt) + psiM(z0m_m*zt/(z + z0m_m));
fH = @(zt) lnh - psiH(zt) + psiH(z0h_m*zt/(z + z0m_m));
% implicit relation Ri_b = zeta * fH / fM^2  [eq. 3.25], solved by bisection
if abs(Ri_b) < 1e-9
    zeta = 0;
elseif Ri_b > 0
    lo = 0; hi = 5;                             % stable; zeta capped at 5 as in the IFS
    if Ri_b >= hi*fH(hi)/fM(hi)^2, zeta = hi;
    else
        for it = 1:60
            mid = 0.5*(lo+hi);
            if mid*fH(mid)/fM(mid)^2 < Ri_b, lo = mid; else, hi = mid; end
        end
        zeta = 0.5*(lo+hi);
    end
else
    lo = -50; hi = 0;                           % unstable; g(zeta) = zeta fH/fM^2 increases with zeta
    for it = 1:60
        mid = 0.5*(lo+hi);
        if mid*fH(mid)/fM(mid)^2 < Ri_b, lo = mid; else, hi = mid; end   % g too negative -> move toward 0
    end
    zeta = 0.5*(lo+hi);
end
C_M = kap^2 / fM(zeta)^2;                        % [eq. 3.16]
C_H = kap^2 / (fM(zeta) * fH(zeta));             % [eq. 3.17]
end

function p = psiM(zt)
if zt < 0
    x = (1 - 16*zt)^0.25;
    p = pi/2 - 2*atan(x) + log((1 + x)^2 * (1 + x^2) / 8);       % [eq. 3.22]
else
    a = 1; b = 2/3; c = 5; d = 0.35;
    p = -b*(zt - c/d)*exp(-d*zt) - a*zt - b*c/d;                  % [eq. 3.24]
end
end

function p = psiH(zt)
if zt < 0
    x = (1 - 16*zt)^0.25;
    p = 2*log((1 + x^2)/2);                                       % [eq. 3.22]
else
    a = 1; b = 2/3; c = 5; d = 0.35;
    p = -b*(zt - c/d)*exp(-d*zt) - (1 + 2*a*zt/3)^1.5 - b*c/d + 1; % [eq. 3.24]
end
end
