function LW = lw_grey_slabs(T_slabs_K, emis_slabs, T_sfc_K, emis_sfc, sigma)
% LW_GREY_SLABS  Longwave fluxes through a stack of grey (non-scattering) slabs.
%
%   LW = lw_grey_slabs(T_slabs_K, emis_slabs, T_sfc_K, emis_sfc, sigma)
%
%   T_slabs_K   [n x 1] slab temperatures, ordered BOTTOM (1) to TOP (n)
%   emis_slabs  [n x 1] slab broadband emissivities (= absorptivities, Kirchhoff)
%   T_sfc_K     surface (skin) temperature
%   emis_sfc    surface emissivity
%   sigma       Stefan-Boltzmann constant
%
%   LW.D(j)      downward flux at interface j   [W m-2], j = 1 is the surface,
%                j = k+1 is the top of slab k, so D(n+1) is the flux from space (0)
%   LW.U(j)      upward flux at interface j     [W m-2]; U(n+1) is the OLR
%   LW.LWD_sfc   = D(1), downwelling longwave at the surface
%   LW.LWU_sfc   = U(1), upwelling longwave at the surface (emission + reflection)
%   LW.OLR       = U(n+1)
%   LW.Q_slab(k) net longwave heating of slab k [W m-2] (absorbed - emitted)
%   LW.Q_sfc     net longwave into the surface, LWD - LWU
%   LW.residual  conservation check: Q_sfc + sum(Q_slab) + OLR  (should be 0)
%
% THE MODEL.  Each slab emits eps*sigma*T^4 up and down, absorbs the fraction eps
% of what falls on it, and transmits the rest (1 - eps); there is no scattering,
% which is an excellent approximation in the thermal infrared for gases and a
% good one for clouds [Z23; St78]. This is the classic grey-slab atmosphere
% (Pierrehumbert 2010, sec. 4.3 [P10]) extended to n slabs. The recursion is:
%
%     downward:  D(n+1) = 0;            D(k)   = eps_k sigma T_k^4 + (1-eps_k) D(k+1)
%     upward:    U(1) = eps_s sigma T_s^4 + (1-eps_s) D(1);
%                U(k+1) = eps_k sigma T_k^4 + (1-eps_k) U(k)
%     heating:   Q_k = eps_k [ U(k) + D(k+1) ] - 2 eps_k sigma T_k^4
%
% The surface reflects (1 - eps_s) of the downwelling flux, so LWU is not just
% emission; this matters at the 1-2 W m-2 level for eps_s = 0.985 [S25].
%
% Energy is conserved exactly by construction: everything emitted is either
% absorbed by another layer or leaves at the top. LW.residual verifies it.

T_slabs_K  = T_slabs_K(:);  emis_slabs = emis_slabs(:);
n = numel(T_slabs_K);
B_slab = sigma .* T_slabs_K.^4;             % blackbody emission of each slab [W m-2]
B_sfc  = sigma .* T_sfc_K.^4;

% --- downward sweep from space to the surface ---------------------------
D = zeros(n+1, 1);
D(n+1) = 0;                                 % no downwelling LW from space
for k = n:-1:1
    D(k) = emis_slabs(k)*B_slab(k) + (1 - emis_slabs(k))*D(k+1);
end

% --- upward sweep from the surface to space -----------------------------
U = zeros(n+1, 1);
U(1) = emis_sfc*B_sfc + (1 - emis_sfc)*D(1);   % emission + reflected downwelling
for k = 1:n
    U(k+1) = emis_slabs(k)*B_slab(k) + (1 - emis_slabs(k))*U(k);
end

% --- net heating of each slab and of the surface -------------------------
Q_slab = emis_slabs .* (U(1:n) + D(2:n+1)) - 2*emis_slabs.*B_slab;
Q_sfc  = D(1) - U(1);

LW.D = D;  LW.U = U;
LW.LWD_sfc = D(1);  LW.LWU_sfc = U(1);  LW.OLR = U(n+1);
LW.Q_slab = Q_slab;  LW.Q_sfc = Q_sfc;
LW.residual = Q_sfc + sum(Q_slab) + LW.OLR;   % = 0 to round-off
end
