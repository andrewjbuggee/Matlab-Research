function [E_W_m2, dEdT_W_m2_K] = planck_emit(T_K, P, T_ref_K)
% PLANCK_EMIT  Blackbody emission, full or linearised -- the 'lw_emission' option.
%
%   [E, dEdT] = planck_emit(T_K, P)              reference = P.opt.planck_T_ref_fallback_K
%   [E, dEdT] = planck_emit(T_K, P, T_ref_K)     explicit reference temperature
%
% P.opt.lw_emission = 'planck'      E = sigma T^4,            dE/dT = 4 sigma T^3
%                   = 'linearized'  E = sigma T_ref^4 + 4 sigma T_ref^3 (T - T_ref)
%
% THE PHYSICS.  Every layer radiates as a grey body: emission = emissivity x
% sigma T^4 (Stefan-Boltzmann).  The energy-balance climate models of Budyko
% (1969), Sellers (1969) and North et al. (1981) replace the quartic by a straight
% line in temperature so the balance can be solved by hand: North et al. (1981),
% eq. (7), write the infrared loss to space as I = A + B T with A = 203.3 W m-2
% and B = 2.09 W m-2 degC-1 fitted to observations (their sec. 2); the same
% section notes that a Taylor expansion of sigma (273 + T)^4 about 0 degC gives
% the 'black radiator' coefficients A_b = 314.9 W m-2 and B_b = 4.61 W m-2 degC-1,
% which is exactly the linearisation implemented here, about an arbitrary T_ref:
%
%     sigma T^4  ~  sigma T_ref^4 + 4 sigma T_ref^3 (T - T_ref)
%
% so B = 4 sigma T_ref^3 is the 'Planck feedback' of the layer: 4.61 W m-2 K-1 at
% 273 K, 3.54 at 250 K.  The IFS does the same thing to solve its skin balance
% implicitly (IFS Cy50r1, Part IV, Section 3.6: "linearization of the long wave
% radiation", eq. 3.84; and Section 8.3 after eq. 8.25).  The linearised option is
% here so the whole column can be read as a LINEAR system C dT/dt = -Lambda T + f,
% the form of the PRISM three-slab EBM; its error against the full form is a
% fraction (T - T_ref)^2 ... of the emission, i.e. ~1 % for a 10 K excursion.
%
% WHY T_ref MATTERS.  The linearisation is exact at T_ref and its slope is
% right there; 20 K away the emission is 2 % low and the slope 20 % low.  The
% callers pass the layer's own background temperature as T_ref (the atmosphere
% slabs) or the initial skin temperature (the surface), so no layer is far from
% its reference; P.opt.planck_T_ref_fallback_K is used only when none is given.

sig = P.const.sigma_SB_W_m2_K4;
switch lower(P.opt.lw_emission)
    case 'planck'
        E_W_m2      = sig .* T_K.^4;
        dEdT_W_m2_K = 4 .* sig .* T_K.^3;
    case 'linearized'
        if nargin < 3 || isempty(T_ref_K) || ischar(T_ref_K) || isstring(T_ref_K)
            T_ref_K = P.opt.planck_T_ref_fallback_K;
        end
        B = 4 .* sig .* T_ref_K.^3;                       % the Planck slope at T_ref
        E_W_m2      = sig .* T_ref_K.^4 + B .* (T_K - T_ref_K);
        dEdT_W_m2_K = B .* ones(size(T_K));
    otherwise
        error('planck_emit: unknown lw_emission option "%s"', P.opt.lw_emission);
end
end
