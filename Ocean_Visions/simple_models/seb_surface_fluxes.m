function F = seb_surface_fluxes(T_sk_K, T_air_K, P)
% SEB_SURFACE_FLUXES  Sensible and latent heat fluxes between the skin and the
%                     air at the reference height -- the 'turb_scheme' option.
%                     Both POSITIVE DOWNWARD (into the surface), the opposite of
%                     [S25] eq. (1) and the same as ERA5 / [LSM].
%
%   F = seb_surface_fluxes(T_sk_K, T_air_K, P)
%
%   T_air_K   air temperature at the reference height z_ref (the caller takes it
%             from the lowest atmospheric layer: its mean for an isothermal slab,
%             the profile value at z_ref for a lapse-rate slab)
%
%   F.SH_down_W_m2, F.LH_down_W_m2   the fluxes
%   F.lambda_SH, F.lambda_LH          d(SH_down)/d(T_air - T_sk), -d(LH_down)/dT_sk  [W m-2 K-1]
%   F.C_H, F.zeta, F.Ri_b             the exchange coefficient and stability diagnostics
%   F.rho_air_kg_m3, F.q_air_kg_kg, F.q_sk_kg_kg
%
% FOUR SETTINGS (P.opt.turb_scheme), from least to most physics:
%
%  'none'        no turbulent exchange (the atmosphere of Thorndike 1992 talks to
%                the ice by radiation only); for analytic checks
%
%  'linear'      SH = lambda_SH (T_air - T_sk) with a PRESCRIBED coefficient
%                P.opt.lambda_SH_lin_W_m2_K, and LH = lambda_LH (T_air - T_sk)
%                with lambda_LH = lambda_SH (L_s/c_p) dq_sat/dT [LSM]: the
%                regression coupling of the ERA5 analysis (9.9 W m-2 K-1 over
%                pack ice, 2-m air); the column is then linear.  NOTE: this LH
%                ignores RH_ice_la (saturated air assumed), unlike 'bulk_const'.
%
%  'bulk_const'  bulk aerodynamic formulae with CONSTANT coefficients
%                   SH_down = rho c_p C_H U (T_air - T_sk)              [IFS eq. 8.8 without the g z/c_p term]
%                   LH_down = rho L_s C_E U (q_air - q_sat,ice(T_sk))   [IFS eq. 8.9]
%                q_air = RH_ice q_sat,ice(T_air) [A02].  C_H = P.atm.C_H, an
%                EFFECTIVE value for exchange with a thick slab (see the
%                parameter file).  The default.
%
%  'most'        the same bulk formulae with C_H from Monin-Obukhov similarity
%                exactly as the IFS computes it (most_exchange_coefficient.m: IFS
%                Cy50r1 Part IV eqs. 3.15-3.27; Dyer & Hicks unstable, Holtslag &
%                de Bruin 1988 stable functions, bulk-Richardson closure), with
%                the potential-temperature correction g z_ref/c_p of IFS eq. 8.8.
%                The coefficient now FALLS as the surface layer stabilises --
%                the physics the constant scheme lacks -- and refers to the
%                height z_ref (10 m), so it belongs with a lapse-rate profile
%                that supplies T(z_ref), not with a 400-m slab mean.
%
% P.opt.latent_scheme = 'bulk' (as above) | 'off' (LH = 0, for the bare
% three-slab model of the 22 Oct cartoon).
%
% Saturation over ice from Murphy & Koop (2005) [MK05]; the latent damping
% coefficient lambda_LH = rho L_s C_E U dq_sat/dT with dq_sat/dT = q L_s/(R_v T^2)
% (Clausius-Clapeyron).

c = P.const;  a = P.atm;  o = P.opt;

rho_air = a.p_sfc_Pa / (c.R_dry_J_kg_K * T_air_K);          % ideal gas at the surface pressure
[~, q_sat_air] = sat_vapor_ice(T_air_K, a.p_sfc_Pa, c.eps_molar);
[~, q_sat_sk]  = sat_vapor_ice(T_sk_K,  a.p_sfc_Pa, c.eps_molar);
q_air = a.RH_ice_la * q_sat_air;
R_v   = c.R_dry_J_kg_K / c.eps_molar;
dq_dT = q_sat_sk * c.L_sub_J_kg / (R_v * T_sk_K^2);          % Clausius-Clapeyron over ice

C_H = NaN;  zeta = NaN;  Ri_b = NaN;
switch lower(o.turb_scheme)
    case 'none'                     % no turbulent exchange at all (Thorndike 1992 neglects it)
        lambda_SH = 0;  lambda_LH = 0;  SH_down = 0;  LH_down = 0;
    case 'linear'
        lambda_SH = o.lambda_SH_lin_W_m2_K;
        if ischar(o.lambda_LH_lin) || isstring(o.lambda_LH_lin)      % 'clausius'
            lambda_LH = lambda_SH * (c.L_sub_J_kg / c.cp_air_J_kg_K) * dq_dT;
        else
            lambda_LH = o.lambda_LH_lin;
        end
        SH_down = lambda_SH * (T_air_K - T_sk_K);
        LH_down = lambda_LH * (T_air_K - T_sk_K);
    case 'bulk_const'
        C_H = a.C_H;  C_E = a.C_E;
        lambda_SH = rho_air * c.cp_air_J_kg_K * C_H * a.U_m_s;
        SH_down   = lambda_SH * (T_air_K - T_sk_K);
        LH_down   = rho_air * c.L_sub_J_kg * C_E * a.U_m_s * (q_air - q_sat_sk);
        lambda_LH = rho_air * c.L_sub_J_kg * C_E * a.U_m_s * dq_dT;
    case 'most'
        % potential temperature of the air at z_ref, and virtual temperatures
        theta_air = T_air_K + c.g_m_s2 * o.z_ref_m / c.cp_air_J_kg_K;   % [IFS eq. 8.8]
        Tv_air = theta_air * (1 + 0.61*q_air);
        Tv_sfc = T_sk_K    * (1 + 0.61*q_sat_sk);
        [C_H, ~, zeta, Ri_b] = most_exchange_coefficient(a.U_m_s, Tv_air, Tv_sfc, ...
                                   o.z_ref_m, o.z0m_m, o.z0h_m, P);
        C_E = C_H;                                                 % z0q = z0h over ice
        lambda_SH = rho_air * c.cp_air_J_kg_K * C_H * a.U_m_s;
        SH_down   = lambda_SH * (theta_air - T_sk_K);
        LH_down   = rho_air * c.L_sub_J_kg * C_E * a.U_m_s * (q_air - q_sat_sk);
        lambda_LH = rho_air * c.L_sub_J_kg * C_E * a.U_m_s * dq_dT;
    otherwise
        error('seb_surface_fluxes: unknown turb_scheme "%s"', o.turb_scheme);
end
if strcmpi(o.latent_scheme, 'off')
    LH_down = 0;  lambda_LH = 0;
end

F.SH_down_W_m2 = SH_down;   F.LH_down_W_m2 = LH_down;
F.lambda_SH = lambda_SH;    F.lambda_LH = lambda_LH;
F.C_H = C_H;  F.zeta = zeta;  F.Ri_b = Ri_b;
F.rho_air_kg_m3 = rho_air;  F.q_air_kg_kg = q_air;  F.q_sk_kg_kg = q_sat_sk;
end
