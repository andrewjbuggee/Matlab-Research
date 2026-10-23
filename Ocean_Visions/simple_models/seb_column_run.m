function R = seb_column_run(P, varargin)
% SEB_COLUMN_RUN  Integrate the coupled ocean / ice / atmosphere / cloud column.
%
%   R = seb_column_run(P)                          start from the default initial state
%   R = seb_column_run(P, 'init_state', S0)        continue from a saved state
%                                                  (S0 = R_previous.final_state)
%   R = seb_column_run(P, ..., 'verbose', false)
%
% Returns R with time series (one column per step) of every temperature, the ice
% layer temperatures, the thickness, every flux, and an energy-conservation
% diagnostic.  R.final_state can be handed to another call -- that is how the
% glaciation experiments branch a control and a perturbed run from the same
% spun-up state.  R.t_days is the time since the start of THIS run (so a
% branched run starts at 0); R.t_abs_days is the time since the original
% initialisation, carried in the state.
%
% =========================================================================
% THE EQUATIONS, layer by layer (all heat capacities C in J m-2 K-1, fluxes in
% W m-2; the slab bookkeeping follows the three-slab EBM of [PRISM])
% =========================================================================
%
%   upper atmosphere   C_ua dT_ua/dt = Q_LW,ua + C_ua (T_bg,ua - T_ua)/tau
%   cloud layer        C_cl dT_cl/dt = Q_LW,cl - gamma_mix (T_cl - T_la)
%                                      + C_cl (T_bg,cl - T_cl)/tau
%   lower atmosphere   C_la dT_la/dt = Q_LW,la + gamma_mix (T_cl - T_la)
%                                      - SH_down - LH_down + C_la (T_bg,la - T_la)/tau
%   skin (no C)                    0 = eps_s (LWD - sigma T_sk^4) + SH_down + LH_down
%                                      + F_cond_top_up                   [S25 eq. 5]
%   ice layers         rho c dT/dt = d/dz (k dT/dz),  T(top) = T_sk, T(base) = T_f
%   ice thickness      rho L dh/dt = F_cond_bot_up - F_ocean_up        [Stefan]
%   mixed layer        C_ml dT_ml/dt = gamma_deep (T_deep - T_ml) - F_ocean_up
%   deep ocean         C_deep dT_deep/dt = -gamma_deep (T_deep - T_ml)
%
% with   Q_LW,k  the net longwave heating of slab k from lw_grey_slabs,
%        SH_down, LH_down  the bulk turbulent fluxes (seb_surface_fluxes),
%        F_ocean_up = rho_w c_w c_h u* (T_ml - T_f)                   [McP92].
%
% WHY THE RELAXATION TERM EXISTS.  In polar night the atmosphere over sea ice
% has no energy source but horizontal advection; without it every slab would
% cool without limit.  Newtonian relaxation toward a background on a timescale
% tau_adv is the standard single-column stand-in (e.g. the Larcform SCM setups
% of Pithan et al. 2016, JAMES 8, 1345).  It is the ONE unavoidable
% parameterization here, and tau_adv is the knob that sets how tightly the air
% follows the surface -- the alpha of [LSM].
%
% NUMERICS.  Operator-split, one pass per step:
%   (1) longwave from the current slab temperatures -> LWD at the surface
%   (2) ocean heat flux from the current T_ml
%   (3) skin + ice: implicit conduction, Newton skin solve, Stefan growth, remap
%   (4) longwave again with the NEW skin temperature -> slab heating, OLR
%   (5) atmosphere and ocean slabs: forward Euler with the fluxes just computed
% Every flux is applied with equal and opposite sign to the two layers it
% connects, so the scheme conserves energy to the Newton tolerance; the
% diagnostic at the end checks that.  Forward Euler is safe because the slab
% relaxation times (C / sum of coupling coefficients) are hours to weeks, far
% longer than dt = 10 min.

opts.init_state = [];
opts.verbose    = true;
for i = 1:2:numel(varargin), opts.(varargin{i}) = varargin{i+1}; end

c = P.const;  a = P.atm;  o = P.ocean;  ice = P.ice;
sig = c.sigma_SB_W_m2_K4;
dt  = P.num.dt_s;
nt  = round(P.num.t_end_days * 86400 / dt);
N   = ice.N_layers;
T_f = o.T_f_K;

% ------------------------------------------------------------------------
% Initial state
% ------------------------------------------------------------------------
if isempty(opts.init_state)
    S.T_la_K = a.T_bg_la_K;  S.T_cl_K = a.T_bg_cl_K;  S.T_ua_K = a.T_bg_ua_K;
    if isfinite(ice.T_sk_prescribed_K), S.T_sk_K = ice.T_sk_prescribed_K;
    else,                                S.T_sk_K = ice.T_sk_init_K; end
    S.h_m = ice.h_init_m;
    % steady-conduction initial profile through snow (if any) and ice: the
    % temperature drop across each layer is proportional to its thermal
    % resistance (thickness / conductivity)
    R_snow = ice.h_snow_m / ice.k_snow_W_m_K;
    R_ice  = ice.h_init_m / ice.k_W_m_K;
    dT_tot = T_f - S.T_sk_K;
    T_snow_ice_interface = S.T_sk_K + dT_tot * R_snow / (R_snow + R_ice);
    S.T_snow_K = 0.5*(S.T_sk_K + T_snow_ice_interface);          % snow-cell centre
    depth_frac = ice.centre_fractions(:);            % ice layer centres, fraction of h
    S.T_ice_K  = T_snow_ice_interface + (T_f - T_snow_ice_interface) * depth_frac;
    S.T_ml_K   = o.T_ml_init_K;  S.T_deep_K = o.T_deep_init_K;
    S.t_s = 0;
else
    S = opts.init_state;
    if numel(S.T_ice_K) ~= N
        error('seb_column_run: init_state has %d ice layers, P has %d', numel(S.T_ice_K), N);
    end
end
% Slab emissivities, bottom to top. The cloud slab combines its gases with the
% cloud particles as independent absorbers (transmissions multiply), so a
% glaciated or vanished cloud still leaves the layer's clear-air emissivity.
eps_cloud_slab = 1 - (1 - a.emis_gas_cl) * (1 - P.cloud.emis);
eps_slabs = [a.emis_la; eps_cloud_slab; a.emis_ua];

% ------------------------------------------------------------------------
% Storage
% ------------------------------------------------------------------------
z = @() zeros(1, nt);
R.t_days = z();  R.t_abs_days = z();
R.T_sk_K = z(); R.T_la_K = z(); R.T_cl_K = z(); R.T_ua_K = z();
R.T_ml_K = z(); R.T_deep_K = z(); R.h_m = z(); R.T_snow_K = z();
R.T_ice_K = zeros(N, nt);
R.LWD_sfc_W_m2 = z(); R.LWU_sfc_W_m2 = z(); R.OLR_W_m2 = z();
R.SH_down_W_m2 = z(); R.LH_down_W_m2 = z();
R.F_cond_top_up_W_m2 = z(); R.F_cond_bot_up_W_m2 = z(); R.F_ocean_up_W_m2 = z();
R.dhdt_bot_m_s = z(); R.melt_W_m2 = z();
R.Q_LW_la_W_m2 = z(); R.Q_LW_cl_W_m2 = z(); R.Q_LW_ua_W_m2 = z();
R.Q_relax_W_m2 = z(); R.Q_mix_W_m2 = z();
R.lambda_SH_W_m2_K = z(); R.lambda_LH_W_m2_K = z(); R.lambda_LW_W_m2_K = z();
R.newton_iter = z(); R.E_ice_sens_J_m2 = z();

% ------------------------------------------------------------------------
% Time loop
% ------------------------------------------------------------------------
for n = 1:nt
    % (1) longwave with the current atmospheric state -> LWD at the surface.
    %     (LWD does not depend on T_sk, so this is exactly what the skin sees.)
    LW1 = lw_grey_slabs([S.T_la_K; S.T_cl_K; S.T_ua_K], eps_slabs, S.T_sk_K, ice.emis_sfc, sig);

    % (2) ocean heat flux into the ice base [McP92], or a fixed value
    if isfinite(o.F_ocean_fixed_W_m2)
        F_ocean_up = o.F_ocean_fixed_W_m2;
    else
        F_ocean_up = o.coeff_iceocean_W_m2_K * (S.T_ml_K - T_f);
    end

    % (3) skin balance + ice conduction + Stefan growth
    forcing.LWD_sfc_W_m2      = LW1.LWD_sfc;
    forcing.F_extra_down_W_m2 = a.F_extra_down_W_m2;
    forcing.T_la_K            = S.T_la_K;
    forcing.F_ocean_up_W_m2 = F_ocean_up;
    [S, F] = ice_column_step(S, forcing, P, dt);

    % (4) longwave again with the new skin temperature -> slab heating, OLR.
    %     LWD is unchanged; only the upward stream (and hence slab absorption)
    %     changes, so steps (3) and (4) see one consistent radiation field.
    LW2 = lw_grey_slabs([S.T_la_K; S.T_cl_K; S.T_ua_K], eps_slabs, S.T_sk_K, ice.emis_sfc, sig);
    Q_la = LW2.Q_slab(1);  Q_cl = LW2.Q_slab(2);  Q_ua = LW2.Q_slab(3);

    % (5) atmosphere slabs (forward Euler) ---------------------------------
    Q_mix = a.gamma_mix_W_m2_K * (S.T_cl_K - S.T_la_K);   % cloud -> sub-cloud
    rel_la = a.C_la_J_m2_K * (a.T_bg_la_K - S.T_la_K) / a.tau_adv_s;
    rel_cl = a.C_cl_J_m2_K * (a.T_bg_cl_K - S.T_cl_K) / a.tau_adv_s;
    rel_ua = a.C_ua_J_m2_K * (a.T_bg_ua_K - S.T_ua_K) / a.tau_adv_s;
    if strcmpi(a.mode, 'prognostic')
        S.T_la_K = S.T_la_K + dt/a.C_la_J_m2_K * (Q_la + Q_mix - F.SH_down_W_m2 - F.LH_down_W_m2 + rel_la);
        S.T_cl_K = S.T_cl_K + dt/a.C_cl_J_m2_K * (Q_cl - Q_mix + rel_cl);
        S.T_ua_K = S.T_ua_K + dt/a.C_ua_J_m2_K * (Q_ua + rel_ua);
        Q_relax  = rel_la + rel_cl + rel_ua;
    else
        % prescribed atmosphere: temperatures held; the implied external source
        % is whatever balances the budget (recorded for the energy diagnostic)
        Q_relax = -(Q_la + Q_cl + Q_ua - F.SH_down_W_m2 - F.LH_down_W_m2);
    end

    %     ocean slabs -----------------------------------------------------
    F_deep_up = o.gamma_deep_W_m2_K * (S.T_deep_K - S.T_ml_K);
    S.T_ml_K  = S.T_ml_K + dt/o.C_ml_J_m2_K * (F_deep_up - F_ocean_up);
    if o.deep_prognostic
        S.T_deep_K = S.T_deep_K - dt/o.C_deep_J_m2_K * F_deep_up;
    end
    S.t_s = S.t_s + dt;

    % store ------------------------------------------------------------------
    R.t_days(n)     = n * dt / 86400;          % time since the start of this run
    R.t_abs_days(n) = S.t_s / 86400;           % time since the original initial state
    R.T_sk_K(n) = S.T_sk_K;  R.T_la_K(n) = S.T_la_K;  R.T_cl_K(n) = S.T_cl_K;  R.T_ua_K(n) = S.T_ua_K;
    R.T_ml_K(n) = S.T_ml_K;  R.T_deep_K(n) = S.T_deep_K;  R.h_m(n) = S.h_m;
    if ice.h_snow_m > 0, R.T_snow_K(n) = S.T_snow_K; else, R.T_snow_K(n) = NaN; end
    R.T_ice_K(:, n) = S.T_ice_K;
    R.LWD_sfc_W_m2(n) = LW2.LWD_sfc;  R.LWU_sfc_W_m2(n) = LW2.LWU_sfc;  R.OLR_W_m2(n) = LW2.OLR;
    R.SH_down_W_m2(n) = F.SH_down_W_m2;  R.LH_down_W_m2(n) = F.LH_down_W_m2;
    R.F_cond_top_up_W_m2(n) = F.F_cond_top_up_W_m2;  R.F_cond_bot_up_W_m2(n) = F.F_cond_bot_up_W_m2;
    R.F_ocean_up_W_m2(n) = F_ocean_up;
    R.dhdt_bot_m_s(n) = F.dhdt_bot_m_s;  R.melt_W_m2(n) = F.melt_W_m2;
    R.Q_LW_la_W_m2(n) = Q_la;  R.Q_LW_cl_W_m2(n) = Q_cl;  R.Q_LW_ua_W_m2(n) = Q_ua;
    R.Q_relax_W_m2(n) = Q_relax;  R.Q_mix_W_m2(n) = Q_mix;
    R.lambda_SH_W_m2_K(n) = F.lambda_SH_W_m2_K;  R.lambda_LH_W_m2_K(n) = F.lambda_LH_W_m2_K;
    R.lambda_LW_W_m2_K(n) = F.lambda_LW_W_m2_K;
    R.newton_iter(n) = F.newton_iter;  R.E_ice_sens_J_m2(n) = F.E_ice_sens_J_m2;
end

R.P = P;
R.final_state = S;
R.eps_cloud = P.cloud.emis;  R.eps_cloud_slab = eps_cloud_slab;
R.L_eff_J_kg = F.L_eff_J_kg;

% ------------------------------------------------------------------------
% Energy-conservation diagnostic.  Total column energy (relative to fixed
% reference temperatures, which drop out of differences):
%   E = sum_atm C_k T_k + E_ice,sens(rel. to T_f) - rho L h + C_ml T_ml + C_deep T_deep
% Its change over the run must equal the time integral of the only external
% exchanges: -OLR at the top and the relaxation (advective) heating.  The
% ice-base latent term is in -rho*L*h.  Residual reported as a mean flux.
% ------------------------------------------------------------------------
E = a.C_la_J_m2_K*R.T_la_K + a.C_cl_J_m2_K*R.T_cl_K + a.C_ua_J_m2_K*R.T_ua_K ...
    + R.E_ice_sens_J_m2 - ice.rho_kg_m3*R.L_eff_J_kg*R.h_m ...
    + o.C_ml_J_m2_K*R.T_ml_K + o.C_deep_J_m2_K*R.T_deep_K;
if ~o.deep_prognostic
    % the deep slab is a fixed reservoir: its supply is an external source
    src_deep = o.gamma_deep_W_m2_K*(R.T_deep_K - R.T_ml_K);
else
    src_deep = 0*R.t_days;
end
if ~isempty(opts.init_state)
    E0 = energy_of_state(opts.init_state, P);
else
    E0 = NaN;      % unknown reference for the very first step; use step 1 onward
end
sources = -R.OLR_W_m2 + R.Q_relax_W_m2 + src_deep + a.F_extra_down_W_m2;   % F_extra is external
if isfinite(E0)
    dE_model = E - E0;                    % includes the first step
    dE_flux  = cumsum(sources) * dt;
else
    dE_model = E(2:end) - E(1);
    dE_flux  = cumsum(sources(2:end)) * dt;
end
R.energy.dE_model_J_m2 = dE_model;
R.energy.dE_flux_J_m2  = dE_flux;
R.energy.residual_mean_W_m2 = (dE_model(end) - dE_flux(end)) / (numel(dE_model)*dt);
R.energy.residual_max_W_m2  = max(abs(dE_model - dE_flux)) / (numel(dE_model)*dt);

if opts.verbose
    fprintf('seb_column_run: %d steps of %g s (%.1f days) | cloud eps = %.3f (%s, %g g/m2) | atm %s\n', ...
        nt, dt, P.num.t_end_days, P.cloud.emis, P.cloud.phase, P.cloud.WP_g_m2, a.mode);
    fprintf('   final: T_sk = %.2f K, T_la = %.2f K, T_cl = %.2f K, h = %.3f m, LWD = %.1f W/m2\n', ...
        S.T_sk_K, S.T_la_K, S.T_cl_K, S.h_m, LW2.LWD_sfc);
    fprintf('   energy residual: %.2e W/m2 (mean), %.2e W/m2 (max) | Newton iterations: max %d\n', ...
        R.energy.residual_mean_W_m2, R.energy.residual_max_W_m2, max(R.newton_iter));
end
end

% ------------------------------------------------------------------------
function E = energy_of_state(S, P)
% Column energy of a state, with the same references as above.
a = P.atm;  o = P.ocean;  ice = P.ice;
N = ice.N_layers;
[~, rc, L_eff] = ice_properties(S.T_ice_K, P);
E_ice = sum(rc .* (S.T_ice_K(:) - o.T_f_K) .* (S.h_m * ice.layer_fractions(:)));
if ice.h_snow_m > 0
    E_ice = E_ice + ice.rho_c_snow_J_m3_K*ice.h_snow_m*(S.T_snow_K - o.T_f_K);
end
E = a.C_la_J_m2_K*S.T_la_K + a.C_cl_J_m2_K*S.T_cl_K + a.C_ua_J_m2_K*S.T_ua_K ...
    + E_ice - ice.rho_kg_m3*L_eff*S.h_m + o.C_ml_J_m2_K*S.T_ml_K + o.C_deep_J_m2_K*S.T_deep_K;
end
