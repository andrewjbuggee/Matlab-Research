function R = seb_column_run(P, varargin)
% SEB_COLUMN_RUN  Integrate the coupled ocean / ice / atmosphere / cloud column,
%                 at whatever complexity level and with whatever physics options
%                 P carries (model_config.m, seb_column_params.m).
%
%   R = seb_column_run(P)                          start from the default initial state
%   R = seb_column_run(P, 'init_state', S0)        continue from a saved state
%                                                  (S0 = R_previous.final_state)
%   R = seb_column_run(P, ..., 'verbose', false)
%
% Returns R with time series (one column per step) of every temperature, the ice
% layer temperatures, the thickness, every flux, and an energy-conservation
% diagnostic.  R.final_state can be handed to another call -- that is how the
% experiments branch a control and a perturbed run from the same spun-up state.
% R.t_days is the time since the start of THIS run; R.t_abs_days since the
% original initialisation.
%
% =========================================================================
% THE EQUATIONS, layer by layer.  All heat capacities C in J m-2 K-1, fluxes in
% W m-2; the slab bookkeeping follows the three-slab EBM of [PRISM] and the
% zero-dimensional models of North et al. (1981, sec. 2): C dT/dt = sum of fluxes.
% =========================================================================
%
%   atmosphere layer k (bottom = 1, top = n), prognostic mean temperature T_k:
%     C_k dT_k/dt = Q_LW,k                         net longwave (lw_grey_slabs)
%                   - gamma_k (T_k - T_k+1) + gamma_k-1 (T_k-1 - T_k)   turbulent mixing
%                   - (SH_down + LH_down) [k = 1 only]      the surface takes them
%                   + C_k (T_bg,k - T_k) / tau_adv          advection (relaxation)
%   skin (no C)   0 = eps_s (LWD - B(T_sk)) + F_extra + SH_down + LH_down + F_cond_top_up
%                                                                [S25 eq. 5; IFS eq. 3.76]
%   ice           'skin_profile': rho c dT/dt = d/dz (k dT/dz), T(top) = T_sk, T(base) = T_f
%                 'slab':         rho c h dT_i/dt = surface fluxes + k (T_f - T_i)/(h/2)
%   thickness     rho L dh/dt = F_cond_bot_up - F_ocean_up                         [Stefan]
%   mixed layer   C_ml dT_ml/dt = gamma_deep (T_deep - T_ml) - F_ocean_up
%   deep ocean    C_deep dT_deep/dt = -gamma_deep (T_deep - T_ml)    (or held fixed)
%
% with B the Planck emission (full or linearised, planck_emit.m), the layer
% emissions from slab_emission.m (isothermal or lapse-rate profile), SH/LH from
% seb_surface_fluxes.m (three schemes), F_ocean_up = rho_w c_w c_h u* (T_ml - T_f)
% [McP92].
%
% WHY THE ADVECTION TERM EXISTS.  In polar night the atmosphere over sea ice
% has no energy source but horizontal advection; without it every slab would
% cool without limit.  Two ways to represent it (P.opt.atm_forcing):
%   'relaxation'  Newtonian relaxation toward a background on a timescale
%                 tau_adv, the standard single-column stand-in (e.g. the
%                 Larcform SCM setups of Pithan et al. 2016, JAMES 8, 1345);
%                 tau_adv sets how tightly the air follows the surface.
%   'heating'     a prescribed heat convergence D (W m-2) shared among the
%                 layers, the closure of Thorndike (1992, sec. 2): the layer
%                 temperatures are then whatever radiative balance makes them,
%                 and the one external number is observable (Nakamura & Oort
%                 1988).  The atmosphere is free to follow the surface.
% Either way this is the ONE unavoidable parameterization of a column model.
%
% NUMERICS.  Operator-split, one pass per step:
%   (1) each layer's up/down emission, then the longwave sweep -> LWD at the surface
%   (2) the air temperature at the reference height, from the lowest layer
%   (3) ocean heat flux from the current T_ml
%   (4) skin + ice: implicit conduction, Newton skin solve, Stefan growth, remap
%       (or the linearly implicit slab step)
%   (5) longwave again with the NEW surface emission -> layer heating, OLR
%   (6) atmosphere and ocean slabs: forward Euler with the fluxes just computed
% Every flux is applied with equal and opposite sign to the two layers it
% connects, so the scheme conserves energy to the Newton tolerance; the
% diagnostic at the end checks that.  Forward Euler is safe because the slab
% relaxation times (C / sum of coupling coefficients) are hours to weeks, far
% longer than dt = 10 min.

opts.init_state = [];
opts.verbose    = true;
for i = 1:2:numel(varargin), opts.(varargin{i}) = varargin{i+1}; end

c = P.const;  a = P.atm;  o = P.opt;  oc = P.ocean;  ice = P.ice;
sig = c.sigma_SB_W_m2_K4;
dt  = P.num.dt_s;
nt  = round(P.num.t_end_days * 86400 / dt);
T_f = oc.T_f_K;
nL  = a.n_layers;  lay = a.layers;
slab_ice = strcmpi(o.ice_model, 'slab');
N_ice = ice.N_layers;  if slab_ice, N_ice = 1; end
has_snow = ice.h_snow_m > 0 && ~slab_ice;

% emissivity of each slab: gases, plus the cloud in the layer that holds it
eps_slabs = [lay.emis_gas]';
if ~isempty(a.cloud_layer), [eps_c_now, eps_slabs(a.cloud_layer)] = cloud_slab_emissivity(P); else, eps_c_now = 0; end
gam_up = [lay.gamma_mix_above_W_m2_K]';          % coupling of layer k to k+1
C_lay  = [lay.C_J_m2_K]';
T_bg   = [lay.T_bg_K]';

% ------------------------------------------------------------------------
% Initial state
% ------------------------------------------------------------------------
if isempty(opts.init_state)
    S.T_atm_K = T_bg;
    if isfinite(ice.T_sk_prescribed_K), S.T_sk_K = ice.T_sk_prescribed_K;
    else,                                S.T_sk_K = ice.T_sk_init_K; end
    S.h_m = ice.h_init_m;
    if slab_ice
        S.T_ice_K = 0.5*(S.T_sk_K + T_f);         % mean of the linear profile
        S.T_sk_K  = S.T_ice_K;                    % the slab is the surface
    else
        % steady-conduction profile through snow (if any) and ice: the drop
        % across each layer is proportional to its thermal resistance
        R_snow = ice.h_snow_m / ice.k_snow_W_m_K;   R_ice = ice.h_init_m / ice.k_W_m_K;
        T_si = S.T_sk_K + (T_f - S.T_sk_K) * R_snow / (R_snow + R_ice);
        S.T_snow_K = 0.5*(S.T_sk_K + T_si);
        S.T_ice_K  = T_si + (T_f - T_si) * ice.centre_fractions(:);
    end
    S.T_ml_K = oc.T_ml_init_K;  S.T_deep_K = oc.T_deep_init_K;
    S.t_s = 0;
else
    S = opts.init_state;
    if numel(S.T_atm_K) ~= nL
        error('seb_column_run: init_state has %d atmosphere layers, P has %d', numel(S.T_atm_K), nL);
    end
    if numel(S.T_ice_K) ~= N_ice
        error('seb_column_run: init_state has %d ice layers, P has %d', numel(S.T_ice_K), N_ice);
    end
end

% ------------------------------------------------------------------------
% Storage
% ------------------------------------------------------------------------
z = @() zeros(1, nt);
R.t_days = z();  R.t_abs_days = z();
R.T_atm_K = zeros(nL, nt);  R.T_ice_K = zeros(N_ice, nt);
R.T_sk_K = z(); R.T_snow_K = z(); R.T_air_ref_K = z(); R.T_ml_K = z(); R.T_deep_K = z(); R.h_m = z();
R.LWD_sfc_W_m2 = z(); R.LWU_sfc_W_m2 = z(); R.OLR_W_m2 = z();
R.SH_down_W_m2 = z(); R.LH_down_W_m2 = z();
R.F_cond_top_up_W_m2 = z(); R.F_cond_bot_up_W_m2 = z(); R.F_ocean_up_W_m2 = z();
R.dhdt_bot_m_s = z(); R.melt_W_m2 = z();
R.Q_LW_atm_W_m2 = zeros(nL, nt); R.Q_relax_W_m2 = z(); R.Q_mix_W_m2 = z();
R.lambda_SH_W_m2_K = z(); R.lambda_LH_W_m2_K = z(); R.lambda_LW_W_m2_K = z();
R.C_H = z(); R.zeta = z(); R.Ri_b = z(); R.Gamma_lowest_K_m = z();
R.newton_iter = z(); R.E_ice_sens_J_m2 = z(); R.E_runoff_J_m2 = z();

% ------------------------------------------------------------------------
% Time loop
% ------------------------------------------------------------------------
for n = 1:nt
    % (1) emission of each slab (isothermal or profile; full or linearised Planck)
    E_dn = zeros(nL,1);  E_up = zeros(nL,1);
    for k = 1:nL
        [E_dn(k), E_up(k), prof] = slab_emission(S.T_atm_K(k), eps_slabs(k), lay(k), P);
        if k == 1, prof1 = prof; end
    end
    E_sfc = ice.emis_sfc * planck_emit(S.T_sk_K, P, ice.T_sk_init_K);
    LW1 = lw_grey_slabs(S.T_atm_K, eps_slabs, S.T_sk_K, ice.emis_sfc, sig, E_dn, E_up, E_sfc);

    % (2) air temperature at the reference height, from the lowest layer
    if strcmpi(o.atm_profile, 'lapse')
        T_air = prof1.T_base_K - prof1.Gamma_K_m * o.z_ref_m;
    else
        T_air = S.T_atm_K(1);
    end

    % (3) ocean heat flux into the ice base [McP92], or a fixed value
    if isfinite(oc.F_ocean_fixed_W_m2), F_ocean_up = oc.F_ocean_fixed_W_m2;
    else,                               F_ocean_up = oc.coeff_iceocean_W_m2_K * (S.T_ml_K - T_f); end

    % (4) the ice
    forcing.LWD_sfc_W_m2      = LW1.LWD_sfc;
    forcing.F_extra_down_W_m2 = a.F_extra_down_W_m2;
    forcing.T_air_K           = T_air;
    forcing.T_la_K            = T_air;           % older name, same quantity
    forcing.F_ocean_up_W_m2   = F_ocean_up;
    if slab_ice, [S, F] = ice_slab_step(S, forcing, P, dt);
    else,        [S, F] = ice_column_step(S, forcing, P, dt); end

    % (5) longwave again with the new surface emission -> layer heating, OLR
    E_sfc = F.LWU_sfc_W_m2 - (1 - ice.emis_sfc)*LW1.LWD_sfc;
    LW2 = lw_grey_slabs(S.T_atm_K, eps_slabs, S.T_sk_K, ice.emis_sfc, sig, E_dn, E_up, E_sfc);
    Q_LW = LW2.Q_slab;

    % (6) atmosphere layers (forward Euler) ---------------------------------
    mix_up = gam_up .* (S.T_atm_K - [S.T_atm_K(2:end); 0]);   % flux from k to k+1
    mix_up(end) = 0;
    mix_in = [0; mix_up(1:end-1)];                            % what k receives from k-1
    if strcmpi(o.atm_forcing, 'heating')
        relax = a.D_adv_W_m2 * a.D_weight_vec(:);        % prescribed convergence [Thorndike 1992]
    else
        relax = C_lay .* (T_bg - S.T_atm_K) / a.tau_adv_s;
    end
    turb   = zeros(nL,1);  turb(1) = F.SH_down_W_m2 + F.LH_down_W_m2;
    if strcmpi(a.mode, 'prognostic')
        S.T_atm_K = S.T_atm_K + dt ./ C_lay .* (Q_LW - mix_up + mix_in - turb + relax);
        Q_relax = sum(relax);
    else
        Q_relax = -(sum(Q_LW) - turb(1));                     % implied external source
    end

    %     ocean slabs -----------------------------------------------------
    F_deep_up = oc.gamma_deep_W_m2_K * (S.T_deep_K - S.T_ml_K);
    S.T_ml_K  = S.T_ml_K + dt/oc.C_ml_J_m2_K * (F_deep_up - F_ocean_up);
    if oc.deep_prognostic, S.T_deep_K = S.T_deep_K - dt/oc.C_deep_J_m2_K * F_deep_up; end
    S.t_s = S.t_s + dt;

    % store ------------------------------------------------------------------
    R.t_days(n) = n*dt/86400;  R.t_abs_days(n) = S.t_s/86400;
    R.T_atm_K(:,n) = S.T_atm_K;  R.T_ice_K(:,n) = S.T_ice_K;
    R.T_sk_K(n) = S.T_sk_K;  R.T_air_ref_K(n) = T_air;
    if has_snow, R.T_snow_K(n) = S.T_snow_K; else, R.T_snow_K(n) = NaN; end
    R.T_ml_K(n) = S.T_ml_K;  R.T_deep_K(n) = S.T_deep_K;  R.h_m(n) = S.h_m;
    R.LWD_sfc_W_m2(n) = LW2.LWD_sfc;  R.LWU_sfc_W_m2(n) = LW2.LWU_sfc;  R.OLR_W_m2(n) = LW2.OLR;
    R.SH_down_W_m2(n) = F.SH_down_W_m2;  R.LH_down_W_m2(n) = F.LH_down_W_m2;
    R.F_cond_top_up_W_m2(n) = F.F_cond_top_up_W_m2;  R.F_cond_bot_up_W_m2(n) = F.F_cond_bot_up_W_m2;
    R.F_ocean_up_W_m2(n) = F_ocean_up;
    R.dhdt_bot_m_s(n) = F.dhdt_bot_m_s;  R.melt_W_m2(n) = F.melt_W_m2;
    R.Q_LW_atm_W_m2(:,n) = Q_LW;  R.Q_relax_W_m2(n) = Q_relax;  R.Q_mix_W_m2(n) = mix_up(1);
    R.lambda_SH_W_m2_K(n) = F.lambda_SH_W_m2_K;  R.lambda_LH_W_m2_K(n) = F.lambda_LH_W_m2_K;
    R.lambda_LW_W_m2_K(n) = F.lambda_LW_W_m2_K;
    R.C_H(n) = F.C_H;  R.zeta(n) = F.zeta;  R.Ri_b(n) = F.Ri_b;  R.Gamma_lowest_K_m(n) = prof1.Gamma_K_m;
    R.newton_iter(n) = F.newton_iter;  R.E_ice_sens_J_m2(n) = F.E_ice_sens_J_m2;  R.E_runoff_J_m2(n) = F.E_runoff_J_m2;
end

R.P = P;  R.final_state = S;
R.eps_cloud = eps_c_now;  R.eps_cloud_slab = eps_slabs(max([a.cloud_layer 1]));  R.eps_slabs = eps_slabs;
R.L_eff_J_kg = F.L_eff_J_kg;
% level-3 names kept for older scripts and figures
R.T_la_K = R.T_atm_K(1,:);  R.T_ua_K = R.T_atm_K(end,:);
if ~isempty(a.cloud_layer), R.T_cl_K = R.T_atm_K(a.cloud_layer,:); else, R.T_cl_K = nan(1,nt); end
R.Q_LW_la_W_m2 = R.Q_LW_atm_W_m2(1,:);  R.Q_LW_ua_W_m2 = R.Q_LW_atm_W_m2(end,:);
if ~isempty(a.cloud_layer), R.Q_LW_cl_W_m2 = R.Q_LW_atm_W_m2(a.cloud_layer,:); else, R.Q_LW_cl_W_m2 = nan(1,nt); end

% ------------------------------------------------------------------------
% Energy-conservation diagnostic.  Column energy (relative to fixed references):
%   E = sum_k C_k T_k + E_ice,sens(rel. T_f) - rho L h + C_ml T_ml + C_deep T_deep
% Its change must equal the integral of the external exchanges: -OLR at the top,
% the relaxation heating, F_extra, a fixed deep reservoir's supply, and with
% ice_growth = 'fixed' the latent heat the fixed base would have released (ERA5).
% With use_brine the ice heat capacity depends on T, E_ice,sens is a linearisation, and the residual is an approximation check, not round-off.
% ------------------------------------------------------------------------
E = (C_lay' * R.T_atm_K) + R.E_ice_sens_J_m2 - ice.rho_kg_m3*R.L_eff_J_kg*R.h_m ...
    + oc.C_ml_J_m2_K*R.T_ml_K + oc.C_deep_J_m2_K*R.T_deep_K;
src = -R.OLR_W_m2 + R.Q_relax_W_m2 + a.F_extra_down_W_m2 - R.E_runoff_J_m2/dt;   % meltwater runoff removes its sensible heat (zero in winter)
if ~oc.deep_prognostic, src = src + oc.gamma_deep_W_m2_K*(R.T_deep_K - R.T_ml_K); end
if strcmpi(o.ice_growth, 'fixed'), src = src + (R.F_cond_bot_up_W_m2 - R.F_ocean_up_W_m2) - R.melt_W_m2; end
if ~isempty(opts.init_state), E0 = energy_of_state(opts.init_state, P); else, E0 = NaN; end
if isfinite(E0), dE_model = E - E0;        dE_flux = cumsum(src) * dt;
else,            dE_model = E(2:end) - E(1); dE_flux = cumsum(src(2:end)) * dt; end
R.energy.dE_model_J_m2 = dE_model;  R.energy.dE_flux_J_m2 = dE_flux;
R.energy.residual_mean_W_m2 = (dE_model(end) - dE_flux(end)) / (numel(dE_model)*dt);
R.energy.residual_max_W_m2  = max(abs(dE_model - dE_flux)) / (numel(dE_model)*dt);

if opts.verbose
    fprintf('seb_column_run: level %d | %d steps of %g s (%.1f d) | %d atm layers, ice %s (%d), ocean %s | LW %s, turb %s, profile %s | cloud %s %g g/m2 (eps %.3f)\n', ...
        o.complexity, nt, dt, P.num.t_end_days, nL, o.ice_model, N_ice, o.ocean_model, o.lw_emission, o.turb_scheme, o.atm_profile, ...
        P.cloud.phase, P.cloud.WP_g_m2, eps_c_now);
    fprintf('   final: T_sk = %.2f K, T_atm(1) = %.2f K, h = %.3f m, LWD = %.1f W/m2 | energy residual %.1e (mean) %.1e (max) W/m2 | Newton max %d\n', ...
        S.T_sk_K, S.T_atm_K(1), S.h_m, LW2.LWD_sfc, R.energy.residual_mean_W_m2, R.energy.residual_max_W_m2, max(R.newton_iter));
end
end

% ------------------------------------------------------------------------
function E = energy_of_state(S, P)
a = P.atm;  oc = P.ocean;  ice = P.ice;
C_lay = [a.layers.C_J_m2_K]';
[~, rc, L_eff] = ice_properties(S.T_ice_K, P);
if strcmpi(P.opt.ice_model, 'slab')
    E_ice = rc(1) * S.h_m * (S.T_ice_K(1) - oc.T_f_K);
elseif strcmpi(P.opt.ice_model, 'zero_layer')
    E_ice = 0;                                   % no heat capacity by construction
else
    E_ice = sum(rc .* (S.T_ice_K(:) - oc.T_f_K) .* (S.h_m * ice.layer_fractions(:)));
    if ice.h_snow_m > 0, E_ice = E_ice + ice.rho_c_snow_J_m3_K*ice.h_snow_m*(S.T_snow_K - oc.T_f_K); end
end
E = C_lay' * S.T_atm_K(:) + E_ice - ice.rho_kg_m3*L_eff*S.h_m + oc.C_ml_J_m2_K*S.T_ml_K + oc.C_deep_J_m2_K*S.T_deep_K;
end
