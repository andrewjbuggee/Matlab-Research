function P = seb_column_finalize(P)
% SEB_COLUMN_FINALIZE  Derived quantities, and the LAYER STRUCTURE for the chosen
%                      complexity level.  Called by seb_column_params; call it
%                      again after editing P by hand (e.g. P.atm.layers(2).emis_gas).
%
% THE ATMOSPHERE AS A LIST OF LAYERS.  P.atm.layers is a struct array, BOTTOM to
% TOP; the integrator loops over it, so one code path serves every level.  Each
% layer carries:
%   name        label for the diagram
%   dp_hPa      pressure thickness -> heat capacity C = c_p dp/g [PRISM], and the
%               geometric depth dz by the hypsometric equation
%   emis_gas    clear-air (water vapour + CO2) broadband emissivity of the layer
%   has_cloud   true for the one layer that holds the cloud; its emissivity is
%               1 - (1 - emis_gas)(1 - eps_cloud), independent absorbers
%   T_bg_K      advective background temperature it relaxes toward
%   gamma_mix_above_W_m2_K   turbulent coupling to the layer above (0 = none)
%   lapse       lapse-rate kind used when atm_profile = 'lapse' ('auto' = P.opt.lapse_kind)
%
% THE LEVELS.  Level 3 is the 22 Oct 2026 model and its numbers come from the
% scalar fields of seb_column_params (so overrides such as 'T_bg_la_K' still
% work).  Levels 1 and 2 merge those layers; their gas emissivities are set so
% that the clear-sky LWD at the surface stays near 165 W m-2 (the tuned target):
%   level 2  BL emissivity 0.49 = 1 - (1 - 0.40)(1 - 0.15), the sub-cloud and
%            cloud-layer gases of level 3 combined;
%   level 1  one slab of 790 hPa with emissivity 0.75 relaxing to 252 K: at its
%            equilibrium (~250 K) it emits 0.75 sigma 250^4 = 166 W m-2 toward
%            the surface in clear sky and sigma 250^4 = 221 with a black cloud.
% Test T4 prints the equilibrium states for every level.

c = P.const;  o = P.opt;

% ------------------------------------------------------------------------
% Freezing point of seawater [Mil78] and ocean/ice derived quantities
% ------------------------------------------------------------------------
S = P.ocean.S_psu;
P.ocean.T_f_K = c.T0_K + (-0.0575*S + 1.710523e-3*S^1.5 - 2.154996e-4*S^2);
P.ocean.C_ml_J_m2_K   = P.ocean.rho_w_kg_m3 * P.ocean.c_w_J_kg_K * P.ocean.h_ml_m;
P.ocean.C_deep_J_m2_K = P.ocean.rho_w_kg_m3 * P.ocean.c_w_J_kg_K * P.ocean.h_deep_m;
P.ocean.coeff_iceocean_W_m2_K = P.ocean.rho_w_kg_m3 * P.ocean.c_w_J_kg_K ...
                                * P.ocean.c_h_mcphee * P.ocean.ustar_m_s;      % [McP92]
P.ocean.T_ml_init_K   = P.ocean.T_f_K + P.ocean.dT_ml_init_K;
P.ocean.T_deep_init_K = P.ocean.T_f_K + P.ocean.dT_deep_init_K;
%  (P.ocean.deep_prognostic is set below, once the 'auto' options are resolved)

P.ice.rho_c_J_m3_K      = P.ice.rho_kg_m3 * P.ice.c_J_kg_K;             % 1.93e6 ([IFS]: 1.88e6)
P.ice.rho_c_snow_J_m3_K = P.ice.rho_snow_kg_m3 * P.ice.c_snow_J_kg_K;   % 6.3e5
P.ice.kappa_m2_s        = P.ice.k_W_m_K / P.ice.rho_c_J_m3_K;           % ~1.05e-6 m2/s
P.ice.h_eff_init_m      = P.ice.h_init_m + P.ice.h_snow_m * P.ice.k_W_m_K / P.ice.k_snow_W_m_K;  % [S25 eq. 10]

% ------------------------------------------------------------------------
% Structure by complexity level (only if no layer list exists yet)
% ------------------------------------------------------------------------
% 'auto' options take the level's default; anything set explicitly is kept at
% every level (so e.g. 'zero_layer' ice or a lapse-rate profile can be asked for
% at level 1).
auto = @(v) strcmpi(v, 'auto');
a = P.atm;
if ~isfield(a, 'layers') || isempty(a.layers)
    switch o.complexity
        case 1
            layers = mk_layer('atmosphere', 790, 0.75, 252.0, true, 0, 'auto');
            if auto(o.ice_model),   o.ice_model   = 'slab'; end
            if auto(o.ocean_model), o.ocean_model = 'ml_fixed_deep'; end
            if auto(o.atm_profile), o.atm_profile = 'isothermal'; end
            %  ice layers: resolved below ('slab'/'zero_layer' -> 1; 'auto' at level 1 -> 3 [0.1 0.8 0.1])
        case 2
            layers = [mk_layer('boundary layer + cloud', 90, 1 - (1 - a.emis_la)*(1 - a.emis_gas_cl), a.T_bg_la_K, true, 0, 'auto'), ...
                      mk_layer('free troposphere', a.dp_ua_hPa, a.emis_ua, a.T_bg_ua_K, false, 0, 'auto')];
            if auto(o.ice_model),   o.ice_model   = 'skin_profile'; end
            if auto(o.ocean_model), o.ocean_model = 'ml_deep'; end
            if auto(o.atm_profile), o.atm_profile = 'isothermal'; end
            %  ice layers: 'auto' -> 2 equal layers (resolved below); an explicit N_layers is kept
        case {3, 4}
            layers = [mk_layer('sub-cloud air', a.dp_la_hPa, a.emis_la, a.T_bg_la_K, false, a.gamma_mix_W_m2_K, 'auto'), ...
                      mk_layer('cloud layer',   a.dp_cl_hPa, a.emis_gas_cl, a.T_bg_cl_K, true, 0, 'auto'), ...
                      mk_layer('free troposphere', a.dp_ua_hPa, a.emis_ua, a.T_bg_ua_K, false, 0, 'auto')];
            if auto(o.ice_model),   o.ice_model   = 'skin_profile'; end
            if auto(o.ocean_model), o.ocean_model = 'ml_deep'; end
            if auto(o.atm_profile), if o.complexity == 4, o.atm_profile = 'lapse'; else, o.atm_profile = 'isothermal'; end; end
        otherwise
            error('seb_column_finalize: complexity must be 1, 2, 3 or 4');
    end
    a.layers = layers;
end
% a caller may edit P and re-finalize: resolve any 'auto' still standing
if auto(o.ice_model),   o.ice_model   = 'skin_profile'; end
if auto(o.ocean_model), o.ocean_model = 'ml_deep'; end
if auto(o.atm_profile), o.atm_profile = 'isothermal'; end
if any(strcmpi(o.ice_model, {'zero_layer', 'slab'}))     % one layer: a slab, or a diagnostic linear profile
    P.ice.N_layers = 1;  P.ice.layer_fractions = 1;
elseif auto(P.ice.N_layers)                              % the level's default layering (see seb_column_params)
    if o.complexity == 2, P.ice.N_layers = 2;  fr_default = [0.5 0.5];
    else,                 P.ice.N_layers = 3;  fr_default = [0.10 0.80 0.10]; end
    if auto(P.ice.layer_fractions), P.ice.layer_fractions = fr_default; end
end
P.opt = o;
P.ocean.deep_prognostic = strcmpi(o.ocean_model, 'ml_deep');   % AFTER the 'auto' resolution (bug fixed in the 26 Oct 2026 review:
                                                               % evaluated before it, a fresh config had a fixed reservoir and a
                                                               % re-finalized copy a prognostic one, so control and perturbed runs differed)

% ------------------------------------------------------------------------
% Per-layer derived quantities: heat capacity, pressure bounds, depth
% ------------------------------------------------------------------------
% (rebuilt into a fresh struct array: every element must carry the same fields)
p_bot = a.p_sfc_Pa;
for k = 1:numel(a.layers)
    L = a.layers(k);
    L.C_J_m2_K = c.cp_air_J_kg_K * L.dp_hPa*100 / c.g_m_s2;          % c_p dp / g [PRISM]
    L.p_bot_Pa = p_bot;  L.p_top_Pa = p_bot - L.dp_hPa*100;
    L.p_mid_Pa = 0.5*(L.p_bot_Pa + L.p_top_Pa);
    % hypsometric depth at the background temperature (Wallace & Hobbs 2006, eq. 3.29)
    L.dz_m = c.R_dry_J_kg_K * L.T_bg_K / c.g_m_s2 * log(L.p_bot_Pa / L.p_top_Pa);
    if ~isfield(L, 'lapse') || isempty(L.lapse), L.lapse = 'auto'; end     % 'auto' = P.opt.lapse_kind, read at run time by slab_emission
    if k == 1, layers_out = L; else, layers_out(k) = L; end %#ok<AGROW>
    p_bot = L.p_top_Pa;
end
a.layers = layers_out;
a.n_layers = numel(a.layers);
a.cloud_layer = find([a.layers.has_cloud], 1);
% share of the advective heating D received by each layer (atm_forcing = 'heating')
if ischar(a.D_weights) || isstring(a.D_weights)
    switch lower(a.D_weights)
        case 'mass',  w = [a.layers.dp_hPa];
        case 'equal', w = ones(1, a.n_layers);
        otherwise, error('seb_column_finalize: D_weights must be ''mass'', ''equal'' or a vector');
    end
else
    w = a.D_weights(:)';
    if numel(w) ~= a.n_layers, error('seb_column_finalize: D_weights needs one entry per layer'); end
end
a.D_weight_vec = w / sum(w);
% convenience copies for the level-3 names used by older scripts
a.C_la_J_m2_K = a.layers(1).C_J_m2_K;
a.C_ua_J_m2_K = a.layers(end).C_J_m2_K;
if a.n_layers >= 2, a.C_cl_J_m2_K = a.layers(min(2, a.n_layers)).C_J_m2_K; else, a.C_cl_J_m2_K = NaN; end
P.atm = a;

% ------------------------------------------------------------------------
% Ice layer fractions
% ------------------------------------------------------------------------
if ischar(P.ice.layer_fractions) || isstring(P.ice.layer_fractions), P.ice.layer_fractions = []; end  % 'auto' with an explicit N: equal layers
if isempty(P.ice.layer_fractions) || numel(P.ice.layer_fractions) ~= P.ice.N_layers
    if ~isempty(P.ice.layer_fractions)
        warning('seb_column_finalize:fractions', 'layer_fractions has %d entries but N_layers = %d; using equal layers', ...
            numel(P.ice.layer_fractions), P.ice.N_layers);
    end
    P.ice.layer_fractions = ones(1, P.ice.N_layers) / P.ice.N_layers;
end
P.ice.layer_fractions  = P.ice.layer_fractions(:)' / sum(P.ice.layer_fractions);
P.ice.cum_fractions    = [0, cumsum(P.ice.layer_fractions)];
P.ice.centre_fractions = P.ice.cum_fractions(1:end-1) + P.ice.layer_fractions/2;

% ------------------------------------------------------------------------
% Cloud emissivity and the emissivity of the slab that holds the cloud
% ------------------------------------------------------------------------
[P.cloud.emis, P.cloud.k_abs_m2_g] = cloud_lw_emissivity(P.cloud.WP_g_m2, P.cloud.phase, P);
[~, P.cloud.emis_slab] = cloud_slab_emissivity(P);     % the integrator recomputes this at run time too

% ------------------------------------------------------------------------
% Consistency warnings (26 Oct 2026 review)
% ------------------------------------------------------------------------
if strcmpi(o.ice_model, 'slab') && P.ice.h_snow_m > 0
    warning('seb_column_finalize:snow', ...
        'ice_model = ''slab'' has no skin and ignores the snow layer (h_snow_m = %.2f m); use ''skin_profile'' or ''zero_layer'' for snow', P.ice.h_snow_m);
end
end

% ------------------------------------------------------------------------
function L = mk_layer(name, dp_hPa, emis_gas, T_bg_K, has_cloud, gamma_mix_above, lapse)
L = struct('name', name, 'dp_hPa', dp_hPa, 'emis_gas', emis_gas, 'T_bg_K', T_bg_K, ...
           'has_cloud', has_cloud, 'gamma_mix_above_W_m2_K', gamma_mix_above, 'lapse', lapse);
end
