function P = seb_column_params(varargin)
% SEB_COLUMN_PARAMS  Default parameters for the zero-dimensional ocean / sea-ice /
%                    atmosphere / cloud column model (polar night, no solar).
%
%   P = seb_column_params()                       all defaults
%   P = seb_column_params('LWP_g_m2', 50, ...)    override any leaf field by name
%
% Every physical quantity carries its units in its name.  Every value carries a
% comment saying where it came from, and whether it is (a) a physical constant,
% (b) a measured / published value, (c) a TUNED value chosen here so the model
% reproduces an observed state, or (d) a free knob the user is meant to turn.
% Values marked TUNED are the ones to be sceptical of.
%
% The model structure (bottom to top):
%
%     deep ocean  ->  ocean mixed layer  ->  sea ice (N_ice layers, bottom at
%     the freezing point, growing by the Stefan condition)  ->  skin (zero heat
%     capacity, diagnostic)  ->  lower atmosphere (sub-cloud)  ->  cloud layer
%     ->  upper atmosphere  ->  space
%
% Longwave radiation is exchanged between every layer with a grey-slab model;
% sensible and latent heat are exchanged between the skin and the lower
% atmosphere with bulk formulae; heat is conducted through the ice with the heat
% equation; the ice base is held at the seawater freezing point and grows or
% melts according to the latent-heat balance there (Stefan 1891).
%
% SOURCES used throughout this directory (short tags used in the comments):
%   [S25]   Sledd, Shupe, Solomon & Cox (2025), JGR Atmos 130, e2024JD042578.
%           Surface energy balance responses to radiative forcing, MOSAiC.
%           - eq. (1) SEB with turbulent fluxes positive AWAY from the surface
%           - eq. (5) the "infinitely thin interface" surface energy balance
%           - eq. (7) the 1-D heat equation used for the snow/ice column
%           - eq. (10) effective depth h_eff = h_ice + h_snow*(k_ref/k_snow)
%           - Table 2: IFS/ERA5 uses constant 1.5 m ice with no snow
%           - sec. 2.1: skin temperature from LWU with emissivity 0.985
%   [M17]   Miller et al. (2017), The Cryosphere 11, 497-516. Summit SEB responses.
%   [PRISM] multliayerModel_ajb.ipynb (W. Camp / A. Buggee, PRISM project):
%           three coupled slabs C_i dT_i/dt = sum(fluxes) with bulk exchange
%           coefficients gamma [W m-2 K-1], atmospheric heat capacity written as
%           c_p*dp/g, equilibrium and timescale printouts, stiff ODE integration.
%           The slab bookkeeping here follows that notebook.
%   [V22]   Villanueva et al. (2022), Environ. Res. Lett. 17, 114057. Mixed-phase
%           regime cloud thinning. Table 2: LES glaciation of an Arctic deck at
%           -15 C: total water 0.170 -> 0.009 g/kg and LW CRE 57 -> 3 W m-2 at a
%           1 %/h freezing rate; mechanism = WBF growth -> precipitation.
%   [L26]   Li, Evan, Eisenman & Kok (2026, preprint), Limited climate response to
%           extreme Arctic MCT. Winter sea-ice LWP ~30 g m-2 vs ~114 g m-2 over
%           open ocean; surface LW forcing of complete cloud removal -14 W m-2
%           over ice, -39 over ocean; clouds over ice often WARMER than the
%           surface (embedded in surface-based inversions).
%   [Z23]   Zapalac (2023, preprint). eps_c = 1 - exp(-k_e LWP), k_e = 0.156 m2/g;
%           clear-sky effective emissivity ~0.7 at T_A = 235 K; ice growth upper
%           bound dh = dF*t/(rho*L) -> 0.85 cm per W m-2 per month.
%   [St78]  Stephens (1978), J. Atmos. Sci. 35, 2123-2132. Cloud LW emissivity
%           eps = 1 - exp(-a0 LWP), a0 = 0.130 (downward) / 0.158 (upward) m2 g-1.
%   [BL99]  Bitz & Lipscomb (1999), JGR 104, 15669-15677. Energy-conserving
%           thermodynamic sea-ice model: brine-pocket heat capacity, salinity-
%           dependent conductivity k = k0 + beta*S/T, conservative layer remap.
%   [MU71]  Maykut & Untersteiner (1971), JGR 76, 1550-1575. 1-D sea-ice thermo-
%           dynamics; ocean heat flux F_w = 2 W m-2.
%   [Sem76] Semtner (1976), J. Phys. Oceanogr. 6, 379-389. Three-layer ice model.
%   [W00]   Winton (2000), J. Atmos. Ocean. Tech. 17, 525-531. Two-layer ice with
%           a linearised surface balance (the skin solve here is the same idea).
%   [McP92] McPhee (1992), JGR 97, 5365-5379. Ocean-ice heat flux
%           F_w = rho_w c_w c_h u* (T_ml - T_f), c_h ~ 0.0056.
%   [Mil78] Millero (1978), UNESCO Tech. Pap. Mar. Sci. 28. Seawater freezing point.
%   [MK05]  Murphy & Koop (2005), QJRMS 131, 1539-1565. Vapour pressure over ice.
%   [A10]   Andreas et al. (2010), J. Hydrometeorol. 11, 87-104. Bulk transfer
%           coefficients over winter sea ice (SHEBA), C_H of order 1e-3.
%   [A02]   Andreas et al. (2002), JGR 107, 8033. Near-surface air over polar
%           sea ice is always close to ice saturation.
%   [CJ59]  Carslaw & Jaeger (1959), Conduction of Heat in Solids, 2nd ed.
%           Ch. XI: the Stefan / Neumann moving-boundary solutions.
%   [IFS]   ECMWF IFS Documentation Cy50r1 Part IV (2025): eqs. 3.76/8.25 skin balance, 3.84 LW linearisation,
%           8.8-8.9 bulk H and E, Table 8.3 roughness (ice caps), sec. 8.11 sea ice: eq. 8.189, k = 2.03, rho*c = 1.88e6, 1.5 m, 4 layers.
%   [P10]   Pierrehumbert (2010), Principles of Planetary Climate, sec. 4.3:
%           grey-gas slab atmospheres.
%   [LSM]   linearized_seb_model.ipynb (this project, Sep 2026): ERA5 Barrow
%           regressions giving rho*c_p*C_H*U = 9.9 W m-2 K-1 over pack ice and
%           lambda_LH = 2.9 W m-2 K-1; the frequency-dependent ice admittance.
%
% Andrew Buggee / Claude, Ocean Visions simple models, October 2026.

% =========================================================================
% Physical constants  (CODATA 2018 unless noted)
% =========================================================================
P.const.sigma_SB_W_m2_K4 = 5.670374419e-8;  % Stefan-Boltzmann constant
P.const.cp_air_J_kg_K    = 1004.7;          % specific heat of dry air at const. p [IFS Ch.12]
P.const.R_dry_J_kg_K     = 287.06;          % gas constant of dry air
P.const.g_m_s2           = 9.80665;         % standard gravity
P.const.L_sub_J_kg       = 2.834e6;         % latent heat of sublimation (ice -> vapour) [IFS]
P.const.L_vap_J_kg       = 2.501e6;         % latent heat of vaporisation (moist lapse rate)
P.const.L_fus_J_kg       = 3.34e5;          % latent heat of fusion, FRESH water [BL99 L0; Z23 uses 333.4 kJ/kg]
P.const.eps_molar        = 0.622;           % M_water / M_dry_air
P.const.T0_K             = 273.15;          % 0 degC

% =========================================================================
% MODEL OPTIONS: structure (complexity) and physics, all switchable
% =========================================================================
% The 26 Oct 2026 notes asked for two kinds of modularity: the NUMBER OF LAYERS
% and the PHYSICS used between them.  Both live here; run_simple_model.m is the
% one place to set them, model_config.m applies them.  Every option names its
% physics in the module that implements it.
%
%   complexity   1  atmosphere slab + ice slab + ocean mixed layer (the simplest
%                   cartoon: LW/SH/LH between air and ice, conduction ice-ocean)
%                2  boundary-layer slab (with the cloud) + free troposphere;
%                   skin + 2-layer ice; mixed layer + deep ocean
%                3  sub-cloud slab + cloud layer + free troposphere; skin +
%                   3-layer ice [0.1 0.8 0.1]; mixed layer + deep ocean (the
%                   model of 22 Oct 2026)
%                4  as 3, with a lapse-rate temperature profile inside each
%                   boundary-layer slab (atm_profile = 'lapse')
P.opt.complexity      = 3;
P.opt.lw_emission     = 'planck';      % 'planck' sigma T^4 | 'linearized' tangent at T_ref   [planck_emit.m]
P.opt.planck_T_ref_fallback_K = 250;   % reference when a caller supplies none
P.opt.turb_scheme     = 'bulk_const';  % 'none' | 'linear' | 'bulk_const' | 'most' [seb_surface_fluxes.m]
P.opt.latent_scheme   = 'bulk';        % 'bulk' | 'off'
P.opt.lambda_SH_lin_W_m2_K = 3.3;      % coefficient for turb_scheme = 'linear'
P.opt.lambda_LH_lin   = 'clausius';    % 'clausius' (from lambda_SH) or a number [W m-2 K-1]
P.opt.z_ref_m         = 10;            % reference height for 'most' (IFS lowest level ~10 m)
P.opt.z0m_m           = 1.3e-3;        % roughness for momentum, ice caps/glaciers [IFS Cy50r1 Table 8.3]
P.opt.z0h_m           = 1.3e-4;        % roughness for heat                        [IFS Cy50r1 Table 8.3]
P.opt.ice_model       = 'auto';        % 'auto' = the level's default (level 1: 'slab'; levels 2-4: 'skin_profile')
                                       % | 'slab' (one temperature) | 'skin_profile' (zero-capacity skin + N layers)
                                       % | 'zero_layer' (skin + linear profile with NO heat capacity: Semtner 1976
                                       %   "0-layer", the ice of Thorndike 1992)
P.opt.atm_forcing     = 'relaxation';  % how the advective heat source enters the atmosphere:
                                       % 'relaxation' Newtonian relaxation of each layer toward T_bg on tau_adv
                                       % 'heating'    a prescribed heat convergence D_adv_W_m2 shared among the
                                       %              layers (Thorndike 1992, sec. 2: "northward heat transport D
                                       %              distributed evenly throughout the polar atmosphere"); the
                                       %              layer temperatures then follow from radiative balance and
                                       %              no background temperatures are needed
P.opt.ice_growth      = 'stefan';      % 'stefan' | 'fixed' (constant thickness, as ERA5 [IFS 8.11])
P.opt.atm_profile     = 'auto';        % 'auto' = the level's default (level 4: 'lapse'; else 'isothermal')
                                       % | 'isothermal' | 'lapse' (profile inside each BL slab) [slab_emission.m]
P.opt.lapse_kind      = 'dry';         % 'dry' | 'moist' | a number in K m-1 (negative = inversion) [lapse_rate.m]; per layer: P.atm.layers(k).lapse ('auto' = this)
P.opt.ocean_model     = 'auto';        % 'auto' = the level's default (level 1: 'ml_fixed_deep'; else 'ml_deep')
                                       % | 'ml_deep' | 'ml_fixed_deep'
P.opt.draw_diagram    = true;

% =========================================================================
% Numerics
% =========================================================================
P.num.dt_s          = 600;     % time step [s]. Ice conduction is implicit; the
                               % atmosphere/ocean slabs are explicit and have
                               % relaxation times of hours or longer, so 10 min
                               % is far inside the stability limit (see README).
P.num.t_end_days    = 30;      % run length [days]
P.num.newton_tol_W_m2 = 1e-4;  % convergence of the skin Newton solve: |net flux into the skin| below this [W m-2]
P.num.newton_maxit  = 60;
P.num.T_sk_bounds_K = [150 300];  % bisection safeguard for the skin solve

% =========================================================================
% Sea ice
% =========================================================================
P.ice.N_layers        = 'auto'; % number of ice layers, or 'auto' = the level's default:
                                % 1 for ice_model 'slab'/'zero_layer', 2 at level 2,
                                % else 3 = the user's "ocean-ice interface / intermediate /
                                % ice-atmosphere interface" layers. Any N is accepted
                                % ([Sem76] 2 + snow; [IFS] 4; [S25] Table 2: 7-12).
P.ice.layer_fractions = 'auto'; % thickness of each layer as a fraction of h, TOP to
                                % BOTTOM (sum = 1); 'auto' = [0.10 0.80 0.10] for 3 layers,
                                % [0.5 0.5] for 2: thin INTERFACE layers and a thick interior
                                % (22 Oct 2026 notes); [] = equal layers. WHY IT MATTERS: the skin feels
                                % the ice through the top layer's half-thickness,
                                % K_top = k/(dz_1/2). The response of a conductor to a
                                % forcing of period P penetrates only sqrt(kappa P/pi)
                                % (3.5 cm at 1 h, 17 cm at 1 day for ice), so a layer
                                % cannot resolve responses faster than its own
                                % diffusion time dz^2/kappa (6 h for 15 cm; 2.7 days
                                % for the 50 cm of an equal 3-layer split of 1.5 m).
                                % The IFS puts its top ice layer at 7 cm for this
                                % reason [IFS eq. 8.150; LSM]. Steady states and the
                                % Stefan growth law are exact for ANY layering.
P.ice.h_init_m        = 1.5;    % initial thickness [m]. 1.5 m is the constant ice
                                % thickness of ERA5/IFS [IFS eq. 8.150; S25 Table 2],
                                % so the default column is "what ERA5 has".
% --- optional snow layer on top of the ice (OFF by default: h_snow_m = 0) --
% Snow is the first-order modulator of winter conduction: k_snow/k_ice ~ 0.15,
% so 20 cm of snow adds as much thermal resistance as ~1.3 m of ice [S25 eq. 10,
% "effective depth"]. ERA5 carries no snow on sea ice, which is the source of
% its 5-10 K warm skin bias (Batrak & Mueller 2019). Ian Eisenman's point at the
% OV meeting (22 Oct 2026) -- earlier freeze-up puts snow on ice instead of in
% the ocean, insulating the ice and slowing later growth -- is this term.
P.ice.h_snow_m        = 0.0;    % snow depth [m]; 0.15-0.25 m is typical mid-winter MOSAiC
P.ice.k_snow_W_m_K    = 0.30;   % snow conductivity [Sturm et al. 1997; S25 sec. 4.5, CAFS value]
P.ice.rho_snow_kg_m3  = 300;    % wind-packed Arctic snow
P.ice.c_snow_J_kg_K   = 2090;   % specific heat of snow (ice value at -5 C)
P.ice.rho_kg_m3       = 917;    % density of sea ice [Z23; BL99]
P.ice.c_J_kg_K        = 2106;   % specific heat of pure ice near -10 C [BL99 c0]
P.ice.k_W_m_K         = 2.03;   % thermal conductivity of pure ice [Untersteiner 1961; IFS eq. 8.149]
P.ice.emis_sfc        = 0.985;  % broadband LW emissivity of snow/ice [S25 sec. 2.1]
P.ice.T_melt_sfc_K    = 273.15; % surface melt temperature (snow / fresh ice). Not
                                % reached in polar night; kept so the regime
                                % switch of [S25] (growth vs melt) is in the code.
P.ice.T_sk_init_K     = 250.0;  % initial skin temperature; the ice starts with a
                                % linear (steady-conduction) profile from this to T_f
P.ice.T_sk_prescribed_K = NaN;  % if finite, the skin temperature is PRESCRIBED and
                                % the surface energy balance is bypassed. Used by the
                                % analytic Stefan-problem test; NaN = normal operation.
P.ice.h_min_m         = 0.01;   % thickness floor; the run stops if the ice vanishes
% --- optional brine corrections [BL99], OFF by default -------------------
P.ice.use_brine       = false;  % true: salinity-dependent k and c, effective L
P.ice.S_ice_ppt       = 4.0;    % bulk ice salinity [ppt] for the brine corrections
P.ice.beta_k_W_m_ppt  = 0.13;   % k_i = k0 + beta*S/T(degC)       [BL99 eq. (6), after Untersteiner 1961]
P.ice.mu_K_ppt        = 0.054;  % liquidus slope T_m = -mu*S      [BL99]

% =========================================================================
% Ocean: mixed layer directly under the ice, and a "deep" reservoir below it
% =========================================================================
P.ocean.S_psu            = 34.0;    % mixed-layer salinity -> freezing point via [Mil78]
P.ocean.rho_w_kg_m3      = 1025;    % seawater density
P.ocean.c_w_J_kg_K       = 3985;    % seawater specific heat
P.ocean.h_ml_m           = 30;      % winter mixed-layer depth [m] (20-50 m typical,
                                    % Peralta-Ferriz & Woodgate 2015, Prog. Oceanogr.)
P.ocean.h_deep_m         = 200;     % thickness of the reservoir slab below the ML
P.ocean.c_h_mcphee       = 0.0056;  % ice-ocean heat transfer coefficient [McP92]
P.ocean.ustar_m_s        = 0.005;   % friction velocity under the ice [m/s]; quiet
                                    % winter pack. rho_w*c_w*c_h*u* = 114 W m-2 K-1.
P.ocean.gamma_deep_W_m2_K = 4.0;    % ML <-> deep exchange coefficient. TUNED so the
                                    % steady ocean heat flux to the ice base is
                                    % ~2 W m-2, the classical value of [MU71].
                                    % (Slab-exchange form after [PRISM] gamma_ul.)
P.ocean.F_ocean_fixed_W_m2 = NaN;   % if finite, use this constant ocean heat flux
                                    % into the ice base instead of the bulk formula
                                    % (e.g. 2 W m-2 as in [MU71], or 0 for the
                                    % analytic Stefan test). NaN = bulk formula.
P.ocean.dT_ml_init_K     = 0.02;    % initial T_ml - T_f  [K]
P.ocean.dT_deep_init_K   = 0.50;    % initial T_deep - T_f [K]
%  (deep reservoir prognostic or fixed: chosen by P.opt.ocean_model 'ml_deep' | 'ml_fixed_deep'; finalize derives P.ocean.deep_prognostic)

% =========================================================================
% Atmosphere: three slabs.  Heat capacity C = c_p * dp / g, after [PRISM].
% =========================================================================
P.atm.mode            = 'prognostic'; % 'prognostic' (slab temperatures evolve) or
                                      % 'prescribed' (held at their initial values;
                                      % the "fixed atmosphere" limit, alpha = 0 in
                                      % the notation of [LSM]).
P.atm.dp_la_hPa       = 50;     % sub-cloud layer, ~400 m
P.atm.dp_cl_hPa       = 40;     % cloud layer, ~300 m
P.atm.dp_ua_hPa       = 700;    % the rest of the troposphere above the cloud
P.atm.p_sfc_Pa        = 1.013e5;
% --- grey-slab longwave emissivities (clear-sky gases) -------------------
P.atm.emis_la         = 0.40;   % TUNED (see T_bg below). Together with emis_gas_cl,
P.atm.emis_gas_cl     = 0.15;   % emis_ua and the backgrounds this gives clear-sky LWD
P.atm.emis_ua         = 0.63;   % ~165 W m-2, the clear cluster of the winter LWD
                                % distribution over pack ice ([S25 Fig. 1-2]: ~150-170
                                % at MOSAiC; ERA5 Barrow pack ice a few W m-2 higher),
                                % consistent with the ~0.7 whole-column clear-sky
                                % emissivity of [Z23].  emis_gas_cl is the clear-air
                                % (water vapour + CO2) emissivity of the 300 m cloud
                                % layer itself: when the cloud glaciates or dissipates
                                % that air still radiates. It combines with the cloud as
                                % independent absorbers, transmissions multiplying:
                                %     eps_slab = 1 - (1 - emis_gas_cl)(1 - eps_cloud)
                                % so a black cloud is unaffected (0.98 -> 0.983) and a
                                % vanished cloud leaves 0.15. emis_ua was lowered from
                                % 0.75 to 0.63 to keep the clear-sky LWD unchanged.
% --- background (advective reservoir) temperatures and relaxation ----------
% In polar night the only energy SOURCE for the atmosphere above sea ice is
% horizontal advection; the column model represents it as Newtonian relaxation
% of each slab toward a background temperature on a timescale tau_adv. This is
% THE unavoidable parameterization of the model. tau_adv -> inf recovers the
% "air follows the surface" regime (alpha ~ 1, what ERA5 regressions show
% [LSM]); tau_adv -> 0 recovers the fixed atmosphere (alpha = 0).
P.atm.T_bg_la_K       = 260.0;  % TUNED. One boundary-layer background shared by the
P.atm.T_bg_cl_K       = 260.0;  % sub-cloud and cloud slabs (the cloud sits in the BL
P.atm.T_bg_ua_K       = 232.0;  % air mass) and a colder free troposphere. Chosen by a
                                % bounded grid search (22 Oct 2026) so the ERA5-like
                                % default column -- bare 1.5 m ice, no snow -- sits at
                                % the two observed winter states of the Arctic surface
                                % (Stramler et al. 2011; [S25] Fig. 1 clusters; ERA5
                                % Barrow pack-ice means [LSM]):
                                %     clear sky   : T_sk 246.7 K, LWD 163 W m-2,
                                %                   T_la - T_sk = +4.5 K, SH_down +16 W m-2
                                %     LWP 30 g/m2 : T_sk 254.3 K, LWD 222 W m-2
                                % (targets 245/165 and 255/230; test T4 prints them;
                                % values after emis_gas_cl was added). The backgrounds sit
                                % ~5-8 K above the slab temperatures they produce: that
                                % offset is the radiative cooling the advection has to
                                % replace each day (~20 W m-2 over the BL at tau = 1 d).
                                % A warmer/colder regime: shift all three together.
P.atm.tau_adv_s       = 86400;  % 1 day. L/U ~ 500 km / 5 m/s.  FREE KNOB.
P.atm.D_adv_W_m2      = 100;    % advective heat convergence for atm_forcing = 'heating' [W m-2].
                                % Thorndike (1992) reads ~100 W m-2 for the polar cap from
                                % Nakamura & Oort (1988, JGR 93, 9510-9524, Fig. 7a).
P.atm.D_weights       = 'mass'; % how D is shared among layers: 'mass' (by dp), 'equal', or a vector
P.atm.gamma_mix_W_m2_K = 10.0;  % cloud <-> sub-cloud turbulent coupling. Cloud-top
                                % LW cooling drives mixing that couples the cloud to
                                % the layer below (Morrison et al. 2012, Nat. Geosci.
                                % 5, 11-17; Shupe et al. 2013, ACP 13, 9379). FREE KNOB.
% --- bulk turbulent exchange with the skin --------------------------------
P.atm.U_m_s           = 5.0;    % near-surface wind speed [m/s]
P.atm.C_H             = 0.5e-3; % EFFECTIVE bulk transfer coefficient for heat between
                                % the skin and the ~400 m sub-cloud SLAB. With U = 5
                                % and rho = 1.3 kg/m3: rho*c_p*C_H*U = 3.3 W m-2 K-1.
                                % This is deliberately ~3x smaller than the 2-m-based
                                % coefficient over winter sea ice (order 1e-3 [A10];
                                % 9.9 W m-2 K-1 from ERA5 regressions over Barrow
                                % pack ice [LSM]), because in a stable boundary layer
                                % the 2-m air follows the skin almost one-for-one
                                % ([LSM]: alpha = dT_2m/dT_skin = 1.00) and most of
                                % the skin-to-slab temperature difference sits in the
                                % surface layer. Magnitude chosen (TUNED) so a clear-
                                % sky inversion of 6-8 K gives SH_down = +15-25 W m-2,
                                % the SHEBA/MOSAiC winter range (Persson et al. 2002;
                                % [S25 Fig. 1b]). FREE KNOB -- the SH share of the
                                % response scales with it; the driver runs a 3x case.
P.atm.C_E             = 0.5e-3; % bulk coefficient for moisture (= C_H)
P.atm.F_extra_down_W_m2 = 0.0;  % an extra, prescribed flux absorbed by the skin
                                % [W m-2]. The clean "perturbation and response"
                                % experiment: +10 means 10 W m-2 more into the skin,
                                % with nothing else changed. Counted as an external
                                % energy source in the conservation diagnostic.
P.atm.RH_ice_la       = 0.95;   % relative humidity w.r.t. ice of the sub-cloud air [A02]

% =========================================================================
% Cloud
% =========================================================================
P.cloud.phase         = 'liquid';  % 'liquid' | 'ice' | 'none'
P.cloud.WP_g_m2       = 30;        % condensed water path [g m-2]. 30 g m-2 is the
                                   % winter mean LWP over Arctic sea ice in CESM2
                                   % [L26] and the value above which [Z23] calls
                                   % the cloud black (eps > 0.99).
P.cloud.r_e_liq_um    = 10;        % droplet effective radius [um]; the [St78] anchor
P.cloud.r_e_ice_um    = 30;        % ice effective radius [um] after glaciation
P.cloud.k_abs_liq_ref_m2_g = 0.130;% LW mass absorption coeff. of liquid at r_e = 10 um,
                                   % DOWNWARD emission [St78] (0.158 upward; [Z23] 0.156)
P.cloud.Q_abs_ice     = 1.0;       % absorption efficiency in the large-particle limit
P.cloud.diffusivity   = 1.66;      % Elsasser (1942) diffusivity factor, beam -> flux
P.cloud.rho_ice_kg_m3 = 917;
P.cloud.rho_liq_kg_m3 = 1000;

% =========================================================================
% Derived quantities (filled by finalize)
% =========================================================================
P = apply_overrides(P, varargin{:});
P = seb_column_finalize(P);
end

% -------------------------------------------------------------------------
function P = apply_overrides(P, varargin)
% Set leaf fields by name, e.g. ('WP_g_m2', 50) or ('complexity', 1). The name is
% searched in every sub-struct (const, num, ice, ocean, atm, cloud, opt) so the
% caller does not need to remember which group it lives in.  Layer-by-layer
% atmosphere settings are edited on P.atm.layers after construction, followed by
% P = seb_column_finalize(P).
if mod(numel(varargin), 2) ~= 0
    error('seb_column_params: overrides must be name/value pairs');
end
groups = fieldnames(P);
for i = 1:2:numel(varargin)
    name = varargin{i};  value = varargin{i+1};
    found = false;
    for g = 1:numel(groups)
        if isfield(P.(groups{g}), name)
            P.(groups{g}).(name) = value;
            found = true;
        end
    end
    if ~found
        error('seb_column_params: unknown parameter "%s"', name);
    end
end
end
