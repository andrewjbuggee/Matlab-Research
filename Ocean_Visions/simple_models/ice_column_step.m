function [S, F] = ice_column_step(S, forcing, P, dt_s)
% ICE_COLUMN_STEP  Advance the skin + (snow) + ice column by one time step.
%
%   [S, F] = ice_column_step(S, forcing, P, dt_s)
%
%   S        state:  S.T_sk_K      skin (surface) temperature              [K]
%                    S.T_snow_K    snow-layer temperature (if h_snow_m > 0) [K]
%                    S.T_ice_K     N x 1 ice layer temperatures, TOP to BOTTOM [K]
%                    S.h_m         ice thickness                           [m]
%   forcing  forcing.LWD_sfc_W_m2      downwelling longwave at the surface
%            forcing.F_extra_down_W_m2 extra prescribed flux absorbed by the skin
%            forcing.T_la_K            sub-cloud air temperature (for SH, LH)
%            forcing.F_ocean_up_W_m2   ocean heat flux into the ice base (positive up)
%   P        parameters (seb_column_params); P.ice.T_sk_prescribed_K finite ->
%            skip the surface balance and impose that skin temperature
%   dt_s     time step [s]
%
%   F        fluxes and diagnostics for this step (all W m-2 unless noted):
%            LWU_sfc, SH_down, LH_down, F_cond_top_up, F_cond_bot_up, melt,
%            dhdt_bot_m_s, dhdt_top_m_s, lambda_SH, lambda_LH, lambda_LW,
%            newton_iter, skin_residual, E_ice_sens_J_m2 (snow + ice, rel. to T_f)
%
% =========================================================================
% WHAT THIS SOLVES, AND WHERE EACH PIECE COMES FROM
% =========================================================================
%
% 1. THE SKIN (surface) ENERGY BALANCE -- an interface with no heat capacity
%    [S25 eq. (5), "an infinitely thin layer"; IFS eq. 8.22; W00]:
%
%        eps_s*(LWD - sigma*T_sk^4) + F_extra + SH_down + LH_down + F_cond_top_up = 0
%
%    Every term is a flux INTO the skin.  F_cond_top_up is the conductive flux
%    from the layer just below (snow if present, else the first ice layer) up
%    into the skin; in winter the interior is warmer than the surface, so it is
%    positive and is what keeps the surface from cooling faster.  The equation
%    is nonlinear in T_sk (sigma*T^4, the saturation humidity, the implicit
%    conduction) and is solved by a safeguarded Newton iteration.  Because every
%    term decreases monotonically with T_sk, the root is unique.
%
% 2. HEAT CONDUCTION IN SNOW AND ICE -- Fourier's law / the heat equation
%    [S25 eq. (7) without the solar term; MU71; Sem76; BL99]:
%
%        rho*c * dT/dt = d/dz ( k dT/dz )
%
%    discretised on one snow cell (fixed depth, optional) plus N ice layers of
%    equal thickness h/N, advanced with backward Euler (fully implicit) so the
%    step is stable for any dt and any layer thickness.  Boundary conditions:
%    T = T_sk at the top (from step 1) and T = T_f at the bottom -- the ice-
%    ocean interface IS at the freezing point by definition; that is a physical
%    statement, not a parameterization.  Snow matters because k_snow/k_ice ~ 0.15:
%    20 cm of snow adds the thermal resistance of ~1.3 m of ice [S25 eq. 10].
%
% 3. ICE GROWTH -- THE STEFAN CONDITION [Stefan 1891; CJ59 ch. XI; MU71 eq. (5)]
%
%        rho_i * L * dh/dt = F_cond_bot_up - F_ocean_up
%
%    This is the answer to "how do I account for phase change": energy
%    conservation at the moving interface.  Freezing releases latent heat L per
%    kg; the interface sits at T_f; so the rate at which seawater freezes onto
%    the base is set by how fast conduction carries heat UP into the ice
%    (F_cond_bot_up = k dT/dz at the base) minus how fast the ocean supplies heat
%    to the interface from below (F_ocean_up).  No parameterization is involved
%    beyond the ocean heat flux itself and the choice of L (fresh water,
%    3.34e5 J/kg, or the brine-adjusted value of [BL99] if use_brine).  Surface
%    melt, when T_sk would exceed the melting point, is the same condition at
%    the top (Sledd's "M", the melt regime of [S25]) -- never reached in winter.
%
% 4. REGRIDDING.  After the thickness changes, the N ice layers are re-cut to
%    equal thickness and the temperature profile is remapped CONSERVING
%    ENTHALPY (thickness x rho*c x T), after [BL99 sec. 2e] and [Sem76].  New
%    ice forms at T_f.  The snow cell is not regridded (fixed depth, no
%    accumulation: this model has no precipitation).
%
% Sign conventions (say them out loud):  SH_down, LH_down, F_extra positive INTO
% the surface; F_cond_top_up, F_cond_bot_up, F_ocean_up positive UPWARD; dh/dt
% positive for growth.  [S25] eq. (1) has turbulent fluxes positive AWAY from
% the surface, so SH here = -SH there.

c  = P.const;  ice = P.ice;
N  = ice.N_layers;
T_f = P.ocean.T_f_K;
sig = c.sigma_SB_W_m2_K4;
h_old  = S.h_m;
T_old  = S.T_ice_K(:);
frac   = ice.layer_fractions(:);               % N x 1, sums to 1
dz_old = h_old * frac;                          % N x 1 layer thicknesses [m]
has_snow = ice.h_snow_m > 0;
if has_snow && ~isfield(S, 'T_snow_K'), S.T_snow_K = S.T_sk_K; end

% ------------------------------------------------------------------------
% Layer thermal properties (pure ice by default; brine corrections of [BL99])
% ------------------------------------------------------------------------
[k_l, rho_c_l, L_eff] = ice_properties(T_old, P);
% Conductances [W m-2 K-1] between adjacent layer CENTRES (harmonic mean so a
% poorly conducting layer controls the flux) and between the outer centres and
% the boundary surfaces a half-layer away.
% centre-to-centre: two half-layers in series, each with its own k
K_int  = 1 ./ ( (dz_old(1:end-1)/2)./k_l(1:end-1) + (dz_old(2:end)/2)./k_l(2:end) );  % (N-1) x 1
K_bot  = k_l(end) / (dz_old(end)/2);            % bottom ice centre -> base (T_f)
cap_ice = rho_c_l .* dz_old / dt_s;             % N x 1
if has_snow
    hs = ice.h_snow_m;  ks = ice.k_snow_W_m_K;
    K_top   = ks / (hs/2);                       % skin -> snow centre
    K_si    = 1 / ( (hs/2)/ks + (dz_old(1)/2)/k_l(1) );   % snow centre -> ice-1 centre
    cap_snow = ice.rho_c_snow_J_m3_K * hs / dt_s;
    T_snow_old = S.T_snow_K;
else
    K_top = k_l(1) / (dz_old(1)/2);              % skin -> ice-1 centre
end

% ------------------------------------------------------------------------
% 1 + 2.  Solve the skin balance with the implicit conduction underneath
% ------------------------------------------------------------------------
if isfinite(ice.T_sk_prescribed_K)
    % --- prescribed skin temperature: no surface balance (analytic tests) --
    T_sk = ice.T_sk_prescribed_K;
    [T_new, T_snow_new, F_top_up, F_bot_up] = conduction_solve(T_sk);
    SF = seb_surface_fluxes(T_sk, forcing.T_la_K, P);
    n_iter = 0;  resid = NaN;  melt = 0;
else
    % --- safeguarded Newton on the scalar residual R(T_sk) ----------------
    lo = P.num.T_sk_bounds_K(1);  hi = P.num.T_sk_bounds_K(2);
    T_sk = min(max(S.T_sk_K, lo + 1), hi - 1);      % start from last step
    n_iter = 0;  resid = skin_residual(T_sk);
    dT_fd = 0.01;                                   % finite-difference step [K]
    while abs(resid) > P.num.newton_tol_K * 10 && n_iter < P.num.newton_maxit
        % R is monotone decreasing in T_sk: a positive residual means the skin
        % is too cold, so the root lies above; shrink the bracket accordingly.
        if resid > 0, lo = T_sk; else, hi = T_sk; end
        dR = (skin_residual(T_sk + dT_fd) - skin_residual(T_sk - dT_fd)) / (2*dT_fd);
        T_try = T_sk - resid / dR;                     % Newton step
        if ~(T_try > lo && T_try < hi)                 % left the bracket -> bisect
            T_try = 0.5*(lo + hi);
        end
        T_sk  = T_try;
        resid = skin_residual(T_sk);
        n_iter = n_iter + 1;
    end
    % --- melt regime [S25]: cap the skin at the melting point ---------------
    if T_sk > ice.T_melt_sfc_K
        T_sk  = ice.T_melt_sfc_K;
        resid = skin_residual(T_sk);        % > 0: surplus energy melts the surface
        melt  = resid;                       % Sledd's M [W m-2]
    else
        melt  = 0;
    end
    [T_new, T_snow_new, F_top_up, F_bot_up] = conduction_solve(T_sk);
    SF = seb_surface_fluxes(T_sk, forcing.T_la_K, P);
end

% ------------------------------------------------------------------------
% 3.  Stefan condition at the base (and surface melt at the top)
% ------------------------------------------------------------------------
dhdt_bot = (F_bot_up - forcing.F_ocean_up_W_m2) / (ice.rho_kg_m3 * L_eff);  % m/s, >0 grows
dhdt_top = -melt / (ice.rho_kg_m3 * L_eff);                                 % m/s, <=0
dh_bot = dhdt_bot * dt_s;
dh_top = dhdt_top * dt_s;
if dh_bot < -dz_old(end) || dh_top < -dz_old(1)
    error('ice_column_step:bigmelt', ...
        'melt in one step (%.3g m) exceeds the boundary layer thickness (%.3g m); reduce dt_s', ...
        max(-dh_bot, -dh_top), min(dz_old([1 end])));
end
if dh_bot > dz_old(end)
    warning('ice_column_step:biggrowth', ...
        'growth in one step (%.3g m) exceeds the bottom layer (%.3g m); reduce dt_s for accuracy', ...
        dh_bot, dz_old(end));
end
h_new = h_old + dh_bot + dh_top;
if h_new < ice.h_min_m
    error('ice_column_step:vanished', 'ice thickness fell below h_min (%.3f m)', ice.h_min_m);
end

% ------------------------------------------------------------------------
% 4.  Enthalpy-conserving remap of the ICE onto N equal layers [BL99 sec. 2e]
% ------------------------------------------------------------------------
% Old cells in depth-from-old-top coordinates. Surface melt removes |dh_top|
% from the top cell; basal growth adds a cell of thickness dh_bot at T_f;
% basal melt removes |dh_bot| from the bottom cell.
edges_old = h_old * ice.cum_fractions(:);  % N+1 edges, depth from the old top
T_cells   = T_new;
rc_cells  = rho_c_l;
if dh_bot > 0
    edges_old = [edges_old; h_old + dh_bot];
    T_cells   = [T_cells; T_f];
    [~, rc_new] = ice_properties(T_f, P);
    rc_cells  = [rc_cells; rc_new];
elseif dh_bot < 0
    edges_old(end) = h_old + dh_bot;      % |dh_bot| < dz guaranteed above
end
if dh_top < 0
    edges_old(1) = -dh_top;               % the melted sliver is gone
end
edges_new = edges_old(1) + h_new * ice.cum_fractions(:);   % new cells, same fractions
T_remap = zeros(N, 1);
for j = 1:N
    H = 0;  Cap = 0;                       % enthalpy and heat capacity collected
    for i = 1:numel(T_cells)
        overlap = min(edges_new(j+1), edges_old(i+1)) - max(edges_new(j), edges_old(i));
        if overlap > 0
            H   = H   + rc_cells(i) * T_cells(i) * overlap;
            Cap = Cap + rc_cells(i) * overlap;
        end
    end
    T_remap(j) = H / Cap;
end

% ------------------------------------------------------------------------
% Pack the state and the diagnostics
% ------------------------------------------------------------------------
S.T_sk_K  = T_sk;
S.T_ice_K = T_remap;
S.h_m     = h_new;
if has_snow, S.T_snow_K = T_snow_new; end

F.LWU_sfc_W_m2       = ice.emis_sfc*sig*T_sk^4 + (1 - ice.emis_sfc)*forcing.LWD_sfc_W_m2;
F.SH_down_W_m2       = SF.SH_down_W_m2;
F.LH_down_W_m2       = SF.LH_down_W_m2;
F.F_cond_top_up_W_m2 = F_top_up;
F.F_cond_bot_up_W_m2 = F_bot_up;
F.melt_W_m2          = melt;
F.dhdt_bot_m_s       = dhdt_bot;
F.dhdt_top_m_s       = dhdt_top;
F.lambda_SH_W_m2_K   = SF.lambda_SH;
F.lambda_LH_W_m2_K   = SF.lambda_LH;
F.lambda_LW_W_m2_K   = 4*ice.emis_sfc*sig*T_sk^3;   % the linearised Planck damping
F.L_eff_J_kg         = L_eff;
F.newton_iter        = n_iter;
F.skin_residual_W_m2 = resid - melt;      % 0 at convergence (melt removed)
[~, rc_now] = ice_properties(T_remap, P);
F.E_ice_sens_J_m2    = sum(rc_now .* (T_remap - T_f) .* (h_new * frac));
if has_snow
    F.E_ice_sens_J_m2 = F.E_ice_sens_J_m2 + ice.rho_c_snow_J_m3_K*ice.h_snow_m*(T_snow_new - T_f);
end

% ========================================================================
% nested functions
% ========================================================================
    function [T_ice_new_, T_snow_new_, F_top_up_, F_bot_up_] = conduction_solve(T_sk_)
        % Backward-Euler step of the heat equation with Dirichlet boundaries
        % T_sk_ (top) and T_f (bottom).  Unknowns: [snow; ice(1..N)] or ice only.
        if has_snow
            main = [cap_snow + K_top + K_si; ...
                    cap_ice + [K_si; K_int] + [K_int; K_bot]];
            off  = [K_si; K_int];                        % N couplings
            b    = [cap_snow*T_snow_old; cap_ice.*T_old];
        else
            main = cap_ice + [K_top; K_int] + [K_int; K_bot];
            off  = K_int;
            b    = cap_ice .* T_old;
        end
        A = diag(main) - diag(off, 1) - diag(off, -1);
        b(1)   = b(1)   + K_top * T_sk_;
        b(end) = b(end) + K_bot * T_f;
        T_all = A \ b;
        if has_snow
            T_snow_new_ = T_all(1);  T_ice_new_ = T_all(2:end);
        else
            T_snow_new_ = NaN;       T_ice_new_ = T_all;
        end
        F_top_up_ = K_top * (T_all(1) - T_sk_);        % up into the skin if interior warmer
        F_bot_up_ = K_bot * (T_f - T_ice_new_(end));   % up into the ice from the base
    end

    function R = skin_residual(T_sk_)
        % Net flux into the zero-heat-capacity skin; zero at the solution.
        [~, ~, F_top_up_] = conduction_solve(T_sk_);
        SF_ = seb_surface_fluxes(T_sk_, forcing.T_la_K, P);
        R = ice.emis_sfc * (forcing.LWD_sfc_W_m2 - sig*T_sk_^4) + forcing.F_extra_down_W_m2 ...
            + SF_.SH_down_W_m2 + SF_.LH_down_W_m2 + F_top_up_;
    end
end
