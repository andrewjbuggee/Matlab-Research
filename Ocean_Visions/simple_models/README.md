# A zero-dimensional ocean / sea-ice / atmosphere / cloud column for polar night

MATLAB, R2024b. Andrew Buggee with Claude, Ocean Visions, October 2026.
Rebuilt 26 Oct 2026 as a modular hierarchy: start from `run_simple_model.m`.
The report "Modular ocean–ice–atmosphere column model: how the code works" (Claude Doc,
26 Oct 2026) documents every equation with its source and file:line.

The purpose is the one in the 22 Oct 2026 meeting notes: *all of the physics is well
known, but the pieces have not been put together in the simplest model that can show
how the ice responds when the cloud above it changes.* Everything here is standard
physics with a citation attached; the one unavoidable parameterization (advection,
see §4) is labelled as such; and the model is checked against the exact analytic
solutions it should reproduce before it is used for anything.

## 1. Structure

```
            space            OLR
   ┌──────────────────────┐
   │   upper atmosphere   │  T_ua   grey slab, ε_ua, heat capacity c_p Δp/g
   ├──────────────────────┤
   │     cloud layer      │  T_cl   grey slab, ε_c(LWP or IWP), c_p Δp/g
   ├──────────────────────┤
   │   lower atmosphere   │  T_la   grey slab, ε_la;  exchanges SH, LH with the skin
   ╞══════════════════════╡
   │   skin (no heat cap.)│  T_sk   DIAGNOSTIC: fluxes into it sum to zero
   ├ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─┤
   │   snow (optional)    │  T_snow fixed depth, k = 0.3 W m⁻¹ K⁻¹
   ├──────────────────────┤
   │  ice layer 1 (top)   │  T_ice(1)   ┐
   │  ice layer 2         │  T_ice(2)   │ heat equation, implicit; layers are
   │  ice layer 3 (base)  │  T_ice(3)   ┘ fixed FRACTIONS of h (default 0.1/0.8/0.1)
   ╞══════════════════════╡  T_f      base held at the seawater freezing point;
   │  ocean mixed layer   │  T_ml     h grows by the STEFAN CONDITION
   ├──────────────────────┤
   │  deep ocean          │  T_deep
   └──────────────────────┘
```

State vector: `T_ua, T_cl, T_la, T_sk, (T_snow), T_ice(1..N), h, T_ml, T_deep`.
Winter only — no solar terms anywhere.

## 2. The equations (where each comes from)

All fluxes W m⁻², heat capacities J m⁻² K⁻¹. Sign convention — stated once, used
everywhere: **turbulent fluxes positive downward into the surface; conductive and ocean
fluxes positive upward; dh/dt positive for growth.** Sledd et al. (2025) eq. (1) has
turbulent fluxes positive *away* from the surface, so their SH is minus ours.

| process | equation | source |
|---|---|---|
| longwave between slabs | grey, non-scattering slabs: `D_k = ε_k σT_k⁴ + (1−ε_k) D_{k+1}` downward, `U_{k+1} = ε_k σT_k⁴ + (1−ε_k) U_k` upward, `U_1 = ε_s σT_sk⁴ + (1−ε_s) D_1`; slab heating `Q_k = ε_k(U_k + D_{k+1}) − 2ε_kσT_k⁴` | Pierrehumbert (2010) §4.3; Zapalac (2023) for "no LW scattering by clouds" |
| cloud emissivity | `ε_c = 1 − exp(−k_abs·WP)`, liquid `k_abs = 0.130 m² g⁻¹ × (10 µm / r_e)`, ice `k_abs = 1.66·3Q_abs/(4ρ_i r_e)`; the cloud *slab* combines cloud and its own gases as independent absorbers, `ε_slab = 1 − (1−ε_gas)(1−ε_c)` | Stephens (1978) eq. 13 (liquid anchor); geometric optics + Elsasser (1942) diffusivity factor for the r_e scaling and for ice — an **estimate**, see `cloud_lw_emissivity.m` |
| sensible heat | `SH↓ = ρ c_p C_H U (T_la − T_sk)` | bulk aerodynamic; Andreas et al. (2010) for winter sea ice |
| latent heat | `LH↓ = ρ L_s C_E U (q_la − q_sat,ice(T_sk))`, `q_la = RH_ice·q_sat,ice(T_la)` | Murphy & Koop (2005) eq. 7 for `e_sat,ice`; Andreas et al. (2002) for RH ≈ ice saturation |
| skin balance | `ε_s(LWD − σT_sk⁴) + F_extra + SH↓ + LH↓ + F_cond,top↑ = 0` | Sledd et al. (2025) eq. (5), the "infinitely thin interface"; IFS Cy41r2 eq. 8.22; Winton (2000) |
| conduction | `ρc ∂T/∂t = ∂/∂z (k ∂T/∂z)`, `T(top) = T_sk`, `T(base) = T_f`; backward Euler | Sledd et al. (2025) eq. (7) (no solar term); Maykut & Untersteiner (1971); Semtner (1976) |
| **ice growth** | **`ρ_i L dh/dt = F_cond,bot↑ − F_ocean↑`** | **Stefan (1891); Carslaw & Jaeger (1959) ch. XI; Maykut & Untersteiner (1971) eq. 5** |
| surface melt | if `T_sk` would exceed `T_melt`: cap it, surplus `M` melts the top | Sledd et al. (2025) "melt regime"; not reached in winter |
| regridding | enthalpy-conserving remap onto the fixed layer fractions after `h` changes | Bitz & Lipscomb (1999) §2e |
| ocean → ice | `F_ocean↑ = ρ_w c_w c_h u*(T_ml − T_f)`, `c_h = 0.0056` | McPhee (1992); gives ≈ 2 W m⁻², the Maykut & Untersteiner (1971) value |
| freezing point | `T_f(S) = −0.0575S + 1.710523×10⁻³S^1.5 − 2.154996×10⁻⁴S²` °C → 271.28 K at 34 psu | Millero (1978), UNESCO |
| ocean slabs | `C_ml dT_ml/dt = γ_deep(T_deep − T_ml) − F_ocean↑`, `C_deep dT_deep/dt = −γ_deep(T_deep − T_ml)` | slab exchange after the PRISM three-layer EBM (`multliayerModel_ajb.ipynb`) |
| atmosphere slabs | `C_k dT_k/dt = Q_LW,k ± γ_mix(T_cl − T_la) − (SH↓+LH↓ for la) + C_k(T_bg,k − T_k)/τ_adv` | slab bookkeeping and `C = c_p Δp/g` after the PRISM notebook; relaxation: Pithan et al. (2016) SCM practice |

Optional brine corrections to k, c and L (Bitz & Lipscomb 1999 eqs. 3, 4, 6) are in
`ice_properties.m` behind `use_brine` (off by default).

### How phase change is handled — the answer to the first question in the notes

Ice growth is **not** parameterized. At the ice–ocean interface the temperature is the
freezing point of seawater (a thermodynamic property, not a choice), and energy
conservation across that moving interface — the Stefan condition — says that whatever
heat conduction carries up into the ice faster than the ocean can supply it from below
must come from freezing seawater, releasing `L` per kilogram:

```
ρ_i L (dh/dt) = k (∂T/∂z)|base − F_ocean
```

The only parameterized inputs are the ocean heat flux (a bulk coefficient, McPhee 1992)
and the choice of `L` (fresh, or brine-adjusted). Surface melt is the same statement at
the top. That is why `test_column_model.m` can check the model against the *exact*
analytic solution of the Stefan problem (§5).

The mechanism asked about in the notes — "glaciate the cloud, the top of the ice cools,
how does that propagate down to the base where freezing happens" — is then: the skin
cools within an hour (zero heat capacity), the conductive flux into the top layer rises,
the cooling diffuses downward on the timescale `dz²/κ` (κ = k/ρc ≈ 1.05×10⁻⁶ m² s⁻¹:
3 h for 10 cm, 1.2 d for 30 cm, 2.4 d for the whole of 1.5 m), the temperature gradient
at the base steepens, and the basal growth rate `k ∂T/∂z|base / (ρ_i L)` increases.
Figure 4 of the driver shows exactly this.

## 3. Files

| file | what it is |
|---|---|
| `run_simple_model.m` | **master script**: choose the complexity level and the physics in one block, draw the configuration, run the perturbation experiment, compare all levels |
| `model_config.m` | builds the parameter/option struct for a level (`model_config(3, 'turb_scheme', 'most', ...)`) and prints a summary |
| `seb_column_params.m` | every parameter and option, with units in the name and a source in the comment |
| `seb_column_finalize.m` | derived quantities and the layer structure of each complexity level (edit `P.atm.layers` by hand, then re-run it) |
| `seb_column_run.m` | the integrator: any number of air layers, either ice representation, either advective closure; energy diagnostic |
| `draw_model_diagram.m` | renders the configured column (layers, fluxes, options) as a PNG |
| `ice_column_step.m` | skin balance (Newton) + implicit conduction (+ snow) + Stefan growth + conservative regrid; also the zero-layer ice |
| `ice_slab_step.m` | the single-temperature ice slab of level 1 |
| `ice_properties.m` | k, ρc, L for pure ice, or the Bitz & Lipscomb (1999) brine corrections |
| `lw_grey_slabs.m` | the grey-slab longwave solver (any number of slabs; emission supplied per slab) |
| `planck_emit.m` | σT⁴ or its linearisation (North et al. 1981; IFS eq. 3.84) |
| `slab_emission.m` | up/down emission of a slab: isothermal, or with a lapse-rate profile inside it |
| `lapse_rate.m` | dry, saturated-adiabatic, or prescribed lapse rate |
| `cloud_lw_emissivity.m`, `cloud_slab_emissivity.m` | Stephens (1978) + geometric-optics scaling; cloud + gases combined per slab |
| `seb_surface_fluxes.m`, `most_exchange_coefficient.m` | SH and LH: none / linear / bulk constant / IFS Monin–Obukhov |
| `sat_vapor_ice.m` | Murphy & Koop (2005) over ice |
| `response_partition.m`, `print_column_state.m` | the partition bookkeeping and a state printout |
| `new_figure_clean.m` | invisible, factory-default figure (R2024b export-stall workaround; resets `startup.m`'s dark defaults) |
| `test_column_model.m` | the validation suite, T1–T9 (§5) |
| `run_glaciation_experiment.m` | the 22 Oct experiment at level 3: glaciation scenarios, thickness dependence, profiles |
| `plot_dlr_vs_lwp.m` | side question from the notes: what moves the LWP at which DLR saturates |

### The hierarchy (`complexity` 1–4) and the physics switches

| level | atmosphere | ice | ocean |
|---|---|---|---|
| 1 | one slab (790 hPa) holding the cloud as an emissivity | single slab, one temperature | mixed layer over a fixed reservoir |
| 2 | boundary layer with the cloud (90 hPa) + free troposphere | skin + 2 layers | mixed layer + reservoir |
| 3 | sub-cloud air + cloud layer + free troposphere | skin + 3 layers [0.1 0.8 0.1] | mixed layer + reservoir |
| 4 | as 3, with a lapse-rate profile inside each boundary-layer slab | as 3 | as 3 |

Independently of the level: `lw_emission` (planck | linearized), `turb_scheme` (bulk_const | most | linear | none),
`latent_scheme` (bulk | off), `atm_profile` (isothermal | lapse) with `lapse_kind` (dry | moist | K m⁻¹),
`ice_model` (skin_profile | slab | zero_layer), `ice_growth` (stefan | fixed, ERA5-like), `atm_forcing`
(relaxation | heating, Thorndike 1992's prescribed convergence D), `ocean_model`, `h_snow_m`, `use_brine`,
`layer_fractions`. Anything set explicitly is honoured at every level; `'auto'` takes the level's default.

Run from the directory (`addpath` is handled):

```matlab
run_simple_model             % ~3 min: diagram, experiment at the chosen level, all levels compared
test_column_model            % ~2 min; prints PASS/FAIL for T1-T9
run_glaciation_experiment    % ~1-2 min; writes figures/fig1..fig4 (level 3)
plot_dlr_vs_lwp
```

Minimal programmatic use:

```matlab
P  = seb_column_params('WP_g_m2', 30, 'h_init_m', 1.5, 'h_snow_m', 0.2, 't_end_days', 30);
R  = seb_column_run(P);                               % spin up
P2 = P;  P2.cloud.phase = 'ice';  P2.cloud.WP_g_m2 = 1.5;
[P2.cloud.emis, P2.cloud.k_abs_m2_g] = cloud_lw_emissivity(P2.cloud.WP_g_m2, 'ice', P2);
R2 = seb_column_run(P2, 'init_state', R.final_state);  % branch from the spun-up state
plot(R2.t_days, 100*(R2.h_m - R.h_m(end)))             % extra ice, cm
```

## 4. What is tuned, and what is a free knob

In polar night the atmosphere above sea ice has no energy source except horizontal
advection. A column model must represent that somehow; here each atmospheric slab
relaxes toward a background temperature on a timescale `τ_adv` (default 1 day ≈ 500 km /
5 m s⁻¹). **This is the one unavoidable parameterization**, and it is also the knob that
decides how tightly the air follows the surface — the α = dT_air/dT_skin of the ERA5
analysis (`linearized_seb_model.ipynb`): `τ_adv → 0` pins the air (α = 0, the
fixed-atmosphere limit, run with `mode = 'prescribed'`); `τ_adv → ∞` lets the air cool
with the surface (α → 1, what the ERA5 regressions show over pack ice).

Values marked **TUNED** in `seb_column_params.m`:

* `T_bg_la_K = T_bg_cl_K = 260 K`, `T_bg_ua_K = 232 K`, `emis_la = 0.40`, `emis_gas_cl = 0.15`,
  `emis_ua = 0.63` — the backgrounds and `emis_la` were chosen by a
  bounded grid search (one shared boundary-layer background, colder free troposphere)
  so that the default ERA5-like column (bare 1.5 m ice) sits at the two observed winter
  states: clear sky `T_sk = 247.8 K, LWD = 165 W m⁻²` (4-K inversion, SH↓ ≈ +16 W m⁻²);
  LWP = 30 g m⁻² `T_sk = 255.4 K, LWD = 227 W m⁻²` (targets 245/165 and 255/230; the
  ERA5 Barrow pack-ice means of `linearized_seb_model.ipynb`). An unconstrained search
  found a lower cost at a 301-K cloud background — the optimizer using the cloud slab's
  relaxation as a heater — which is why the constraint is there. The backgrounds sit
  5–8 K above the slab temperatures they produce; that offset is the radiative cooling
  that advection must replace (~20 W m⁻² over the boundary layer at `τ_adv = 1 d`).
  `emis_gas_cl` (the clear-air emissivity of the 300-m cloud layer, so that a vanished
  cloud leaves air that still radiates) was then set to 0.15 and `emis_ua` lowered from
  0.75 to 0.63 to leave the clear-sky LWD unchanged; the cloudy state is unaffected
  (0.98 → 0.983).
* `C_H = 0.5×10⁻³` — the *effective* coefficient between the skin and a 400-m slab.
  Deliberately ~3× below the 2-m-based coefficient (ERA5 gives ρc_pC_HU = 9.9 W m⁻² K⁻¹
  over pack ice), because the 2-m air follows the skin almost one-for-one in a stable
  boundary layer. **The sensible-heat share of any response scales with this number**;
  the driver runs a 3× case so the sensitivity is on the page.
* `gamma_mix`, `gamma_deep` — order-of-magnitude choices, documented in the file.

Everything else is a physical constant or a published value.

## 5. What the tests check (`test_column_model.m`)

| test | what | expected |
|---|---|---|
| T1 | grey-slab LW conserves energy | residual ~10⁻¹³ W m⁻² |
| T2 | ice growth vs the exact one-phase Neumann solution of the Stefan problem (`h = 2λ√(κt)`, `λe^{λ²}erf λ = St/√π`) and the quasi-steady law `h² = h₀² + 2kΔT t/(ρL)` | within 1 % of Neumann for every layering; the quasi-steady law overshoots by ~2 % because it neglects heat capacity (St = 0.13) |
| T3 | coupled run conserves energy | residual ~10⁻⁷ W m⁻² |
| T4 | the two winter states (clear / cloudy) | skin ~245 / 255 K, LWD ~165 / 230 W m⁻² for the ERA5-like column |
| T6 | the options: lapse-rate emission reduces to isothermal at Γ = 0; Monin–Obukhov C_H equals κ²/(ln ln) when neutral and falls when stable; linearised vs full Planck; fixed thickness; linear vs bulk turbulence; the MOST scheme runs | all pass |
| T7 | every level, clear and cloudy, energy closure | residual < 4×10⁻⁶ W m⁻² at all four levels |
| T8 | zero-layer ice + linearised Planck + no turbulence under a fixed atmosphere vs its analytic growth ODE (Thorndike 1992, eqs. 13–15) | 0.00 % after 60 days |
| T9 | `atm_forcing = 'heating'`: the Thorndike closure, D = 100/140/180 W m⁻² | energy closes; D ≈ 180 W m⁻² reproduces the winter states with these emissivities |
| T5 | step response of the skin to +10 W m⁻² vs the frequency-dependent admittance `λ_G(P) = Re[k q coth(qh)]` at `P = 4t` (exact for a semi-infinite conductor) | 30 equal layers: 0.28 / 0.52 / 0.74 / 0.89 K at 1 h / 6 h / 1 d / 3 d against 0.27 / 0.49 / 0.69 / 0.83 from the admittance; **3 equal layers give 0.63 K at 1 h** (50-cm top layer) and converge after ~3 days; the default 0.1/0.8/0.1 split is within 30 % at 1 h, ~15 % high at 1–3 d, converged by 1 week |

The T5 result is the practical limitation to carry: a 3-layer ice model is exact for
steady states and within 0.5 % for the Stefan growth law, and good for responses slower
than the top layer's diffusion time, but it cannot resolve the hourly skin response
unless the top layer is a few centimetres thick (which is why the IFS puts its top ice
layer at 7 cm). The thick middle layer of the 0.1/0.8/0.1 split also makes the 1–3-day
skin response ~15 % too large. Use `N_layers = 12` or more when sub-weekly skin
temperatures matter; the 30-day energy partition and the extra ice are insensitive to
the layering (the driver prints the 12-layer check next to the 3-layer headline).

## 6. What the glaciation experiment does

`run_glaciation_experiment.m` spins up under a liquid cloud (LWP = 30 g m⁻², the
CESM2 winter mean over sea ice in Li et al. 2026) and branches: glaciation at constant
mass (ε_c 0.98 → 0.74), glaciation with the 95 % condensate loss to precipitation found
in the Villanueva et al. (2022) LES (IWP = 1.5 g m⁻², ε_c → 0.07), complete
dissipation, the same with a fixed atmosphere (α = 0 bound), with 3× C_H, and with
20 cm of snow. For each it partitions the LWD forcing — at every instant, by the skin
balance, into less LWU, more downward SH, LH and more conduction out of the ice; and the
conduction, by the ice energy budget, into **cooling the ice** and **making ice** — and
compares the extra ice with the Zapalac (2023) bound `Δh = ΔF·t/(ρL)` in which every
W m⁻² freezes water. Figure 3 repeats this across initial thickness (0.5–2.5 m) with and
without snow.

### Results with the defaults (22 Oct 2026; `run_glaciation_experiment` prints these)

Spun-up cloudy state over 1.67 m of bare ice: `T_sk = 254.4 K`, `LWD = 222 W m⁻²`,
conduction 21 W m⁻² up to the skin, growth 0.51 cm d⁻¹, λ_LW / λ_SH / λ_LH = 3.7 / 3.5 /
0.7 W m⁻² K⁻¹. Glaciation + precipitation drops LWD by **−54 W m⁻²** at t = 0⁺ (ε_c 0.98
→ 0.07); the surface relaxes toward the clear-sky state over ~10 days (skin −7.0 K after
30 d; the model's two steady states differ by 7.6 K and 59 W m⁻²).

Instantaneous shares of the forcing (sum = 1 by the skin balance):

| after | LWU | SH | LH | conduction out of the ice | ΔT_skin |
|---|---|---|---|---|---|
| 1 h | 0.15 | 0.12 | 0.02 | **0.71** | −2.0 K |
| 6 h | 0.24 | 0.16 | 0.03 | 0.57 | −3.4 K |
| 1 d | 0.37 | 0.27 | 0.04 | 0.32 | −5.4 K |
| 3 d | 0.41 | 0.31 | 0.05 | 0.24 | −6.1 K |
| 30 d | 0.47 | 0.34 | 0.05 | 0.14 | −7.0 K |

The conductive share falls from 71 % to 14 % as the ice column cools — the same
frequency dependence of the ice admittance found in the ERA5 analysis (`λ_G` from 58 to
1.3 W m⁻² K⁻¹ between hourly and steady forcing).

Cumulative over 30 days (141 MJ m⁻² of forcing): **LWU 44 %, SH 32 %, LH 5 %, conduction
19 % — of which cooling the ice 9 % and making ice 10 %.** Extra ice **4.5 cm**, against
46 cm if all the forcing had frozen water (Zapalac's bound): a growth efficiency of 10 %.
For comparison, the observed winter responses to LWD are LWU ≈ 0.5, SH ≈ 0.25, G ≈ 0.2 at
MOSAiC (Sledd et al. 2025, Fig. 2a) and LWU 0.65–0.85, SH 0.11, C 0.23 at Summit (Miller
et al. 2017).

Sensitivities: 12 ice layers instead of 3 → shares unchanged, extra ice 4.6 cm. Fixed
atmosphere (α = 0 bound) → SH share 0.40 at 30 d, extra ice 4.0 cm. 3× C_H → SH 0.50,
extra ice 3.2 cm. τ_adv = 5 d → extra ice 3.0 cm (smaller forcing, −38 W m⁻², because
the cloudy equilibrium itself is colder). The 400-m slab follows the skin with
α = ΔT_air/ΔT_skin ≈ 0.24 in every case — the slab is held by mixing with the cloud
layer and relaxation, which is why the *effective* C_H is below the 2-m value.

Thickness (extra ice after 30 d, glaciate + precipitate):

| h₀ | 0.5 m | 1.0 m | 1.5 m | 2.5 m |
|---|---|---|---|---|
| bare ice | 8.5 cm (18 % of bound) | 6.3 (14 %) | 4.5 (10 %) | 2.2 (5 %) |
| + 20 cm snow | 5.0 cm (11 %) | 3.7 (8 %) | 2.6 (6 %) | 1.2 (3 %) |
| conduction share: cooling / growth (bare) | 0.06 / 0.18 | 0.07 / 0.14 | 0.09 / 0.10 | 0.12 / 0.05 |

So Ray's "more energy goes into cooling the ice than making ice" is true for thick ice
(≥ 1.5 m) over a month and reversed for thin ice; Matt Shupe's "growth is fastest for the
thinnest ice" follows from the Stefan law (`dh/dt ∝ 1/h`); and 20 cm of snow roughly
halves the response (Ian Eisenman's point).

## 7. Known limitations (honest list)

* No snow accumulation, no precipitation, no solar, no leads, no ice dynamics.
* Grey emissivities for the gaseous atmosphere are tuned, not computed from water vapour.
* `C_H` constant: no stability dependence, so the SH response is an upper bound for a
  surface layer that becomes very stable after glaciation (Grachev et al. 2007).
* The cloud slab is a radiative-only object: its liquid water path is prescribed, not
  predicted; cloud-top cooling does not feed back on the cloud's own water.
* The ice-side admittance is resolved only to the top layer's thickness (T5).
* Ocean: a bulk mixed layer and a reservoir; no halocline structure.

## 7b. Code review, 26 Oct 2026 (what changed)

* `seb_column_finalize.m`: the deep-ocean switch `P.ocean.deep_prognostic` was evaluated
  *before* the `'auto'` options were resolved, so a fresh configuration ran with a fixed
  reservoir and a re-finalized copy (the perturbed run of every experiment) with a
  prognostic one. Control and perturbed runs therefore used different ocean models
  (≤ 0.01 cm of extra ice at levels 2 and 4; nothing at level 3 at the printed precision).
  Now set after the resolution. `N_layers` / `layer_fractions` default to `'auto'` (the
  level's layering) and an explicit value is honoured at every level — level 2 no longer
  overrides it silently. A layer's `lapse` is `'auto'` = `P.opt.lapse_kind`, read at run
  time, so editing `lapse_kind` after construction now propagates. Warning when snow is
  set with `ice_model = 'slab'` (the slab has no skin and ignored it silently).
* `slab_emission.m`: the 41-point trapezoid rule (0.09 % low at ε = 0.98, NaN at ε = 1)
  is replaced by the exact exponential-kernel integration between nodes (linear-in-τ,
  Clough et al. 1992); the isothermal limit is now exact to round-off and the black
  limit finite.
* `ice_slab_step.m`: the melt surplus lacked the slab heat-capacity term (dormant in
  winter). The energy diagnostic now counts the sensible heat that leaves with melted
  ice (`E_runoff`), so it closes in the melt regime as well (it was 2–3 W m⁻² off
  under a forced melt).
* `P.num.newton_tol_W_m2 = 1e-4` replaces `newton_tol_K` (it is a flux tolerance); the
  skin solve warns if it does not converge.
* Comments and labels: IFS citations updated to Cy50r1 (eqs. 3.76/8.25, 8.8–8.9,
  Table 8.3, §8.11 eq. 8.189); `planck_emit` reference note; test header;
  `plot_dlr_vs_lwp` labels follow `emis_la`; the dead parameter `deep_prognostic`
  removed from `seb_column_params`; the run summary prints the live cloud emissivity.
  The Monin–Obukhov sensible heat flux is *not* monotone in the inversion strength
  (it saturates near 32 W m⁻² for 8–16 K inversions and falls by 0.2 W m⁻² K⁻¹ before
  the ζ = 5 cap); the skin solve keeps a unique root because the Planck and conduction
  slopes dominate — documented in `most_exchange_coefficient.m`.
* Known and not fixed: with `use_brine` the ice heat capacity depends on T and the
  energy diagnostic is a linearisation (residual ≈ 2 W m⁻² over 10 days); a proper
  enthalpy formulation (Bitz & Lipscomb 1999) would be needed.

All nine tests pass after the changes; the results of §6 are unchanged at the printed
precision.

## 8. Sources

* Andreas, E. L., et al. (2002). Near-surface water vapor over polar sea ice is always
  near ice saturation. *JGR* 107, 8033.
* Andreas, E. L., et al. (2010). Parameterizing turbulent exchange over sea ice in
  winter. *J. Hydrometeorol.* 11, 87–104.
* Batrak, Y. & Müller, M. (2019). On the warm bias in atmospheric reanalyses induced by
  the missing snow over Arctic sea-ice. *Nat. Commun.* 10, 4170.
* Bitz, C. M. & Lipscomb, W. H. (1999). An energy-conserving thermodynamic model of
  sea ice. *JGR* 104, 15669–15677.
* Carslaw, H. S. & Jaeger, J. C. (1959). *Conduction of Heat in Solids*, 2nd ed.
* ECMWF (2016). IFS Documentation Cy41r2, Part IV (eq. 8.22, 8.149–8.150).
* Elsasser, W. M. (1942). Heat transfer by infrared radiation in the atmosphere.
  *Harvard Meteorol. Studies* 6 (the 1.66 diffusivity factor).
* Li, X., Evan, A. T., Eisenman, I. & Kok, J. F. (2026, preprint). Limited climate
  response to extreme Arctic mixed-phase cloud thinning. Research Square.
* Maykut, G. A. & Untersteiner, N. (1971). Some results from a time-dependent
  thermodynamic model of sea ice. *JGR* 76, 1550–1575.
* McPhee, M. G. (1992). Turbulent heat flux in the upper ocean under sea ice. *JGR* 97,
  5365–5379.
* Millero, F. J. (1978). Freezing point of seawater. UNESCO Tech. Pap. Mar. Sci. 28.
* Miller, N. B., et al. (2017). Surface energy budget responses to radiative forcing at
  Summit, Greenland. *The Cryosphere* 11, 497–516.
* Murphy, D. M. & Koop, T. (2005). Review of the vapour pressures of ice and
  supercooled water. *QJRMS* 131, 1539–1565.
* Pierrehumbert, R. T. (2010). *Principles of Planetary Climate*, §4.3.
* Pithan, F., et al. (2016). Select strengths and biases of models in representing the
  Arctic winter boundary layer over sea ice: the Larcform 1 SCM intercomparison.
  *JAMES* 8, 1345–1357.
* PRISM project: `multliayerModel_ajb.ipynb` (W. Camp / A. Buggee) — three-slab EBM
  structure, `C = c_p p_s/g`, slab exchange coefficients, equilibrium/timescale
  printouts, stiff integration.
* Semtner, A. J. (1976). A model for the thermodynamic growth of sea ice in numerical
  investigations of climate. *J. Phys. Oceanogr.* 6, 379–389.
* Sledd, A., Shupe, M. D., Solomon, A. & Cox, C. J. (2025). Surface energy balance
  responses to radiative forcing in the central Arctic from MOSAiC and models. *JGR
  Atmos.* 130, e2024JD042578.
* Stefan, J. (1891). Über die Theorie der Eisbildung, insbesondere über die Eisbildung
  im Polarmeere. *Ann. Phys.* 278, 269–286.
* Stephens, G. L. (1978). Radiation profiles in extended water clouds. II. *J. Atmos.
  Sci.* 35, 2123–2132.
* Sturm, M., Holmgren, J., König, M. & Morris, K. (1997). The thermal conductivity of
  seasonal snow. *J. Glaciol.* 43, 26–41.
* Villanueva, D., et al. (2022). Mixed-phase regime cloud thinning could help restore
  sea ice. *Environ. Res. Lett.* 17, 114057.
* Winton, M. (2000). A reformulated three-layer sea ice model. *J. Atmos. Ocean.
  Technol.* 17, 525–531.
* Zapalac, G. (2023, preprint). Glaciation of mixed-phase clouds by dry ice to cool the
  Arctic. Research Square.
* This project: `ERA5/surface_energy_budget/linearized_seb_model.ipynb` (Sep 2026) —
  ERA5 Barrow regressions and the frequency-dependent ice admittance.
