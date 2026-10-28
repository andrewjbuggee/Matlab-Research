# SW transmittance vs. CCN, grouped by cloud supersaturation (EPCAPE)

`sw_transmittance_vs_ccn.m` plots the Pier SW transmittance against the Pier CCN concentration, averaged over
consecutive `averageWindow_min` windows (e.g. 10 min or 1 h), coloured by the cloud's effective supersaturation
(binned, or on a continuous scale). Every setting, including the
filters, is in Section 1 at the top of the script. Press Run. The first run reads every instrument (~1–2 min) and
caches the 1-min inputs in `output/epcape_1min_inputs.mat`; later runs take seconds (`rebuildCache = true` re-reads).

## Main settings
- `caseSet`: `"dedrick24"` (the 24 Mt. Soledad cloud events), `"soledadInCloud"`, `"lowCloudPeriods"`, `"all"`
- filters, each with an on/off switch: daytime (`szaMax_deg`), LWP (`lwpMin_gm2`, `lwpMax_gm2`), no rain
  (`rainWindow_min`, `rainSites`), single layer, liquid (height proxy), persistence (`persistenceSites`,
  `persistHalfWindow_min`, ≥ 90 %), all tested minute by minute
- `averageWindow_min` (10, 15, 20, 30, 60, …) and `minPassFraction`: T_SW and LWP are averaged over the passing
  minutes of each window; a window is plotted if ≥ `minPassFraction` of its minutes pass. Independent of the
  persistence window: 10-min windows give 6 points for an hour in which every minute passes
- `ccnSetpoint_pct` (0.1–1.0 %), or `ccnAtRetrievedS = true` for the CCN spectrum at the retrieved s
- `supersaturationSource`: `"posted"` (Dedrick; `supersaturationMethod` = sNd/sCCN/skappa/sARG00) or `"AJB"`
  (Andrew's s_Nd,FM from `../supersaturation`, not yet verified)
- `binBySupersaturation`, `nSupersaturationBins` (5, 6, 7, …), `binEdgeRule` (`"quantile"` or `"equalWidth"`)

## functions/
| File | What |
|---|---|
| `build_epcape_1min_inputs.m` | RADFLUX, MWRRET v2 LWP, ceilometer, KAZR-ARSCL layers, disdrometers (Pier, Soledad), Mt. Soledad cloud condition → one 1-min grid |
| `analysis_filters_1min.m` | the filters, minute by minute (definitions and reasons in its header) |
| `load_hourly_ccn_and_supersaturation.m` | hourly Pier CCN spectra, Dedrick's posted s (skappa/sCCN label swap fixed), Andrew's s_Nd,FM |
| `case_set_mask_1min.m` | which minutes belong to each case set |
| `average_passing_minutes.m` | T_SW, LWP and SZA of each averaging window, from the passing minutes only |
| `hourly_values_on_windows.m` | the hourly CCN and s put onto the averaging windows (mean over the passing minutes) |
| `hourly_from_passing_minutes.m` | 60-min call of `average_passing_minutes.m`, kept for the data-loss script |

## With the default settings (27 Oct 2026)
First version (hourly, rain ±30 min, persistence ±30 min, LWP < 250 g m⁻²): 13 hours from 6 of the 24 events.
With rain ±15 min, persistence ±10 min, LWP < 300 g m⁻²: 14 hours with 60-min windows; 100 points from 25 hours
with 10-min windows (≥ 5 passing minutes each), the same 6 events. Daytime and the Mt. Soledad half of the persistence test remove
the most; see `../Dedrick_2026_coupled_cases/dedrick2026_data_loss_analysis.m`. Dedrick's low-cloud periods with
`persistenceSites = "Pier"` give 191 hours.

## Caveats
- **"Liquid"** is a proxy (cloud top and base < 3 km). No cloud-top temperature is on disk (only 6 radiosondes, no INTERPSONDE).
- **"Cloud at Mt. Soledad"** can only mean Mt. Soledad *inside* the cloud (visibility < 1 km, FM-120 LWC > 0.01 g m⁻³).
  A deck with its base above 250 m MSL is invisible there with the local data.
- **Single layer** uses ARSCL + ceilometer, not the NDROP `cloud_base_type` that Dedrick used (NDROP is not downloaded).
- **CCN and s are hourly.** The CCN-200 measures each setpoint about once per 1-h cycle, and the posted
  supersaturations are hourly, so windows shorter than an hour share their hour's CCN and s: in panel (a) the
  10-min points of one hour sit in a vertical column. Only T_SW and LWP vary within the hour.
- **LWP dominates T_SW.** Look for a CCN/supersaturation signal at fixed LWP, and remember that hours within one
  event are not independent.

Data and references: see the script header. Dependencies: `../SW_transmittance`, `../Dedrick_2026_coupled_cases`,
`../supersaturation/functions`.
