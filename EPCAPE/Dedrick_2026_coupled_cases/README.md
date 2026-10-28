# Dedrick et al. (2026) Mt. Soledad cloud events: filters and data-loss analysis

Functions that restrict any EPCAPE data to the 24 Mt. Soledad cloud events of Dedrick et al. (2026), rebuild
how those events were selected, and a script that counts the data lost at every selection step and filter.

## Credits
- **Method and event list:** Dedrick, J. L., and coauthors, 2026: Quantitative estimates of upwind and on-site
  particle scavenging and retrieved cloud supersaturations from urban coastal stratocumulus during EPCAPE.
  *J. Atmos. Sci.*, 83, 1061–1082, doi:10.1175/JAS-D-26-0004.1; supplement doi:10.1175/JAS-D-26-0004.s1
  (Text S7, Tables S4–S6). The event algorithm is credited there to Williams et al. (2026, submitted).
- **Data:** Berta, V., S. Han, A. Osawa, J. Dedrick, N. Maneenoi, and L. M. Russell, 2025: Aerosol quantification
  products at Mt. Soledad and Scripps Pier during EPCAPE 2023-24. UC San Diego Library Digital Collections,
  doi:10.6075/J0WM1DSC (GCVI 1-s file, corrected Pier CCN, supersaturations, low-cloud periods). Mt. Soledad
  visibility (UCSD group); FM-120 5-min files (Dalhousie: L. Robinson, R. Chang); ARM EPCAPE products.

## Files
| File | What |
|---|---|
| `dedrick2026_soledad_cloud_events.m` | The 24 events (pre-, in-, post-cloud windows, UTC), copied from supplement Table S5 with `pdftotext` |
| `filter_to_dedrick2026_events.m` | **Filter any data** (datetime vector, timetable or table) to the events; handles averaged samples (`SampleLength`, `TimeLabel`, `MinOverlap`) and can keep only the 5-min cloud intervals |
| `read_soledad_visibility.m` | The ~1-s Mt. Soledad visibility record (both CSVs, overlap removed) |
| `soledad_cloud_condition_5min.m` | Table S4's 5-min cloud / haze / out-of-cloud classes (visibility + FM-120 LWC) |
| `gcvi_sampling_segments.m` | Continuous GCVI sampling periods from the 1-s enhancement-factor file (port of the Python `gcvi_segments`) |
| `find_soledad_cloud_events.m` | Rebuilds the Table S4 selection (154 → 95 → 44 → 24) with every threshold an input |
| `dedrick2026_data_loss_analysis.m` | **The data-loss script** (Parts A–D below); writes `output/dedrick2026_data_loss_tables.xlsx` and a figure |

Example:
```matlab
SS = read_posted_product('.../bb2743661p_34_1_supersaturation.txt');      % hourly
[keep, event_id, SS24] = filter_to_dedrick2026_events(SS, 'SampleLength', hours(1));
```

## Things found while building these (27 Oct 2026)
- **Table S5 "cloud hours" are the window lengths.** For all 24 events they equal the in-cloud window length to
  within ±5 min, and the windows add up to the paper's 214 h. Counting only the 5-min intervals that meet the
  cloud condition gives ~137 h. Table S4's step totals (342 h, 272 h), by contrast, behave like cloud-interval counts.
- **The posted `*_SoledadCloudEvent` supersaturations are not the 24 events:** 278 h from 12 Mar 2023 to 21 Jan
  2024; 141 h fall inside the 24 windows, none in events 23–24 (October). 266 of the 278 h are hours that Abbey
  Williams' 15-min file classes as cloud. Times are UTC (a 0-h shift fits best).
- **Rebuilt selection:** 162 / 96 / 46 events at steps 1–3 (paper 154 / 95 / 44). The final SEMS-coverage step
  cannot be tested (no SEMS data here). 20 of the 24 events are recovered; event 3 has no GCVI data in the
  posted file (30 Apr 2023 00:32–16:32 UTC), and events 15, 16, 22 have ≤ 2 h of cloud intervals by this classification.

## Needs
Base MATLAB R2024b. The data-loss script also needs `../SW_trans_vs_CCN/functions`, `../SW_transmittance`,
`../supersaturation/functions`, and the 1-min cache built by `../SW_trans_vs_CCN/sw_transmittance_vs_ccn.m` (run it once first).
