# SW transmittance at the Scripps Pier (MATLAB translation of Andrew's Python method)

`T_SW = <SW down> / <clear-sky SW down>` from ARM RADFLUX at the Pier (1-min, QC = 0, sun above 10°), averaged
over fixed windows. MATLAB versions of the functions in
`VS_CODE/Python-Research/EPCAPE/comparisons/seasonal_averages/quantities.py` and `kavin.py`.

## Files
| File | Python original | What |
|---|---|---|
| `read_radflux_sw.m` | `qc0_values` | RADFLUX SW, clear-sky SW, cos(SZA); QC = 0; minute START and END times |
| `window_mean.m` | `grid_mean`, `window_mean` | means over fixed windows anchored at midnight |
| `sw_transmittance.m` | `sw_transmittance` | the transmittance; adds an optional `Keep` mask for extra filters |
| `hourly_ratio.m` | `kavin.hourly_ratio` | hourly `<num>/<den>` with a minimum sample count |
| `sw_clearness_index.m` | `sw_clearness_index` | SW / top-of-atmosphere SW (sensitivity test) |
| `read_lubin_sw_transmittance.m` | `kavin.read_lubin_sw`, `lubin_sw_series` | Dan Lubin's hourly tables as one UTC series |
| `compare_sw_transmittance_with_lubin.m` | — | **comparison script** → `output/` (figure + matched hours) |

## Verification (27 Oct 2026)
On the full RADFLUX file, `sw_transmittance` (with `'TimeField','time_end'`, the Python time base) reproduces the
Python output exactly: 45,458 five-minute and 225,868 one-minute values, max difference 2×10⁻⁷ (float32). The
Python unit tests (`test_sw_and_lw_window_minutes`, `test_window_mean_...`) also pass. The MATLAB default bins on
minute START times, so 5-min values differ slightly from the Python ones.

## What the comparison with Dan Lubin's tables shows
- **Alignment:** Dan's column *h* is the PST hour [h, h+1): r = 0.958 at shift 0, versus 0.85 at +1 h and 0.91 for half-hour bins.
- **February file:** rows 1–14 = 15–28 Feb 2023 and rows 15–28 = 1–14 Feb 2024 (r = 0.97). The Python
  `lubin_sw_series` assumption (row = day of month) gives r = 0.03, so that assumption is wrong.
- **Agreement:** 3,741 matched hours, r = 0.970, Dan = 0.86 × mine + 0.03, median ratio Dan/mine = 0.91.
- **Cause = the clear-sky model:** on 149 clear-sky hours mine is 1.000 and Dan's is 0.897 (IQR 0.86–0.93), so his
  clear-sky irradiance (Atwater & Ball 1981, per Kavin's note) is ~11 % brighter than RADFLUX's (Long & Ackerman
  2000). The ratio falls with sun angle, from 0.95 at SZA < 30° to 0.83 at 70–80°.

## Data and references
- ARM RADFLUX (radflux1long c2, EPCAPE M1); cite with the DOI on its ARM Data Discovery page.
- Long, C. N., and T. P. Ackerman, 2000: *J. Geophys. Res.*, 105, 15609–15626, doi:10.1029/2000JD900077.
- Kopp, G., and J. L. Lean, 2011: *Geophys. Res. Lett.*, 38, L01706, doi:10.1029/2010GL045777 (clearness index).
- Dan Lubin's `<MON>_SW_Boxplot.csv` tables: `VS_CODE/Python-Research/EPCAPE/Kavin_data_for_BAMS_paper` (8 months)
  and `EPCAPE/Kavin_code_for_BAMS/Annual Echo Plot` (all 12; APR byte-identical).
