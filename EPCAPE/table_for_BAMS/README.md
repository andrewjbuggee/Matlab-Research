# table_for_BAMS: re-deriving Kavin's rows of the EPCAPE seasonal-averages table

This MATLAB code recomputes every "Andrew Buggee" line that checks a row of
`Table-EPCAPE_SeasonalAverages.xlsx` (sheet *Cloud&MetQuantities*) attributed
to Kavin. It is a translation of Andrew's Python notebook
(`EPCAPE/comparisons/seasonal_averages/check_seasonal_averages.ipynb`), reduced
to the calculations that appear in the table.

## How to run

1. upzip data folder `table_for_BAMS_data_complete.zip` from Google Drive: 
ScrippsAerosolResearchGroup@UCSD/Andrew Buggge/EPCAPE/BAMS paper. (see its `README.txt`).
2. Open `make_table_for_BAMS.m` and edit the folder paths in Section 1:
   * `data_folder`: the `table_for_BAMS_data` folder;
   * `kavin_folder`: Kavin's folder with Lubin's `*_SW_Boxplot.csv`,
     `EPCAPE_15mindata.nc`, `ET_Dcrit_nativetime.csv` and
     `epcmwrret_turnM1_hourly_averages.nc`.
3. Press **Run**. It takes about 1 minute and uses ~4-5 GB of memory (the
   MWRRET v2 and visibility records hold 22 and 25 million samples).
4. Results are written to `output/`:
   * `table_for_BAMS_MATLAB.xlsx`, sheet **Table**: one line per new table
     row, laid out like the sheet (columns A-D, the five season cells,
     author, notes);
   * sheet **Numbers** (and `table_for_BAMS_MATLAB_numbers.csv`): the
     unrounded mean, std and sample count of every cell.

Requires MATLAB R2024b or later, base MATLAB only (no toolboxes).

## Which file computes which sheet row

| Sheet row(s) | Quantity | File in `calculations/` |
|---|---|---|
| 7 | SW transmittance | `calc_row07_SW_transmittance.m` |
| 8 | LW down | `calc_row08_LW_down.m` |
| 10 | Liquid water path | `calc_row10_LWP.m` |
| 15 | Cloud-top height | `calc_row15_cloud_top.m` |
| 17 | PBL height | `calc_row17_PBL_height.m` |
| 20 | LCL | `calc_row20_LCL.m` |
| 21 | Cloud-base height | `calc_row21_cloud_base.m` |
| 22-24 | Rain | `calc_rows22to24_rain.m` |
| 25, 26, 28, 29, 33-42, 54, 55 | Kavin's method on Abbey's 15-min file | `calc_rows_15min_file.m` |
| 25, 26, 28, 29 | FM-120 5-min rows with Abbey's classes | `calc_rows25to29_FM120_5min.m` |
| 27, 30 | Activated fraction | `calc_rows27_30_activated_fraction.m` |
| 31, 32 | FM-120 cloud and haze hours | `calc_rows31_32_FM120_hours.m` |

Rows 18 (PBLH-MPL) and 19 (BL-Depth Miller) are not computed: no data
source for them was found.

Every function starts with a comment block saying which data it reads,
which samples it keeps, and how it averages. The same text, shortened,
is in the Notes column of the output.

## Conventions used everywhere

* **Seasons:** EPCAPE = 2023-02-15 to 2024-02-14; Spring = Apr-Jun 2023;
  Summer = Jul-Sep 2023; Fall = Oct-Dec 2023; Winter = 15 Feb-Mar 2023 plus
  Jan-14 Feb 2024 (UTC), as in the sheet's headers (`helpers/sheet_season_masks.m`).
  Where Kavin's own method is re-run, his calendar quarters are used instead
  (winter = Jan-Mar; `helpers/calendar_quarter_masks.m`). The Notes say which.
* **Statistics:** mean ± sample standard deviation over all samples in the
  season (`helpers/season_stats.m`). The campaign value is therefore not the
  mean of the four seasons.
* **Quality control:** an ARM variable with a `qc_<variable>` field keeps
  only samples with QC = 0 (`readers/read_arm_qc0.m`).
* **Rounding:** the std to 2 significant figures, and the mean to the same
  decimal place. Hours and fractions get 3 significant figures
  (`helpers/format_value.m`).

## Check against the Python version

On 2026-10-08 the output was compared with the Python notebook (run in
MATLAB R2024b):
* all 210 cells (42 lines × 5 seasons) are identical to the pink rows of
  `Table-EPCAPE_SeasonalAverages - AJB verification - ver 2.xlsx`;
* the unrounded means and stds agree to within 1e-6 (relative), and every
  sample count is the same. The tiny differences come from Python averaging
  some single-precision values in single precision.
