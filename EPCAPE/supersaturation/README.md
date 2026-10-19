# EPCAPE cloud supersaturation: Jeramy Dedrick's retrievals in MATLAB

This is the MATLAB version of the Python notebook
`VS_CODE/Python-Research/EPCAPE/comparisons/supersaturation/dedrick2026_supersaturation_retrievals.ipynb`.
It explains and re-derives the effective cloud supersaturations of Dedrick et al. (2026, *J. Atmos. Sci.* 83, 1061,
doi:10.1175/JAS-D-26-0004.1) and recreates their Fig. 10. It gives the same numbers as the Python notebook.
That was checked on 19 Oct 2026 in R2024b: every printed statistic and table is identical.

## Credits: this work is based on Jeramy L. Dedrick's methods, data and code

The retrievals, data and parts of the code here are Jeramy Dedrick's (Scripps Institution of Oceanography). Please
cite his paper, his data and his code in anything that uses this folder:

- **Method:** Dedrick, J. L., and coauthors, 2026: Quantitative estimates of upwind and on-site particle scavenging
  and retrieved cloud supersaturations from urban coastal stratocumulus during EPCAPE. *J. Atmos. Sci.*, 83,
  1061–1082, doi:10.1175/JAS-D-26-0004.1 (supplement doi:10.1175/JAS-D-26-0004.s1).
- **Data:** his posted products, "Aerosol Quantification Products at Mt. Soledad and Scripps Pier during EPCAPE
  2023-24" (Maneenoi, Berta, Osawa, Han, Dedrick and Russell), UC San Diego Library, doi:10.6075/J0WM1DSC.
- **Code (MATLAB), GitHub user jdedrick95:**
  - `kappakohler` (2024), https://github.com/jdedrick95/kappakohler. `kohler_constants.m`,
    `kohler_saturation_ratio.m`, `critical_supersaturation_pct.m` and `critical_diameter_m.m` are based on his
    `kappakohler.m`: the same constants and the same κ-Köhler equations (his options 4, 1 and 2).
  - `ParcelActivation_ARG2000` (2025), https://github.com/jdedrick95/ParcelActivation_ARG2000.
    `smax_posted_arg2000_pct.m` is a line-by-line port of his `ParcelActivation_ARG2000.m`.
    `smax_arg2000_pct.m` follows its structure and constants, with the corrections described in the live script
    (Part 7).
  - `EPCAPE_Pier_CCNclosure_and_correction.m` and `EPCAPE_Pier_kappa4CCN.mat` (in the Library item; local copy in
    `../EPCAPE_Pier_CCNclosure_and_correction/`). This is his CCN closure, which produced the 1.57 correction
    (1.569 exactly) in the posted CCN data. It is not used here yet. Its `.mat` file holds the uncorrected Pier
    CCN, the 30-min Pier κ (mean 0.37 ± 0.14) and the measured SMPS size distributions.
- **Reconstructed here, not his code:** Jeramy has not posted the code for the two CCN-based retrievals
  (`Corrected_PierCCN_ARM_EPCAPE` on GitHub is empty). `hourly_ccn_spectrum.m`, `extend_ccn_spectrum.m`,
  `retrieve_supersaturation_pct.m` and `number_above_dhm_cm3.m` were written from the method description in his
  paper (Sec. 2d, pp. 1065–1066). The details the paper leaves open, such as how missing CCN setpoints are filled,
  were chosen to reproduce his posted retrievals: r = 0.95 for s_Nd and r = 0.94–0.97 for s_CCN. Differences from
  his own code are possible.
- **Licence:** his repositories have no licence file, so by default he keeps all rights. Ask Jeramy before sharing
  code derived from his functions outside the group, or before publishing it.

## Files

| File | What |
|---|---|
| `dedrick2026_supersaturation_retrievals.mlx` | **The live script.** Open it in MATLAB and press Run (~20 s). It is 
|   saved with its outputs, so you can read it without running it.
| `plain_text/dedrick2026_supersaturation_retrievals_plain.m` | The same script as plain text: easy to diff and 
|  run with `-batch`. The `.mlx` was generated from it with `matlab.internal.liveeditor.openAndSave`. Edit one
|  or the other, and regenerate if you want them in sync. 
| `functions/` | One function per step, each with comments on the physics and units 
| `figures/dedrick2026_fig10_recreated.png` | Fig. 10 recreation

### Functions
- `read_posted_product.m`: reads Jeramy's tab-separated Library files (MATLAB datenum times → timetable)
- `read_fm120_daily_files.m`: Mt. Soledad FM-120 5-min files (copied from `table_for_BAMS/readers`)
- `kohler_constants.m`, `kohler_saturation_ratio.m`, `critical_supersaturation_pct.m`, `critical_diameter_m.m`: κ-Köhler theory (Petters & Kreidenweis 2007), constants as in Jeramy's `kappakohler.m`
- `hourly_ccn_spectrum.m`: hourly Pier CCN spectra (hourly mean per setpoint, gaps of ≤ 2 h filled in time)
- `extend_ccn_spectrum.m`, `retrieve_supersaturation_pct.m`: the CCN-spectrum lookup behind s_Nd and s_CCN (spectrum extended to 0 % and 1.2 %)
- `number_above_dhm_cm3.m`: particles larger than the Hoppel minimum, from the lognormal mode fits
- `smax_arg2000_pct.m`: textbook Abdul-Razzak & Ghan (2000) parcel S_max
- `smax_posted_arg2000_pct.m`: line-by-line port of Jeramy's posted `ParcelActivation_ARG2000.m`, issues included (for comparison only)
- `compare_to_posted.m`: r, median difference and agreement with the posted values
- `plot_fig10_recreation.m`: the Fig. 10 recreation
- `new_figure.m`: a figure with MATLAB's factory graphics defaults (see "Figures" below)
- `ccn_spectrum_at_hour.m`, `mean_snd_with_scaled_ccn01.m`: the hourly spectrum at one hour, and the s_Nd mean with a rescaled 0.1 % CCN count (Sourita's correction)

## Physics sources
Physics statements in the live script and in the function headers cite their sources with page (and equation) numbers. The main source is Rogers & Yau (1989), *A Short Course in Cloud Physics*, 3rd ed. (printed page numbers). Others are Petters & Kreidenweis (2007), Dedrick et al. (2026) and its supplement, Abdul-Razzak et al. (1998), and Abdul-Razzak & Ghan (2000). The full bibliography is the last section of the `.mlx`. Two short physics checks use the textbook directly:
- Twomey's peak supersaturation from the measured CCN spectrum (after Part 3).
- The parcel terms α and γ against Rogers & Yau's Fig. 7.2 (Part 7).

## Data (set the two folders in the script's Settings section)
- `products_folder`: Jeramy's "Aerosol Quantification Products" (UC San Diego Library, doi:10.6075/J0WM1DSC). The files used are `bb2743661p_34_1_supersaturation.txt`, `bb2743661p_35_1_CCN.txt`, `bb2743661p_32_1_PSD_fits.txt` and `bb2743661p_30_1_lowCloudPeriods.txt`.
- `fm120_folder`: the FM-120 daily 5-min CSV files (Dalhousie).

**Known issue in the posted data:** in `EPCAPE_Supersaturation.txt` (`bb2743661p_34_1_supersaturation.txt`), the `skappa_*` and `sCCN_*` columns are swapped. The script corrects the labels and shows the evidence (Part 2).

## Figures, and why MATLAB hung
- **Hangs.** Exporting a *visible* figure with `exportgraphics` has hung R2024b for more than 10 min with MATLAB idle. It happened in `-batch` runs (this script, 19 Oct 2026) and in the Live Editor (EPCAPE BAMS Fig. 7, Oct 2026). The same figure drawn invisibly saves in about 1 s. So the Fig. 10 function draws on an **invisible** figure, saves the PNG and returns. The script closes it and shows the saved PNG, which is a single image object.
- **Dark theme.** `startup.m` sets dark root graphics defaults (black background, white 18-pt bold text). These plots use light and dark greys, so `new_figure.m` gives each figure the **factory** defaults at the figure level. Other figures and the root defaults are left untouched.
- **Testing a live script.** A plain `run()` in `-batch` does not test Live Editor behaviour. Run `matlab.internal.liveeditor.executeAndSave` on a *copy* of the `.mlx`, in the background, with a time limit.

Base MATLAB only (R2024b or later); no toolboxes. `boxchart` is used instead of `boxplot` (Statistics Toolbox).
