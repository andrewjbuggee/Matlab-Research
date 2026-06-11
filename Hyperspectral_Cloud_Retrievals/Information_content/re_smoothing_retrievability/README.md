# Retrievability of vertical droplet-size structure (HySICS)

Tests whether HySICS can detect the high-frequency vertical fluctuations in the
in-situ measured cloud droplet effective-radius profiles from VOCALS-REx (73
profiles) and ORACLES (237 profiles). For each profile we compute the 636-band
top-of-atmosphere nadir reflectance spectrum for the **raw** r_e(z) profile and
several **vertically smoothed** versions, holding the cloud optical depth fixed,
then compare the reflectance differences against the HySICS measurement
uncertainty (0.3%). If a difference is below the uncertainty, that vertical
structure is not retrievable.

## The controlled variable

Only the **vertical shape of r_e(z)** changes between versions. We hold the
cloud optical depth `tau_c` fixed (libRadtran `wc_modify tau set`) and keep the
measured LWC(z) shape; its magnitude is rescaled internally to hit `tau_c`, so
the per-level LWC is free to follow the droplet size. We do NOT hold LWP fixed —
when you smooth r_e you can hold tau OR LWP, not both, and LWP is not an
independent radiative-transfer input. Holding tau fixed makes the VIS bands
nearly identical across versions and isolates the r_e-shape effect in the SWIR.

We do **not** use `write_wc_file.m` for this: in its vertically-resolved path it
imposes a linearly-increasing (adiabatic-like) LWC profile and resamples r_e
onto a uniform altitude grid, which would inject two confounds. Instead we use
`write_wc_file_from_in_situ` (keeps the measured LWC(z)) + `wc_modify tau set`.

## Files

| File | Role | Where to run |
|------|------|--------------|
| `define_re_smoothing_windows.m` | Single source of truth for the moving-average windows (meters). | — |
| `smooth_re_profile.m` | Applies the vertical moving average (movmean in meters). Shared by preview and RT. | — |
| `preview_re_smoothing.m` | **Run first.** Visualizes raw vs smoothed r_e(z) and quantifies how much each window removes. Tweak windows, repeat. | Local (Mac) |
| `hysics_refl_smoothingTest_from_insitu.m` | Core: computes the spectra for one profile (raw + smoothed). | Alpine (via batch) |
| `../../Batch_Scripts/bash_scripts/re_smoothing_retrievability/smoothingTest_VR.sh` | SLURM array over the 73 VOCALS profiles. | Alpine |
| `../../Batch_Scripts/bash_scripts/re_smoothing_retrievability/smoothingTest_OR.sh` | SLURM array over the 237 ORACLES profiles. | Alpine |
| `plot_refl_diff_smoothing_aggregate.m` | Aggregate figure: median \|ΔR\| + percentile bands vs wavelength, with the 0.3% uncertainty envelope. | Local (Mac) |

## Workflow

1. `preview_re_smoothing` — inspect the smoothing; edit `define_re_smoothing_windows.m` until happy.
2. Edit the two `.sh` scripts: set `SZA`, the `output_dir` date stamp, and `--array`/`%` concurrency. Then `sbatch smoothingTest_VR.sh` and `sbatch smoothingTest_OR.sh`.
3. Copy the resulting `refl_smoothingTest_*.mat` files locally.
4. `plot_refl_diff_smoothing_aggregate('/path/to/results/')`.

## Geometry

Single fixed geometry: nadir view (vza = 0, vaz = 0), solar azimuth 0, and a
fixed SZA (default 30°, set in the batch scripts). Per profile this is
636 channels × (1 + number of windows) uvspec runs.
