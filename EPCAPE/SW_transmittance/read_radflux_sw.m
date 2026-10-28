%% Read the shortwave fields of the ARM RADFLUX product (Scripps Pier, M1) for the transmittance
%
% DESCRIPTION:
%   Reads the four RADFLUX variables that the SW transmittance needs, keeps
%   only best-quality samples, and returns them as column vectors.
%
%   THE PRODUCT
%   RADFLUX ("Radiative Flux Analysis", ARM product radflux1long, here the
%   c2 file epcradflux1longM1.c2 combined over 15 Feb 2023 - 14 Feb 2024)
%   holds 1-min means of the SKYRAD broadband radiometers at the Pier and,
%   from the Long & Ackerman (2000) method, an estimate of what the
%   irradiance would have been under a clear sky:
%     downwelling_shortwave           measured global SW down (W m-2),
%                                     from the sum of direct + diffuse or the
%                                     global pyranometer ("cell_methods: time:
%                                     mean": a 1-min mean of 1-s samples)
%     clearsky_downwelling_shortwave  clear-sky SW down (W m-2): an empirical
%                                     fit, as a power law of cos(SZA), to this
%                                     site's own detected clear-sky periods
%                                     (Long & Ackerman 2000, Eq. 1), so it
%                                     includes the site's typical aerosol and
%                                     water vapour. No QC field.
%     cosine_zenith                   cosine of the solar zenith angle (point
%                                     value, from the PNNL solarposition code)
%
%   QUALITY CONTROL (the same rule as qc0_values in the Python
%   comparisons/seasonal_averages/quantities.py): a sample is kept only if
%   its qc_<variable> is exactly 0, i.e. it passed every test, including the
%   tests ARM rates "Indeterminate". For downwelling_shortwave bit 3 is
%   "Indeterminate": the value was estimated from the direct and diffuse
%   components because the global measurement failed QCRAD; QC = 0 drops it.
%   Variables without a qc_ field (clear-sky SW, cos SZA) are used as they
%   are. ncread turns the _FillValue (NaN here) into NaN.
%
%   TIME STAMPS
%   The 1-min means are stamped at the END of their minute (ARM RADFLUX
%   time_bounds = [-60 s, 0], per the original b1/c2 files; the Python
%   docstring of sw_transmittance says the same). This function returns
%   both the file's stamps (rad.time_end) and the START of each minute
%   (rad.time = rad.time_end - 60 s). Every other EPCAPE product used here
%   labels intervals by their start, so rad.time is the one to use when
%   averaging into windows or matching other data.
%
% INPUTS:
%   radflux_file - path to the combined RADFLUX netCDF file, e.g.
%                  .../data/processed/radflux_M1_20230215_20240214.nc
%
% OUTPUTS:
%   rad - struct of column vectors (one element per 1-min sample)
%       time            start of each minute (datetime, UTC)
%       time_end        the file's time stamp = end of the minute
%       sw_Wm2          measured SW down, QC = 0 only (NaN otherwise)
%       sw_clear_Wm2    clear-sky SW down (no QC field)
%       cos_sza         cosine of the solar zenith angle
%       qc_info         struct: n samples, n finite, n kept after QC = 0
%
% REFERENCES:
%   Long, C. N., and T. P. Ackerman, 2000: Identification of clear skies
%     from broadband pyranometer measurements and calculation of downwelling
%     shortwave cloud effects. J. Geophys. Res., 105, 15609-15626,
%     doi:10.1029/2000JD900077.
%   Long, C. N., and D. D. Turner, 2008: A method for continuous estimation
%     of clear-sky downwelling longwave radiative flux developed using ARM
%     surface measurements. J. Geophys. Res., 113, D18206,
%     doi:10.1029/2008JD009936.
%   ARM data: Radiative Flux Analysis (radflux1long, c2), EPCAPE M1,
%     Atmospheric Radiation Measurement (ARM) User Facility. Cite it with the
%     DOI given on the product's ARM Data Discovery page (not checked here).
%   Python original: VS_CODE/Python-Research/EPCAPE/comparisons/
%     seasonal_averages/quantities.py (qc0_values, sw_transmittance).
%
% By Andrew John Buggee
%%

function rad = read_radflux_sw(radflux_file)


% -------------------------------------------------------------------
% ----------------------------- Time --------------------------------
% -------------------------------------------------------------------

% The combined file stores "seconds since 1970-01-01" (UTC)
units = string(ncreadatt(radflux_file, 'time', 'units'));
if ~startsWith(units, "seconds since 1970-01-01")
    error('read_radflux_sw: unexpected time units "%s" in %s', units, radflux_file);
end
rad.time_end = datetime(double(ncread(radflux_file, 'time')), 'ConvertFrom', 'posixtime');
rad.time_end = rad.time_end(:);
rad.time = rad.time_end - seconds(60);    % start of each 1-min averaging interval



% -------------------------------------------------------------------
% ---------- Measured SW down, best-quality samples only ------------
% -------------------------------------------------------------------

sw = double(ncread(radflux_file, 'downwelling_shortwave'));
qc = double(ncread(radflux_file, 'qc_downwelling_shortwave'));
n_finite = nnz(isfinite(sw));
sw(qc ~= 0) = NaN;                          % keep QC == 0 only (passed every test)
rad.sw_Wm2 = sw(:);

rad.qc_info = struct('variable', 'downwelling_shortwave', 'n', numel(sw), ...
                     'n_finite', n_finite, 'n_kept', nnz(isfinite(sw)));



% -------------------------------------------------------------------
% -------- Clear-sky SW down and the solar zenith angle -------------
% -------------------------------------------------------------------

% Neither has a qc_ companion in RADFLUX, so they are used as distributed
rad.sw_clear_Wm2 = double(ncread(radflux_file, 'clearsky_downwelling_shortwave'));
rad.sw_clear_Wm2 = rad.sw_clear_Wm2(:);
rad.cos_sza = double(ncread(radflux_file, 'cosine_zenith'));
rad.cos_sza = rad.cos_sza(:);


end
