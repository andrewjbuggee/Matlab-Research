%% Clearness index = measured SW down / top-of-atmosphere SW down, averaged over fixed windows
%
% DESCRIPTION:
%   MATLAB version of sw_clearness_index in the Python module
%   comparisons/seasonal_averages/quantities.py. Not the transmittance used
%   in the EPCAPE table: a sensitivity test. It divides by the irradiance at
%   the TOP of the atmosphere instead of a clear-sky surface value, so it
%   also includes the clear atmosphere's own extinction (Rayleigh
%   scattering, water vapour and aerosol absorption, ~20-30 % at noon). A
%   clear-sky ratio is bounded between this index and 1; it brackets how
%   bright a clear-sky denominator would have to be to reproduce a given
%   transmittance.
%
%   TOA irradiance on a horizontal surface:
%       F_TOA = S0 * (1 + 0.033 cos(2 pi n / 365)) * cos(SZA)
%   with S0 = 1361 W m-2 (Kopp & Lean 2011) and n the day of the year; the
%   cosine term approximates the Earth-Sun distance correction (Duffie &
%   Beckman, Solar Engineering of Thermal Processes, Eq. 1.4.1a).
%
%   Samples: SW down with QC = 0 and cos(SZA) > cos(SzaMax_deg). The same
%   samples go into both window means, as in sw_transmittance.m.
%
% INPUTS:
%   rad - struct from read_radflux_sw.m
%   Optional name-value pairs: 'SzaMax_deg' (80), 'WindowMin' (5),
%   'TimeField' ('time' = minute start, or 'time_end' as in Python)
%
% OUTPUTS:
%   K     - clearness index in each window
%   t_win - window start times (datetime column, UTC)
%
% REFERENCES:
%   Kopp, G., and J. L. Lean, 2011: A new, lower value of total solar
%     irradiance: Evidence and climate significance. Geophys. Res. Lett., 38,
%     L01706, doi:10.1029/2010GL045777.
%   Duffie, J. A., and W. A. Beckman: Solar Engineering of Thermal
%     Processes, Wiley, Eq. 1.4.1a (as cited in the Python original).
%   Python original: quantities.py (sw_clearness_index, SOLAR_CONSTANT_WM2).
%
% By Andrew John Buggee
%%

function [K, t_win] = sw_clearness_index(rad, varargin)


p = inputParser;
addParameter(p, 'SzaMax_deg', 80);
addParameter(p, 'WindowMin', 5);
addParameter(p, 'TimeField', 'time', @(x) any(strcmp(x, {'time', 'time_end'})));
parse(p, varargin{:});
opt = p.Results;

S0_Wm2 = 1361;                                      % total solar irradiance at 1 AU (Kopp & Lean 2011)

t   = rad.(opt.TimeField);
mu0 = rad.cos_sza;
doy = day(rad.time_end, 'dayofyear');              % Python uses the file's own time stamps for n
toa_Wm2 = S0_Wm2 * (1 + 0.033 * cos(2*pi*doy/365)) .* mu0;

% same samples in numerator and denominator: SW QC = 0 and the sun above the cut
ok = isfinite(rad.sw_Wm2) & (mu0 > cosd(opt.SzaMax_deg));
sw = rad.sw_Wm2;  sw(~ok) = NaN;
toa_Wm2(~ok) = NaN;

[t_win, sw_avg]  = window_mean(t, sw, opt.WindowMin);
[~,     toa_avg] = window_mean(t, toa_Wm2, opt.WindowMin);

K = sw_avg ./ toa_avg;
K(~(toa_avg > 0)) = NaN;


end
