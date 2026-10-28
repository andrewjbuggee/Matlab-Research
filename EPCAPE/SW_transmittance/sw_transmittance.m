%% Shortwave transmittance = measured SW down / clear-sky SW down, averaged over fixed windows
%
% DESCRIPTION:
%   MATLAB version of sw_transmittance in the Python module
%   comparisons/seasonal_averages/quantities.py (Andrew's method for the
%   EPCAPE seasonal-averages table, row 7).
%
%   DEFINITION
%       T_SW = < SW_down > / < SW_down,clear >
%   Both irradiances come from ARM RADFLUX at the Pier (read_radflux_sw.m):
%   the measured 1-min global SW down and the Long & Ackerman (2000)
%   clear-sky estimate. Each is averaged over the SAME valid samples in a
%   window of window_min minutes, and the ratio of the two window means is
%   the transmittance. Because the clear-sky value already accounts for the
%   sun's elevation, the Earth-Sun distance and the site's clear-sky
%   atmosphere, T_SW isolates the effect of cloud: ~1 under a clear sky and
%   well below 1 under stratocumulus. RADFLUX itself carries the same
%   1-min quantity as cloud_transmissivity_shortwave ("calculated as ratio
%   of downwelling_shortwave to clearsky_downwelling_shortwave", its file
%   comment); this function adds QC = 0, the sun-elevation cut and the
%   window averaging. Values > 1 (cloud-edge enhancement of the diffuse
%   field) are kept.
%
%   WHICH SAMPLES ARE USED
%     1. measured SW down passed every QC test (QC = 0; read_radflux_sw.m)
%     2. the clear-sky estimate is finite and > 0
%     3. the sun is higher than 90 - SzaMax_deg degrees: cos(SZA) >
%        cos(SzaMax_deg). Near the horizon both irradiances are small and
%        the ratio is noisy, and a horizontal pyranometer's cosine response
%        is least accurate there.
%     4. (MATLAB addition, not in the Python) the optional 'Keep' mask, so a
%        caller can apply further filters (cloud type, LWP, ...) BEFORE the
%        averaging. All four conditions select one set of samples, used in
%        both the numerator and the denominator, so the ratio is not biased
%        by unequal gaps.
%
%   AVERAGING
%   Windows come from window_mean.m (anchored at midnight, labelled by their
%   start). WindowMin = 0 means no averaging: one ratio per valid 1-min
%   sample. The Python default is 5 min, "Kavin's note: ratio of five-
%   minute-averaged downwelling shortwave irradiances".
%
%   A NOTE ON THE CLEAR-SKY MODEL
%   Kavin's table note names an "idealized clear sky analytical model
%   calculation (Atwater and Ball, 1981)", which is what Dan Lubin's tables
%   use. This function uses RADFLUX's empirical Long & Ackerman fit instead;
%   the Python docstring records that the Atwater & Ball equations could not
%   be checked from an accessible source. Different clear-sky models move
%   the ratio by a few percent or more; compare_sw_transmittance_with_lubin.m
%   measures how much for these two.
%
%   DIFFERENCE FROM THE PYTHON TIME BASE
%   The Python function bins the samples on the RADFLUX file's own time
%   stamps, which mark the END of each minute. Here the START of each minute
%   is used by default (rad.time), consistent with every other EPCAPE product
%   in this repository. Set 'TimeField', 'time_end' to reproduce the Python
%   numbers exactly. With 5-min windows the two differ by one 1-min sample
%   per window boundary.
%
% INPUTS:
%   rad - struct from read_radflux_sw.m (time, time_end, sw_Wm2,
%         sw_clear_Wm2, cos_sza)
%
%   Optional name-value pairs:
%   'SzaMax_deg'  largest solar zenith angle used (default 80, as in Python)
%   'WindowMin'   averaging window in minutes (default 5; 0 = no averaging)
%   'MinCount'    fewest valid 1-min samples a window needs (default 1, as
%                 in Python; hourly_ratio.m uses 30 for 1-h windows)
%   'Keep'        logical vector, one per sample; false = do not use
%                 (default: all true)
%   'TimeField'   'time' (default, minute start) or 'time_end' (as Python)
%
% OUTPUTS:
%   T_sw   - transmittance in each window (NaN where no valid samples)
%   t_win  - window start times (datetime column, UTC)
%   info   - struct: SzaMax_deg, WindowMin, MinCount, n_values (finite
%            ratios), n_samples_used (1-min samples entering any window),
%            qc (from read_radflux_sw), sw_mean_Wm2 and sw_clear_mean_Wm2
%            (the window means behind each ratio), n_per_window
%
% REFERENCES:
%   Long, C. N., and T. P. Ackerman, 2000: J. Geophys. Res., 105,
%     15609-15626, doi:10.1029/2000JD900077.
%   Python original: VS_CODE/Python-Research/EPCAPE/comparisons/
%     seasonal_averages/quantities.py (sw_transmittance) and its test
%     test_sw_and_lw_window_minutes in tests/test_seasonal_averages.py.
%
% By Andrew John Buggee
%%

function [T_sw, t_win, info] = sw_transmittance(rad, varargin)


% -------------------------------------------------------------------
% ------------------------ Parse the inputs -------------------------
% -------------------------------------------------------------------

p = inputParser;
p.FunctionName = 'sw_transmittance';
addParameter(p, 'SzaMax_deg', 80, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 90);
addParameter(p, 'WindowMin', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'MinCount', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'Keep', [], @(x) isempty(x) || islogical(x));
addParameter(p, 'TimeField', 'time', @(x) any(strcmp(x, {'time', 'time_end'})));
parse(p, varargin{:});
opt = p.Results;

t   = rad.(opt.TimeField);
sw  = rad.sw_Wm2;
cs  = rad.sw_clear_Wm2;
mu0 = rad.cos_sza;

keep = true(size(sw));
if ~isempty(opt.Keep)
    if numel(opt.Keep) ~= numel(sw)
        error('sw_transmittance: Keep must have one element per RADFLUX sample (%d)', numel(sw));
    end
    keep = opt.Keep(:);
end



% -------------------------------------------------------------------
% ------------- One set of samples for numerator and denominator ----
% -------------------------------------------------------------------

% daytime: sun higher than 90 - SzaMax_deg degrees
day = mu0 > cosd(opt.SzaMax_deg);

% the same samples go into numerator and denominator: SW passed QC = 0 (NaN
% otherwise), clear-sky estimate > 0, the sun above the cut, and the
% caller's mask
both = isfinite(sw) & isfinite(cs) & (cs > 0) & day & keep;

sw_used = sw;  sw_used(~both) = NaN;
cs_used = cs;  cs_used(~both) = NaN;



% -------------------------------------------------------------------
% ---------------- Window means and their ratio ---------------------
% -------------------------------------------------------------------

[t_win, sw_avg, n_sw] = window_mean(t, sw_used, opt.WindowMin, opt.MinCount);
[~,     cs_avg]       = window_mean(t, cs_used, opt.WindowMin, opt.MinCount);

T_sw = sw_avg ./ cs_avg;
T_sw(~(cs_avg > 0)) = NaN;                 % no valid clear-sky mean: no ratio



% -------------------------------------------------------------------
% ----------------------- What was done -----------------------------
% -------------------------------------------------------------------

info.SzaMax_deg        = opt.SzaMax_deg;
info.WindowMin         = opt.WindowMin;
info.MinCount          = opt.MinCount;
info.TimeField         = opt.TimeField;
info.n_values          = nnz(isfinite(T_sw));
info.n_samples_used    = nnz(both);
info.n_per_window      = n_sw;
info.sw_mean_Wm2       = sw_avg;
info.sw_clear_mean_Wm2 = cs_avg;
if isfield(rad, 'qc_info')
    info.qc = rad.qc_info;
end


end
