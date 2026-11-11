%% Read the ARM NDROP (NDROPMFRSR) droplet number concentration and average it onto a 1-min grid
%
% DESCRIPTION:
%   The ARM Droplet Number Concentration value-added product (NDROP; Zhang
%   et al. 2023), datastream epcndropmfrsrM1.c1 at the Scripps Pier,
%   retrieves the cloud-layer droplet number concentration N_d every 20 s
%   from the MFRSR cloud optical depth, the MWR liquid water path and the
%   cloud base and top (adiabatic cloud model). It only exists for daytime
%   overcast, single-layer, warm, liquid cloud (Dedrick et al. 2026,
%   supplement Text S5): the retrieval needs the MFRSR optical depth, which
%   needs sunlight (cos SZA >= 0.2) and cloud fraction >= 0.9.
%
%   This function reads every daily file in the folder and returns N_d on
%   the 1-min grid t_min (minute starts, UTC): the mean of the 20-s values
%   in [t, t + 1 min) that pass QC, in cm^-3.
%
%   QC (qc_drop_number_conc, bit-packed; descriptions from the file):
%     "strict"  keep only qc == 0 (no test failed or indeterminate)
%     "bad"     reject only bits ARM assesses as Bad (1 no optical depth,
%               2 LWP < 0.02 kg m-2, 3 no cloud top, 4 cloud base < 260 K,
%               6/7 cloud-base T or pressure out of range) and bit 9
%               (N_d > valid_max, "not physically reasonable"); keep the
%               Indeterminate bits 5 (default cloud base used), 8 (cloud-
%               base QC indeterminate) and 10 (value < valid_min, reset to 0)
%
% INPUTS:
%   folder  - folder with epcndropmfrsrM1.c1.*.nc
%   t_min   - the 1-min grid (datetime column, minute starts, UTC)
%   qc_mode - "strict" (default) or "bad"
%
% OUTPUTS:
%   nd_cm3 - N_d on t_min (cm^-3; NaN where no QC-passing 20-s value)
%   n_nd   - number of 20-s values behind each minute (0-3)
%
% REFERENCES:
%   Zhang, D., and coauthors, 2023: Evaluation of four ground-based
%     retrievals of cloud droplet number concentration in marine
%     stratocumulus with aircraft in situ measurements. Atmos. Meas. Tech.,
%     16, 5827-5846, doi:10.5194/amt-16-5827-2023 (the NDROP reference
%     cited by Dedrick et al. 2026).
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1, supplement Text S5.
%
% By Andrew John Buggee
%%

function [nd_cm3, n_nd] = read_ndrop_1min(folder, t_min, qc_mode)

if nargin < 3 || isempty(qc_mode)
    qc_mode = "strict";
end
bad_bits = [1 2 3 4 6 7 9];                       % used by qc_mode "bad"

f = dir(fullfile(folder, 'epcndropmfrsrM1.c1.*.nc'));
assert(~isempty(f), 'read_ndrop_1min: no NDROP files in %s', folder);

n_min = numel(t_min);
sum_nd = zeros(n_min, 1);
n_nd   = zeros(n_min, 1);
t0 = t_min(1);

for k = 1:numel(f)
    file = fullfile(f(k).folder, f(k).name);

    % time = base_time (s since 1970, UTC) + time_offset (s)
    base = double(ncread(file, 'base_time'));
    off  = double(ncread(file, 'time_offset'));
    t = datetime(base + off, 'ConvertFrom', 'posixtime', 'TimeZone', '');

    nd = double(ncread(file, 'drop_number_conc')) / 1e6;     % m-3 -> cm-3
    qc = double(ncread(file, 'qc_drop_number_conc'));

    switch qc_mode
        case "strict"
            ok = qc == 0;
        case "bad"
            ok = true(size(qc));
            for b = bad_bits
                ok = ok & ~bitget(qc, b);
            end
        otherwise
            error('read_ndrop_1min: qc_mode must be "strict" or "bad"');
    end
    ok = ok & isfinite(nd) & nd >= 0 & nd < 1e5;             % also drops the -9999 missing value

    % minute of each sample on the grid
    i = floor(minutes(t - t0)) + 1;
    ok = ok & i >= 1 & i <= n_min;
    sum_nd = sum_nd + accumarray(i(ok), nd(ok), [n_min 1]);
    n_nd   = n_nd   + accumarray(i(ok), 1,      [n_min 1]);
end

nd_cm3 = sum_nd ./ n_nd;                                     % NaN where n_nd = 0

end
