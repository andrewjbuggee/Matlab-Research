function C = ceilometer_30min_cbh(file_path)
% CEILOMETER_30MIN_CBH  30-min mean cloud-base height from the ARM ceilometer (CL31) at the Scripps Pier.
%
%   C = ceilometer_30min_cbh(file_path)
%
%   Reads the combined ceilometer file made by the Python tool EPCAPE/download_data
%   (product cbh_ceil_M1: data/processed/cbh_ceil_M1_20230215_20240214.nc, one
%   sample every ~16 s) and averages the lowest cloud base, first_cbh, over 30-min
%   windows [t, t + 30 min), the resolution of s_ARG00 (Dedrick et al. 2026,
%   supplement Table S7). Dedrick et al. (2026, p. 1065) take the cloud-base
%   thermodynamics from INTERPSONDE "at the height level of the ceilometer-
%   determined cloud-base height"; how they averaged the cloud base in time is not
%   stated, so the 30-min mean of the detected bases is my choice.
%
%   A sample counts as a cloud base when detection_status is 1, 2 or 3 (one, two
%   or three cloud bases detected; 0 = no significant backscatter, 4 = full
%   obscuration, 5 = some obscuration, from the file's attributes / ARM's CL31
%   handbook) and first_cbh has qc_first_cbh = 0.
%
%   Heights: first_cbh is above the instrument, not sea level (ARM's ceilometer
%   handbook, as quoted in the Python tool's README); alt_m is returned so the
%   caller can add the instrument altitude.
%
%   Output struct
%     C.time      start of each 30-min window (datetime column, UTC)
%     C.cbh_m     mean first_cbh of the detected samples in the window, m above the instrument (NaN if none)
%     C.frac      fraction of the window's samples with a detected cloud base (0-1)
%     C.alt_m     instrument altitude, m above mean sea level

t = datetime(double(ncread(file_path, 'time')), 'ConvertFrom', 'posixtime');   % UTC
cbh = double(ncread(file_path, 'first_cbh'));                                  % m above the instrument
qc = double(ncread(file_path, 'qc_first_cbh'));
status = double(ncread(file_path, 'detection_status'));
C.alt_m = double(ncread(file_path, 'alt'));

detected = ismember(status, [1 2 3]) & isfinite(cbh) & qc == 0;

% 30-min window of each sample, counted from the first window
t0 = dateshift(t(1), 'start', 'hour');
w = floor(minutes(t - t0) / 30) + 1;
n_w = max(w);
n_all = accumarray(w, 1, [n_w 1]);
n_det = accumarray(w(detected), 1, [n_w 1]);
s_det = accumarray(w(detected), cbh(detected), [n_w 1]);

C.time = t0 + minutes(30 * (0:n_w - 1)');
C.cbh_m = s_det ./ n_det;            % NaN where nothing was detected (0/0)
C.frac = n_det ./ max(n_all, 1);

end
