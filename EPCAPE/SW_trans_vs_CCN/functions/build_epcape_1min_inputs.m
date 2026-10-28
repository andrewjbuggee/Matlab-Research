%% Put every measurement the cloud-scene filters need on one 1-min grid (Scripps Pier and Mt. Soledad)
%
% DESCRIPTION:
%   The filters for the SW transmittance vs. CCN analysis (single-layer,
%   liquid, non-precipitating, persistent, daytime cloud with 20 < LWP <
%   250 g m-2) combine instruments with very different sampling: 1-s
%   visibility, ~1.3-s microwave radiometer, 4-s ARSCL, ~16-s ceilometer,
%   1-min RADFLUX and disdrometers, 5-min FM-120. This function averages
%   (or summarises) each onto the same 1-min grid, [t, t + 1 min) labelled by
%   its start, so every filter can be evaluated minute by minute and the
%   SW transmittance computed from exactly the minutes that pass.
%
%   WHAT EACH COLUMN IS
%   Pier (ARM M1, 8 m MSL)
%     sw_Wm2, sw_clear_Wm2, cos_sza  RADFLUX 1-min means (read_radflux_sw.m:
%                                    SW QC = 0; times moved to minute START)
%     lwp_gm2, n_lwp                 MWRRET v2 (mwrret2turn.c1) phys_lwp,
%                                    mean of the samples in the minute that
%                                    pass qc_phys_lwp = 0, phys_qc_flag = 0
%                                    ("retrievals_good") and are under a
%                                    detected cloud base (detection_status
%                                    1-3 and cbh_detected > 0), i.e. the
%                                    Python cloudy_mwrret2 + lwp_cloudy_gm2
%                                    rules; n_lwp = number of such samples
%     pier_cloud_frac                ceilometer (ceil.b1): fraction of the
%                                    minute's samples with a cloud base
%                                    (detection_status 1-3, first_cbh > 0)
%     pier_cbh_m                     ceilometer: median first_cbh (m above
%                                    ground) of those samples
%     pier_multi_base                ceilometer: 1 if any sample reported a
%                                    2nd or 3rd base (detection_status 2-3)
%     arscl_n, arscl_max_layers,     KAZR-ARSCL (arsclkazr1kollias.c1): number
%     arscl_top_max_m                of 4-s samples, the largest number of
%                                    hydrometeor layers in any of them, and
%                                    the highest layer top (m above ground).
%                                    Empty layer slots hold -1, not NaN; a
%                                    layer is present when its base is a
%                                    real height >= 0.
%     rain_pier                      laser disdrometer at the Pier (ld.b1): 1
%                                    if any sample in the minute had
%                                    precip_rate > 0 (qc = 0), 0 if none,
%                                    NaN if no data
%   Mt. Soledad (250 m MSL)
%     rain_soledad                   laser disdrometer at ARM site S2, as above
%     soledad_cloud                  1 if the minute lies in a 5-min interval
%                                    meeting Dedrick's cloud condition
%                                    (visibility < 1 km throughout and FM-120
%                                    LWC > 0.01 g m-3), 0 if not, NaN if the
%                                    visibility or FM-120 value is missing
%                                    (soledad_cloud_condition_5min.m). This
%                                    is "Mt. Soledad is INSIDE the cloud",
%                                    i.e. cloud base below ~250 m MSL.
%
% INPUTS:
%   files - struct with the paths:
%       radflux, mwrret2, ceil, arscl, rain_pier, rain_soledad  (netCDF)
%       visibility    string array of the Mt. Soledad visibility CSVs
%       fm120_folder  folder of the FM-120 daily 5-min CSVs
%   time_range - (optional) [start end) of the grid, default
%                [datetime(2023,2,15) datetime(2024,2,15)] (UTC)
%
% OUTPUTS:
%   D      - timetable, one row per minute (525,600 for the campaign)
%   cloud5 - the 5-min Mt. Soledad classification (soledad_cloud_condition_5min.m)
%   fm     - the FM-120 5-min record (read_fm120_daily_files.m)
%
% DEPENDENCIES (add to the path first):
%   EPCAPE/SW_transmittance            read_radflux_sw.m
%   EPCAPE/Dedrick_2026_coupled_cases  read_soledad_visibility.m, soledad_cloud_condition_5min.m
%   EPCAPE/supersaturation/functions   read_fm120_daily_files.m
%
% REFERENCES:
%   ARM products (EPCAPE, M1 and S2): RADFLUX (radflux1long c2), MWRRET v2
%     (mwrret2turn c1), ceilometer (ceil b1), KAZR-ARSCL (arsclkazr1kollias
%     c1), laser disdrometer (ld b1). Atmospheric Radiation Measurement
%     (ARM) User Facility; cite each with the DOI on its ARM Data Discovery
%     page.
%   Python rules copied: comparisons/seasonal_averages/quantities.py
%     (cloudy_mwrret2, lwp_cloudy_gm2, lowest_layer_top_m, near_rain).
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082
%     (supplement Table S4 for the Mt. Soledad cloud condition).
%
% By Andrew John Buggee
%%

function [D, cloud5, fm] = build_epcape_1min_inputs(files, time_range)


if nargin < 2 || isempty(time_range)
    time_range = [datetime(2023, 2, 15) datetime(2024, 2, 15)];
end

% the 1-min grid and a helper that gives the grid row of any time stamp
t0 = time_range(1);
n_min = round(minutes(time_range(2) - t0));
time = t0 + minutes(0:n_min - 1)';
row_of = @(t) floor(seconds(t - t0) / 60 + 1e-9) + 1;    % minute [t, t + 1 min) holds t
in_grid = @(k) k >= 1 & k <= n_min;

D = timetable(time);



% -------------------------------------------------------------------
% ------------------- RADFLUX shortwave (1 min) ---------------------
% -------------------------------------------------------------------

fprintf('  RADFLUX ...\n');
rad = read_radflux_sw(files.radflux);
k = row_of(rad.time);               % rad.time = minute START
ok = in_grid(k);
D.sw_Wm2       = nan(n_min, 1);  D.sw_Wm2(k(ok))       = rad.sw_Wm2(ok);
D.sw_clear_Wm2 = nan(n_min, 1);  D.sw_clear_Wm2(k(ok)) = rad.sw_clear_Wm2(ok);
D.cos_sza      = nan(n_min, 1);  D.cos_sza(k(ok))      = rad.cos_sza(ok);
clear rad



% -------------------------------------------------------------------
% ------------- MWRRET v2 liquid water path (cloudy, QC) ------------
% -------------------------------------------------------------------

fprintf('  MWRRET v2 LWP ...\n');
t      = nc_time(files.mwrret2);
lwp    = ncread(files.mwrret2, 'phys_lwp');          % g m-2 (units checked: "g/m^2")
qc     = ncread(files.mwrret2, 'qc_phys_lwp');
physqc = ncread(files.mwrret2, 'phys_qc_flag');      % 0 = retrievals_good
cbh_km = ncread(files.mwrret2, 'cbh_detected');      % km; -0.001 when no cloud
status = ncread(files.mwrret2, 'detection_status');  % ceilometer code copied by the VAP

% keep: QC = 0, retrieval good, and a cloud base detected (status 1-3, cbh > 0)
good = isfinite(lwp) & qc == 0 & physqc == 0 & status >= 1 & status <= 3 & cbh_km > 0;
k = row_of(t(good));
ok = in_grid(k);
v = double(lwp(good));
D.n_lwp   = accumarray(k(ok), 1,     [n_min 1], @sum, 0);
D.lwp_gm2 = accumarray(k(ok), v(ok), [n_min 1], @sum, 0) ./ D.n_lwp;   % 0/0 = NaN where no sample
clear t lwp qc physqc cbh_km status good v



% -------------------------------------------------------------------
% ------------------- Ceilometer at the Pier ------------------------
% -------------------------------------------------------------------

fprintf('  ceilometer ...\n');
t      = nc_time(files.ceil);
cbh_m  = double(ncread(files.ceil, 'first_cbh'));     % m above ground
status = double(ncread(files.ceil, 'detection_status'));

valid  = status >= 0 & status <= 4;                   % a real detection code (-9999 = missing)
cloudy = valid & status >= 1 & status <= 3 & cbh_m > 0;
multi  = valid & (status == 2 | status == 3);

k = row_of(t);
ok = in_grid(k) & valid;
n_valid  = accumarray(k(ok), 1, [n_min 1], @sum, 0);
n_cloudy = accumarray(k(ok), double(cloudy(ok)), [n_min 1], @sum, 0);
D.pier_cloud_frac = n_cloudy ./ n_valid;              % NaN where no ceilometer data
D.pier_multi_base = accumarray(k(ok), double(multi(ok)), [n_min 1], @max, NaN);

okc = in_grid(k) & cloudy;
D.pier_cbh_m = accumarray(k(okc), cbh_m(okc), [n_min 1], @median, NaN);
clear t cbh_m status valid cloudy multi



% -------------------------------------------------------------------
% --------------------- KAZR-ARSCL layers ---------------------------
% -------------------------------------------------------------------

fprintf('  ARSCL (7.9 million 4-s samples x 10 layers) ...\n');
t    = nc_time(files.arscl);
base = ncread(files.arscl, 'cloud_layer_base_height');   % (layer x time) as stored, m above ground
top  = ncread(files.arscl, 'cloud_layer_top_height');
if size(base, 1) ~= numel(t)                              % make it (time x layer)
    base = base';  top = top';
end
present  = isfinite(base) & base >= 0;                    % empty slots are -1 (flag "clear_sky")
n_layers = double(sum(present, 2));
top(~present) = NaN;
top_max = double(max(top, [], 2));                        % highest layer top in each sample (NaN if none)
clear base top present

k = row_of(t);
ok = in_grid(k);
D.arscl_n          = accumarray(k(ok), 1,              [n_min 1], @sum, 0);
D.arscl_max_layers = accumarray(k(ok), n_layers(ok),   [n_min 1], @max, NaN);
okt = ok & isfinite(top_max);
D.arscl_top_max_m  = accumarray(k(okt), top_max(okt),  [n_min 1], @max, NaN);
clear t n_layers top_max



% -------------------------------------------------------------------
% ---------------- Laser disdrometers (Pier and Soledad) ------------
% -------------------------------------------------------------------

fprintf('  disdrometers ...\n');
D.rain_pier    = rain_per_minute(files.rain_pier, t0, n_min);
D.rain_soledad = rain_per_minute(files.rain_soledad, t0, n_min);



% -------------------------------------------------------------------
% -------- Mt. Soledad in-cloud condition (visibility + FM-120) -----
% -------------------------------------------------------------------

fprintf('  Mt. Soledad visibility and FM-120 ...\n');
vis = read_soledad_visibility(files.visibility);
fm  = read_fm120_daily_files(files.fm120_folder);
cloud5 = soledad_cloud_condition_5min(vis, fm, 'TimeRange', time_range);
clear vis

% each minute takes the class of the 5-min interval that holds it
k5 = floor((0:n_min - 1)' / 5) + 1;
known = cloud5.has_vis & cloud5.has_fm120;
sc = double(cloud5.cloud);
sc(~known) = NaN;                       % unknown without both measurements
D.soledad_cloud = sc(k5);


end




% ===================================================================
% ========================= Local functions =========================
% ===================================================================

function r = rain_per_minute(file, t0, n_min)
% 1 if any QC = 0 sample in the minute had precip_rate > 0, 0 if the minute
% had QC = 0 samples but none with rain, NaN if no data
t    = nc_time(file);
rate = double(ncread(file, 'precip_rate'));          % mm h-1
qc   = double(ncread(file, 'qc_precip_rate'));
k    = floor(seconds(t - t0) / 60 + 1e-9) + 1;       % grid row of each sample
ok   = isfinite(rate) & qc == 0 & k >= 1 & k <= n_min;
r = accumarray(k(ok), double(rate(ok) > 0), [n_min 1], @max, NaN);
end


function t = nc_time(file)
% "seconds since 1970-01-01" -> datetime (UTC), as in every combined ARM file here
units = string(ncreadatt(file, 'time', 'units'));
if ~startsWith(units, "seconds since 1970-01-01")
    error('build_epcape_1min_inputs: unexpected time units "%s" in %s', units, file);
end
t = datetime(double(ncread(file, 'time')), 'ConvertFrom', 'posixtime');
t = t(:);
end
