function results = calc_row10_LWP(arm_folder, kavin_folder, rain_window_min)
% CALC_ROW10_LWP  Sheet row 10: liquid water path (LWP) at the Scripps Pier.
%
%   results = calc_row10_LWP(arm_folder, kavin_folder, rain_window_min)
%
%   THREE NEW LINES (all in g m-2, sheet seasons in UTC):
%
%   Line 1, Kavin's own file, corrected: epcmwrret_turnM1_hourly_averages.nc
%     holds hourly means of MWRRET v2 phys_lwp (the mean of the QC = 0
%     samples with LWP > 0 in each hour; all sky, no cloud or rain screen).
%     Kavin's code put these into a nan(500,13) array; months with more than
%     500 hours made MATLAB grow the array and fill the new cells of the
%     other months with ZEROS, which entered his means. Here every hour is
%     used once and there are no padding zeros.
%
%   Line 2, MWRRET v1 (epcmwrret1liljclouM1.c2, 2-channel MWR, ~20-s samples;
%     Turner et al. 2007, IEEE TGRS 45, 3680), cloudy scenes, rain-free:
%       * phys_lwp with qc_phys_lwp = 0
%       * cloudy: qc_cloud_base_height = 0 and cloud_base_height >= 0
%         (that QC field's bit 2 marks clear sky)
%       * rain-free: no rain at the Pier Parsivel disdrometer (LD M1,
%         precip_rate > 0 with QC = 0) within ±rain_window_min (30 min);
%         a wet radome makes the retrieved LWP far too high.
%
%   Line 3, MWRRET v2 (epcmwrret2turnM1.c1, 3-channel MWR, ~1-s samples),
%     cloudy scenes, rain-free:
%       * phys_lwp with qc_phys_lwp = 0 AND phys_qc_flag = 0
%       * cloudy: ceilometer detection_status 1, 2 or 3 (a cloud base was
%         found) and cbh_detected > 0
%       * rain-free as in line 2.
%   No time averaging: every sample that passes the filters counts once.

% ===================== Rain at the Pier (for lines 2 and 3) ============================
ld_file   = fullfile(arm_folder, 'rain_ld_M1_20230215_20240214.nc');
t_ld      = read_netcdf_time(ld_file);
rain_rate = read_arm_qc0(ld_file, 'precip_rate');          % mm/h, NaN where QC ~= 0
rain_times = t_ld(rain_rate > 0);                           % every 1-min sample with rain

% ===================== Line 1: Kavin's hourly file ======================================
kfile  = fullfile(kavin_folder, 'epcmwrret_turnM1_hourly_averages.nc');
t_hour = read_netcdf_time(kfile);                           % hour START (UTC)
lwp_h  = double(ncread(kfile, 'phys_lwp'));                 % g m-2
stats1 = season_stats(lwp_h, sheet_season_masks(t_hour));
results = make_result(10, "Liquid Water Path", "g m-2", "Pier", ...
    "epcmwrret_turnM1_hourly_averages.nc (Kavin's file; MWRRET v2)", stats1, ...
    "Same data as Kavin. filters: his hourly all-sky means (QC = 0 samples with LWP > 0); MATLAB padding " + ...
    "zeros removed; no cloud or rain filter; sheet seasons (UTC).");

% ===================== Line 2: MWRRET v1 =================================================
f1  = fullfile(arm_folder, 'lwp_mwrret1_M1_20230215_20240214.nc');
t1  = read_netcdf_time(f1);
lwp = read_arm_qc0(f1, 'phys_lwp');                          % g m-2, QC = 0 only
cbh_km = double(ncread(f1, 'cloud_base_height'));           % km
qc_cbh = double(ncread(f1, 'qc_cloud_base_height'));
cloudy = (qc_cbh == 0) & isfinite(cbh_km) & (cbh_km >= 0);
lwp(~cloudy) = NaN;                                          % cloudy samples only
lwp(near_rain(rain_times, t1, rain_window_min)) = NaN;       % rain-free only
stats2 = season_stats(lwp, sheet_season_masks(t1));
results(2) = make_result(10, "Liquid Water Path", "g m-2", "Pier", "epcmwrret1liljclouM1.c2 (MWRRET v1)", stats2, ...
    "filters: QC = 0 on phys_lwp; cloudy samples only (qc_cloud_base_height = 0); rain-free: no rain " + ...
    sprintf("(any rate > 0, QC = 0) at the Pier Parsivel (LD M1) within ±%g min; every sample; sheet seasons (UTC).", rain_window_min));
clear t1 lwp cbh_km qc_cbh cloudy

% ===================== Line 3: MWRRET v2 =================================================
f2  = fullfile(arm_folder, 'lwp_mwrret2_M1_20230215_20240214.nc');
t2  = read_netcdf_time(f2);
lwp = read_arm_qc0(f2, 'phys_lwp');                          % g m-2, QC = 0 only
phys_qc_flag = double(ncread(f2, 'phys_qc_flag'));           % 0 = retrieval good
lwp(phys_qc_flag ~= 0) = NaN;
status = double(ncread(f2, 'detection_status'));             % ceilometer: 1-3 = cloud base(s) found
cbh_km = double(ncread(f2, 'cbh_detected'));                 % km; -0.001 when no cloud
cloudy = ismember(status, [1 2 3]) & isfinite(cbh_km) & (cbh_km > 0);
lwp(~cloudy) = NaN;
lwp(near_rain(rain_times, t2, rain_window_min)) = NaN;
stats3 = season_stats(lwp, sheet_season_masks(t2));
results(3) = make_result(10, "Liquid Water Path", "g m-2", "Pier", "epcmwrret2turnM1.c1 (MWRRET v2)", stats3, ...
    "filters: QC = 0 on phys_lwp and phys_qc_flag = 0; cloudy samples only (cbh_detected > 0 and " + ...
    sprintf("detection_status 1-3); rain-free (no LD M1 rain within ±%g min); every sample; sheet seasons.", rain_window_min));

end
