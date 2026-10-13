%% CloudCube (JPL Ka/W/G-band radar) at EPCAPE: data availability and low-cloud fraction
%
% PURPOSE
%   Answers two questions for the EPCAPE BAMS table (DataAvailability sheet, row
%   "JPL-Pier | Ka, W, G band Radar | https://zenodo.org/records/14811996"):
%     Q1. What percentage of the time is the JPL CloudCube radar data available?
%     Q2. Can a liquid (low) cloud fraction be estimated from it, and from what data?
%
% SHORT ANSWERS (details and numbers are printed when the script runs)
%   * The Zenodo record (Socuellamos et al., doi:10.5281/zenodo.14811996) holds
%     PROCESSED radar data, not raw I/Q data, but only Level-1-type products:
%       - CloudCube_EPCAPE_Gband_Moments.zip     (0.20 GB) G-band reflectivity,
%             mean Doppler velocity and spectrum width, 10 m range gates
%       - CloudCube_EPCAPE_Multifrequency.zip    (0.04 GB) Ka-, W-, G-band
%             reflectivity and dual-frequency ratios, 60 m gates
%       - CloudCube_EPCAPE_Gband_Spectra.zip     (1.5 GB)  G-band Doppler spectra
%             with the noise masked out
%       - CloudCube_EPCAPE_Gband_Spectra_Noise.zip (47 GB) the same spectra WITH
%             noise and SNR
%     There is NO cloud mask, cloud fraction, LWC, LWP, r_e or phase product.
%     Any cloud fraction has to be derived from the reflectivity yourself (done
%     below as a "hydrometeor occurrence" with a height cut).
%   * The 47 GB Spectra_Noise file is NOT needed: availability needs only the
%     time stamps, and echo occurrence needs only reflectivity. Both come from
%     the 0.24 GB Moments + Multifrequency zips. (The noise-masked Spectra zip
%     is not needed either.)
%   * CloudCube only ran for short, operator-selected periods between
%     2023-03-23 and 2023-04-27 (~51 h of G-band data in 56 files). That is far
%     too little, and too non-random, for a Feb-2023-to-Feb-2024 seasonal table
%     entry of cloud fraction. For cloud fraction, ARM value-added products at
%     the Pier (KAZR-ARSCL, ceilometer) cover the whole campaign; Section 5
%     computes those over the same minutes as CloudCube and over the full
%     campaign so the two can be compared directly.
%
% WHAT THE SCRIPT DOES
%   2. Reads the time stamps of every G-band Moments and Multifrequency file and
%      counts the 1-min bins holding >= 1 radar profile -> hours available per
%      month, % of each calendar month, % of the deployment window.
%   3. G-band hydrometeor occurrence: a profile is "cloudy" if any range gate
%      between min_range_m and low_cloud_max_m has a finite reflectivity. The
%      noise is already masked out in the Moments files (NaN = no signal), so
%      "finite" = detected echo. Averaged into 1-min bins, then over bins.
%   4. Ka-band (and W-band) echo occurrence from the Multifrequency files, for
%      comparison with the more sensitive G band.
%   5. ARM comparison at the Pier (files already used for the BAMS table):
%        KAZR-ARSCL (epcarsclkazr1kolliasM1.c1) lowest hydrometeor layer, and
%        Vaisala ceilometer (epcceilM1.b1) cloud-base detection,
%      plus ARSCL's lidar-only cloud_base_best_estimate, evaluated (a) only
%      during CloudCube minutes, (b) for full Mar and Apr 2023, (c) for the
%      table's seasons.
%   6. Saves a summary text file and a figure in ./output.
%
% WHY THESE ARE "HYDROMETEOR" FRACTIONS, NOT STRICTLY "LIQUID CLOUD" FRACTIONS
%   - Radar reflectivity scales with D^6, so drizzle/rain below cloud base is
%     detected as "cloud", and the radar cannot tell liquid from ice. With the
%     height cut (<= 3.5 km, same as the BAMS table's low-cloud definition)
%     the echoes are mostly the marine stratocumulus deck, but frontal or
%     mixed-phase cloud below 3.5 km in March-April would also be counted.
%   - KAZR-ARSCL radar layers miss much of the thin stratocumulus (see 5a).
%   - The ceilometer is the conventional liquid-cloud proxy: lidar is
%     extinguished within tens of meters inside liquid cloud, so a detected
%     base <= 3.5 km is very likely liquid cloud base.
%   - ARSCL combines KAZR + lidars (Clothiaux et al., 2000, J. Appl. Meteor.,
%     39, 645-665, doi:10.1175/1520-0450(2000)039<0645:ODOCHA>2.0.CO;2).
%
% DATA (edit paths in Section 1)
%   CloudCube: https://zenodo.org/records/14811996 . Unzip
%     CloudCube_EPCAPE_Gband_Moments.zip and CloudCube_EPCAPE_Multifrequency.zip.
%     Product description: Socuellamos, J. M., Rodriguez Monje, R., Lebsock,
%     M. D., Cooper, K. B., Beauchamp, R. M., and Umeyama, A.: Multifrequency
%     radar observations of marine clouds during the EPCAPE campaign, Earth
%     Syst. Sci. Data, 16, 2701-2715, doi:10.5194/essd-16-2701-2024, 2024.
%   ARM: the combined files in table_for_BAMS_data/ARM (see
%     EPCAPE/table_for_BAMS/README.md).
%
% NOTE ON THE NUMBERS ALREADY IN THE TABLE
%   The DataAvailability sheet lists 0.1075 (Mar) and 0.375 (Apr) for this
%   radar. Section 2 prints the fraction of TIME with data (from the file time
%   stamps, like the other rows of the sheet). Section 2b shows that March
%   = 10/93 = band-days with data in Socuellamos et al. (2024) Fig. 3 divided
%   by 3 bands x 31 days; April's 0.375 is not reproduced by the same rule.
%   For the G band the Zenodo files are the complete record (Fig. 3 shows no
%   unreleased G-band days), so 51.8 h is all the G-band data that exist;
%   more Ka/W data exist (W 23-24 Mar, Ka/W 12 and 17-27 Apr) but are unreleased.
%
% Written for MATLAB R2024b, base MATLAB only. Run time ~1-2 min (the ARSCL
% file holds 7.9 million 4-s samples).
%
% Andrew Buggee, Scripps Institution of Oceanography, Oct 2026

clear; clc;

%% 1. Settings: EDIT THESE PATHS --------------------------------------------------------

% Folder holding the unzipped CloudCube folders
cloudcube_folder = "/Users/andrewbuggee/Documents/Scripps - UCSD/EPCAPE/Lebsock_radar_data";
moments_folder   = fullfile(cloudcube_folder, "CloudCube_EPCAPE_Gband_Moments");
multifreq_folder = fullfile(cloudcube_folder, "CloudCube_EPCAPE_Multifrequency");

% Folder with the ARM files used for the BAMS table (set to "" to skip Section 5)
arm_folder = "/Users/andrewbuggee/Documents/Scripps - UCSD/EPCAPE/table_for_BAMS_data/ARM";

% Where the summary and figure are written
output_folder = fullfile(fileparts(mfilename('fullpath')), 'output');

% The table's season helper (sheet_season_masks.m) lives in the BAMS-table code
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'table_for_BAMS', 'helpers'));

% ---- thresholds ----
% Lowest usable gate. In the G-band Moments files the gates at 0-100 m are
% filled in ~100% of profiles with values that do not change with the weather
% (antenna/transmitter leakage in the FMCW near field), so they are not cloud.
% Gates >= 150 m are used. (Checked on 2023-04-11 and 2023-04-26 files.)
settings.min_range_m     = 150;      % m above the radar
% Low-cloud top of the detection window: same 3.5 km cut used for the BAMS
% table's low-cloud rows (cloud top / cloud base <= 3500 m)
settings.low_cloud_max_m = 3500;     % m above the radar / ground
% Time bin used to define "available" and to average the cloud flags
settings.bin_min         = 1;        % minutes
% CloudCube deployment window stated in the Zenodo record ("March 23 to April 27")
settings.deploy_start    = datetime(2023, 3, 23);
settings.deploy_end      = datetime(2023, 4, 28);   % exclusive (end of Apr 27)

bin_s = 60 * settings.bin_min;      % s per bin


%% 2. CloudCube data availability -------------------------------------------------------
% A 1-min bin counts as "available" if at least one profile falls in it. The
% files are short bursts (2-380 min) with 0.1-1.9 s profile spacing, so 1-min
% bins lose < 1 min per file at the edges.

fprintf('Reading CloudCube G-band Moments files ...\n');
G = read_cloudcube_folder(moments_folder, "reflectivity", settings);
fprintf('Reading CloudCube Multifrequency files ...\n');
M = read_cloudcube_folder(multifreq_folder, ["reflectivity_ka", "reflectivity_w", "reflectivity_g"], settings);

% Unique 1-min bins (integer minutes since 1970-01-01) with data
G_bins = unique(floor(G.t_epoch_s / bin_s));
M_bins = unique(floor(M.t_epoch_s / bin_s));

% Availability per calendar month and over the deployment window
months = [datetime(2023,3,1), datetime(2023,4,1)];
avail  = table('Size', [0 7], ...
    'VariableTypes', ["string","double","double","double","double","double","double"], ...
    'VariableNames', ["period","hours_in_period","Gband_hours","Gband_pct","Multifreq_hours","Multifreq_pct","days_with_Gband_data"]);

for m = 1:numel(months)
    t1 = months(m);  t2 = t1 + calmonths(1);
    avail = [avail; availability_row(string(datetime(t1, 'Format', 'MMM yyyy')), t1, t2, G_bins, M_bins, bin_s)]; %#ok<AGROW>
end
avail = [avail; availability_row("Deployment 23 Mar - 27 Apr 2023", settings.deploy_start, settings.deploy_end, G_bins, M_bins, bin_s)];
avail = [avail; availability_row("EPCAPE 15 Feb 2023 - 14 Feb 2024", datetime(2023,2,15), datetime(2024,2,15), G_bins, M_bins, bin_s)];

fprintf('\n===== Q1: CloudCube data availability (1-min bins with >= 1 profile) =====\n');
disp(avail);
fprintf(['Values currently in the BAMS DataAvailability sheet: Mar 2023 = 10.75%%, Apr 2023 = 37.5%%.\n' ...
         'From the file time stamps:                           Mar 2023 = %.2f%%, Apr 2023 = %.2f%% (G band).\n'], ...
         avail.Gband_pct(1), avail.Gband_pct(2));
fprintf('Days with G-band data: %s\n', strjoin(string(unique(dateshift(epoch_to_datetime(G_bins*bin_s), 'start', 'day')), 'dd-MMM'), ', '));


%% 2b. Where the table's 10.75% (Mar) and 37.5% (Apr) probably come from ------------------
% Socuellamos et al. (2024, ESSD 16, 2701-2715, doi:10.5194/essd-16-2701-2024),
% Sect. 2.3: CloudCube ran on weekdays, ~06:00-18:00 PT, 23 Mar - 27 Apr 2023,
% and "only when cloud targets were present"; its three modules ran
% independently. Their Fig. 3 marks, for each band and day, whether data exist
% (including days with data "available but not provided" on Zenodo).
% Fig. 3 transcribed below: 1 = data exist (provided or not), 0 = no data.
%            Ka  W   G
fig3_days = datetime(2023,3,[23 24 30 31]');                 % March days in Fig. 3
fig3_mar  = [0 1 1;  0 1 1;  1 1 1;  1 1 1];
fig3_apr_days = datetime(2023,4,[1 3 11 12 13 17 18 25 26 27]');
fig3_apr  = [1 1 1;  1 1 1;  1 1 1;  1 1 0;  1 1 1;  1 1 1;  1 1 1;  1 1 1;  1 1 1;  1 1 1];
% "Band-days" with data / (3 bands x days in month):
%   March: 10 / 93 = 0.10753 -> reproduces the table's 0.1075268817 EXACTLY.
%   April: 29 / 90 = 0.322   -> does NOT reproduce 0.375 (0.375 x 90 = 33.75,
%          not an integer), so April was computed some other way.
fprintf('\nFig. 3 band-days: Mar %d/93 = %.4f (table 0.1075), Apr %d/90 = %.4f (table 0.375)\n', ...
    sum(fig3_mar(:)), sum(fig3_mar(:))/93, sum(fig3_apr(:)), sum(fig3_apr(:))/90);
fprintf(['A day counts in full even if the radar ran for minutes, so this metric is not\n' ...
         'comparable to the other rows of the sheet (fraction of TIME with data).\n']);


%% 3. G-band hydrometeor (low-cloud) occurrence -----------------------------------------
% Per profile: true if any gate in [min_range_m, low_cloud_max_m] has a
% finite (= above noise) reflectivity. Also flag echoes ABOVE low_cloud_max_m
% (G-band unambiguous range is 6.3 km) to show how much of the sampled time
% was mid/high (likely ice) cloud.
% Profiles are averaged within each 1-min bin first so that periods sampled
% at 0.37 s and at 1.9 s spacing get the same weight per minute.

[G_bin_id, G_low_frac]  = bin_average(G.t_epoch_s, G.echo_low(:, 1),  bin_s);
[~,        G_high_frac] = bin_average(G.t_epoch_s, G.echo_high(:, 1), bin_s);
G_bin_time = epoch_to_datetime(G_bin_id * bin_s);

fprintf('\n===== Q2: CloudCube G-band hydrometeor occurrence (%g-%g m) =====\n', settings.min_range_m, settings.low_cloud_max_m);
cc = table('Size', [0 5], 'VariableTypes', ["string","double","double","double","double"], ...
    'VariableNames', ["period","Gband_hours","low_echo_pct","echo_above_low_max_pct","n_minutes"]);
for m = 1:numel(months)
    in_m = G_bin_time >= months(m) & G_bin_time < months(m) + calmonths(1);
    cc = [cc; {string(datetime(months(m), 'Format', 'MMM yyyy')), sum(in_m)/60*settings.bin_min, ...
        100*mean(G_low_frac(in_m)), 100*mean(G_high_frac(in_m)), sum(in_m)}]; %#ok<AGROW>
end
cc = [cc; {"All CloudCube minutes", numel(G_bin_id)/60*settings.bin_min, 100*mean(G_low_frac), 100*mean(G_high_frac), numel(G_bin_id)}];
disp(cc);


%% 4. Ka- and W-band echo occurrence (Multifrequency files) -----------------------------
% Same test as Section 3 on each band of the Multifrequency product (60 m gates,
% so the first usable gate is the first one >= min_range_m). Ka/W are less
% sensitive than G band to small cloud droplets, but less attenuated by rain.

band_names = ["Ka", "W", "G"];
fprintf('\n===== Multifrequency product: echo occurrence %g-%g m during its %.1f h =====\n', ...
    settings.min_range_m, settings.low_cloud_max_m, numel(M_bins)/60*settings.bin_min);
for b = 1:numel(band_names)
    [~, f_low] = bin_average(M.t_epoch_s, M.echo_low(:, b), bin_s);
    fprintf('  %-2s band: low-echo occurrence = %5.1f %%\n', band_names(b), 100*mean(f_low));
end


%% 5. ARM comparison: ARSCL and ceilometer ----------------------------------------------
arm_summary = table();
if strlength(arm_folder) > 0 && isfolder(arm_folder)

    % ---- 5a. KAZR-ARSCL ----
    % Two ARSCL variables, which answer different questions:
    %  * cloud_layer_base_height(layer, time), layer 1 = lowest: hydrometeor
    %    layers from the radar (KAZR) + lidars. -1 = no layer. Read in MATLAB
    %    the -1 flag is kept (ncread does not apply valid_range).
    %  * cloud_base_best_estimate(time): cloud base from the ceilometer and
    %    micropulse lidar only. -1 clear, -2 possibly clear, >= 0 base (m AGL).
    % In this file the 4-s time grid is complete and neither variable is ever
    % NaN, so ARSCL itself cannot tell us when the instruments were down.
    % FINDING (checked when this script was written): when the lidar best
    % estimate has a base <= 3.5 km, ARSCL layer 1 is EMPTY 33% of the time and
    % above 3.5 km another 5%. KAZR misses a large share of the thin, weakly
    % reflecting marine stratocumulus that the lidars see, so the radar-based
    % layers UNDER-estimate liquid cloud occurrence (31% vs 49% campaign-wide).
    arscl_file = fullfile(arm_folder, 'cloud_arscl_M1_20230215_20240214.nc');
    fprintf('\nReading ARSCL (7.9 million samples) ...\n');
    t_arscl  = ncread_epoch_time(arscl_file);
    base1_m  = double(ncread(arscl_file, 'cloud_layer_base_height', [1 1], [1 Inf])).';   % m AGL
    top1_m   = double(ncread(arscl_file, 'cloud_layer_top_height',  [1 1], [1 Inf])).';   % m AGL
    cbbe_m   = double(ncread(arscl_file, 'cloud_base_best_estimate'));  cbbe_m = cbbe_m(:); % m AGL
    arscl_valid = isfinite(base1_m) & isfinite(top1_m) & isfinite(cbbe_m);
    % Radar+lidar low hydrometeor layer: lowest layer exists and its base <= 3.5 km
    arscl_low   = arscl_valid & base1_m >= 0 & top1_m >= 0 & base1_m <= settings.low_cloud_max_m;
    % Lidar-only low cloud base (the liquid-cloud proxy)
    cbbe_low    = arscl_valid & cbbe_m >= 0 & cbbe_m <= settings.low_cloud_max_m;
    clear base1_m top1_m cbbe_m

    % ---- 5b. Ceilometer cloud base ----
    % detection_status: 0 clear; 1-3 one to three bases; 4 obscured (fog/heavy
    % precip, no base). first_cbh is missing when status is 0 or 4, so the
    % qc_first_cbh "missing value" bit is set on every clear sample: the valid
    % sample set is therefore defined by detection_status, not by qc_first_cbh.
    % (status 4 never occurs in this EPCAPE file, so fog is not split out.)
    % Unlike ARSCL, the ceilometer file has gaps (1.84 M samples at 16 s =
    % ~8200 h of 8760 h), so its denominator excludes downtime.
    ceil_file = fullfile(arm_folder, 'cbh_ceil_M1_20230215_20240214.nc');
    t_ceil    = ncread_epoch_time(ceil_file);
    status    = double(ncread(ceil_file, 'detection_status'));  status = status(:);
    cbh_m     = double(ncread(ceil_file, 'first_cbh'));         cbh_m  = cbh_m(:);   % m above instrument
    ceil_valid   = ismember(status, 0:4);
    ceil_liquid  = ceil_valid & ismember(status, 1:3) & isfinite(cbh_m) & cbh_m <= settings.low_cloud_max_m;

    % ---- 5c. Evaluate over three sets of times ----
    % (i) only minutes when CloudCube G band has data
    in_cc_arscl = ismember(floor(t_arscl / bin_s), G_bin_id);
    in_cc_ceil  = ismember(floor(t_ceil  / bin_s), G_bin_id);
    dt_arscl = epoch_to_datetime(t_arscl);
    dt_ceil  = epoch_to_datetime(t_ceil);

    periods = {"CloudCube minutes only", in_cc_arscl, in_cc_ceil};
    % (ii) the full calendar months CloudCube sampled
    for m = 1:numel(months)
        ma = dt_arscl >= months(m) & dt_arscl < months(m) + calmonths(1);
        mc = dt_ceil  >= months(m) & dt_ceil  < months(m) + calmonths(1);
        periods(end+1, :) = {"Full " + string(datetime(months(m), 'Format', 'MMM yyyy')), ma, mc}; %#ok<AGROW>
    end
    % (iii) the table's seasons
    season_names = ["EPCAPE (F23-F24)", "Spring (AMJn23)", "Summer (JlAS23)", "Fall (OND23)", "Winter (FM23+JF24)"];
    Sa = sheet_season_masks(dt_arscl);
    Sc = sheet_season_masks(dt_ceil);
    for s = 1:5
        periods(end+1, :) = {season_names(s), Sa(:, s), Sc(:, s)}; %#ok<AGROW>
    end

    arm_summary = table('Size', [0 5], 'VariableTypes', ["string","double","double","double","double"], ...
        'VariableNames', ["period","ARSCL_radar_lidar_layer_pct","ARSCL_lidar_base_pct", ...
                          "ceilometer_base_pct","ceilometer_hours"]);
    for p = 1:size(periods, 1)
        ka = periods{p, 2} & arscl_valid;
        kc = periods{p, 3} & ceil_valid;
        arm_summary = [arm_summary; {periods{p, 1}, 100*mean(arscl_low(ka)), 100*mean(cbbe_low(ka)), ...
            100*mean(ceil_liquid(kc)), sum(kc)*16/3600}]; %#ok<AGROW>
    end

    fprintf('\n===== ARM at the Pier: low-cloud occurrence (%%, base <= %g m) =====\n', settings.low_cloud_max_m);
    disp(arm_summary);
    fprintf(['Compare row 1 (CloudCube minutes) with CloudCube''s own %.1f %% above: same times,\n' ...
             'different instruments. Compare rows 2-3 (full months) with row 1 to see how\n' ...
             'unrepresentative CloudCube''s sampling is (it was run on cloudy days).\n' ...
             'For a liquid cloud fraction use the lidar columns, not the radar-layer column.\n'], 100*mean(G_low_frac));
else
    fprintf('\nARM folder not found; Section 5 skipped.\n');
    dt_arscl = []; cbbe_low = []; arscl_valid = [];
end


%% 6. Save summary and figure -----------------------------------------------------------
if ~isfolder(output_folder)
    mkdir(output_folder);
end
writetable(avail, fullfile(output_folder, 'cloudcube_availability.csv'));
writetable(cc,    fullfile(output_folder, 'cloudcube_Gband_low_echo_occurrence.csv'));
if ~isempty(arm_summary)
    writetable(arm_summary, fullfile(output_folder, 'ARM_low_cloud_occurrence.csv'));
end

% Figure: (top) when CloudCube ran, coloured by whether the 1-min bin had a low
% echo, with the hourly ARSCL low-cloud fraction behind it; (bottom) daily
% hours of CloudCube data. Factory-style colours are set on the figure so the
% root defaults in startup.m (dark theme) do not change the result.
fig = figure('Color', 'w', 'Position', [100 100 1100 650], ...
    'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', 'DefaultAxesYColor', 'k', ...
    'DefaultTextColor', 'k', 'DefaultAxesFontSize', 11, 'DefaultLineLineWidth', 1);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
xl = [datetime(2023,3,20), datetime(2023,4,30)];

ax1 = nexttile(tl); hold(ax1, 'on');
if ~isempty(cbbe_low)
    % Hourly ARSCL lidar (ceilometer + MPL) low-cloud-base fraction for context
    in_win = dt_arscl >= xl(1) & dt_arscl < xl(2) & arscl_valid;
    hr = dateshift(dt_arscl(in_win), 'start', 'hour');
    [uh, ~, ih] = unique(hr);
    frac_hr = accumarray(ih, double(cbbe_low(in_win)), [], @mean);
    area(ax1, uh, 100*frac_hr, 'FaceColor', [0.80 0.85 0.95], 'EdgeColor', 'none', 'DisplayName', 'ARSCL lidar hourly low-cloud-base %');
end
scatter(ax1, G_bin_time, 100*G_low_frac, 8, [0.85 0.33 0.10], 'filled', 'DisplayName', 'CloudCube G-band 1-min low-echo %');
ylim(ax1, [0 105]); xlim(ax1, xl); ylabel(ax1, 'Occurrence (%)');
legend(ax1, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
title(ax1, sprintf('Low-cloud (%g-%g m) occurrence', settings.min_range_m, settings.low_cloud_max_m));
grid(ax1, 'on');

ax2 = nexttile(tl);
day_of_bin = dateshift(G_bin_time, 'start', 'day');
[ud, ~, id] = unique(day_of_bin);
bar(ax2, ud, accumarray(id, 1)/60*settings.bin_min, 1, 'FaceColor', [0.30 0.45 0.70], 'EdgeColor', 'none');
xlim(ax2, xl); ylabel(ax2, 'CloudCube G-band data (h day^{-1})'); grid(ax2, 'on');
title(ax2, sprintf('CloudCube availability: %.1f h total (%.1f%% of Mar, %.1f%% of Apr 2023)', ...
    numel(G_bin_id)/60*settings.bin_min, avail.Gband_pct(1), avail.Gband_pct(2)));

% Same weekly date ticks on both panels (bar() otherwise puts a tick on every
% bar), factory-like font scaling (startup.m sets large title/label
% multipliers), and no axes toolbar in the exported image
for ax = [ax1, ax2]
    xticks(ax, xl(1):days(7):xl(2));
    xtickformat(ax, 'dd-MMM');
    set(ax, 'FontWeight', 'normal', 'TitleFontSizeMultiplier', 1.1, 'LabelFontSizeMultiplier', 1.0, 'TitleFontWeight', 'bold');
    ax.Toolbar.Visible = 'off';
end

exportgraphics(fig, fullfile(output_folder, 'cloudcube_availability_and_low_cloud.png'), 'Resolution', 200);
fprintf('\nOutputs written to %s\n', output_folder);


%% ======================================================================================
%  Local functions
%  ======================================================================================

function S = read_cloudcube_folder(folder, z_vars, settings)
% READ_CLOUDCUBE_FOLDER  Read time and echo flags from every CloudCube .nc in a folder.
%
%   S.t_epoch_s   [n_profiles x 1] profile time, s since 1970-01-01 UTC
%   S.echo_low    [n_profiles x n_vars] logical: any finite reflectivity in
%                 [settings.min_range_m, settings.low_cloud_max_m]
%   S.echo_high   [n_profiles x n_vars] logical: any finite reflectivity above
%                 settings.low_cloud_max_m
%   z_vars        reflectivity variable name(s); one column per name
%
% The files store reflectivity(time, range) (netCDF/C order), which ncread
% returns as [range x time]. Noise is masked to NaN by the data producer.

files = dir(fullfile(folder, '*.nc'));
if isempty(files)
    error('No .nc files in %s. Unzip the CloudCube zip files there first.', folder);
end

t_all = cell(numel(files), 1);
low_all = cell(numel(files), 1);
high_all = cell(numel(files), 1);
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    t_all{k} = ncread_epoch_time(f);
    r_m = double(ncread(f, 'range'));                                 % m above radar
    in_low  = r_m >= settings.min_range_m & r_m <= settings.low_cloud_max_m;
    in_high = r_m >  settings.low_cloud_max_m;
    low_k  = false(numel(t_all{k}), numel(z_vars));
    high_k = false(numel(t_all{k}), numel(z_vars));
    for v = 1:numel(z_vars)
        Z_dBZ = ncread(f, z_vars(v));                                 % [range x time], dBZ
        low_k(:, v)  = any(isfinite(Z_dBZ(in_low, :)),  1).';
        high_k(:, v) = any(isfinite(Z_dBZ(in_high, :)), 1).';
    end
    low_all{k} = low_k;
    high_all{k} = high_k;
end
S.t_epoch_s = vertcat(t_all{:});
S.echo_low  = vertcat(low_all{:});
S.echo_high = vertcat(high_all{:});
S.n_files   = numel(files);
end


function t_epoch_s = ncread_epoch_time(file)
% NCREAD_EPOCH_TIME  Read 'time' and return seconds since 1970-01-01 00:00 UTC.
% Both CloudCube products and the combined ARM files use an epoch reference;
% this is checked so a file with a different reference fails loudly.
units = string(ncreadatt(file, 'time', 'units'));
if ~startsWith(units, "seconds since 1970-01-01")
    error('Unexpected time units "%s" in %s', units, file);
end
t_epoch_s = double(ncread(file, 'time'));
t_epoch_s = t_epoch_s(:);
end


function t = epoch_to_datetime(t_epoch_s)
% EPOCH_TO_DATETIME  Seconds since 1970-01-01 UTC -> datetime (no time zone).
t = datetime(1970, 1, 1) + seconds(t_epoch_s(:));
end


function [bin_id, frac] = bin_average(t_epoch_s, flag, bin_s)
% BIN_AVERAGE  Mean of a logical flag inside each time bin that holds data.
%   bin_id  integer bin index = floor(t / bin_s) (bin start = bin_id*bin_s)
%   frac    fraction of the bin's profiles with flag = true (0-1)
b = floor(t_epoch_s / bin_s);
[bin_id, ~, ib] = unique(b);
frac = accumarray(ib, double(flag(:)), [], @mean);
end


function row = availability_row(name, t1, t2, G_bins, M_bins, bin_s)
% AVAILABILITY_ROW  Hours and % of [t1, t2) covered by 1-min bins with data.
b1 = floor(seconds(t1 - datetime(1970,1,1)) / bin_s);
b2 = floor(seconds(t2 - datetime(1970,1,1)) / bin_s);
nG = sum(G_bins >= b1 & G_bins < b2);
nM = sum(M_bins >= b1 & M_bins < b2);
hours_total = hours(t2 - t1);
h_per_bin = bin_s / 3600;
days_G = numel(unique(floor(G_bins(G_bins >= b1 & G_bins < b2) * bin_s / 86400)));
row = table(name, hours_total, nG*h_per_bin, 100*nG*h_per_bin/hours_total, ...
    nM*h_per_bin, 100*nM*h_per_bin/hours_total, days_G, ...
    'VariableNames', ["period","hours_in_period","Gband_hours","Gband_pct","Multifreq_hours","Multifreq_pct","days_with_Gband_data"]);
end
