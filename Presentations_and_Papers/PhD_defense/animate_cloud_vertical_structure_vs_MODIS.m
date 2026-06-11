%% Animated figure: complex vertical structure of a cloud (in-situ) vs. a coincident MODIS retrieval
%
%   Built for the PhD-defense introduction to motivate why a single satellite
%   column-effective-radius retrieval cannot capture the vertical structure of
%   a stratocumulus cloud.
%
%   The animation has three panels, side by side, and plays in two acts:
%
%   ACT 1  -- "complex vertical structure" (the C-130 flies a vertical sounding
%             from cloud base to cloud top through the SE-Pacific Sc deck):
%       Panel 1 (left)   : true-colour MODIS image of the scene. Aircraft marker
%                          + trailing flight track build as the plane climbs.
%       Panel 2 (middle) : altitude (m) vs. along-track distance (km). The slant
%                          ascent draws progressively with a marker at the
%                          current aircraft position.
%       Panel 3 (right)  : droplet size distribution log10(N_c) over (radius,
%                          altitude). Rows fill in from cloud base upward as the
%                          plane samples each level; the in-situ r_e(z) line
%                          draws on top.
%
%   ACT 2  -- "the coincident MODIS measurement" (revealed, then held):
%       Panel 1 : the coincident MODIS pixel footprints appear at true ~1 km
%                 scale next to the thin in-situ track -> footprint vs.
%                 point-sampling contrast.
%       Panel 3 : the single MODIS column-effective-radius r_e (bands 1+7)
%                 appears as a vertical line against the full in-situ r_e(z).
%
%   DATA: each example is regenerated from the raw VOCALS-REx .nc flight file and
%   its coincident MODIS granule, so the full 2-D droplet size distribution
%   Nc(radius, altitude) AND the MODIS coincidence index are available for ANY
%   flight (not just a single pre-saved product). Results are cached to
%   anim_cache/ so switching/re-running an example is fast after the first build.
%
%   To animate a different example, change SELECT below (see the `ex` registry).
%
%   Output: an MPEG-4 movie written next to this script.
%
% By Andrew John Buggee
%%

clear variables; close all;

% =========================================================================
% ---------------------------- CONFIG -------------------------------------
% =========================================================================

repo = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research';

% --- which example to animate (field name in the `ex` registry below) ---
SELECT = 'RF11_0911';

% --- panel-1 backdrop source toggle ---
%     'geotiff'  : georeferenced NASA Worldview true-colour GeoTIFF (smooth, ~250 m)
%     'pansharp' : pan-sharpened 250 m true colour from the L1B files (crisp, true 250 m)
%     'auto'     : fall back to whatever the example provides (GeoTIFF > HKM > QKM > 1 km)
PANEL1_BG = 'geotiff';

% --- set true to rebuild an example's cache from raw data even if it exists ---
FORCE_REBUILD = false;

% --- Set QUICK_TEST = true to render a coarse, fast movie (for checking layout) ---
QUICK_TEST = false;

% -------------------------------------------------------------------------
% Example registry. Each example points at:
%   .gn_file     - a saved (modis, vocalsRex) product holding the single
%                  coincident vertical profile + its MODIS coincidence index.
%   .vocals_nc   - the raw VOCALS-REx C-130 flight, used ONLY to recover the
%                  full 2-D droplet size distribution Nc(radius, altitude) for
%                  products that stored Nc collapsed to a vector.
%   .granule_dir - the coincident MODIS granule folder (trailing slash required)
%                  for the true-colour background, footprints, and retrieved r_e.
% Optional manual background image (e.g. a NASA Worldview true-colour export):
%   .bg_image  = '/path/to/worldview.png';
%   .bg_lonlim = [lon_left, lon_right];   .bg_latlim = [lat_bottom, lat_top];
% If omitted, the true-colour background is generated from the granule itself.
% (These three are the Fig. 3 examples from the paper; add more as needed.)
% -------------------------------------------------------------------------
nc_dir    = [repo '/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];
modis_dir = [repo '/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/'];

ex.RF12_1850 = struct( ...
    'label',       '2008-11-11, Aqua 18:50 UTC', ...
    'gn_file',     [modis_dir '2008_11_11_1850/Retrieval_outputs_05-Apr-2025/GN_inputs_outputs_withAdvection_rt-cov_8.76_rb-cov_100_tc-cov_100_05-Apr-2025_rev1.mat'], ...
    'vocals_nc',   [nc_dir 'RF12.20081111.125000_214500.PNI.nc'], ...
    'granule_dir', [modis_dir '2008_11_11_1850/']);

ex.RF12_1430 = struct( ...
    'label',       '2008-11-11, Terra 14:30 UTC', ...
    'gn_file',     [modis_dir '2008_11_11_1430/Retrieval_outputs_24-Sep-2024/GN_inputs_outputs_withAdvection_rt-cov_6.52_rb-cov_100_tc-cov_100_25-Sep-2024_rev1.mat'], ...
    'vocals_nc',   [nc_dir 'RF12.20081111.125000_214500.PNI.nc'], ...
    'granule_dir', [modis_dir '2008_11_11_1430/']);

ex.RF11_0911 = struct( ...
    'label',       '2008-11-09, Terra 14:30 UTC', ...
    'gn_file',     [modis_dir '2008_11_09/Retrieval_outputs_16-Jan-2025/GN_inputs_outputs_withAdvection_rt-cov_6.74_rb-cov_100_tc-cov_100_16-Jan-2025_rev1.mat'], ...
    'vocals_nc',   [nc_dir 'RF11.20081109.125700_213600.PNI.nc'], ...
    'granule_dir', [modis_dir '2008_11_09/'], ...
    'geotiff_file',[modis_dir '2008_11_09/snapshot-2008-11-09T00_00_00Z_60m/snapshot-2008-11-09T00_00_00Z.tif'], ... % 60 m Worldview true colour
    'qkm_file',    [modis_dir '2008_11_09/MOD02QKM.A2008314.1440.007.2025073193155.nc'], ...   % 250 m band 1 (pan-sharpen luminance)
    'hkm_file',    [modis_dir '2008_11_09/MOD02HKM.A2008314.1440.007.2025073193155.nc']);      % 500 m colour (pan-sharpen / fallback)

% --- playback speed. SLOWDOWN scales the number of frames at a fixed frame rate,
%     so the motion stays smooth (Act 1 interpolates the aircraft position per
%     frame -> more frames = smoother slow-motion, not just frames held longer).
%     1 = normal (~12 s), 2 = half speed (~24 s), 3 = third speed, etc. ---
SLOWDOWN = 2;

% --- frame counts and pacing ---
if QUICK_TEST
    Nf_act1    = 24;     % frames to build the in-situ profile
    Nf_reveal  = 14;     % frames to reveal the MODIS overlay
    Nf_hold    = 8;      % frames to hold the final comparison
    frame_rate = 12;
else
    Nf_act1    = round(165 * SLOWDOWN);   % frames to build the in-situ profile
    Nf_reveal  = round(60  * SLOWDOWN);   % frames to reveal the MODIS overlay
    Nf_hold    = round(72  * SLOWDOWN);   % frames to hold the final comparison
    frame_rate = 24;
end

% --- droplet-radius display cap (microns). The cloud mode and r_e (~6-9 um)
%     live below this; the sparse drizzle tail above it is suppressed so the
%     evolving cloud distribution and r_e are clearly visible. ---
r_display_max_um = 20;

% --- panel-1 zoom: degrees of margin added around the track + footprints.
%     Larger = more zoomed out (more cloud structure visible). ~0.012 is tight. ---
% MAP_PAD_DEG = 0.045;
MAP_PAD_DEG = 0.15;
% MAP_PAD_DEG = 0.25;      % Good cloud contrast, but pixels are perhaps too small

% --- fonts / styling ---
ax_fnt  = 20;
ttl_fnt = 23;
cb_fnt  = 20;
ln_sz   = 2;
geoLineSz = 2;
tck_label_sz = 18;

% --- paths ---
addpath(genpath([repo '/Generally_useful_functions']));
addpath(genpath([repo '/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval']));
addpath(genpath([repo '/Hyperspectral_Cloud_Retrievals/VOCALS_REx']));

out_dir   = [repo '/Presentations_and_Papers/PhD_defense/'];
cache_dir = [out_dir 'anim_cache/'];
out_name  = ['cloud_vertical_structure_vs_MODIS_' SELECT '.mp4'];


% =========================================================================
% ---------------------- LOAD / BUILD EXAMPLE DATA ------------------------
% =========================================================================
cfg = ex.(SELECT);
% the resolved backdrop source is part of the cache name so toggling PANEL1_BG or
% changing the zoom auto-rebuilds rather than loading a stale cache.
bg_tag = resolve_bg(cfg, PANEL1_BG);
cache_file = [cache_dir sprintf('prep_%s_pad%.3g_%s.mat', SELECT, MAP_PAD_DEG, bg_tag)];

dat = prep_example(cfg, cache_file, FORCE_REBUILD, MAP_PAD_DEG, bg_tag);
vr  = dat.vr;

% --- in-situ profile, ordered cloud base -> cloud top (altitude ascending) ---
[altitude_m, sort_idx] = sort(vr.altitude(:));            % m, ascending
n_alt   = numel(altitude_m);
re_um   = vr.re(sort_idx);   re_um   = re_um(:);          % effective radius, um
lat_deg = vr.latitude(sort_idx);   lat_deg = lat_deg(:);  % deg
lon_deg = vr.longitude(sort_idx);  lon_deg = lon_deg(:);  % deg
horz_dist_km = vr.horz_dist(sort_idx);  horz_dist_km = horz_dist_km(:) / 1000;   % m -> km

% --- droplet size distribution: Nc is [n_radius x n_altitude] ---
r_bins_um = double(vr.drop_radius_bin_center(:));         % um
Nc        = double(vr.Nc(:, sort_idx));                   % cm^-3, [n_radius x n_alt]
r_mask    = r_bins_um <= r_display_max_um;
r_plot_um = r_bins_um(r_mask);
Nc_plot   = Nc(r_mask, :);

% log10 for display. Empty (zero-count) and not-yet-revealed cells map to a
% below-range "sentinel" that renders WHITE (renderer-independent; transparency
% would composite to black under headless capture).
Nc_log = log10(Nc_plot);
Nc_log(~isfinite(Nc_log)) = NaN;
Nc_disp_full = Nc_log';                                   % [n_alt x n_radius]
clim_log   = [max(-1, min(Nc_disp_full(:))), max(Nc_disp_full(:))];
white_val  = clim_log(1) - 1;                             % sentinel -> white
Nc_disp_full(isnan(Nc_disp_full)) = white_val;            % empty bins -> white

% --- MODIS coincident pixels (precomputed in prep) ---
re_modis_mean_um = dat.re_modis_mean;
re_modis_lo_um   = dat.re_modis_lo;
re_modis_hi_um   = dat.re_modis_hi;
modis_dt_min     = dat.dt_min;                            % time offset to overpass (min)
foot_lat = dat.foot_lat;
foot_lon = dat.foot_lon;
latlim   = dat.latlim;
lonlim   = dat.lonlim;

fprintf('In-situ: %d levels, %.0f-%.0f m, r_e %.1f-%.1f um, track %.2f km\n', ...
    n_alt, altitude_m(1), altitude_m(end), min(re_um), max(re_um), horz_dist_km(end));
fprintf('MODIS  : %d coincident pixels, r_e = %.2f um (%.2f-%.2f), dt = %.2f min\n', ...
    numel(foot_lat), re_modis_mean_um, re_modis_lo_um, re_modis_hi_um, modis_dt_min);


% =========================================================================
% ----------------------------- COLOURS -----------------------------------
% =========================================================================
clr_track   = [0.85 0.10 0.10];     % in-situ flight track (red)
clr_plane   = [0.85 0.10 0.10];
clr_re      = [0 0 0];              % in-situ r_e line (black)
clr_modis   = [0.10 0.30 0.85];     % MODIS accent (blue)
turbo_cmap  = turbo(256);


% =========================================================================
% --------------------------- FIGURE SETUP --------------------------------
% =========================================================================
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [60 60 1860 620]);
set(fig, 'InvertHardcopy', 'off');     % keep white backgrounds when captured

% --- Panel 1: true-colour MODIS scene with the flight track (left) ---
% A regular lon/lat axes (NOT geoaxes): geoaxes cannot drape an RGB image, and a
% regular axes gives full control over tick-label fonts. DataAspectRatio is set
% so 1 km looks the same horizontally and vertically (footprint scale is honest).
% Expand the map limits so the panel box is filled exactly while 1 km reads the
% same horizontally and vertically (honest footprint scale). This replaces
% daspect, whose letterboxing made the axes placement unpredictable and pushed
% the y-label off-canvas.
ax1_pos = [0.07 0.13 0.225 0.74];
box_w = ax1_pos(3) * 1860;  box_h = ax1_pos(4) * 620;     % panel size in pixels
mlat  = mean(latlim);
kmpp  = max(diff(lonlim)*111*cosd(mlat)/box_w, diff(latlim)*111/box_h);   % km per pixel
lonlim = mean(lonlim) + [-0.5 0.5] * (kmpp * box_w / (111*cosd(mlat)));
latlim = mean(latlim) + [-0.5 0.5] * (kmpp * box_h / 111);

ax1 = axes('Parent', fig, 'Position', ax1_pos, 'Color', 'k', ...
    'XColor', 'k', 'YColor', 'k');
hold(ax1, 'on'); box(ax1, 'on');
% drape the true-colour image below everything (z = -1) so track/footprints sit on top
surface(ax1, dat.rlon, dat.rlat, -ones(size(dat.rlat)), dat.rgb, ...
    'EdgeColor', 'none', 'FaceColor', 'flat');   % flat -> crisp per-pixel (no smoothing)
view(ax1, 2);
xlim(ax1, lonlim); ylim(ax1, latlim);
xt = get(ax1, 'XTick');  set(ax1, 'XTickLabel', deg_labels(xt, 'W'));
yt = get(ax1, 'YTick');  set(ax1, 'YTickLabel', deg_labels(yt, 'S'));
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', tck_label_sz, 'Layer', 'top');
xlabel(ax1, 'Longitude', 'Interpreter', 'latex', 'FontSize', ax_fnt);
ylabel(ax1, 'Latitude',  'Interpreter', 'latex', 'FontSize', ax_fnt);
title(ax1,  'Aircraft track \& MODIS footprint', 'Interpreter', 'latex', 'FontSize', ttl_fnt);

% --- Panel 2: altitude vs along-track distance (middle) ---
ax2 = axes('Parent', fig, 'Position', [0.375 0.13 0.245 0.74], 'Color', 'w', ...
    'XColor', 'k', 'YColor', 'k');
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on'); grid(ax2, 'minor');
xlim(ax2, [0, max(horz_dist_km)*1.02]);
ylim(ax2, [altitude_m(1)-10, altitude_m(end)+10]);
xlabel(ax2, 'Along-track distance (km)', 'Interpreter', 'latex', 'FontSize', ax_fnt);
ylabel(ax2, 'Altitude (m)',              'Interpreter', 'latex', 'FontSize', ax_fnt);
title(ax2,  'Aircraft vertical sounding', 'Interpreter', 'latex', 'FontSize', ttl_fnt);
set(ax2, 'FontSize', tck_label_sz, 'TickLabelInterpreter', 'latex');
yline(altitude_m(end), 'LineWidth', 2, 'Color', 'k', 'LabelHorizontalAlignment',...
    'center', 'Label', 'Cloud Top', 'FontSize', 20, 'Interpreter','latex',...
    'LabelVerticalAlignment','middle','LabelColor', 'k');
yline(altitude_m(1), 'LineWidth', 2, 'Color', 'k', 'LabelHorizontalAlignment',...
    'center', 'Label', 'Cloud Base', 'FontSize', 20, 'Interpreter','latex',...
    'LabelVerticalAlignment','middle','LabelColor', 'k');

% --- Panel 3: droplet size distribution (right) ---
ax3 = axes('Parent', fig, 'Position', [0.70 0.13 0.20 0.74], 'Color', 'w', ...
    'XColor', 'k', 'YColor', 'k');
% Truecolor (explicit RGB per cell) rendering -> independent of the figure colormap.
RGB_full = nc_to_rgb(Nc_disp_full, turbo_cmap, clim_log);          % [n_alt x n_r x 3]
img = image(ax3, r_plot_um, altitude_m, ones([size(Nc_disp_full) 3]));   % start white
set(ax3, 'YDir', 'normal');
hold(ax3, 'on');
xlim(ax3, [min(r_plot_um), max(r_plot_um)]);
ylim(ax3, [altitude_m(1), altitude_m(end)]);

% manual truecolor colorbar strip (also RGB -> independent of the figure colormap)
cax = axes('Parent', fig, 'Position', [0.915 0.13 0.013 0.74], 'Color', 'w',...
    'FontSize', 18, 'XColor', 'k', 'YColor', 'k');
n_strip   = 256;
strip_rgb = nc_to_rgb(linspace(clim_log(1), clim_log(2), n_strip)', turbo_cmap, clim_log);
image(cax, [0 1], clim_log, strip_rgb);
set(cax, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', ...
    'TickLabelInterpreter', 'latex', 'FontSize', 13, 'XColor', 'k', 'YColor', 'k');
ylabel(cax, '$\log_{10}(N_c)\;\;(\mathrm{cm}^{-3})$', 'Interpreter', 'latex', 'FontSize', cb_fnt, 'Color', 'k');
xlabel(ax3, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt,...
    'Color','k');
ylabel(ax3, 'Altitude (m)',            'Interpreter', 'latex', 'FontSize', ax_fnt,...
    'Color','k');
title(ax3,  'Droplet size distribution', 'Interpreter', 'latex', 'FontSize', ttl_fnt,...
    'Color','k');
set(ax3, 'FontSize', tck_label_sz, 'TickLabelInterpreter', 'latex',...
    'XColor', 'k', 'YColor', 'k', 'Color', 'w');

% animated graphics handles -------------------------------------------------
track_h = plot(ax1, lon_deg(1), lat_deg(1), '-', 'Color', clr_track, 'LineWidth', geoLineSz);
plane_h = plot(ax1, lon_deg(1), lat_deg(1), '^', 'MarkerFaceColor', clr_plane, ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 12, 'LineWidth', 1);

alt_line_h   = plot(ax2, horz_dist_km(1), altitude_m(1), '-', 'Color', clr_track, 'LineWidth', ln_sz);
alt_marker_h = plot(ax2, horz_dist_km(1), altitude_m(1), '^', 'MarkerFaceColor', clr_plane, ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 12);

% -- plot effective radius profile --
% re_line_h = plot(ax3, re_um(1), altitude_m(1), '-', 'Color', clr_re, 'LineWidth', ln_sz, ...
%     'DisplayName', 'in-situ $r_e(z)$');

% MODIS overlay handles (created hidden, revealed in Act 2) -----------------
foot_h = gobjects(numel(foot_lat),1);
for k = 1:numel(foot_lat)
    foot_h(k) = plot(ax1, foot_lon{k}([1:end 1]), foot_lat{k}([1:end 1]), '-', ...
        'Color', clr_modis, 'LineWidth', 2, 'Visible', 'off');
end

% MODIS r_e on the distribution panel (hidden until Act 2)
yb = [altitude_m(1), altitude_m(end)];
modis_re_h = plot(ax3, [re_modis_mean_um re_modis_mean_um], yb, ':', 'Color', 'k', ...
    'LineWidth', ln_sz+1, 'Visible', 'off', 'DisplayName', 'MODIS $r_{2.1}$');

% annotation handle (created in Act 2)
ann_foot  = gobjects(0);

drawnow; pause(0.5);


% =========================================================================
% ----------------------------- VIDEO -------------------------------------
% =========================================================================
v = VideoWriter([out_dir, out_name], 'MPEG-4');
v.FrameRate = frame_rate;
v.Quality   = 100;
open(v);

% -------------------------------------------------------------------------
% ACT 1 -- build the in-situ vertical profile
% -------------------------------------------------------------------------
fprintf('Rendering Act 1 (%d frames)...\n', Nf_act1);
for f = 1:Nf_act1

    frac = f / Nf_act1;                 % 0 -> 1
    pos  = 1 + frac * (n_alt - 1);      % continuous level for smooth marker
    k    = floor(pos);                  % integer level revealed
    k    = min(max(k, 1), n_alt);
    fr   = pos - k;                     % fractional part for interpolation

    % interpolated aircraft position (smooth motion)
    if k < n_alt
        lat_i = lat_deg(k) + fr*(lat_deg(k+1)-lat_deg(k));
        lon_i = lon_deg(k) + fr*(lon_deg(k+1)-lon_deg(k));
        hd_i  = horz_dist_km(k) + fr*(horz_dist_km(k+1)-horz_dist_km(k));
        alt_i = altitude_m(k) + fr*(altitude_m(k+1)-altitude_m(k));
    else
        lat_i = lat_deg(end); lon_i = lon_deg(end);
        hd_i  = horz_dist_km(end); alt_i = altitude_m(end);
    end

    % geomap: track + plane
    set(track_h, 'XData', lon_deg(1:k), 'YData', lat_deg(1:k));
    set(plane_h, 'XData', lon_i, 'YData', lat_i);

    % altitude panel: slant ascent + marker
    set(alt_line_h,   'XData', [horz_dist_km(1:k); hd_i], 'YData', [altitude_m(1:k); alt_i]);
    set(alt_marker_h, 'XData', hd_i, 'YData', alt_i);

    % distribution panel: reveal rows up to current level, draw r_e(z)
    rgb = RGB_full;
    rgb(k+1:end, :, :) = 1;                 % not-yet-sampled levels stay white
    set(img, 'CData', rgb);
    % set(re_line_h, 'XData', re_um(1:k), 'YData', altitude_m(1:k));

    drawnow;
    writeVideo(v, getframe(fig));
end

% ensure the full profile is shown
set(img, 'CData', RGB_full);
% set(re_line_h, 'XData', re_um, 'YData', altitude_m);
% legend(ax3, re_line_h, 'Location', 'northeast', 'Interpreter', 'latex', ...
%     'FontSize', 18, 'Color', 'w', 'TextColor', 'k');

% -------------------------------------------------------------------------
% ACT 2 -- reveal the coincident MODIS measurement
% -------------------------------------------------------------------------
fprintf('Rendering Act 2 reveal (%d frames)...\n', Nf_reveal);

for k = 1:numel(foot_h); set(foot_h(k), 'Visible', 'on'); end
set(modis_re_h, 'Visible', 'on');

% legend(ax3, [re_line_h, modis_re_h], 'Location', 'southeast', 'Interpreter', 'latex', ...
%     'FontSize', 18, 'Color', 'w', 'TextColor', 'k');
legend(ax3, modis_re_h, 'Location', 'southeast', 'Interpreter', 'latex', ...
    'FontSize', 18, 'Color', 'w', 'TextColor', 'k');

str_foot = sprintf('MODIS pixel $\\approx$ 1 km\nin-situ track = %.1f km', horz_dist_km(end));
ann_foot = annotation(fig, 'textbox', [0.0375 0.4 0.27 0.10], 'String', str_foot, ...
    'Interpreter', 'latex', 'FontSize', 18, 'EdgeColor', clr_modis, 'LineWidth', 1.2, ...
    'BackgroundColor', 'w', 'Color', clr_modis, 'FitBoxToText', 'on', 'Visible', 'off');

for f = 1:Nf_reveal
    a = f / Nf_reveal;                          % 0 -> 1 emphasis ramp
    for k = 1:numel(foot_h)
        set(foot_h(k), 'LineWidth', 2 + 1.5*a);
    end
    if f == 1
        set(ann_foot,  'Visible', 'on');
    end
    drawnow;
    writeVideo(v, getframe(fig));
end

% -------------------------------------------------------------------------
% Hold the final comparison
% -------------------------------------------------------------------------
fprintf('Holding final frame (%d frames)...\n', Nf_hold);
final_frame = getframe(fig);
for f = 1:Nf_hold
    writeVideo(v, final_frame);
end

close(v);
fprintf('\nDone. Movie written to:\n  %s%s\n', out_dir, out_name);


% =========================================================================
% ----------------------------- HELPERS -----------------------------------
% =========================================================================
function dat = prep_example(cfg, cache_file, force_rebuild, map_pad, bg_src)
% Build (or load from cache) everything the animation needs for one example:
% a single coincident VOCALS-REx vertical profile WITH the full 2-D droplet size
% distribution Nc(radius, altitude) and the MODIS coincidence, plus a true-colour
% background image of the scene.
%
% The saved product already holds the coincident profile + modisIndex_minDist.
% If that product stored Nc collapsed to a vector, the full 2-D Nc(radius,
% altitude) is recovered by reading the raw flight and slicing the size-
% distribution columns whose UTC timestamps match the profile (no need to rerun
% the heavy find_verticalProfiles / Mie retrieval machinery).

    if ~force_rebuild && isfile(cache_file)
        fprintf('Loading cached example: %s\n', cache_file);
        S = load(cache_file);  dat = S.dat;  return
    end

    fprintf('Building example "%s"...\n', cfg.label);

    % 1) saved coincident vertical profile (carries modisIndex_minDist)
    S  = load(cfg.gn_file, 'vocalsRex');
    vr = S.vocalsRex;

    % 2) ensure a full 2-D droplet size distribution Nc(radius, altitude).
    %    If Nc was stored collapsed, recover it from the raw flight by matching
    %    UTC timestamps to the size-distribution columns.
    if size(vr.Nc, 1) ~= numel(vr.drop_radius_bin_center)
        fprintf('  recovering full Nc(r,z) from raw flight %s...\n', cfg.vocals_nc);
        vr_full = readVocalsRex(cfg.vocals_nc);
        idx_t = nearest_time_idx(vr.time_utc(:), vr_full.time_utc(:));
        vr.Nc                    = vr_full.Nc(:, idx_t);
        vr.drop_radius_bin_center = vr_full.drop_radius_bin_center;
        vr.drop_radius_bin_edges  = vr_full.drop_radius_bin_edges;
    end

    % 3) coincident MODIS granule (true-colour background, footprints, r_e)
    [modis, ~] = retrieveMODIS_data(cfg.granule_dir);

    % 4) MODIS coincident-pixel summary
    idx = unique(vr.modisIndex_minDist(:));
    re_pix = modis.cloud.effRadius17(idx);
    dat.re_modis_lo   = min(re_pix);
    dat.re_modis_hi   = max(re_pix);
    dat.re_modis_mean = mean(re_pix);
    if isfield(modis, 'EV1km') && isfield(modis.EV1km, 'pixel_time_decimal')
        tpix = modis.EV1km.pixel_time_decimal(idx);
        dat.dt_min = abs(mean(tpix(:)) - mean(vr.time_utc(:))) * 60;   % min
    elseif isfield(vr, 'timeDiff_withMODIS')
        dat.dt_min = mean(vr.timeDiff_withMODIS);
    else
        dat.dt_min = NaN;
    end

    % 6) MODIS footprint polygons (true near-nadir ~1 km cells)
    n_idx = numel(idx);
    fl = cell(n_idx,1);  fo = cell(n_idx,1);
    for k = 1:n_idx
        [rr, cc] = ind2sub(size(modis.geo.lat), idx(k));
        [fl{k}, fo{k}] = modis_footprint_corners(modis, rr, cc);
    end
    dat.foot_lat = fl;  dat.foot_lon = fo;

    % 7) map extent enclosing the flight track and the footprints
    all_lat = [vr.latitude(:);  cell2mat(fl)];
    all_lon = [vr.longitude(:); cell2mat(fo)];
    pad = map_pad;                                  % larger pad -> more zoomed out
    dat.latlim = [min(all_lat)-pad, max(all_lat)+pad];
    dat.lonlim = [min(all_lon)-pad, max(all_lon)+pad];

    % 8) scene background. bg_src is resolved by resolve_bg() from the PANEL1_BG
    %    toggle and what the example provides.
    switch bg_src
        case 'gtiff'
            fprintf('  reading GeoTIFF background %s...\n', cfg.geotiff_file);
            [dat.rgb, dat.rlat, dat.rlon] = read_geotiff_bg(cfg.geotiff_file, dat.latlim, dat.lonlim);
        case 'pansharp'
            fprintf('  building pan-sharpened 250 m true colour (QKM x HKM)...\n');
            [dat.rgb, dat.rlat, dat.rlon] = read_pansharp_250m(cfg.qkm_file, cfg.hkm_file, dat.latlim, dat.lonlim);
        case 'manual'
            im = double(imread(cfg.bg_image));
            if max(im(:)) > 1; im = im / 255; end
            [ny, nx, ~] = size(im);
            [gx, gy] = meshgrid(linspace(cfg.bg_lonlim(1), cfg.bg_lonlim(2), nx), ...
                                linspace(cfg.bg_latlim(2), cfg.bg_latlim(1), ny));
            dat.rgb = im;  dat.rlon = gx;  dat.rlat = gy;
        case 'hkm500'
            fprintf('  building 500 m true-colour background from %s...\n', cfg.hkm_file);
            [dat.rgb, dat.rlat, dat.rlon] = read_hkm_truecolor(cfg.hkm_file, dat.latlim, dat.lonlim);
        case 'qkm250'
            fprintf('  building 250 m grayscale background from %s...\n', cfg.qkm_file);
            [dat.rgb, dat.rlat, dat.rlon] = read_qkm_250m(cfg.qkm_file, dat.latlim, dat.lonlim);
        otherwise   % 'g1km'
            [rgbF, latF, lonF] = create_modis_true_color(modis);
            [dat.rgb, dat.rlat, dat.rlon] = crop_region(rgbF, latF, lonF, dat.latlim, dat.lonlim);
    end

    dat.vr    = vr;
    dat.label = cfg.label;

    if ~isfolder(fileparts(cache_file)); mkdir(fileparts(cache_file)); end
    save(cache_file, 'dat', '-v7.3');
    fprintf('Cached example to: %s\n', cache_file);
end


function idx = nearest_time_idx(t_query, t_ref)
% For each query time, return the index of the nearest reference time.
    idx = zeros(numel(t_query), 1);
    for i = 1:numel(t_query)
        [~, idx(i)] = min(abs(t_ref - t_query(i)));
    end
end


function src = resolve_bg(cfg, mode)
% Resolve the panel-1 backdrop source string from the PANEL1_BG toggle and what
% files the example provides. Honoured by both the cache key and prep_example.
    has = @(f) isfield(cfg, f) && ~isempty(cfg.(f));
    switch lower(mode)
        case 'geotiff'
            if has('geotiff_file'); src = 'gtiff'; return; end
        case 'pansharp'
            if has('qkm_file') && has('hkm_file'); src = 'pansharp'; return; end
    end
    % 'auto' or requested source unavailable -> best available
    if      has('geotiff_file'); src = 'gtiff';
    elseif  has('bg_image');     src = 'manual';
    elseif  has('hkm_file');     src = 'hkm500';
    elseif  has('qkm_file');     src = 'qkm250';
    else;                        src = 'g1km';
    end
end


function [rgb, rlat, rlon] = read_pansharp_250m(qkm_f, hkm_f, latlim, lonlim)
% Pan-sharpened 250 m TRUE colour: 500 m colour (HKM bands 1/4/3) modulated by the
% 250 m band-1 detail (QKM) via a Brovey-style ratio. Gives real per-pixel 250 m
% structure with true colour -- no resample blur or JPEG blocks. No toolbox needed.

    G = '/HDFEOS/SWATHS/MODIS_SWATH_Type_L1B/';
    lat1 = double(ncread(qkm_f, [G 'Geolocation Fields/Latitude']));
    lon1 = double(ncread(qkm_f, [G 'Geolocation Fields/Longitude']));
    m = lat1 >= latlim(1)-0.05 & lat1 <= latlim(2)+0.05 & ...
        lon1 >= lonlim(1)-0.05 & lon1 <= lonlim(2)+0.05;
    [rr, cc] = find(m);
    r0 = min(rr); r1 = max(rr);  c0 = min(cc); c1 = max(cc);
    nr = r1-r0+1;  nc = c1-c0+1;

    % 250 m band-1 (high-res luminance)
    R0 = 4*r0-3;  C0 = 4*c0-3;  nR = 4*nr;  nC = 4*nc;
    pan = double(ncread(qkm_f, [G 'Data Fields/EV_250_RefSB'], [R0 C0 1], [nR nC 1]));
    pan(pan > 32767 | pan == 65535) = NaN;
    pan = 5.1927e-05 * pan;

    % 500 m colour bands (band 1 red, band 4 green, band 3 blue)
    S0 = 2*r0-1;  T0 = 2*c0-1;  nS = 2*nr;  nT = 2*nc;
    red5 = 5.1927e-05 * double(ncread(hkm_f, [G 'Data Fields/EV_250_Aggr500_RefSB'], [S0 T0 1], [nS nT 1]));
    grn5 = 3.6262e-05 * double(ncread(hkm_f, [G 'Data Fields/EV_500_RefSB'],        [S0 T0 2], [nS nT 1]));
    blu5 = 4.7776e-05 * double(ncread(hkm_f, [G 'Data Fields/EV_500_RefSB'],        [S0 T0 1], [nS nT 1]));

    % upsample 500 m colour to the 250 m grid (bilinear, via interp2)
    [Xs, Ys] = meshgrid(((1:nT)-0.5)/nT, ((1:nS)-0.5)/nS);
    [Xq, Yq] = meshgrid(((1:nC)-0.5)/nC, ((1:nR)-0.5)/nR);
    up = @(z) interp2(Xs, Ys, z, Xq, Yq, 'linear');
    red = up(red5);  grn = up(grn5);  blu = up(blu5);  panlo = up(red5);

    % Brovey-style detail injection: colour x (pan_hi / pan_lo), ratio clamped
    ratio = pan ./ panlo;  ratio(~isfinite(ratio)) = 1;  ratio = min(max(ratio, 0.5), 2);
    red = red .* ratio;  grn = grn .* ratio;  blu = blu .* ratio;

    % gentle joint stretch (shared black/white point -> no false colour) + mild gamma
    v  = [red(:); grn(:); blu(:)];  v = v(isfinite(v));
    lo = prctile(v, 1);  hi = prctile(v, 99);
    f  = @(x) min(max((x - lo) ./ (hi - lo), 0), 1);
    rgb = cat(3, f(red), f(grn), f(blu)) .^ (1/1.3);
    rgb(~isfinite(rgb)) = 0;

    % 250 m geolocation (interpolate the 1 km lat/lon up to the 250 m grid)
    p1 = 4*(r0:r1) - 1.5;   q1 = 4*(c0:c1) - 1.5;
    [Q1, P1] = meshgrid(q1, p1);
    [Qq, Pp] = meshgrid(C0:C0+nC-1, R0:R0+nR-1);
    rlat = interp2(Q1, P1, lat1(r0:r1, c0:c1), Qq, Pp, 'linear');
    rlon = interp2(Q1, P1, lon1(r0:r1, c0:c1), Qq, Pp, 'linear');
end


function [rgb, rlat, rlon] = read_geotiff_bg(tif, latlim, lonlim)
% Read a georeferenced true-colour GeoTIFF (geographic WGS84 lat/lon, e.g. a NASA
% Worldview "Corrected Reflectance" snapshot) and crop it to the region. Returns
% an RGB array (0-1) plus its lat/lon grids. Needs the Mapping Toolbox.

    [A, R] = readgeoraster(tif);
    A = double(A(:, :, 1:3)) / 255;                 % drop alpha if present
    ny = R.RasterSize(1);  nx = R.RasterSize(2);
    dLat = R.CellExtentInLatitude;  dLon = R.CellExtentInLongitude;

    % cell-centre coordinate vectors, honouring the raster's row/col orientation
    if strcmpi(R.ColumnsStartFrom, 'north')
        latv = R.LatitudeLimits(2) - ((1:ny)' - 0.5) * dLat;   % north -> south
    else
        latv = R.LatitudeLimits(1) + ((1:ny)' - 0.5) * dLat;
    end
    if strcmpi(R.RowsStartFrom, 'west')
        lonv = R.LongitudeLimits(1) + ((1:nx) - 0.5) * dLon;   % west -> east
    else
        lonv = R.LongitudeLimits(2) - ((1:nx) - 0.5) * dLon;
    end

    % crop to the region (with margin) for lean plotting
    mar = 0.05;
    ri = latv >= latlim(1)-mar & latv <= latlim(2)+mar;
    ci = lonv >= lonlim(1)-mar & lonv <= lonlim(2)+mar;
    A = A(ri, ci, :);  latv = latv(ri);  lonv = lonv(ci);
    [rlon, rlat] = meshgrid(lonv, latv);

    % Worldview Corrected Reflectance is already colour-balanced, so apply only a
    % gentle JOINT contrast stretch -- one black/white point shared by all three
    % channels. A per-channel stretch would invent hue shifts (purples) and
    % amplify the JPEG block edges of the resampled product; a joint stretch lifts
    % the low-contrast export while preserving the true colour.
    v  = A(isfinite(A));
    lo = prctile(v, 0.5);  hi = prctile(v, 99.5);
    rgb = min(max((A - lo) ./ (hi - lo), 0), 1);
end


function [rgb, rlat, rlon] = read_hkm_truecolor(f, latlim, lonlim)
% Build a 500 m TRUE-COLOUR background from a MOD02HKM/MYD02HKM L1B granule.
% True colour uses band 1 (red), band 4 (green), band 3 (blue); in the HKM file
% band 1 lives in EV_250_Aggr500_RefSB (index 1) and bands 3,4 in EV_500_RefSB
% (Band_500M = [3 4 5 6 7] -> blue=index 1, green=index 2). All are 500 m, so no
% pan-sharpening is needed. The 1 km geolocation is interpolated up to 500 m.

    G    = '/HDFEOS/SWATHS/MODIS_SWATH_Type_L1B/';
    lat1 = double(ncread(f, [G 'Geolocation Fields/Latitude']));    % [cross x along], 1 km
    lon1 = double(ncread(f, [G 'Geolocation Fields/Longitude']));

    % 1 km bounding box of the region (generous margin so the zoomed/aspect-
    % expanded view stays fully covered by image data)
    mar = 0.05;
    m = lat1 >= latlim(1)-mar & lat1 <= latlim(2)+mar & ...
        lon1 >= lonlim(1)-mar & lon1 <= lonlim(2)+mar;
    [rr, cc] = find(m);
    r0 = min(rr); r1 = max(rr);  c0 = min(cc); c1 = max(cc);

    % 500 m index ranges (each 1 km pixel k spans 500 m 2k-1..2k)
    R0 = 2*r0-1;  R1 = 2*r1;   C0 = 2*c0-1;  C1 = 2*c1;
    nR = R1-R0+1; nC = C1-C0+1;
    rd = double(ncread(f, [G 'Data Fields/EV_250_Aggr500_RefSB'], [R0 C0 1], [nR nC 1])); % band 1
    gn = double(ncread(f, [G 'Data Fields/EV_500_RefSB'],         [R0 C0 2], [nR nC 1])); % band 4
    bl = double(ncread(f, [G 'Data Fields/EV_500_RefSB'],         [R0 C0 1], [nR nC 1])); % band 3

    inval = @(x) (x == 65535) | (x > 32767);
    R = 5.1927e-05 * rd;   R(inval(rd)) = NaN;     % band-1 reflectance scale
    Gc = 3.6262e-05 * gn;  Gc(inval(gn)) = NaN;    % band-4
    B = 4.7776e-05 * bl;   B(inval(bl)) = NaN;     % band-3

    % per-channel contrast stretch (white-balances out the Rayleigh blue cast so
    % clouds read neutral white) + mild gamma to brighten, like Worldview.
    s = @(x) min(max((x - prctile(x(isfinite(x)),2)) ./ ...
        (prctile(x(isfinite(x)),98) - prctile(x(isfinite(x)),2)), 0), 1);
    rgb = cat(3, s(R), s(Gc), s(B)) .^ (1/1.5);
    rgb(~isfinite(rgb)) = 0;

    % interpolate the 1 km lat/lon up to the 500 m subset grid (1 km centre at 2k-0.5)
    p1 = 2*(r0:r1) - 0.5;   q1 = 2*(c0:c1) - 0.5;
    [Q1, P1] = meshgrid(q1, p1);
    [Qq, Pp] = meshgrid(C0:C1, R0:R1);
    rlat = interp2(Q1, P1, lat1(r0:r1, c0:c1), Qq, Pp, 'linear');
    rlon = interp2(Q1, P1, lon1(r0:r1, c0:c1), Qq, Pp, 'linear');
end


function [rgb, rlat, rlon] = read_qkm_250m(f, latlim, lonlim)
% Build a grayscale 250 m background from a MOD02QKM/MYD02QKM L1B granule
% (band 1, 645 nm reflectance) cropped to the region of interest. The file
% carries 1 km geolocation, which is interpolated up to the 250 m grid (the
% 250 m bands are nested 4x4 within each 1 km pixel). Returns an RGB array
% (grey replicated to 3 channels) plus its 250 m lat/lon grids.

    G  = '/HDFEOS/SWATHS/MODIS_SWATH_Type_L1B/';
    lat1 = double(ncread(f, [G 'Geolocation Fields/Latitude']));    % [cross x along], 1 km
    lon1 = double(ncread(f, [G 'Geolocation Fields/Longitude']));

    % 1 km bounding box of the region (with margin)
    m = lat1 >= latlim(1)-0.05 & lat1 <= latlim(2)+0.05 & ...
        lon1 >= lonlim(1)-0.05 & lon1 <= lonlim(2)+0.05;
    [rr, cc] = find(m);
    r0 = min(rr); r1 = max(rr);  c0 = min(cc); c1 = max(cc);

    % corresponding 250 m index ranges (each 1 km pixel k spans 250 m 4k-3..4k)
    R0 = 4*r0-3;  R1 = 4*r1;   C0 = 4*c0-3;  C1 = 4*c1;
    raw = double(ncread(f, [G 'Data Fields/EV_250_RefSB'], ...
        [R0 C0 1], [R1-R0+1, C1-C0+1, 1]));               % band 1 only
    bad  = (raw == 65535) | (raw > 32767);
    refl = 5.1927e-05 * raw;                              % band-1 reflectance scale
    refl(bad) = NaN;

    % interpolate the 1 km lat/lon (subset) up to the 250 m subset grid.
    % 1 km pixel k centres at 250 m index 4k-1.5.
    p1 = 4*(r0:r1) - 1.5;   q1 = 4*(c0:c1) - 1.5;
    [Q1, P1] = meshgrid(q1, p1);
    [Qq, Pp] = meshgrid(C0:C1, R0:R1);
    rlat = interp2(Q1, P1, lat1(r0:r1, c0:c1), Qq, Pp, 'linear');
    rlon = interp2(Q1, P1, lon1(r0:r1, c0:c1), Qq, Pp, 'linear');

    % contrast stretch to the patch's own 1-99 percentile range
    v  = refl(isfinite(refl));
    lo = prctile(v, 1);  hi = prctile(v, 99);
    g  = min(max((refl - lo) ./ (hi - lo), 0), 1);
    g(~isfinite(g)) = 0;
    rgb = repmat(g, 1, 1, 3);
end


function [rgbC, latC, lonC] = crop_region(rgb, lat, lon, latlim, lonlim)
% Crop a full MODIS granule's RGB + lat/lon grids to a small bounding box around
% the region of interest (with a small margin) for fast/lean plotting.
    m = lat >= latlim(1)-0.05 & lat <= latlim(2)+0.05 & ...
        lon >= lonlim(1)-0.05 & lon <= lonlim(2)+0.05;
    if ~any(m(:)); rgbC = rgb; latC = lat; lonC = lon; return; end
    [r, c] = find(m);
    rr = min(r):max(r);  cc = min(c):max(c);
    latC = lat(rr, cc);  lonC = lon(rr, cc);  rgbC = rgb(rr, cc, :);

    % Local contrast stretch: create_modis_true_color stretches over the WHOLE
    % granule (here spanning ~20 deg of latitude), so a small all-cloud patch
    % comes out dull grey. Re-stretch to the patch's own 1-99 percentile range.
    v  = rgbC(isfinite(rgbC));
    lo = prctile(v, 1);  hi = prctile(v, 99);
    if hi > lo
        rgbC = min(max((rgbC - lo) ./ (hi - lo), 0), 1);
    end
end


function lbls = deg_labels(ticks, hemi)
% Format numeric lon/lat ticks as e.g. '75.00$^\circ$W' / '24.05$^\circ$S'.
% hemi is the sign label for negative values ('W' for longitude, 'S' for lat).
    lbls = arrayfun(@(t) sprintf('%.2f$^\\circ$%s', abs(t), hemi), ticks, ...
        'UniformOutput', false);
end


function RGB = nc_to_rgb(D, cmap, clim_log)
% Map a matrix of log10(N_c) values to an explicit truecolor (RGB) array using
% the supplied colormap over [clim_log(1), clim_log(2)]. Values below the lower
% limit (empty / not-yet-sampled bins) are rendered WHITE. Returns size [M N 3].
    lo = clim_log(1);  hi = clim_log(2);  n = size(cmap, 1);
    idx = round((D - lo) ./ (hi - lo) .* (n - 1)) + 1;
    idx = min(max(idx, 1), n);
    R = cmap(:,1);  G = cmap(:,2);  B = cmap(:,3);
    rr = R(idx);  gg = G(idx);  bb = B(idx);          % each [M x N]
    empty = D < lo;                                    % below range -> white
    rr(empty) = 1;  gg(empty) = 1;  bb(empty) = 1;
    RGB = cat(3, rr, gg, bb);
end


function [lat_corners, lon_corners] = modis_footprint_corners(modis, row, col)
% True near-nadir MODIS 1 km footprint corners, scaled by sensor zenith.
% (Replicates the logic of get_modis_footprint in plot_instrument_footprints_*.m)

    center_lat = modis.geo.lat(row, col);
    center_lon = modis.geo.long(row, col);

    if isfield(modis, 'sensor') && isfield(modis.sensor, 'zenith')
        sensor_zenith = double(modis.sensor.zenith(row, col));   % degrees
    else
        sensor_zenith = 0;
    end

    % MODIS pixel grows with scan angle: ~1 km x 1 km at nadir.
    along_track_size = 1.0 * (1 + 0.02 * (sensor_zenith/10)^2);  % km
    along_scan_size  = 1.0 * (1 / cosd(sensor_zenith));          % km

    lat_per_km = 1.0 / 111.0;
    lon_per_km = 1.0 / (111.0 * cosd(center_lat));

    half_at = (along_track_size / 2) * lat_per_km;
    half_as = (along_scan_size  / 2) * lon_per_km;

    lat_corners = [center_lat+half_at; center_lat+half_at; ...
                   center_lat-half_at; center_lat-half_at];
    lon_corners = [center_lon-half_as; center_lon+half_as; ...
                   center_lon+half_as; center_lon-half_as];
end
