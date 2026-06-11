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
%       Panel 1 (left)   : geographic map. Aircraft marker + trailing flight
%                          track build as the plane climbs.
%       Panel 2 (middle) : altitude (m) vs. along-track distance (km). The slant
%                          ascent draws progressively with a marker at the
%                          current aircraft position.
%       Panel 3 (right)  : droplet size distribution log10(N_c) over (radius,
%                          altitude). Rows fill in from cloud base upward as the
%                          plane samples each level; the in-situ r_e(z) line
%                          draws on top.
%
%   ACT 2  -- "the coincident MODIS measurement" (fades in, then holds):
%       Panel 1 : the coincident MODIS pixel footprints fade in at true 1 km
%                 scale, coloured by retrieved r_e (bands 1+7), next to the thin
%                 in-situ track -> footprint vs. point-sampling contrast.
%       Panel 3 : the single MODIS column-averaged r_e fades in as a vertical
%                 line (with the spread across coincident pixels as a band)
%                 against the full in-situ r_e(z) curve.
%
%   Data product (single saved file, holds BOTH the full droplet distribution
%   and the MODIS coincidence index):
%       VOCALS-REx vertical profile, 2008-11-11, Aqua overpass 18:50 UTC.
%       Spatial match: all coincident MODIS pixels within ~0.7 km (sub-pixel).
%       Time offset to MODIS: 8.91 min.
%
%   Output: an MPEG-4 movie written next to this script.
%
% By Andrew John Buggee
%%

clear variables; close all;

% =========================================================================
% ---------------------------- CONFIG -------------------------------------
% =========================================================================

% Set QUICK_TEST = true to render a coarse, fast movie (for checking layout);
% set to false for the final high-frame-count render.
QUICK_TEST = false;

% --- frame counts and pacing ---
if QUICK_TEST
    Nf_act1    = 24;     % frames to build the in-situ profile
    Nf_reveal  = 14;     % frames to fade in the MODIS overlay
    Nf_hold    = 8;      % frames to hold the final comparison
    frame_rate = 12;
else
    Nf_act1    = 165;    % frames to build the in-situ profile (~7 s)
    Nf_reveal  = 60;     % frames to fade in the MODIS overlay (~2.5 s)
    Nf_hold    = 72;     % frames to hold the final comparison (~3 s)
    frame_rate = 24;
end

% --- droplet-radius display cap (microns). The cloud mode and r_e (~6-9 um)
%     live below this; the sparse drizzle tail above it is suppressed so the
%     evolving cloud distribution and r_e are clearly visible. ---
r_display_max_um = 20;

% --- fonts / styling ---
ax_fnt  = 20;
ttl_fnt = 23;
cb_fnt  = 20;
ln_sz   = 2;
geoLineSz = 2;
tck_label_sz = 18;

% --- paths ---
repo = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research';
addpath(genpath([repo '/Generally_useful_functions']));
addpath(genpath([repo '/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval']));

vr_file = [repo '/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/', ...
    '2008_11_11_1850/Retrieval_outputs_05-Apr-2025/', ...
    'GN_inputs_outputs_withAdvection_rt-cov_8.76_rb-cov_100_tc-cov_100_05-Apr-2025_rev1.mat'];

% NOTE: retrieveMODIS_data globs [folderName,'*.hdf'] -> the trailing slash is required.
granule_dir = [repo '/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/'];

out_dir = [repo '/Presentations_and_Papers/PhD_defense/'];
out_name = 'cloud_vertical_structure_vs_MODIS.mp4';


% =========================================================================
% ---------------------------- LOAD DATA ----------------------------------
% =========================================================================

fprintf('Loading VOCALS-REx profile...\n');
S  = load(vr_file, 'vocalsRex');
vr = S.vocalsRex;

fprintf('Reloading coincident MODIS (Aqua) granule...\n');
modis = retrieveMODIS_data(granule_dir);

% --- in-situ profile, ordered cloud base -> cloud top (altitude ascending) ---
[altitude_m, sort_idx] = sort(vr.altitude(:));            % m, ascending
n_alt   = numel(altitude_m);
re_um   = vr.re(sort_idx);                                % effective radius, um
re_um   = re_um(:);
lat_deg = vr.latitude(sort_idx);   lat_deg = lat_deg(:);  % deg
lon_deg = vr.longitude(sort_idx);  lon_deg = lon_deg(:);  % deg

% along-track horizontal distance (km), monotonic with the climb
horz_dist_km = vr.horz_dist(sort_idx);  horz_dist_km = horz_dist_km(:) / 1000;   % m -> km

% --- droplet size distribution: Nc is [n_radius x n_altitude] ---
r_bins_um = double(vr.drop_radius_bin_center(:));         % um, 91 bins
Nc        = double(vr.Nc(:, sort_idx));                   % cm^-3, [n_radius x n_alt]

% cap displayed radius range
r_mask        = r_bins_um <= r_display_max_um;
r_plot_um     = r_bins_um(r_mask);
Nc_plot       = Nc(r_mask, :);

% log10 for display. Empty (zero-count) and not-yet-revealed cells are mapped
% to a below-range "sentinel" value that renders WHITE via a white-topped
% colormap. This is renderer-independent (transparency/AlphaData composites
% to black under headless capture), so it looks identical in batch and GUI.
Nc_log = log10(Nc_plot);
Nc_log(~isfinite(Nc_log)) = NaN;
Nc_disp_full = Nc_log';                                   % [n_alt x n_radius]
clim_log   = [max(-1, min(Nc_disp_full(:))), max(Nc_disp_full(:))];
white_val  = clim_log(1) - 1;                             % sentinel -> white
Nc_disp_full(isnan(Nc_disp_full)) = white_val;            % empty bins -> white

% --- MODIS coincident pixels ---
modis_idx   = unique(vr.modisIndex_minDist(:));
re_modis_um = modis.cloud.effRadius17(modis_idx);         % um, one per coincident pixel
re_modis_mean_um = mean(re_modis_um);
re_modis_lo_um   = min(re_modis_um);
re_modis_hi_um   = max(re_modis_um);
modis_dt_min     = 8.91;                                  % time offset to overpass (min)

fprintf('In-situ: %d levels, %.0f-%.0f m, r_e %.1f-%.1f um, track %.2f km\n', ...
    n_alt, altitude_m(1), altitude_m(end), min(re_um), max(re_um), horz_dist_km(end));
fprintf('MODIS  : %d coincident pixels, r_e = %.2f um (%.2f-%.2f)\n', ...
    numel(modis_idx), re_modis_mean_um, re_modis_lo_um, re_modis_hi_um);

% --- MODIS footprint polygons (true near-nadir ~1 km cells) ---
foot_lat = cell(numel(modis_idx),1);
foot_lon = cell(numel(modis_idx),1);
for k = 1:numel(modis_idx)
    [rr, cc] = ind2sub(size(modis.geo.lat), modis_idx(k));
    [foot_lat{k}, foot_lon{k}] = modis_footprint_corners(modis, rr, cc);
end

% --- map extent: enclose the flight track and the MODIS footprints ---
all_lat = [lat_deg; cell2mat(foot_lat)];
all_lon = [lon_deg; cell2mat(foot_lon)];
pad_lat = 0.012;  pad_lon = 0.012;
latlim  = [min(all_lat)-pad_lat, max(all_lat)+pad_lat];
lonlim  = [min(all_lon)-pad_lon, max(all_lon)+pad_lon];


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

% --- Panel 1: geographic map (left) ---
have_mapping = exist('geoaxes', 'file') == 2 || exist('geoaxes', 'builtin') == 5;
if have_mapping
    gx = geoaxes('Parent', fig, 'Position', [0.035 0.13 0.275 0.74]);
    try
        geobasemap(gx, 'satellite');
    catch
        geobasemap(gx, 'grayland');   % offline fallback
    end
    geolimits(gx, latlim, lonlim);
    hold(gx, 'on');
    gx.LatitudeLabel.String  = 'Latitude';
    gx.LongitudeLabel.String = 'Longitude';
    gx.LatitudeLabel.Interpreter  = 'latex';
    gx.LongitudeLabel.Interpreter = 'latex';
    gx.LatitudeLabel.FontSize  = ax_fnt;
    gx.LongitudeLabel.FontSize = ax_fnt;
    % gx.FontSize = 13;
    gx.Title.Interpreter = 'latex';
    gx.Title.FontSize = ttl_fnt;
    title(gx, 'Aircraft track \& MODIS footprint');
else
    % Plain lat/lon axes fallback (no Mapping Toolbox)
    gx = axes('Parent', fig, 'Position', [0.05 0.13 0.255 0.74]);
    hold(gx, 'on'); box(gx, 'on'); grid(gx, 'on');
    xlim(gx, lonlim); ylim(gx, latlim);
    xlabel(gx, 'Longitude', 'Interpreter', 'latex', 'FontSize', ax_fnt);
    ylabel(gx, 'Latitude',  'Interpreter', 'latex', 'FontSize', ax_fnt);
    title(gx, 'Aircraft track \& MODIS footprint', 'Interpreter', 'latex', 'FontSize', ttl_fnt);
end

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
% Truecolor (explicit RGB per cell) rendering. A colormap-based image (imagesc or
% pcolor) is corrupted by the geoaxes/basemap sharing the figure (black edge bands
% + a mangled colorbar under getframe); mapping to RGB ourselves sidesteps the
% shared figure colormap entirely and renders cleanly next to the map.
RGB_full = nc_to_rgb(Nc_disp_full, turbo_cmap, clim_log);          % [n_alt x n_r x 3]
img = image(ax3, r_plot_um, altitude_m, ones([size(Nc_disp_full) 3]));   % start white
set(ax3, 'YDir', 'normal');
hold(ax3, 'on');
% Limits must hug the image extent exactly: any axes area outside the image
% renders black (not white) when a geoaxes/basemap shares the figure.
xlim(ax3, [min(r_plot_um), max(r_plot_um)]);
ylim(ax3, [altitude_m(1), altitude_m(end)]);

% manual truecolor colorbar strip (also RGB -> independent of the figure colormap)
cax = axes('Parent', fig, 'Position', [0.915 0.13 0.013 0.74], 'Color', 'w',...
    'FontSize', 18);
n_strip   = 256;
strip_rgb = nc_to_rgb(linspace(clim_log(1), clim_log(2), n_strip)', turbo_cmap, clim_log);
image(cax, [0 1], clim_log, strip_rgb);
set(cax, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', ...
    'TickLabelInterpreter', 'latex', 'FontSize', 13);
ylabel(cax, '$\log_{10}(N_c)\;\;(\mathrm{cm}^{-3})$', 'Interpreter', 'latex', 'FontSize', cb_fnt);
xlabel(ax3, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt);
ylabel(ax3, 'Altitude (m)',            'Interpreter', 'latex', 'FontSize', ax_fnt);
title(ax3,  'Droplet size distribution', 'Interpreter', 'latex', 'FontSize', ttl_fnt);
set(ax3, 'FontSize', tck_label_sz, 'TickLabelInterpreter', 'latex');

% animated graphics handles -------------------------------------------------
if have_mapping
    track_h = geoplot(gx, lat_deg(1), lon_deg(1), '-', 'Color', clr_track, 'LineWidth', geoLineSz);
    plane_h = geoplot(gx, lat_deg(1), lon_deg(1), '^', 'MarkerFaceColor', clr_plane, ...
        'MarkerEdgeColor', 'w', 'MarkerSize', 12, 'LineWidth', 1);
else
    track_h = plot(gx, lon_deg(1), lat_deg(1), '-', 'Color', clr_track, 'LineWidth', geoLineSz);
    plane_h = plot(gx, lon_deg(1), lat_deg(1), '^', 'MarkerFaceColor', clr_plane, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 12);
end

alt_line_h   = plot(ax2, horz_dist_km(1), altitude_m(1), '-', 'Color', clr_track, 'LineWidth', ln_sz);
alt_marker_h = plot(ax2, horz_dist_km(1), altitude_m(1), '^', 'MarkerFaceColor', clr_plane, ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 12);

re_line_h = plot(ax3, re_um(1), altitude_m(1), '-', 'Color', clr_re, 'LineWidth', ln_sz, ...
    'DisplayName', 'in-situ $r_e(z)$');

% MODIS overlay handles (created hidden, revealed in Act 2) -----------------
foot_h = gobjects(numel(modis_idx),1);
for k = 1:numel(modis_idx)
    if have_mapping
        foot_h(k) = geoplot(gx, foot_lat{k}([1:end 1]), foot_lon{k}([1:end 1]), '-', ...
            'Color', clr_modis, 'LineWidth', 2, 'Visible', 'off');
    else
        foot_h(k) = plot(gx, foot_lon{k}([1:end 1]), foot_lat{k}([1:end 1]), '-', ...
            'Color', clr_modis, 'LineWidth', 2, 'Visible', 'off');
    end
end

% MODIS r_e on the distribution panel: mean as a dashed line, with the spread
% across coincident pixels as thin dotted lines. No transparency -> robust under
% getframe alongside the map. All hidden until Act 2.
yb = [altitude_m(1), altitude_m(end)];
% modis_lo_h = plot(ax3, [re_modis_lo_um re_modis_lo_um], yb, ':', 'Color', clr_modis, ...
%     'LineWidth', 1.5, 'Visible', 'off');
% modis_hi_h = plot(ax3, [re_modis_hi_um re_modis_hi_um], yb, ':', 'Color', clr_modis, ...
%     'LineWidth', 1.5, 'Visible', 'off');
modis_re_h = plot(ax3, [re_modis_mean_um re_modis_mean_um], yb, ':', 'Color', 'k', ...
    'LineWidth', ln_sz+1, 'Visible', 'off', 'DisplayName', 'MODIS $r_{2.1}$');

% annotation handles (created in Act 2)
ann_foot  = gobjects(0);
ann_redif = gobjects(0);

% let the basemap tiles load before capturing frames
drawnow; pause(2.5);


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
    if have_mapping
        set(track_h, 'LatitudeData', lat_deg(1:k), 'LongitudeData', lon_deg(1:k));
        set(plane_h, 'LatitudeData', lat_i, 'LongitudeData', lon_i);
    else
        set(track_h, 'XData', lon_deg(1:k), 'YData', lat_deg(1:k));
        set(plane_h, 'XData', lon_i, 'YData', lat_i);
    end

    % altitude panel: slant ascent + marker
    set(alt_line_h,   'XData', [horz_dist_km(1:k); hd_i], 'YData', [altitude_m(1:k); alt_i]);
    set(alt_marker_h, 'XData', hd_i, 'YData', alt_i);

    % distribution panel: reveal rows up to current level, draw r_e(z)
    rgb = RGB_full;
    rgb(k+1:end, :, :) = 1;                 % not-yet-sampled levels stay white
    set(img, 'CData', rgb);
    set(re_line_h, 'XData', re_um(1:k), 'YData', altitude_m(1:k));

    drawnow;
    writeVideo(v, getframe(fig));
end

% ensure the full profile is shown
set(img, 'CData', RGB_full);
set(re_line_h, 'XData', re_um, 'YData', altitude_m);
legend(ax3, re_line_h, 'Location', 'northeast', 'Interpreter', 'latex', ...
    'FontSize', 18, 'Color', 'w', 'TextColor', 'k');

% -------------------------------------------------------------------------
% ACT 2 -- fade in the coincident MODIS measurement
% -------------------------------------------------------------------------
fprintf('Rendering Act 2 reveal (%d frames)...\n', Nf_reveal);

for k = 1:numel(foot_h); set(foot_h(k), 'Visible', 'on'); end
set(modis_re_h, 'Visible', 'on');

% Act-2 annotations

legend(ax3, [re_line_h, modis_re_h], 'Location', 'southeast', 'Interpreter', 'latex', ...
    'FontSize', 18, 'Color', 'w', 'TextColor', 'k');


str_foot = sprintf('MODIS pixel $\\approx$ 1 km\nin-situ track = %.1f km', horz_dist_km(end));
ann_foot = annotation(fig, 'textbox', [0.0375 0.4 0.27 0.10], 'String', str_foot, ...
    'Interpreter', 'latex', 'FontSize', 18, 'EdgeColor', clr_modis, 'LineWidth', 1.2, ...
    'BackgroundColor', 'w', 'Color', clr_modis, 'FitBoxToText', 'on', 'Visible', 'off');

str_re = sprintf('MODIS $r_{2.1} = %.1f\\;\\mu$m\n($\\Delta t = %.1f$ min)', re_modis_mean_um, modis_dt_min);
% ann_redif = annotation(fig, 'textbox', [0.775 0.16 0.15 0.10], 'String', str_re, ...
%     'Interpreter', 'latex', 'FontSize', 18, 'EdgeColor', clr_modis, 'LineWidth', 1.2, ...
%     'BackgroundColor', 'w', 'Color', clr_modis, 'FitBoxToText', 'on', 'Visible', 'off');

% set(modis_lo_h, 'Visible', 'on');
% set(modis_hi_h, 'Visible', 'on');

for f = 1:Nf_reveal
    a = f / Nf_reveal;                          % 0 -> 1 emphasis ramp
    for k = 1:numel(foot_h)
        set(foot_h(k), 'LineWidth', 2 + 1.5*a);
    end
    if f == 1
        set(ann_foot,  'Visible', 'on');
        set(ann_redif, 'Visible', 'on');
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
