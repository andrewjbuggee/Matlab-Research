%% Solar geometry polar density plot for VOCALS-REx and ORACLES in-situ profiles
%
% Steps through every stored ensemble profile from VOCALS-REx and ORACLES,
% computes the solar zenith and azimuth angle at the mid-point of each
% profile using solar_position_from_lat_lon_time, and renders a polar
% density plot of the joint (azimuth, zenith) distribution.
%
% Conventions in the figure:
%   - Radial axis: solar zenith angle (deg). 0 = sun overhead at center.
%   - Angular axis: compass azimuth (deg). 0 = North at top, increases
%     clockwise. The libRadTran azimuth (clockwise from south) returned
%     by solar_position_from_lat_lon_time is converted here.
%
% Day-of-year is captured implicitly: solar_position_from_lat_lon_time
% computes Julian-day-since-J2000 from the full datetime, which encodes
% year/month/day (and therefore DOY) for the solar declination, mean
% longitude, and mean anomaly terms.
%
% Profiles with the sun below the horizon (SZA >= 90 deg) are excluded
% from the polar density plot — they are pre-sunrise / post-sunset
% samples (e.g. several VOCALS-REx flights take off at 06:00 UTC,
% which is ~01:00 local Chile time, well before sunrise) and are not
% useful for solar-geometry context. The number dropped is reported.

clear variables
close all

%% File locations

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    foldername_oracles = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_vocals  = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    folderpath_figs    = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/saved_figures/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_oracles = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_vocals  = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    folderpath_figs    = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/saved_figures/'];

else
    error('Unrecognized computer in whatComputer.')
end

if ~exist(folderpath_figs, 'dir')
    mkdir(folderpath_figs);
end

oracles_file = ['ensemble_profiles_with_precip_from_33_files_', ...
    'LWC-threshold_0.05_Nc-threshold_10_16-Mar-2026.mat'];

vocals_file  = ['ensemble_profiles_with_precip_from_14_files_', ...
    'LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat'];


%% Compute solar geometry for every profile in both campaigns

sza_all = [];
phi0_lr_all = [];   % libRadTran azimuth (clockwise from south)
src_all = {};       % campaign label per profile

% --- VOCALS-REx ---
S = load(fullfile(foldername_vocals, vocals_file), 'ensemble_profiles');
vocals_profiles = S.ensemble_profiles;
N_v = numel(vocals_profiles);

[sza_v, phi0_v] = compute_geometry_from_profiles(vocals_profiles);

sza_all     = [sza_all,     sza_v];
phi0_lr_all = [phi0_lr_all, phi0_v];
src_all     = [src_all,     repmat({'VOCALS-REx'}, 1, numel(sza_v))];

disp(['VOCALS-REx profiles processed: ', num2str(numel(sza_v)), ' / ', num2str(N_v)])

% --- ORACLES ---
S = load(fullfile(foldername_oracles, oracles_file), 'ensemble_profiles');
oracles_profiles = S.ensemble_profiles;
N_o = numel(oracles_profiles);

[sza_o, phi0_o] = compute_geometry_from_profiles(oracles_profiles);

sza_all     = [sza_all,     sza_o];
phi0_lr_all = [phi0_lr_all, phi0_o];
src_all     = [src_all,     repmat({'ORACLES'}, 1, numel(sza_o))];

disp(['ORACLES profiles processed:    ', num2str(numel(sza_o)), ' / ', num2str(N_o)])
disp(['Total profiles:                ', num2str(numel(sza_all))])


%% Convert azimuth from libRadTran convention to compass convention
% libRadTran: clockwise from south (0=S, 90=W, 180=N, 270=E)
% Compass:    clockwise from north (0=N, 90=E, 180=S, 270=W)

az_compass = mod(phi0_lr_all + 180, 360);   % degrees


%% Drop sun-below-horizon profiles (SZA >= 90)
% These are pre-sunrise / post-sunset samples (mostly the early profiles
% of VOCALS-REx flights that start at 06:00 UTC, which is ~01:00 local
% in Arica, hours before sunrise). They are valid in-situ measurements
% but not meaningful for a solar-geometry density plot.

below_horizon = sza_all >= 90;
N_below = sum(below_horizon);
if N_below > 0
    disp(['Dropping ', num2str(N_below), ' profiles with SZA >= 90 deg ', ...
        '(sun below horizon).'])
    % Per-campaign breakdown
    src_below = src_all(below_horizon);
    disp(['  VOCALS-REx dropped: ', num2str(sum(strcmp(src_below, 'VOCALS-REx')))])
    disp(['  ORACLES dropped:    ', num2str(sum(strcmp(src_below, 'ORACLES')))])
end

sza_all     = sza_all(~below_horizon);
phi0_lr_all = phi0_lr_all(~below_horizon);
az_compass  = az_compass(~below_horizon);
src_all     = src_all(~below_horizon);


%% 2D histogram in polar bins

n_az_bins  = 36;                     % 10 deg per azimuth bin
n_sza_bins = 30;                     % 3 deg per zenith bin
sza_max    = 90;

az_edges  = linspace(0, 360, n_az_bins + 1);
sza_edges = linspace(0, sza_max, n_sza_bins + 1);

az_centers  = 0.5 * (az_edges(1:end-1) + az_edges(2:end));
sza_centers = 0.5 * (sza_edges(1:end-1) + sza_edges(2:end));

counts = histcounts2(az_compass(:), sza_all(:), az_edges, sza_edges);

% Build a Cartesian mesh from the polar grid edges so pcolor can render the
% density on a regular axes (polaraxes does not support pcolor).
% Compass azimuth -> math angle: theta_math = pi/2 - deg2rad(az)
[AzE, RE]  = meshgrid(az_edges, sza_edges);
ThetaE     = pi/2 - deg2rad(AzE);
X_edges    = RE .* cos(ThetaE);
Y_edges    = RE .* sin(ThetaE);


%% Paper-professional polar density figure

fnt_sz   = 18;
ttl_fnt  = 22;
lbl_fnt  = 16;

figure1 = figure('Color', 'w');

ax = axes('Parent', figure1);
hold(ax, 'on');

% pcolor expects Z to be on the same grid as X,Y — counts are bin-valued
% and have size [n_sza_bins x n_az_bins]; pad to match the edge mesh.
Z = zeros(size(X_edges));
Z(1:end-1, 1:end-1) = counts';     % rows = sza bins, cols = az bins

h = pcolor(ax, X_edges, Y_edges, Z);
set(h, 'EdgeColor', 'none');
shading(ax, 'flat');

% Mask exterior to a clean disk by setting the axes background to white
axis(ax, 'equal');
axis(ax, 'off');
xlim(ax, [-sza_max sza_max] * 1.08);
ylim(ax, [-sza_max sza_max] * 1.08);

% Colormap: perceptually uniform with strong contrast for density
colormap(ax, turbo);

% Logarithmic-like scaling so low-density tails are visible without
% washing out the peak. Linear cap at zero for empty bins.
caxis_max = max(counts(:));
if caxis_max > 0
    set(ax, 'CLim', [0, caxis_max]);
end

% --- Radial reference circles (constant SZA) ---
sza_rings = [15, 30, 45, 60, 75, 90];
theta_ring = linspace(0, 2*pi, 361);
for r = sza_rings
    plot(ax, r*cos(theta_ring), r*sin(theta_ring), '-', ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 0.6);
end
% Outer boundary
plot(ax, sza_max*cos(theta_ring), sza_max*sin(theta_ring), '-', ...
    'Color', [0 0 0], 'LineWidth', 1.2);

% --- Angular reference spokes (every 30 deg compass az) ---
spoke_az_deg = 0:30:330;
for a = spoke_az_deg
    th = pi/2 - deg2rad(a);
    plot(ax, [0 sza_max*cos(th)], [0 sza_max*sin(th)], '-', ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 0.5);
end

% --- Radial axis labels (SZA values) along the 60 deg compass spoke
%     for legibility ---
label_az_deg = 120;
th_lbl = pi/2 - deg2rad(label_az_deg);
for r = sza_rings
    text(ax, (r+1.5)*cos(th_lbl), (r+1.5)*sin(th_lbl), ...
        sprintf('%d^{\\circ}', r), ...
        'FontSize', lbl_fnt-3, 'Color', 'w', ...
        'HorizontalAlignment', 'left', 'BackgroundColor', 'none');
end

% --- Cardinal direction labels ---
r_lbl = sza_max * 1.05;
text(ax,  0,        r_lbl,  'N',  'FontSize', lbl_fnt, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax,  r_lbl,    0,      'E',  'FontSize', lbl_fnt, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left',   'VerticalAlignment', 'middle');
text(ax,  0,       -r_lbl,  'S',  'FontSize', lbl_fnt, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(ax, -r_lbl,    0,      'W',  'FontSize', lbl_fnt, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right',  'VerticalAlignment', 'middle');

% --- Title and colorbar ---
title(ax, 'Solar geometry of in-situ profiles (VOCALS-REx + ORACLES)', ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt);

cb = colorbar(ax);
cb.Label.String      = 'Profile count';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize    = fnt_sz;
cb.Label.Color       = 'k';
cb.TickLabelInterpreter = 'latex';
cb.Color             = 'k';
cb.FontSize          = fnt_sz - 4;


% --- Summary text box (counts and ranges) ---
N_total = numel(sza_all);
txt = {sprintf('$N_{\\mathrm{profiles}} = %d$', N_total), ...
       sprintf('VOCALS-REx: $%d$', sum(strcmp(src_all, 'VOCALS-REx'))), ...
       sprintf('ORACLES:   $%d$', sum(strcmp(src_all, 'ORACLES'))), ...
       sprintf('SZA range: $%.1f^{\\circ}$ - $%.1f^{\\circ}$', ...
            min(sza_all), max(sza_all))};
annotation(figure1, 'textbox', [0.02 0.18 0.22 0.16], ...
    'String', txt, 'Interpreter', 'latex', ...
    'FontSize', lbl_fnt-2, 'EdgeColor', 'none', ...
    'BackgroundColor', 'none', 'FitBoxToText', 'on',...
    'Color', 'k');

% Figure size: square, paper-friendly
w = 6.5;   % inches
h = 6.5;   % inches
figure1.Units    = 'inches';
figure1.Position = [1, 1, 1.4*w, 1.4*h];


%% Save figure

base_name = 'solar_geometry_polar_density_VOCALS_and_ORACLES';

saveas(figure1, fullfile(folderpath_figs, [base_name, '.fig']));
exportgraphics(figure1, fullfile(folderpath_figs, [base_name, '.png']), ...
    'Resolution', 500);

disp(['Saved figure to: ', fullfile(folderpath_figs, [base_name, '.png'])])


%% --- helper -----------------------------------------------------------

function [sza_vec, phi0_vec] = compute_geometry_from_profiles(profiles)
%COMPUTE_GEOMETRY_FROM_PROFILES Compute (sza, phi0) for each ensemble profile.
%   Uses the mid-point sample of each profile for date / time / lat / lon.

N = numel(profiles);
sza_vec  = nan(1, N);
phi0_vec = nan(1, N);

for nn = 1:N

    p = profiles{nn};

    if ~isfield(p, 'latitude') || ~isfield(p, 'longitude') || ...
            ~isfield(p, 'dateOfFlight')
        continue
    end

    lat_vec = p.latitude;
    lon_vec = p.longitude;

    if isempty(lat_vec) || isempty(lon_vec)
        continue
    end

    mid_idx = max(1, round(numel(lat_vec)/2));

    % Decimal-hour UTC time at the profile mid-point. Prefer the explicit
    % time_utc field (VOCALS-REx); otherwise derive from the seconds-since-
    % midnight time field (ORACLES).
    if isfield(p, 'time_utc') && ~isempty(p.time_utc)
        t_utc = p.time_utc(min(mid_idx, numel(p.time_utc)));
    elseif isfield(p, 'time') && ~isempty(p.time)
        t_utc = p.time(min(mid_idx, numel(p.time))) / 3600;   % s -> hr
    else
        continue
    end

    lat = lat_vec(min(mid_idx, numel(lat_vec)));
    lon = lon_vec(min(mid_idx, numel(lon_vec)));

    [sza_vec(nn), phi0_vec(nn)] = solar_position_from_lat_lon_time( ...
        p.dateOfFlight, t_utc, lat, lon);
end

% Drop any profiles where computation failed
ok = ~isnan(sza_vec) & ~isnan(phi0_vec);
sza_vec  = sza_vec(ok);
phi0_vec = phi0_vec(ok);

end
