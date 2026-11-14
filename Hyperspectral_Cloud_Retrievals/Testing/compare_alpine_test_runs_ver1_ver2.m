%% Compare HySICS spectra from the ver1 and ver2 synthetic drivers on Alpine
%
% Matches clouds that were computed by both Alpine test runs and compares
% their reflectance spectra band by band:
%
%   ver1  hysics_refl_from_synthetic_NN_inputs.m       (13 Nov 2026, job 35667061)
%         one .mat per cloud, 636 per-band uvspec runs
%   ver2  run_synthetic_spectra_chunk.m                (14 Nov 2026, job 35700478)
%         one chunk .mat for many clouds, 4 uvspec runs per cloud
%
% Both runs read the same input file (training_inputs_jointMVN_N100000_L7_
% 2026-11-12_ver2.nc, MD5 f619a48c...), so a cloud index means the same cloud,
% geometry and atmosphere in both. The script also checks that the geometry
% stored by each driver agrees before comparing spectra.
%
% Outputs (written to OUT_DIR):
%   alpine_ver1_vs_ver2_spectra.png      one panel per shared cloud, spectra overlaid
%   alpine_ver1_vs_ver2_rel_diff.png     relative difference vs wavelength, all clouds
% and a table of per-cloud differences printed to the Command Window.
%
% Figures are drawn invisibly with MATLAB's factory graphics defaults
% (new_figure, at the end of this file) and exported, because startup.m sets
% dark root defaults and exporting visible figures can stall R2024b. When run
% from the desktop, the exported PNGs are opened afterwards.
%
% By Andrew John Buggee; written with Claude (assistant), 14 Nov 2026
%%

clear; close all

% ---- Configuration ------------------------------------------------------
% Folders holding each run's output. If a folder is missing but a .zip of the
% same name exists next to it (as downloaded from CURC OnDemand), it is
% unzipped first.
DIR_VER1 = fullfile(getenv('HOME'), 'Downloads', 'synthetic_dataSet_TEST_12_Nov_2026');
DIR_VER2 = fullfile(getenv('HOME'), 'Downloads', 'synthetic_dataSet_TEST_14_Nov_2026');
OUT_DIR  = fullfile(getenv('HOME'), 'Downloads');      % where the PNGs are written

% Okabe-Ito colour-blind-safe pair, the same as the laptop comparison
C_VER1 = [0 114 178] / 255;      % blue
C_VER2 = [213 94 0]  / 255;      % vermillion


%% Make sure both folders exist (unzip a downloaded archive if needed)
for d = {DIR_VER1, DIR_VER2}
    folder = d{1};
    if ~isfolder(folder)
        zip_file = [folder, '.zip'];
        if isfile(zip_file)
            fprintf('Unzipping %s\n', zip_file);
            unzip(zip_file, folder);
        else
            error('Neither %s nor %s exists.', folder, zip_file);
        end
    end
end


%% Load ver1: one file per cloud
% The cloud index is in the file name ('..._cloud553_sza20_...'). Each file
% holds the clean reflectance (636 x 1) and changing_variables, whose first
% four columns are [sza vza vaz saz] in degrees (identical on every row).
f1 = dir(fullfile(DIR_VER1, 'simulated_spectra_HySICS_reflectance_636bands_*_cloud*.mat'));
n1 = numel(f1);
ver1.idx  = zeros(n1, 1);
ver1.refl = zeros(n1, 636);
ver1.geom = zeros(n1, 4);            % [sza vza vaz saz], deg
for k = 1:n1
    tok = regexp(f1(k).name, '_cloud(\d+)_', 'tokens', 'once');
    ver1.idx(k) = str2double(tok{1});
    S = load(fullfile(f1(k).folder, f1(k).name), ...
        'Refl_model_allStateVectors', 'changing_variables_allStateVectors');
    ver1.refl(k, :) = S.Refl_model_allStateVectors(:)';
    ver1.geom(k, :) = S.changing_variables_allStateVectors(1, 1:4);
end
fprintf('ver1: %d clouds in %s\n', n1, DIR_VER1);


%% Load ver2: chunk file(s), many clouds each
% A '_partial.mat' file is an interrupted checkpoint write and is skipped.
f2 = dir(fullfile(DIR_VER2, 'synthetic_spectra_HySICS_636bands_clouds*.mat'));
f2 = f2(~endsWith({f2.name}, '_partial.mat'));
ver2.idx = []; ver2.refl = []; ver2.geom = [];
for k = 1:numel(f2)
    C = load(fullfile(f2(k).folder, f2(k).name));
    ok = logical(C.done(:));                         % clouds whose RT finished
    ver2.idx  = [ver2.idx;  C.cloud_idx(ok)];
    ver2.refl = [ver2.refl; C.Refl_model(ok, :)];
    ver2.geom = [ver2.geom; [C.sza(ok), C.vza(ok), C.vaz(ok), C.saz(ok)]];
    wl_nm = C.band_center_wavelength_nm(:);          % (636 x 1) nm, band centres
end
fprintf('ver2: %d finished clouds in %d chunk file(s) in %s\n', numel(ver2.idx), numel(f2), DIR_VER2);


%% Match clouds computed by both runs
[shared, i1, i2] = intersect(ver1.idx, ver2.idx, 'stable');
if isempty(shared)
    error('No cloud index appears in both runs.');
end
R1 = ver1.refl(i1, :);           % (n_shared x 636)
R2 = ver2.refl(i2, :);
n  = numel(shared);
fprintf('Clouds in both runs: %s\n\n', mat2str(shared'));

% Same input file, so the geometry must agree; a mismatch would mean the two
% runs did not use the same .nc (e.g. the pre-azimuth-fix file).
geom_diff_deg = max(abs(ver1.geom(i1, :) - ver2.geom(i2, :)), [], 2);
if any(geom_diff_deg > 1e-3)
    warning('Geometry differs between runs for clouds %s (max %.3g deg)', ...
        mat2str(shared(geom_diff_deg > 1e-3)'), max(geom_diff_deg));
end


%% Per-cloud differences
rel_diff = (R2 - R1) ./ abs(R1);                     % (n x 636), dimensionless
max_abs_rel_diff = max(abs(rel_diff), [], 2);
rms_rel_diff     = sqrt(mean(rel_diff.^2, 2));
n_bands_differ   = sum(R2 ~= R1, 2);

T = table(shared(:), geom_diff_deg, n_bands_differ, max_abs_rel_diff, rms_rel_diff, ...
    'VariableNames', {'cloud', 'max_geometry_diff_deg', 'bands_differing', ...
    'max_abs_rel_diff', 'rms_rel_diff'});
disp(T)
fprintf('Largest relative difference over all shared clouds and bands: %.3g\n', ...
    max(max_abs_rel_diff));


%% Figure 1: spectra, one panel per shared cloud
fig1 = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 9 2.6 * n + 0.6]);
tl = tiledlayout(fig1, n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n
    ax = nexttile(tl);
    plot(ax, wl_nm, R1(k, :), '-', 'Color', C_VER1, 'LineWidth', 3, ...
        'DisplayName', 'ver1 (636 uvspec runs), 13 Nov');
    hold(ax, 'on');
    plot(ax, wl_nm, R2(k, :), '--', 'Color', C_VER2, 'LineWidth', 1.2, ...
        'DisplayName', 'ver2 (4 uvspec runs), 14 Nov');
    hold(ax, 'off');
    grid(ax, 'on');
    xlim(ax, [min(wl_nm) - 10, max(wl_nm) + 10]);
    ylabel(ax, 'Reflectance');
    g = ver2.geom(i2(k), :);
    title(ax, {sprintf('Cloud %d:  SZA = %.1f\\circ,  VZA = %.1f\\circ', shared(k), g(1), g(2)), ...
        sprintf('largest relative difference over 636 bands: %.1e  (%d of 636 bands differ)', ...
        max_abs_rel_diff(k), n_bands_differ(k))}, 'FontWeight', 'normal');
    if k == 1
        legend(ax, 'Location', 'northeast', 'Box', 'off');
    end
    if k == n
        xlabel(ax, 'Wavelength (nm)');
    end
end
title(tl, 'Alpine test runs: HySICS reflectance from the ver1 and ver2 drivers', 'FontSize', 13);
png1 = fullfile(OUT_DIR, 'alpine_ver1_vs_ver2_spectra.png');
exportgraphics(fig1, png1, 'Resolution', 400);
close(fig1);


%% Figure 2: relative difference vs wavelength, all shared clouds
% Exact zeros cannot be drawn on a log axis, so the absolute relative
% difference is plotted with zeros replaced by 1e-18 and the axis floor
% marked; anything on the floor is bit-for-bit identical.
FLOOR = 1e-18;
fig2 = new_figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 9 4]);
ax = axes(fig2);
cols = lines(n);
hold(ax, 'on');
for k = 1:n
    y = abs(rel_diff(k, :));
    y(y == 0) = FLOOR;
    semilogy(ax, wl_nm, y, '.', 'Color', cols(k, :), 'MarkerSize', 8, ...
        'DisplayName', sprintf('cloud %d', shared(k)));
end
yline(ax, eps, ':', 'machine epsilon (2.2e-16)', 'Color', [0.4 0.4 0.4], ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
hold(ax, 'off');
set(ax, 'YScale', 'log');
ylim(ax, [FLOOR / 2, max(1e-12, 10 * max(max_abs_rel_diff))]);
xlim(ax, [min(wl_nm) - 10, max(wl_nm) + 10]);
grid(ax, 'on');
xlabel(ax, 'Wavelength (nm)');
ylabel(ax, '|ver2 - ver1| / ver1');
title(ax, sprintf('Relative difference per band (points at %.0e are identical)', FLOOR), ...
    'FontWeight', 'normal');
legend(ax, 'Location', 'northeastoutside', 'Box', 'off');
png2 = fullfile(OUT_DIR, 'alpine_ver1_vs_ver2_rel_diff.png');
exportgraphics(fig2, png2, 'Resolution', 400);
close(fig2);

fprintf('\nFigures written:\n  %s\n  %s\n', png1, png2);

% From the desktop, open the two PNGs for viewing (no export of visible figures)
if usejava('desktop')
    for p = {png1, png2}
        figure('Name', p{1});
        image(imread(p{1}));
        axis image off
    end
end


%% Local function
function fig = new_figure(varargin)
% NEW_FIGURE  A figure using MATLAB's FACTORY graphics defaults, whatever
% startup.m sets: each root default is replaced, on this figure only, by its
% factory value. Copied from EPCAPE_KP_AJB/supersaturation/functions/new_figure.m.
fig = figure(varargin{:});
defs = fieldnames(get(groot, 'default'));
for k = 1:numel(defs)
    name = defs{k};
    if startsWith(name, 'defaultFigure', 'IgnoreCase', true)
        continue
    end
    set(fig, name, get(groot, regexprep(name, '^default', 'factory', 'ignorecase')));
end
fig.Color = 'w';
end
