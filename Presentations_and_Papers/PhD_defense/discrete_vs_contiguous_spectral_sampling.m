%% Discrete vs. contiguous spectral sampling of a reflected-solar cloud spectrum
%
% PhD defense graphic. Two stacked panels share the SAME reflected solar
% spectrum over a cloudy scene (350 - 2300 nm). Each panel overlays its
% spectral sampling scheme as translucent vertical channels:
%
%   TOP    - discrete spectral sampling: ~19 wide channels with gaps
%            (multispectral imager, e.g. MODIS -- not named on the figure)
%   BOTTOM - contiguous spectral sampling: 636 narrow adjacent channels
%            that fill the range (hyperspectral imager, e.g. HySICS /
%            CLARREO Pathfinder -- not named on the figure)
%
% The reflectance spectrum (636 HySICS bands, units 1/sr) is loaded from a
% simulated measurement over an inhomogeneous droplet profile.
%
% Andrew Buggee

clear variables

% --- Make sure the helper functions and band table are on the path ---
addpath('/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/startup/');
addpath('/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/');

% --- Locate the simulated HySICS reflectance file for this machine ---
which_computer = whatComputer();

if strcmp(which_computer, 'anbu8374') == true

    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];
    filename = ['simulated_measurement_HySICS_reflectance_inhomogeneous_', ...
        'droplet_profile_sim-ran-on-02-Jun-2025_ALL_BANDS.mat'];

else   % 'andrewbuggee' macbook (and fallback)

    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];
    filename = ['simulated_measurement_HySICS_reflectance_inhomogeneous_', ...
        'droplet_profile_sim-ran-on-02-Jun-2025_rev1.mat'];

end

ds = load([folderpath, filename]);

% Band-center wavelength (nm) for each of the 636 HySICS channels, and the
% modeled reflectance (1/sr).
wl_nm   = mean(ds.spec_response.wavelength, 2);   % [636 x 1]
refl    = ds.Refl_model;                          % [636 x 1]
n_bands = numel(wl_nm);                            % 636

% --- Plot limits (shared by both panels) ---
wl_lims_nm   = [350, 2300];    % HySICS spectral sensitivity range
refl_lims    = [0, 0.6];       % reflectance, 1/sr

% --- Channel colors ---
color_discrete   = mySavedColors(65, 'fixed');   % salmon (discrete channels)
color_contiguous = [0.10, 0.52, 0.62];           % teal   (contiguous channels)
color_spectrum   = [0, 0, 0];                    % black  (reflected spectrum)

% The user's startup.m applies a dark plotting theme. Force a clean,
% projector-friendly light theme for this figure so the spectrum (black)
% and the translucent channels read correctly on a white background.
light_axis = {'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.15, 0.15, 0.15], 'GridAlpha', 0.25, ...
    'MinorGridColor', [0.55, 0.55, 0.55], 'MinorGridAlpha', 0.15};

% =====================================================================
%                          BUILD THE FIGURE
% =====================================================================
figure;
set(gcf, 'Position', [0, 0, 1250, 850], 'Color', 'w')

tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');


% ---------------------------------------------------------------------
%  TOP PANEL  -  DISCRETE spectral sampling
% ---------------------------------------------------------------------
ax1 = nexttile;

% Discrete channels: the reflective-solar multispectral bands (1-19), drawn
% as translucent vertical patches at their true [lower, upper] widths.
mb = modisBands(1:19);          % columns: [center, lower, upper] in nm
lower_nm = mb(:, 2);
upper_nm = mb(:, 3);

% Vectorized patch (one graphics object, four vertices per channel).
Xd = [lower_nm, upper_nm, upper_nm, lower_nm]';                 % 4 x Nd
Yd = repmat([refl_lims(1); refl_lims(1); refl_lims(2); refl_lims(2)], 1, numel(lower_nm));
hd = patch(Xd, Yd, color_discrete, 'EdgeColor', 'none', 'FaceAlpha', 0.55);

hold on
hs1 = plot(wl_nm, refl, '-', 'Color', color_spectrum, 'LineWidth', 1.5);

ylim(refl_lims)
xlim(wl_lims_nm)
ylabel('Reflectance ($\mathrm{sr}^{-1}$)', 'Interpreter', 'latex', 'FontSize', 22, 'Color', 'k')
set(ax1, 'FontSize', 16, 'TickLabelInterpreter', 'latex', light_axis{:})
grid on; grid minor
box on

% Panel title + channel count
title('Discrete spectral sampling \,\textemdash\, 19 channels', ...
    'Interpreter', 'latex', 'FontSize', 26, 'Color', 'k')

legend([hs1, hd], {'Reflected solar spectrum', 'Spectral channels'}, ...
    'Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast', 'Box', 'on', ...
    'TextColor', 'k', 'Color', 'w', 'EdgeColor', [0.6, 0.6, 0.6])

% Hide x tick labels on the top panel (axis is shared with the bottom)
set(ax1, 'XTickLabel', [])


% ---------------------------------------------------------------------
%  BOTTOM PANEL  -  CONTIGUOUS spectral sampling
% ---------------------------------------------------------------------
ax2 = nexttile;

% Contiguous channels: all 636 narrow bands, drawn as translucent vertical
% patches. Each patch spans roughly half the local band spacing so that
% individual channels read as fine adjacent stripes (contiguous coverage)
% rather than collapsing into one solid block.
dwl = gradient(wl_nm);                 % local band spacing (nm), ~3.1 nm
halfw = 0.42 * dwl;                    % patch half-width (tiny gaps between bands)
Xc = [wl_nm - halfw, wl_nm + halfw, wl_nm + halfw, wl_nm - halfw]';   % 4 x Nc
Yc = repmat([refl_lims(1); refl_lims(1); refl_lims(2); refl_lims(2)], 1, n_bands);
hc = patch(Xc, Yc, color_contiguous, 'EdgeColor', 'none', 'FaceAlpha', 0.30);

hold on
hs2 = plot(wl_nm, refl, '-', 'Color', color_spectrum, 'LineWidth', 1.5);

ylim(refl_lims)
xlim(wl_lims_nm)
xlabel('Wavelength ($nm$)', 'Interpreter', 'latex', 'FontSize', 22, 'Color', 'k')
ylabel('Reflectance ($\mathrm{sr}^{-1}$)', 'Interpreter', 'latex', 'FontSize', 22, 'Color', 'k')
set(ax2, 'FontSize', 16, 'TickLabelInterpreter', 'latex', light_axis{:})
grid on; grid minor
box on

title(sprintf('Contiguous spectral sampling \\,\\textemdash\\, %d channels', n_bands), ...
    'Interpreter', 'latex', 'FontSize', 26, 'Color', 'k')

legend([hs2, hc], {'Reflected solar spectrum', 'Spectral channels'}, ...
    'Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast', 'Box', 'on', ...
    'TextColor', 'k', 'Color', 'w', 'EdgeColor', [0.6, 0.6, 0.6])

% Link the two x-axes so they pan/zoom together
linkaxes([ax1, ax2], 'x')


% =====================================================================
%                          SAVE THE FIGURE
% =====================================================================
save_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/PhD_defense/';

f = gcf;
saveas(f, [save_folder, 'discrete_vs_contiguous_spectral_sampling.fig']);
exportgraphics(f, [save_folder, 'discrete_vs_contiguous_spectral_sampling.png'], ...
    'Resolution', 400);
