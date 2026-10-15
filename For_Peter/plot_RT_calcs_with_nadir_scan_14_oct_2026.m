%% Plot simulated HySICS TOA radiance and reflectance with the CPF nadir scan
% Loads the TOA radiance and TOA reflectance calcs from csv files
% The CPF nadir scan is added to both figures as a dotted black line.
%
% Input files:
%   (1) RT_calcs_over_CPF_spectralRange_tauC_<tau>.csv - one file per cloud
%       optical thickness. Columns:
%           Wavelength (nm), Reflectance, Radiance (W/m^2/sr/nm)
%       The optical thickness is read from the number in the file name.
%   (2) nadir_scan_2026_162.csv - the CPF nadir scan, with the same
%       three columns.
%
% To run on another computer, only edit the "USER SETTINGS" section.
%
% By Andrew John Buggee

clear variables


%% ======================== USER SETTINGS ========================

% The only files needed are the 5 libRadtran csv files and the CPF nadir
% scan csv file. Set the two locations below to where they are saved on
% your computer. Nothing else needs to be changed.

% (1) Folder holding the 5 libRadtran csv files:
%       RT_calcs_over_CPF_spectralRange_tauC_20.csv
%       RT_calcs_over_CPF_spectralRange_tauC_40.csv
%       RT_calcs_over_CPF_spectralRange_tauC_60.csv
%       RT_calcs_over_CPF_spectralRange_tauC_80.csv
%       RT_calcs_over_CPF_spectralRange_tauC_100.csv
rt_csv_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/For_Peter/saved_calcs/';

% (2) Full path to the CPF nadir scan csv file
cpf_csv_file = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/For_Peter/CPF_scans/nadir_scan_2026_162.csv';


% Optional: save the figures as .fig and .jpg files. They are written to
% the same folder as the libRadtran csv files.
save_figures = true;
fig_folder = rt_csv_folder;


% ----- Scene used in the libRadtran calculations (for the subtitle) -----
% These values are not stored in the csv files. They are the settings of
% the 14 Oct 2026 run (saved in HySICS_spectra/simulated_HySICS_TOA_rad_refl_
% homogeneous_cloud_rEff_10um_5-tauC_sza_28_vza_0_sim-ran-on-14-Oct-2026_rev1.mat).
r_e_um = 10;                % microns - cloud effective radius
cloud_top_km = 1.5;         % km - cloud top altitude
cloud_base_km = 1;          % km - cloud base altitude
sza_deg = 28;               % degrees - solar zenith angle
vza_deg = 0;                % degrees - viewing zenith angle
cwv_mm = 14.3;              % mm - column integrated water vapor. The water vapor
                            % profile of the AFGL US standard atmosphere
                            % (afglus.dat) was used unchanged: 14.295 mm.

% ===============================================================



%% Load the simulated spectra (one csv file per optical thickness)

rt_files = dir(fullfile(rt_csv_folder, 'RT_calcs_over_CPF_spectralRange_tauC_*.csv'));

if isempty(rt_files)
    error([newline, 'No RT_calcs_over_CPF_spectralRange_tauC_*.csv files found in: ',...
        rt_csv_folder, newline, 'Update rt_csv_folder in the USER SETTINGS section.', newline])
end

% Read the optical thickness from each file name, e.g. "..._tauC_20.csv" -> 20.
% dir lists files alphabetically (100 before 20), so sort by optical thickness.
tau_c = zeros(length(rt_files), 1);
for ff = 1:length(rt_files)
    tau_str = regexp(rt_files(ff).name, 'tauC_([\d.]+)\.csv$', 'tokens', 'once');
    tau_c(ff) = str2double(tau_str{1});
end

[tau_c, idx_sort] = sort(tau_c);
rt_files = rt_files(idx_sort);

num_tauC = length(tau_c);


% Read each file. Columns are read by position (1 = wavelength,
% 2 = reflectance, 3 = radiance), so the header text can't break this.
wl_nm = cell(num_tauC, 1);              % nm
refl_TOA = cell(num_tauC, 1);           % reflectance
rad_TOA_W = cell(num_tauC, 1);          % W/(m^2 nm sr)

for tt = 1:num_tauC

    T = readtable(fullfile(rt_csv_folder, rt_files(tt).name), 'VariableNamingRule', 'preserve');

    wl_nm{tt} = T{:, 1};
    refl_TOA{tt} = T{:, 2};
    rad_TOA_W{tt} = T{:, 3};

end



%% Load the CPF nadir scan

if ~isfile(cpf_csv_file)
    error([newline, 'CPF file not found: ', cpf_csv_file, newline,...
        'Update cpf_csv_file in the USER SETTINGS section.', newline])
end

T_cpf = readtable(cpf_csv_file, 'VariableNamingRule', 'preserve');

cpf_wl_nm = T_cpf{:, 1};          % nm
cpf_refl = T_cpf{:, 2};           % reflectance
cpf_rad_W = T_cpf{:, 3};          % W/(m^2 nm sr)



%% Plot settings

% Line colors for the simulated spectra, one per optical thickness. These
% are the RGB values of mySavedColors(62:66, 'fixed') used in the original
% figures, written out so no personal function is needed.
line_colors = [0.4471    0.8863    0.4510;
               0.3882    0         0.0980;
               0.2157    0.2235    0.6627;
               0.9333    0.6745    0.1843;
               0.7804    0.6706    1.0000];

% If there are more optical thicknesses than saved colors, use MATLAB's
% default line colors instead
if num_tauC > size(line_colors, 1)
    line_colors = lines(num_tauC);
end

% Font sizes
ttl_fnt = 25;
sub_ttl_fnt = 20;
ax_fnt = 25;
lgnd_fnt = 20;
tick_fnt = 14;           % numbers on the axes (tick labels); MATLAB's default is 10

% Legend entries: one per optical thickness, then the CPF scan
legend_str = cell(num_tauC + 1, 1);
for tt = 1:num_tauC
    legend_str{tt} = ['$\tau_c = $', num2str(tau_c(tt))];
end
legend_str{end} = 'CPF nadir scan';

% Text describing the scene, used as the subtitle of both figures
scene_str = ['$r_e = $', num2str(r_e_um), ' $\mu m$, ',...
    'cloud top $= $', num2str(cloud_top_km), ' km, ',...
    'cloud base $= $', num2str(cloud_base_km), ' km, ',...
    '$sza = $', num2str(sza_deg), '$^{\circ}$, ',...
    '$vza = $', num2str(vza_deg), '$^{\circ}$'];

% Grey note on each figure stating the assumed column water vapor. It is
% placed in normalized axes units (0 to 1), right-aligned just below the
% legend, so it stays put whatever the data limits are.
cwv_str = ['Column integrated water vapor assumed: ', num2str(cwv_mm), ' mm'];
cwv_position = [0.664946255200931 0.959282498523616 0];            % normalized axes units [x, y] of the text's right edge
cwv_color = [0.35, 0.35, 0.35];            % grey

% x-axis limits span both the simulations and the CPF scan
wl_limits = [min([cell2mat(wl_nm); cpf_wl_nm]), max([cell2mat(wl_nm); cpf_wl_nm])];     % nm



%% Figure 1: TOA radiance

fig1 = figure;
use_factory_graphics_defaults(fig1)       % same look on every computer (see local function)

% Simulated spectra, one line per optical thickness
for tt = 1:num_tauC

    plot(wl_nm{tt}, rad_TOA_W{tt}, '-', 'LineWidth', 2, 'Color', line_colors(tt,:))
    hold on

end

% CPF nadir scan as a dotted black line
plot(cpf_wl_nm, cpf_rad_W, 'k:', 'LineWidth', 2)

% Font size of the numbers on the axes. Set before the labels and title,
% which have their own font sizes.
set(gca, 'FontSize', tick_fnt)

grid on; grid minor
xlim(wl_limits)
xlabel('Wavelength ($nm$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel('TOA Radiance ($W \, m^{-2} \, nm^{-1} \, sr^{-1}$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
title('TOA Radiance over HySICS spectral range - homogeneous liquid water cloud', 'Interpreter', 'latex',...
    'FontSize', ttl_fnt)
subtitle(scene_str, 'Interpreter', 'latex', 'FontSize', sub_ttl_fnt)
legend(legend_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', lgnd_fnt)
text(cwv_position(1), cwv_position(2), cwv_str, 'Units', 'normalized', 'Color', cwv_color,...
    'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', lgnd_fnt-5)
set(fig1, 'Position', [0 0 1200 700])



%% Figure 2: TOA reflectance

fig2 = figure;
use_factory_graphics_defaults(fig2)       % same look on every computer (see local function)

% Simulated spectra, one line per optical thickness
for tt = 1:num_tauC

    plot(wl_nm{tt}, refl_TOA{tt}, '-', 'LineWidth', 2, 'Color', line_colors(tt,:))
    hold on

end

% CPF nadir scan as a dotted black line
plot(cpf_wl_nm, cpf_refl, 'k:', 'LineWidth', 2)

% Font size of the numbers on the axes. Set before the labels and title,
% which have their own font sizes.
set(gca, 'FontSize', tick_fnt)

grid on; grid minor
xlim(wl_limits)
xlabel('Wavelength ($nm$)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
% R = pi*L / (mu0*F0) is unitless
ylabel('TOA Reflectance', 'Interpreter', 'latex', 'FontSize', ax_fnt)
title('TOA Reflectance over HySICS spectral range - homogeneous liquid water cloud', 'Interpreter', 'latex',...
    'FontSize', ttl_fnt)
subtitle(scene_str, 'Interpreter', 'latex', 'FontSize', sub_ttl_fnt)
legend(legend_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', lgnd_fnt)
text(cwv_position(1), cwv_position(2), cwv_str, 'Units', 'normalized', 'Color', cwv_color,...
    'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', lgnd_fnt-5)
set(fig2, 'Position', [0 0 1200 700])



%% Save the figures

if save_figures

    if ~exist(fig_folder, 'dir')
        mkdir(fig_folder)
    end

    fig1_name = 'TOA Radiance over HySICS spectral range with CPF nadir scan';
    fig2_name = 'TOA Reflectance over HySICS spectral range with CPF nadir scan';

    % save .fig files
    saveas(fig1, fullfile(fig_folder, [fig1_name, '.fig']));
    saveas(fig2, fullfile(fig_folder, [fig2_name, '.fig']));

    % save .jpg files with 500 DPI resolution
    exportgraphics(fig1, fullfile(fig_folder, [fig1_name, '.jpg']), 'Resolution', 500);
    exportgraphics(fig2, fullfile(fig_folder, [fig2_name, '.jpg']), 'Resolution', 500);

    disp([newline, 'Saved figures to: ', fig_folder, newline])

end


%% Local function

function use_factory_graphics_defaults(fig)
% Give this figure MATLAB's factory graphics settings (white background,
% black text and axes, default font sizes and line widths). This makes the
% figure look the same on every computer, even when a startup.m file has
% changed the root (groot) graphics defaults, e.g. to a dark theme. Without
% this, the dotted black CPF line is invisible on a black background.
% On a computer with no custom defaults this does nothing.

custom_defaults = fieldnames(get(groot, 'default'));

for nn = 1:length(custom_defaults)

    name = custom_defaults{nn};

    % Figure-level defaults only apply when a figure is created, so skip them
    if startsWith(name, 'defaultFigure', 'IgnoreCase', true)
        continue
    end

    % e.g. 'defaultAxesColor' -> 'factoryAxesColor'
    factory_name = regexprep(name, '^default', 'factory', 'ignorecase');
    set(fig, name, get(groot, factory_name));

end

% white figure background
set(fig, 'Color', 'w');

end
