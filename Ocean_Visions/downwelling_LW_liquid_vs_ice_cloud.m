%% Downwelling longwave irradiance at the Arctic Ocean surface: liquid vs. ice cloud
%
% Compares the surface downwelling longwave (LW) irradiance beneath a
% single-layer, vertically homogeneous liquid cloud and beneath an ice cloud
% holding the SAME column condensed water mass. The condensed water path is
% swept from 2 to 100 g/m^2. At each point on the sweep the liquid water path
% (LWP) and the ice water path (IWP) are numerically identical -- all the
% liquid is "converted" to ice -- so the only thing that changes between the
% two curves is the single-scattering and absorption physics of the condensate.
%
% Scene:
%   - Subarctic winter atmosphere (afglsw.dat), representative of 22 December
%     over the Arctic Ocean.
%   - Ocean surface, albedo 0.04 (thermal emissivity 0.96).
%   - No solar source. On 22 December the sun never rises poleward of the
%     Arctic circle, so a thermal-only calculation is the correct idealization
%     rather than an approximation.
%   - Liquid cloud: r_e = 10 micron.  Ice cloud: r_e = 27.5 micron, solid columns.
%   - Cloud between 0.5 and 1.0 km (a typical Arctic low stratus deck). Fixed
%     for both phases, so both clouds sit at the same temperature.
%
% Radiative transfer:
%   - two-stream solver (rte_solver twostr)
%   - REPTRAN coarse band parameterization, band-integrated irradiance summed
%     to broadband. The spectral range follows the ice parameterization toggle
%     (ice_parameterization, set below): 'yang' uses all 260 bands out to
%     99808.6 nm, 'yang2013' drops the topmost band because its tables stop at
%     99000 nm. See the ICE CLOUD OPTICS block for the full argument.
%
% *** A NOTE ON THE CLOUD FILES ***
% The clouds are specified by column water path, not optical depth, and the
% .DAT files are written directly by write_cloud_file_from_waterPath.m rather
% than through write_wc_file.m / write_ic_file.m in the LibRadTran-wrapper
% folder. Reasons:
%   (1) Water path -> water content is exact arithmetic (WC = WP / H). Routing
%       it through optical depth would make the condensed mass depend on a Mie
%       table and a reference wavelength that play no role in this experiment.
%   (2) write_ic_file.m hard-codes a libRadtran-2.0.4 data path for this
%       machine, which is stale after the 2.0.6 migration.
%   (3) The '2limit'/'mono' branch of write_ic_file.m omits the ice density in
%       the IWC calculation (line ~1256), so IWC comes out wrong by a factor of
%       rho_ice.
% Using one small, symmetric helper for both phases keeps the comparison clean.
%
% Requires: runUVSPEC_ver2.m and whatComputer.m to be on the MATLAB path.
%
%
% By Andrew John Buggee
%%

clear variables
close all


%% ------------------------------------------------------------------------
%  ----------------------------- SETTINGS ---------------------------------
%  ------------------------------------------------------------------------

computer_name = whatComputer;

if strcmp(computer_name, 'andrewbuggee')==true

    libRadtran_path      = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/';
    project_folder_path  = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Ocean_Visions/';

elseif strcmp(computer_name, 'anbu8374')==true

    libRadtran_path      = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/';
    project_folder_path  = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Ocean_Visions/';

elseif strcmp(computer_name, 'curc')==true

    libRadtran_path      = '/projects/anbu8374/software/libRadtran-2.0.5/';
    project_folder_path  = '/projects/anbu8374/Matlab-Research/Ocean_Visions/';

else

    error([newline, 'Unrecognized computer. Add its libRadtran path here.', newline])

end

libRadtran_data_path = [libRadtran_path, 'data/'];

% folders for the files this script generates
cloud_folder_path = [project_folder_path, 'cloud_files/'];
INP_folder_path   = [project_folder_path, 'INP_OUT/'];

if exist(cloud_folder_path, 'dir')==0, mkdir(cloud_folder_path), end
if exist(INP_folder_path, 'dir')==0,   mkdir(INP_folder_path),   end


% ----------------------------- THE SWEEP ---------------------------------

% column condensed water path (g/m^2). Log spacing, because the cloud LW
% emissivity saturates steeply: essentially all of the structure lives below
% ~30 g/m^2 and a linear grid wastes points on the flat part of the curve.
waterPath_gPerMeterSquared = logspace(log10(2), log10(150), 30)';        % g/m^2 - column vector

n_runs = length(waterPath_gPerMeterSquared);


% --------------------------- CLOUD PROPERTIES ----------------------------

r_e_liquid_um = 10;                     % micron - liquid droplet effective radius
r_e_ice_um    = 27.5;                   % micron - ice particle effective radius

% z_topBottom_km = [1.0, 0.5];            % km - [cloud top, cloud bottom] above sea level
z_topBottom_km = [0.5, 0.25];            % km - [cloud top, cloud bottom] above sea level


% ------------------------- RADIATIVE TRANSFER ----------------------------

inputs.libRadtran_data_path = libRadtran_data_path;

% subarctic winter -- the 22 December Arctic condition enters here
inputs.atmosphere_file = [libRadtran_data_path, 'atmmod/afglsw.dat'];

inputs.rte_solver = 'twostr';                       % two-stream
inputs.rte_solver_string = 'two-stream';

% inputs.rte_solver = 'disort';                       % Disort
% inputs.n_streams = 16;
% inputs.rte_solver_string = '16-stream';


inputs.mol_abs_param = 'reptran coarse';            % REPTRAN, broadband flux resolution

% *** The spectral selection is set further down, by the ice parameterization
% toggle -- the two ice tables have different upper wavelength limits. ***

% ----- define surface albedo -------
% inputs.albedo_string = 'Ocean albedo';
% inputs.albedo = 0.04;                               % ocean; thermal emissivity = 0.96

inputs.albedo_string = 'Wintertime Tundra albedo';
inputs.albedo = 0.003;                               % From Hori et al. 2006; thermal emissivity = 0.997
% -------------------------------------



inputs.zout_km = 0;                                 % km above the surface

% liquid cloud optics: Hu and Stamnes (1993), valid 0.29 - 150 micron
inputs.wc_properties = 'hu';

% ------------------------- ICE CLOUD OPTICS ------------------------------
%
% Toggle between the two habit-aware parameterizations that cover the thermal.
% The spectral selection is tied to this choice: the two tables have different
% upper wavelength limits, and uvspec validates an ice table against the FULL
% band set of the chosen mol_abs_param -- not against the subset asked for by
% 'wavelength'.
%
%   'yang'      Key et al. (2002)/Yang. Covers 0.2 - 100 micron, so the whole
%               260-band REPTRAN coarse thermal set (out to 99808.6 nm) is
%               usable and the range is set the ordinary way, with 'wavelength'.
%               Habits valid in the thermal at r_e = 27.5 micron:
%                    solid-column   (r_eff 5.96 - 84.22 micron)
%                    hollow-column  (r_eff 4.97 - 70.24 micron)
%                    rough-aggregate(r_eff 3.55 - 108.10 micron)
%                    rosette-6      (r_eff 2.85 - 46.01 micron)
%                    plate          (r_eff 4.87 - 48.18 micron)
%                    droxtal        (r_eff 9.48 - 293.32 micron)
%               NOT usable here: rosette-4 (no data above 3.4 micron),
%               dendrite (r_eff capped at 1.88 micron), and spheroid (rejected
%               by uvspec for this parameterization).
%
%   'yang2013'  Yang et al. (2013). Tables stop at 99000 nm, so the topmost
%               REPTRAN band (index 260 of 260, centred at 93478 nm and running
%               out to 99808.6 nm) has to be dropped with wavelength_index.
%               That band carries 0.85 W/m^2 of the 176.96 W/m^2 clear-sky
%               downwelling flux -- 0.48% of the broadband total -- and the
%               identical 259 bands are used for the clear sky and the liquid
%               cloud, so it cancels exactly in every difference plotted.
%               All 9 habits are valid over r_eff 5 - 90 micron, and each also
%               needs a roughness (smooth | moderate | severe):
%                    solid_column, hollow_column, column_8elements,
%                    plate, plate_5elements, plate_10elements,
%                    solid_bullet_rosette, hollow_bullet_rosette, droxtal
%
% 'fu' (Fu 1996; Fu et al. 1998) is a third option, also randomly oriented
% hexagonal columns; it agrees with yang/solid-column to better than 0.1 W/m^2.

ice_parameterization = 'yang2013';                  % 'yang' or 'yang2013'


if strcmp(ice_parameterization, 'yang')==true

    inputs.ic_properties = 'yang';
    inputs.ic_habit      = 'solid-column';

    % full REPTRAN thermal range. uvspec internally uses 2501.6 - 99808.6 nm;
    % these bounds simply ask for all of it.
    inputs.wavelength_bounds_nm = [2500, 100000];   % nm
    inputs.wavelength_index     = [];               % empty -> use wavelength

    ice_legend_str = inputs.ic_habit;


elseif strcmp(ice_parameterization, 'yang2013')==true

    inputs.ic_properties     = 'yang2013 interpolate';
    inputs.ic_habit_yang2013 = 'solid_column';
    inputs.ic_roughness      = 'severe';            % smooth | moderate | severe

    % bands 1-259 of 260, i.e. everything below 99000 nm. wavelength_index
    % replaces wavelength -- uvspec aborts if both are present.
    inputs.wavelength_index     = [1, 259];
    inputs.wavelength_bounds_nm = [2500, 99000];    % nm - documentation only, ignored

    ice_legend_str = [strrep(inputs.ic_habit_yang2013, '_', '\ '), ', ',...
        inputs.ic_roughness];


else

    error([newline, 'ice_parameterization must be either ''yang'' or ''yang2013''.', newline])

end


inputs.err_msg = 'quiet';


%% ------------------------------------------------------------------------
%  --------------------- CLEAR SKY REFERENCE RUN --------------------------
%  ------------------------------------------------------------------------

inputs_clear = inputs;
inputs_clear.cloud_phase = 'none';

inputFileName_clear  = 'thermal_clearSky.INP';
outputFileName_clear = 'thermal_clearSky';

write_INP_file_thermal_cloud(INP_folder_path, inputFileName_clear, inputs_clear);

runUVSPEC_ver2(INP_folder_path, inputFileName_clear, outputFileName_clear, computer_name);

clearSky = read_uvspec_thermal_perBand(INP_folder_path, [outputFileName_clear, '.OUT']);

LW_down_clearSky_WPerMeterSquared = clearSky.LW_down_WPerMeterSq;       % W/m^2

disp([newline, 'Clear sky downwelling LW at the surface: ',...
    num2str(LW_down_clearSky_WPerMeterSquared, '%.2f'), ' W/m^2', newline])


%% ------------------------------------------------------------------------
%  ------------------------ SWEEP OVER WATER PATH -------------------------
%  ------------------------------------------------------------------------

LW_down_liquid_WPerMeterSquared = zeros(n_runs, 1);         % W/m^2
LW_down_ice_WPerMeterSquared    = zeros(n_runs, 1);         % W/m^2

lwc_gPerMeterCubed = zeros(n_runs, 1);                      % g/m^3
iwc_gPerMeterCubed = zeros(n_runs, 1);                      % g/m^3

tic
for nn = 1:n_runs


    % ------------------------- LIQUID CLOUD -------------------------------

    cloudFileName_liquid = ['wc_LWP', num2str(round(waterPath_gPerMeterSquared(nn), 3)),...
        '_re', num2str(r_e_liquid_um), '_nn', num2str(nn), '.DAT'];

    [~, lwc_gPerMeterCubed(nn)] = write_cloud_file_from_waterPath(waterPath_gPerMeterSquared(nn),...
        r_e_liquid_um, z_topBottom_km, 'liquid', cloud_folder_path, cloudFileName_liquid);

    inputs_liquid = inputs;
    inputs_liquid.cloud_phase     = 'liquid';
    inputs_liquid.cloud_file_path = [cloud_folder_path, cloudFileName_liquid];

    inputFileName_liquid  = ['thermal_liquid_nn', num2str(nn), '.INP'];
    outputFileName_liquid = ['thermal_liquid_nn', num2str(nn)];

    write_INP_file_thermal_cloud(INP_folder_path, inputFileName_liquid, inputs_liquid);

    runUVSPEC_ver2(INP_folder_path, inputFileName_liquid, outputFileName_liquid, computer_name);

    liquid_out = read_uvspec_thermal_perBand(INP_folder_path, [outputFileName_liquid, '.OUT']);

    LW_down_liquid_WPerMeterSquared(nn) = liquid_out.LW_down_WPerMeterSq;           % W/m^2


    % --------------------------- ICE CLOUD --------------------------------

    cloudFileName_ice = ['ic_IWP', num2str(round(waterPath_gPerMeterSquared(nn), 3)),...
        '_re', num2str(r_e_ice_um), '_nn', num2str(nn), '.DAT'];

    [~, iwc_gPerMeterCubed(nn)] = write_cloud_file_from_waterPath(waterPath_gPerMeterSquared(nn),...
        r_e_ice_um, z_topBottom_km, 'ice', cloud_folder_path, cloudFileName_ice);

    inputs_ice = inputs;
    inputs_ice.cloud_phase     = 'ice';
    inputs_ice.cloud_file_path = [cloud_folder_path, cloudFileName_ice];

    inputFileName_ice  = ['thermal_ice_nn', num2str(nn), '.INP'];
    outputFileName_ice = ['thermal_ice_nn', num2str(nn)];

    write_INP_file_thermal_cloud(INP_folder_path, inputFileName_ice, inputs_ice);

    runUVSPEC_ver2(INP_folder_path, inputFileName_ice, outputFileName_ice, computer_name);

    ice_out = read_uvspec_thermal_perBand(INP_folder_path, [outputFileName_ice, '.OUT']);

    LW_down_ice_WPerMeterSquared(nn) = ice_out.LW_down_WPerMeterSq;                 % W/m^2


    disp(['Run ', num2str(nn), '/', num2str(n_runs), ' -- water path = ',...
        num2str(waterPath_gPerMeterSquared(nn), '%.2f'), ' g/m^2 | liquid = ',...
        num2str(LW_down_liquid_WPerMeterSquared(nn), '%.2f'), ' W/m^2 | ice = ',...
        num2str(LW_down_ice_WPerMeterSquared(nn), '%.2f'), ' W/m^2'])


end
run_time_sec = toc;

disp([newline, 'Finished ', num2str(2*n_runs + 1), ' uvspec runs in ',...
    num2str(run_time_sec, '%.1f'), ' seconds.', newline])


%% ------------------------------------------------------------------------
%  ------------------------------- PLOT -----------------------------------
%  ------------------------------------------------------------------------

figure('Position', [0 0 1200 500])

% --- panel 1: the two curves ---
subplot(1,2,1)

plot(waterPath_gPerMeterSquared, LW_down_liquid_WPerMeterSquared, '.-',...
    'LineWidth', 2, 'MarkerSize', 18, 'Color', mySavedColors(1, 'fixed'))
hold on
plot(waterPath_gPerMeterSquared, LW_down_ice_WPerMeterSquared, '.-',...
    'LineWidth', 2, 'MarkerSize', 18, 'Color', mySavedColors(2, 'fixed'))

yline(LW_down_clearSky_WPerMeterSquared, ':', 'Clear sky', 'LineWidth', 2,...
    'FontSize', 14, 'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left',...
    'Color', 'white')

grid on; grid minor
set(gca, 'XScale', 'log')

xlabel('Condensed water path ($g/m^{2}$)', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$F_{LW}^{\downarrow}$ at the surface ($W/m^{2}$)', 'Interpreter', 'latex', 'FontSize', 20)

legend({['Liquid, $r_e = ', num2str(r_e_liquid_um), ' \mu m$'],...
    ['Ice (', ice_legend_str, '), $r_e = ', num2str(r_e_ice_um), ' \mu m$']},...
    'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 16)

title('Surface downwelling longwave', 'Interpreter', 'latex', 'FontSize', 18)


% --- panel 2: the liquid-minus-ice difference ---
subplot(1,2,2)

plot(waterPath_gPerMeterSquared,...
    LW_down_liquid_WPerMeterSquared - LW_down_ice_WPerMeterSquared, '.-',...
    'LineWidth', 2, 'MarkerSize', 18, 'Color', mySavedColors(3, 'fixed'))

yline(0, ':', 'LineWidth', 2)

grid on; grid minor
set(gca, 'XScale', 'log')

xlabel('Condensed water path ($g/m^{2}$)', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$\Delta F_{LW}^{\downarrow}$ (liquid $-$ ice) ($W/m^{2}$)',...
    'Interpreter', 'latex', 'FontSize', 20)

title('Effect of Cloud Phase for fixed column water mass', 'Interpreter', 'latex', 'FontSize', 18)

sgtitle(['Subarctic winter, ', inputs.albedo_string, ' = ', num2str(inputs.albedo),...
    ', cloud ', num2str(z_topBottom_km(2)), '--', num2str(z_topBottom_km(1)),...
    ' km, ', inputs.rte_solver_string, ', ice optics: ', ice_parameterization],...
    'Interpreter', 'latex', 'FontSize', 18)


%% ------------------------------------------------------------------------
%  ------------------------------- SAVE -----------------------------------
%  ------------------------------------------------------------------------

results.waterPath_gPerMeterSquared      = waterPath_gPerMeterSquared;
results.LW_down_liquid_WPerMeterSquared = LW_down_liquid_WPerMeterSquared;
results.LW_down_ice_WPerMeterSquared    = LW_down_ice_WPerMeterSquared;
results.LW_down_clearSky_WPerMeterSquared = LW_down_clearSky_WPerMeterSquared;
results.lwc_gPerMeterCubed              = lwc_gPerMeterCubed;
results.iwc_gPerMeterCubed              = iwc_gPerMeterCubed;
results.r_e_liquid_um                   = r_e_liquid_um;
results.r_e_ice_um                      = r_e_ice_um;
results.z_topBottom_km                  = z_topBottom_km;
results.ice_parameterization            = ice_parameterization;
results.inputs                          = inputs;

save([project_folder_path, 'downwelling_LW_liquid_vs_ice_', char(datetime("today")), '.mat'],...
    'results')
