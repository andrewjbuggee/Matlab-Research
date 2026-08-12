%% Effect of ice particle habit on surface downwelling longwave irradiance
%
% Companion to downwelling_LW_liquid_vs_ice_cloud.m. Same scene, same sweep in
% column condensed water path, but the ice cloud is now run for EVERY habit
% available in the Yang et al. (2013) parameterization, holding effective
% radius, water path, cloud geometry and roughness fixed. A liquid cloud and a
% clear sky are carried along as references.
%
% The 9 yang2013 habits fall into four families:
%       columns   - solid_column, hollow_column, column_8elements
%       plates    - plate, plate_5elements, plate_10elements
%       rosettes  - solid_bullet_rosette, hollow_bullet_rosette
%       droxtal   - droxtal
% All 9 share the same effective radius grid (5 - 90 micron) and the same
% spectral coverage (0.2 - 99 micron), so the comparison is clean: nothing
% varies across the curves except particle shape.
%
% Two figures are produced. Both share the same left panel (absolute
% downwelling longwave for the clear sky, the liquid cloud and all 9 habits);
% they differ in the right panel:
%   figure 1 - habit minus the reference habit, isolating particle shape
%   figure 2 - liquid minus each habit, i.e. the phase effect resolved by habit
% The reference-habit curve in figure 2 reproduces the difference panel of
% downwelling_LW_liquid_vs_ice_cloud.m exactly, as long as that script is set
% to the same ice parameterization, habit and roughness (its
% ice_parameterization toggle set to 'yang2013').
%
%
% *** WHY wavelength_index INSTEAD OF wavelength ***
% uvspec validates an ice optical property table against the FULL band set of
% the chosen mol_abs_param, not against the subset asked for by 'wavelength'.
% With 'reptran coarse' that full set runs to 99808.6 nm, while the yang2013
% tables stop at 99000 nm, so yang2013 is rejected regardless of the bounds
% given to 'wavelength'. Selecting bands by index instead drops the topmost
% band (index 260 of 260, centred at 93478 nm) and the check passes.
%
% That band carries 0.85 W/m^2 of the 176.96 W/m^2 clear-sky downwelling flux,
% i.e. 0.48% of the broadband total. Every run in this script -- clear sky,
% liquid, and all 9 ice habits -- uses the identical 259-band set, so the
% omission is a common offset that cancels exactly in every difference plotted
% here. It does mean the absolute fluxes are ~0.5% lower than those from
% downwelling_LW_liquid_vs_ice_cloud.m, which uses all 260 bands.
%
%
% *** ROUGHNESS ***
% Held fixed across the habit comparison. Surface roughness reshapes the phase
% function but barely touches the absorption and extinction cross sections, so
% its effect on a broadband longwave FLUX is small: at 50 g/m^2 the spread
% between smooth, moderate and severe is under 0.2 W/m^2 for every habit (worst
% case 0.19 W/m^2, the plates; best case 0.004 W/m^2, solid columns), against a
% 2.2 W/m^2 spread across habits at the same water path. Set ic_roughness below
% to check.
%
%
% Requires: runUVSPEC_ver2.m, whatComputer.m and mySavedColors.m on the path.
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
waterPath_gPerMeterSquared = logspace(log10(2), log10(150), 30)';       % g/m^2 - column vector

n_runs = length(waterPath_gPerMeterSquared);


% --------------------------- CLOUD PROPERTIES ----------------------------

r_e_liquid_um = 10;                     % micron - liquid droplet effective radius
r_e_ice_um    = 27.5;                   % micron - ice particle effective radius

z_topBottom_km = [0.5, 0.25];           % km - [cloud top, cloud bottom] above sea level


% ------------------------- THE ICE HABITS --------------------------------

% All 9 habits in the Yang et al. (2013) parameterization, grouped by family
% so that related shapes sit next to each other in the legend.
ice_habits = {'solid_column', 'hollow_column', 'column_8elements',...
    'plate', 'plate_5elements', 'plate_10elements',...
    'solid_bullet_rosette', 'hollow_bullet_rosette', 'droxtal'};

ice_habits_display = {'solid column', 'hollow column', '8-element column agg.',...
    'plate', '5-element plate agg.', '10-element plate agg.',...
    'solid bullet rosette', 'hollow bullet rosette', 'droxtal'};

n_habits = length(ice_habits);

% habit that the difference panel is referenced against
reference_habit = 'solid_column';

idx_reference_habit = find(strcmp(ice_habits, reference_habit));
if isempty(idx_reference_habit)==true
    error([newline, 'The reference habit is not in the list of habits run.', newline])
end


% ------------------------- RADIATIVE TRANSFER ----------------------------

inputs.libRadtran_data_path = libRadtran_data_path;

% subarctic winter -- the 22 December Arctic condition enters here
inputs.atmosphere_file = [libRadtran_data_path, 'atmmod/afglsw.dat'];

inputs.rte_solver = 'twostr';                       % two-stream
inputs.rte_solver_string = 'two-stream';

% inputs.rte_solver = 'disort';                     % Disort
% inputs.n_streams = 16;
% inputs.rte_solver_string = '16-stream';

inputs.mol_abs_param = 'reptran coarse';            % REPTRAN, broadband flux resolution

% Band indices rather than wavelength bounds -- see the header. Bands 1-259 of
% the 260-band reptran coarse thermal set, i.e. everything below 99000 nm.
inputs.wavelength_index = [1, 259];

% kept only so the field exists; ignored while wavelength_index is set
inputs.wavelength_bounds_nm = [2500, 100000];       % nm


% ----- define surface albedo -------
% inputs.albedo_string = 'Ocean albedo';
% inputs.albedo = 0.04;                             % ocean; thermal emissivity = 0.96

inputs.albedo_string = 'Wintertime Tundra albedo';
inputs.albedo = 0.003;                              % From Hori et al. 2006; thermal emissivity = 0.997
% -------------------------------------

inputs.zout_km = 0;                                 % km above the surface

% liquid cloud optics: Hu and Stamnes (1993), valid 0.29 - 150 micron
inputs.wc_properties = 'hu';

% ice cloud optics: Yang et al. (2013), 0.2 - 99 micron, r_eff 5 - 90 micron
inputs.ic_properties = 'yang2013 interpolate';
inputs.ic_roughness  = 'severe';                    % smooth | moderate | severe

inputs.err_msg = 'quiet';


%% ------------------------------------------------------------------------
%  --------------------- CLEAR SKY REFERENCE RUN --------------------------
%  ------------------------------------------------------------------------

inputs_clear = inputs;
inputs_clear.cloud_phase = 'none';

inputFileName_clear  = 'habit_thermal_clearSky.INP';
outputFileName_clear = 'habit_thermal_clearSky';

write_INP_file_thermal_cloud(INP_folder_path, inputFileName_clear, inputs_clear);

runUVSPEC_ver2(INP_folder_path, inputFileName_clear, outputFileName_clear, computer_name);

clearSky = read_uvspec_thermal_perBand(INP_folder_path, [outputFileName_clear, '.OUT']);

LW_down_clearSky_WPerMeterSquared = clearSky.LW_down_WPerMeterSq;       % W/m^2

disp([newline, 'Clear sky downwelling LW at the surface: ',...
    num2str(LW_down_clearSky_WPerMeterSquared, '%.2f'), ' W/m^2',...
    '  (', num2str(length(clearSky.wavelength_nm)), ' bands)', newline])


%% ------------------------------------------------------------------------
%  ------------------------ SWEEP OVER WATER PATH -------------------------
%  ------------------------------------------------------------------------

LW_down_liquid_WPerMeterSquared = zeros(n_runs, 1);                 % W/m^2
LW_down_ice_WPerMeterSquared    = zeros(n_runs, n_habits);          % W/m^2 - one column per habit

lwc_gPerMeterCubed = zeros(n_runs, 1);                              % g/m^3
iwc_gPerMeterCubed = zeros(n_runs, 1);                              % g/m^3

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

    inputFileName_liquid  = ['habit_thermal_liquid_nn', num2str(nn), '.INP'];
    outputFileName_liquid = ['habit_thermal_liquid_nn', num2str(nn)];

    write_INP_file_thermal_cloud(INP_folder_path, inputFileName_liquid, inputs_liquid);

    runUVSPEC_ver2(INP_folder_path, inputFileName_liquid, outputFileName_liquid, computer_name);

    liquid_out = read_uvspec_thermal_perBand(INP_folder_path, [outputFileName_liquid, '.OUT']);

    LW_down_liquid_WPerMeterSquared(nn) = liquid_out.LW_down_WPerMeterSq;           % W/m^2


    % ------------------ ICE CLOUD, ONE RUN PER HABIT ----------------------

    % the ice cloud file depends only on water path and r_e, so it is written
    % once and reused by every habit
    cloudFileName_ice = ['ic_IWP', num2str(round(waterPath_gPerMeterSquared(nn), 3)),...
        '_re', num2str(r_e_ice_um), '_nn', num2str(nn), '.DAT'];

    [~, iwc_gPerMeterCubed(nn)] = write_cloud_file_from_waterPath(waterPath_gPerMeterSquared(nn),...
        r_e_ice_um, z_topBottom_km, 'ice', cloud_folder_path, cloudFileName_ice);

    for hh = 1:n_habits

        inputs_ice = inputs;
        inputs_ice.cloud_phase       = 'ice';
        inputs_ice.cloud_file_path   = [cloud_folder_path, cloudFileName_ice];
        inputs_ice.ic_habit_yang2013 = ice_habits{hh};

        inputFileName_ice  = ['habit_thermal_ice_', ice_habits{hh}, '_nn', num2str(nn), '.INP'];
        outputFileName_ice = ['habit_thermal_ice_', ice_habits{hh}, '_nn', num2str(nn)];

        write_INP_file_thermal_cloud(INP_folder_path, inputFileName_ice, inputs_ice);

        runUVSPEC_ver2(INP_folder_path, inputFileName_ice, outputFileName_ice, computer_name);

        ice_out = read_uvspec_thermal_perBand(INP_folder_path, [outputFileName_ice, '.OUT']);

        LW_down_ice_WPerMeterSquared(nn, hh) = ice_out.LW_down_WPerMeterSq;         % W/m^2

    end


    disp(['Run ', num2str(nn), '/', num2str(n_runs), ' -- water path = ',...
        num2str(waterPath_gPerMeterSquared(nn), '%.2f'), ' g/m^2 | liquid = ',...
        num2str(LW_down_liquid_WPerMeterSquared(nn), '%.2f'), ' W/m^2 | ice spread = ',...
        num2str(min(LW_down_ice_WPerMeterSquared(nn,:)), '%.2f'), ' - ',...
        num2str(max(LW_down_ice_WPerMeterSquared(nn,:)), '%.2f'), ' W/m^2'])


end
run_time_sec = toc;

disp([newline, 'Finished ', num2str(n_runs*(1 + n_habits) + 1), ' uvspec runs in ',...
    num2str(run_time_sec, '%.1f'), ' seconds.', newline])


% largest habit-to-habit spread anywhere on the sweep
[max_spread_WPerMeterSquared, idx_max_spread] = max(max(LW_down_ice_WPerMeterSquared, [], 2) -...
    min(LW_down_ice_WPerMeterSquared, [], 2));

disp(['Largest habit-to-habit spread: ', num2str(max_spread_WPerMeterSquared, '%.2f'),...
    ' W/m^2 at a water path of ',...
    num2str(waterPath_gPerMeterSquared(idx_max_spread), '%.1f'), ' g/m^2'])


% peak phase effect for the reference habit. This is the same quantity plotted
% by downwelling_LW_liquid_vs_ice_cloud.m, and is printed here so the two
% scripts can be checked against each other directly.
[max_phase_WPerMeterSquared, idx_max_phase] = max(LW_down_liquid_WPerMeterSquared -...
    LW_down_ice_WPerMeterSquared(:,idx_reference_habit));

disp(['Peak liquid - ice difference for ', strrep(reference_habit, '_', ' '), ': ',...
    num2str(max_phase_WPerMeterSquared, '%.2f'), ' W/m^2 at a water path of ',...
    num2str(waterPath_gPerMeterSquared(idx_max_phase), '%.1f'), ' g/m^2', newline])


%% ------------------------------------------------------------------------
%  ------------------------------- PLOT -----------------------------------
%  ------------------------------------------------------------------------

% one colour per habit, plus a visually distinct colour for the liquid cloud so
% that it reads as a different family rather than a tenth habit.
% turbo is trimmed at both ends -- its extremes are a very dark navy and a very
% dark maroon, which disappear against a dark figure background.
turbo_full   = turbo(n_habits + 4);
habit_colors = turbo_full(3:(n_habits+2), :);

liquid_color = mySavedColors(1, 'fixed');

% Two figures sharing an identical left panel. They differ only in what the
% right panel differences against:
%   figure 1 - habit minus the reference habit, which isolates particle shape
%   figure 2 - liquid minus each habit, i.e. the phase effect resolved by habit
% The curve for the reference habit in figure 2 reproduces exactly the
% difference panel of downwelling_LW_liquid_vs_ice_cloud.m, provided that
% script is set to the same ice parameterization, habit and roughness.

for ff = 1:2

    figure('Position', [0 0 1300 550])


    % --- panel 1: absolute downwelling longwave (identical in both figures) ---
    subplot(1,2,1)

    hold on

    for hh = 1:n_habits
        plot(waterPath_gPerMeterSquared, LW_down_ice_WPerMeterSquared(:,hh), '.-',...
            'LineWidth', 1.5, 'MarkerSize', 12, 'Color', habit_colors(hh,:))
    end

    plot(waterPath_gPerMeterSquared, LW_down_liquid_WPerMeterSquared, '.-',...
        'LineWidth', 3, 'MarkerSize', 18, 'Color', liquid_color)

    yline(LW_down_clearSky_WPerMeterSquared, ':', 'Clear sky', 'LineWidth', 2,...
        'FontSize', 14, 'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left',...
        'Color', 'white')

    grid on; grid minor
    set(gca, 'XScale', 'log')

    xlabel('Condensed water path ($g/m^{2}$)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('$F_{LW}^{\downarrow}$ at the surface ($W/m^{2}$)', 'Interpreter', 'latex', 'FontSize', 20)

    legend([ice_habits_display, {['Liquid, $r_e = ', num2str(r_e_liquid_um), ' \mu m$']}],...
        'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 12)

    title(['Surface downwelling longwave, ice $r_e = ', num2str(r_e_ice_um), ' \mu m$'],...
        'Interpreter', 'latex', 'FontSize', 18)


    % --- panel 2: the difference ---
    subplot(1,2,2)

    hold on

    if ff==1

        % ---- habit minus the reference habit: the shape effect alone ----
        for hh = 1:n_habits

            if hh==idx_reference_habit
                continue        % identically zero by construction
            end

            plot(waterPath_gPerMeterSquared,...
                LW_down_ice_WPerMeterSquared(:,hh) -...
                LW_down_ice_WPerMeterSquared(:,idx_reference_habit),...
                '.-', 'LineWidth', 1.5, 'MarkerSize', 12, 'Color', habit_colors(hh,:))

        end

        panel2_ylabel = ['$\Delta F_{LW}^{\downarrow}$ (habit $-$ ',...
            strrep(reference_habit, '_', '\ '), ') ($W/m^{2}$)'];

        panel2_title  = 'Effect of ice habit at fixed column water mass';

        panel2_legend = ice_habits_display(setdiff(1:n_habits, idx_reference_habit));

        % the curves converge on zero at large water path, leaving the lower
        % right corner empty
        panel2_legend_location = 'southeast';

    else

        % ---- liquid minus each habit: the phase effect, resolved by habit ----
        for hh = 1:n_habits

            plot(waterPath_gPerMeterSquared,...
                LW_down_liquid_WPerMeterSquared - LW_down_ice_WPerMeterSquared(:,hh),...
                '.-', 'LineWidth', 1.5, 'MarkerSize', 12, 'Color', habit_colors(hh,:))

        end

        panel2_ylabel = '$\Delta F_{LW}^{\downarrow}$ (liquid $-$ ice) ($W/m^{2}$)';

        panel2_title  = 'Effect of cloud phase, resolved by ice habit';

        panel2_legend = ice_habits_display;

        % these curves rise out of the lower left, peak in the middle and fall
        % back to zero on the right, leaving the lower left corner empty
        panel2_legend_location = 'southwest';

    end

    yline(0, ':', 'LineWidth', 2, 'Color', 'white')

    grid on; grid minor
    set(gca, 'XScale', 'log')

    xlabel('Condensed water path ($g/m^{2}$)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel(panel2_ylabel, 'Interpreter', 'latex', 'FontSize', 20)

    legend(panel2_legend, 'Interpreter', 'latex', 'Location', panel2_legend_location,...
        'FontSize', 12)

    title(panel2_title, 'Interpreter', 'latex', 'FontSize', 18)


    sgtitle(['Subarctic winter, ', inputs.albedo_string, ' = ', num2str(inputs.albedo),...
        ', cloud ', num2str(z_topBottom_km(2)), '--', num2str(z_topBottom_km(1)),...
        ' km, ', inputs.rte_solver_string, ', Yang et al. (2013) ', inputs.ic_roughness,...
        ' roughness'], 'Interpreter', 'latex', 'FontSize', 18)

end


%% ------------------------------------------------------------------------
%  ------------------------------- SAVE -----------------------------------
%  ------------------------------------------------------------------------

results.waterPath_gPerMeterSquared        = waterPath_gPerMeterSquared;
results.ice_habits                        = ice_habits;
results.ice_habits_display                = ice_habits_display;
results.reference_habit                   = reference_habit;
results.LW_down_liquid_WPerMeterSquared   = LW_down_liquid_WPerMeterSquared;
results.LW_down_ice_WPerMeterSquared      = LW_down_ice_WPerMeterSquared;
results.LW_down_clearSky_WPerMeterSquared = LW_down_clearSky_WPerMeterSquared;
results.lwc_gPerMeterCubed                = lwc_gPerMeterCubed;
results.iwc_gPerMeterCubed                = iwc_gPerMeterCubed;
results.r_e_liquid_um                     = r_e_liquid_um;
results.r_e_ice_um                        = r_e_ice_um;
results.z_topBottom_km                    = z_topBottom_km;
results.inputs                            = inputs;

save([project_folder_path, 'downwelling_LW_ice_habit_', char(datetime("today")), '.mat'],...
    'results')
