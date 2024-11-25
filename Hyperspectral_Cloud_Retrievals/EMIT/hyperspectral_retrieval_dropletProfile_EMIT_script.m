%% Optimal estimation of a vertical droplet profile using EMIT data


% By Andrew John Buggee

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;

%% Define the EMIT data file to use by defining the data folder

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitDataFolder = '17_Jan_2024_coast/';

% 27 january has overlap with MODIS observations
%emitDataFolder = '27_Jan_2024/';

% -------------------------------------


%% Load EMIT data and define folders 

[emitDataPath, folder2save] = define_EMIT_dataPath_and_saveFolders();

[emit,L1B_fileName] = retrieveEMIT_data([emitDataPath, emitDataFolder]);


%% Define the pixels to use for the retrieval


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% pixels2use.row = 1112;
% pixels2use.col = 974;

% 17_Jan_2024_coast - small optical depth
% pixels2use.row = [912, 913];
% pixels2use.col = [929, 929];

% 17_Jan_2024_coast - clear sky over ocean
% pixels2use.row = 1021;
% pixels2use.col = 536;

% 17_Jan_2024_coast - my TBLUT algorithm found an optical depth of 6.57 and
% an effective radius of 10.79
% pixels2use.row = 932;
% pixels2use.col = 960;

% 17_Jan_2024_coast - optical depth of 3.2 and 3.8
% pixels2use.row = [932, 932];
% pixels2use.col = [970, 969];

% 17_Jan_2024_coast - optical depth of 3.2
% pixels2use.row = 932;
% pixels2use.col = 970;

% 17_Jan_2024_coast - optical depth of 3.8
% pixels2use.row = 932;
% pixels2use.col = 969;

% 17_Jan_2024_coast - optical depth of 8.7
% pixels2use.row = 969;
% pixels2use.col = 991;

% 17_Jan_2024_coast - optical depths of 6.6, 11.6
% pixels2use(1).row = 932; pixels2use(1).col = 960;
% pixels2use(2).row = 969; pixels2use(2).col = 984;

% 17_Jan_2024_coast - optical depths of 9.22, 11.6, 14.53 19.8
% pixels2use.row = [969, 969, 969, 969];
% pixels2use.col = [989, 984, 980, 957];

% 17_Jan_2024_coast - optical depths of 8.7, 9.22, 9.68, 10.3, 11.6, 12.54,
% 13.61, 14.53, 16.79, 19.8
% pixels2use.row = [969, 969, 969, 969, 969, 969, 969, 969, 969, 969];
% pixels2use.col = [991, 989, 987, 986, 984, 980, 976, 974, 966, 957];


% 17_Jan_2024_coast - optical depths of 3.2, 6.6, 10.3, 12.54, 14.53, 19.8
% pixels2use.row = [932, 932, 969, 969, 969, 969];
% pixels2use.col = [970, 960, 986, 980, 974, 987];

% 17_Jan_2024_coast - optical depths of 10.3
% pixels2use.row = 969;
% pixels2use.col = 986;


% 17_Jan_2024_coast - optical depths of 0(clear sky), 6.6, 9.22, 11.6, 14.53 19.8
pixels2use.row = [1021, 932, 969, 969, 969, 969];
pixels2use.col = [536, 960, 989, 984, 980, 957];


% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, [size(emit.radiance.measurements,1),...
    size(emit.radiance.measurements, 2)]);



%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_emit_inputs_hyperspectral_top_bottom(emitDataFolder, folder2save, L1B_fileName, emit);
%inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

% *** Check Inputs ***

%% Define the spectral response function of EMIT for the desired Bands

% create the spectral response functions
emit.spec_response = create_EMIT_specResponse(emit, inputs);


%% Define the solar source file name and read in the solar source data

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);


%% Convert radiance measurements to TOA reflectance for the desired pixels

emit = convert_EMIT_radiance_2_reflectance(emit, inputs);


%% Compute the radiance measurement uncertainty 

emit.radiance.uncertainty = compute_EMIT_radiance_uncertainty(emit);


%% Compute the reflectance uncertainty

emit.reflectance.uncertainty = compute_EMIT_reflectance_uncertainty(emit, inputs);


%% Check the thermodynamic phase of the defined pixels

inputs = determine_cloud_phase_emit(emit, pixels2use);

%% Compute the TBLUT retrieval estimate

tic
tblut_retrieval = TBLUT_forEMIT(emit, emitDataFolder, folder2save, pixels2use);
disp([newline, 'TBLUT retrieval took ', num2str(toc), 'seconds to run', newline])

%% Create the Model and Measurement prior
 

inputs = create_model_prior_covariance_EMIT_top_bottom(inputs, pixels2use, tblut_retrieval, true);
%inputs = create_model_prior_covariance_EMIT_top_middle(inputs, pixels2use, tblut_retrieval, true);

inputs = create_EMIT_measurement_covariance(inputs, emit, pixels2use);


%% Use the tblut retrieval as the initial guess for the hyperspectral retrieval

% Compute the retrieval variables
tic
[retrieval, inputs] = calc_retrieval_gauss_newton_4EMIT_top_bottom(inputs, emit ,pixels2use);
disp([newline, 'Multispectral retrieval took ', num2str(toc), 'seconds to run', newline])
%[retrieval, inputs] = calc_retrieval_gauss_newton_4EMIT_top_middle(inputs, emit ,pixels2use);

% --- save the output ---
 save([inputs.folder2save.reflectance_calcs, inputs.reflectance_calculations_fileName],...
        "inputs", "pixels2use", "retrieval", "tblut_retrieval"); % save inputSettings to the same folder as the input and output file


%% Make plot of the retrieved profile

plot_EMIT_retrieved_vertProf(emit, retrieval, tblut_retrieval)

%% Make one-2-one plot comparing the estimate of reflectance from the retrieved variables and the true measurements
% plot rms as a function of reflectance
figure;

% plot a one-to-one line
plot(linspace(0,1,100), linspace(0,1,100), 'k', "LineWidth", 1)
hold on

% plot the emit reflectance versus the forward model computed reflectance
errorbar(emit.reflectance.value(inputs.bands2run,1), retrieval.computed_reflectance,...
    emit.reflectance.uncertainty(inputs.bands2run,1), 'horizontal', '.', 'markersize', 25,...
    'Color', mySavedColors(1, 'fixed'))

xlim([0.95 * min([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance]),...
    1.05 * max([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance])])
ylim([0.95 * min([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance]),...
    1.05 * max([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance])])
grid on; grid minor

xlabel('EMIT Reflectance ($1/sr$)', 'Interpreter', 'latex', 'Fontsize', 35);
ylabel('Calculated Reflectance ($1/sr$)', 'Interpreter', 'latex', 'Fontsize', 35);
title(['RMS = ', num2str(retrieval.rms_residual{1}(end))], 'Interpreter', 'latex', 'Fontsize', 35);

% set figure size
set(gcf, 'Position', [0 0 700 700])


%% Make one-2-one plot comparing the estimate of reflectance from the retrieved variables and the true measurements
% plot the rms of the percent difference


figure;

% plot a one-to-one line
plot(linspace(0,1,100), linspace(0,1,100), 'k', "LineWidth", 1)
hold on

% plot the emit reflectance versus the forward model computed reflectance
errorbar(emit.reflectance.value(inputs.bands2run,1), retrieval.computed_reflectance,...
    emit.reflectance.uncertainty(inputs.bands2run,1), 'horizontal', '.', 'markersize', 25,...
    'Color', mySavedColors(1, 'fixed'))

xlim([0.95 * min([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance]),...
    1.05 * max([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance])])
ylim([0.95 * min([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance]),...
    1.05 * max([emit.reflectance.value(inputs.bands2run,1); retrieval.computed_reflectance])])
grid on; grid minor

xlabel('EMIT Reflectance ($1/sr$)', 'Interpreter', 'latex', 'Fontsize', 35);
ylabel('Retrieved Reflectance ($1/sr$)', 'Interpreter', 'latex', 'Fontsize', 35);

% compute the rms percent difference
rms_percent_diff = sqrt(mean((100*(1 - retrieval.computed_reflectance./emit.reflectance.value(inputs.bands2run))).^2));


title(['RMS Percent Difference = ', num2str(round(rms_percent_diff,2)), '\%'],...
    'Interpreter', 'latex', 'Fontsize', 35);

% set figure size
set(gcf, 'Position', [0 0 700 700])


%% Make spectral plot comparing the estimate of reflectance from the retrieved variables and the true measurements

figure;

% plot the emit reflectance versus wavelength
errorbar(emit.radiance.wavelength(inputs.bands2run), emit.reflectance.value(inputs.bands2run,1),...
    emit.reflectance.uncertainty(inputs.bands2run), 'vertical', '.-', 'MarkerSize', 25,...
    'LineWidth', 1, 'Color', mySavedColors(1,'fixed'))

hold on;
plot(emit.radiance.wavelength(inputs.bands2run), retrieval.computed_reflectance,...
    '.-', 'MarkerSize', 25,...
    'LineWidth', 1, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor

xlabel('Wavelength ($nm$)', 'Interpreter', 'latex', 'Fontsize', 35);
ylabel('Reflectance ($1/sr$)', 'Interpreter', 'latex', 'Fontsize', 35);
title(['Comparison between Calculated and EMIT Reflectance'],...
    'Interpreter', 'latex', 'Fontsize', 35);
grid on; grid minor

% set legend
legend('EMIT', 'Retrieved', 'location', 'best', 'Interpreter', 'latex', 'FontSize', 30)


% set figure size
set(gcf, 'Position', [0 0 1200 675])
