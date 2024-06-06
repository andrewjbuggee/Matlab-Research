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

% -------------------------------------


%% Load EMIT data and define folders 

[emitDataPath, folder2save] = define_EMIT_dataPath_and_saveFolders();

[emit,L1B_fileName] = retrieveEMIT_data([emitDataPath, emitDataFolder]);


%% Define the pixels to use for the retrieval


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% row = 1112;
% col = 974;

% 17_Jan_2024_coast - small optical depth
% pixels2use.row = [912, 913];
% pixels2use.col = [929, 929];

% 17_Jan_2024_coast - optical depth of 6.6
pixels2use.row = [932];
pixels2use.col = [960];

% 17_Jan_2024_coast - optical depth of 3.2 and 3.8
% pixels2use.row = [932, 932];
% pixels2use.col = [970, 969];

% 17_Jan_2024_coast - optical depth of 3.2
% pixels2use.row = [932];
% pixels2use.col = [970];


% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));



%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
%inputs = create_emit_inputs_hyperspectral_top_bottom(emitDataFolder, folder2save, L1B_fileName, emit);
inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

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

tic
inputs = check_EMIT_therodynamic_phase(emit, inputs);
toc




%% Compute the TBLUT retrieval estimate

tblut_retrieval = TBLUT_forEMIT(emit, emitDataFolder, folder2save, pixels2use);


%% Create the Model and Measurement prior
 

%inputs = create_model_prior_covariance_EMIT_top_bottom(inputs, pixels2use, tblut_retrieval, true);
inputs = create_model_prior_covariance_EMIT_top_middle(inputs, pixels2use, tblut_retrieval, true);

inputs = create_EMIT_measurement_covariance(inputs, emit, pixels2use);


%% Use the tblut retrieval as the initial guess for the hyperspectral retrieval

% Compute the retrieval variables
%[retrieval, inputs] = calc_retrieval_gauss_newton_4EMIT_top_bottom(inputs, emit ,pixels2use);
[retrieval, inputs] = calc_retrieval_gauss_newton_4EMIT_top_middle(inputs, emit ,pixels2use);

% --- save the output ---
 save([inputs.folder2save.reflectance_calcs, inputs.reflectance_calculations_fileName],...
        "inputs", "pixels2use", "retrieval", "tblut_retrieval"); % save inputSettings to the same folder as the input and output file


%% Make plot of the retrieved profile

plot_EMIT_retrieved_vertProf(emit, retrieval, tblut_retrieval)

%% Make plot comparing the estimate of reflectance from the retrieved variables and the true measurements

figure;

% plot a one-to-one line
plot(linspace(0,1,100), linspace(0,1,100), 'k', "LineWidth", 1)
hold on

% plot the emit reflectance versus the forward model computed reflectance
errorbar(emit.reflectance.value(inputs.bands2run,1), retrieval.calculated_reflectance,...
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
