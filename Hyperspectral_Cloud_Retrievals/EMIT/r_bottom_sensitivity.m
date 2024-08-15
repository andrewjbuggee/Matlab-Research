%% Sensitivity to detecting a change in the radius at cloud bottom

% By Andrew John Buggee

%% Compute when the change in reflectivity is greater than the measurement uncertainty for different cloud optical depths

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
% pixels2use.row = 1112;
% pixels2use.col = 974;

% 17_Jan_2024_coast - small optical depth
% pixels2use.row = [912, 913];
% pixels2use.col = [929, 929];

% 17_Jan_2024_coast - optical depth of 6.6
% pixels2use.row = [932];
% pixels2use.col = [960];

% 17_Jan_2024_coast - optical depth of 3.2 and 3.8
% pixels2use.row = [932, 932];
% pixels2use.col = [970, 969];

% 17_Jan_2024_coast - optical depth of 3.2
% pixels2use.row = [932];
% pixels2use.col = [970];

% 17_Jan_2024_coast - optical depths of 3.2, 4.8, 5.4, 5.9, 6.6, 7.2, 7.5
pixels2use.row = [932, 932, 932, 932, 932, 932, 932];
pixels2use.col = [970, 968, 967, 966, 960, 959, 957];

% 17_Jan_2024_coast - optical depths of 8.7, 9.22, 9.68, 10.3, 11.6, 12.54,
% 13.61, 14.53, 16.79, 19.8
% pixels2use.row = [969, 969, 969, 969, 969, 969, 969, 969, 969, 969];
% pixels2use.col = [991, 989, 987, 986, 984, 980, 976, 974, 966, 957];

% 17_Jan_2024_coast - optical depths of 10.3
% pixels2use.row = 969;
% pixels2use.col = 986;



% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));



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


%% Compute the TBLUT retrieval estimate

tblut_retrieval = TBLUT_forEMIT(emit, emitDataFolder, folder2save, pixels2use);


%% Define the state vector

% define the state vector [r_top, r_bottom, tau_c]
%state_vector = [13.275, 6.6, 5.8];

state_vector = [minVals.minRe(4)*1.2, minVals.minRe(4)*0.7, minVals.minTau(4)];

% r_top_bottom = [12, 6];
% tau_vector = 3:2:11;

% r_top_bottom = [8, 4];
% tau_vector = 3:2:19;


%% compute the refltances for the state vector and compute the change due to a changing r_bottom

% Compute the Measurement Change for a single pixel

% Compute the reflectances for the above state vector
measurement_estimate = compute_forward_model_4EMIT_top_bottom(emit, state_vector, inputs, pixels2use, 1)';

[measurement_change, change_in_r_bottom] = compute_reflectanceChange_due_to_rBottom_change(emit, state_vector,...
    measurement_estimate, inputs,pixels2use, 1);






% % Step through different optical depths and compute the measurement change
% measurement_estimate = zeros(length(inputs.bands2run), length(tau_vector));
% measurement_change = cell(1, length(tau_vector));
% 
% for tt = 1:length(tau_vector)
% 
%     disp([newline, 'Iteration: Tau_c = ',num2str(tt), '/', num2str(length(tau_vector)), newline])
% 
%     % define the current state vector [r_top, r_bottom, tau_c]
%     state_vector = [r_top_bottom, tau_vector(tt)];
% 
%     % Compute the reflectances for the above state vector
%     measurement_estimate(:,tt) = compute_forward_model_4EMIT_top_bottom(emit, state_vector, inputs, pixels2use, 1)';
% 
%     [measurement_change{tt}, change_in_r_bottom] = compute_reflectanceChange_due_to_rBottom_change(emit, state_vector,...
%         measurement_estimate(:,tt), inputs,pixels2use, 1);
% 
% 
% 
% end

%% --- save the output ---
 save([inputs.folder2save.reflectance_calcs, 'r_bottom_analysis_', char(datetime("today")),'.mat'],...
        "inputs", "pixels2use", "measurement_change", "change_in_r_bottom", "r_top_bottom", "tau_vector"); % save inputSettings to the same folder as the input and output file


%% Make Plot of single optical depth with changing r_bot

% which tau would you like to plot?
tau_idx = 1;


figure;


for nn = 1:size(measurement_change{tau_idx},2)

    plot(emit.radiance.wavelength(inputs.bands2run), abs(measurement_change{tau_idx}(:,nn)), 'Color', ...
        mySavedColors(nn, 'fixed'))
    hold on
end

grid on; grid minor

% plot the EMIT noise
plot(emit.radiance.wavelength(inputs.bands2run), emit.reflectance.uncertainty(inputs.bands2run),...
    'Color', 'black', 'LineStyle', ':', 'LineWidth', 3)

% make legend
legend([strcat('$\triangle r_{bot} = $', string(change_in_r_bottom(1:size(measurement_change{tau_idx},2))),...
    ' $\mu m$'), 'EMIT Uncertainty'],...
    'Location', 'best', 'Interpreter','latex','FontSize', 25)

xlabel('Wavelength $(nm)$', 'Interpreter','latex')
ylabel('$\triangle$ Reflectance $(1/sr)$', 'Interpreter','latex')
title(['$\tau_c = $ ', num2str(tau_vector(tau_idx))], 'Interpreter', 'latex')
set(gcf, 'Position', [0 0 1200 700])


%% Plot measurement change for a single optical depth with a shaded region for the uncertainty


% which tau would you like to plot?
tau_idx = 1;


figure;

% plot the EMIT noise as an transparent area centered around 0
x = [emit.radiance.wavelength(inputs.bands2run); flipud(emit.radiance.wavelength(inputs.bands2run))];
y = [emit.reflectance.uncertainty(inputs.bands2run); -flipud(emit.reflectance.uncertainty(inputs.bands2run))];
fill(x,y, [0 0.4470 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', 0.4)

hold on

% --- Plot th measurement change for each change in r_bottom ---

for nn = 1:size(measurement_change{tau_idx},2)

    plot(emit.radiance.wavelength(inputs.bands2run), measurement_change{tau_idx}(:,nn), 'Color', ...
        mySavedColors(nn, 'fixed'))
end

grid on; grid minor




% make legend
legend(['EMIT Uncertainty', strcat('$\triangle r_{bot} = $',...
    string(change_in_r_bottom(1:size(measurement_change{tau_idx},2))),' $\mu m$')],...
    'Location', 'best', 'Interpreter','latex','FontSize', 25)

xlabel('Wavelength $(nm)$', 'Interpreter','latex')
ylabel('$\triangle$ Reflectance $(1/sr)$', 'Interpreter','latex')
title(['$\tau_c = $ ', num2str(tau_vector(tau_idx))], 'Interpreter', 'latex')
set(gcf, 'Position', [0 0 1200 700])




%% Find when half of the channels have a relfectance that exceeds the measurement uncertainty

% step through each optical depth and find the number of wavelengths who's
% change in reflectance exceeds the measurement uncertainty

num_exceed_uncertainty = zeros(length(tau_vector), length(change_in_r_bottom));

figure;

for tt = 1:length(tau_vector)

    num_exceed_uncertainty(tt,:) = sum(abs(measurement_change{tt}) > ...
        repmat(emit.reflectance.uncertainty(inputs.bands2run), 1, length(change_in_r_bottom)), 1);

    plot(change_in_r_bottom, num_exceed_uncertainty(tt, :)./length(inputs.bands2run),...
        'Color', mySavedColors(tt, 'fixed'))
    hold on

end

grid on; grid minor

% make legend
legend(strcat('$\tau_{c} = $', string(tau_vector)),...
    'Location', 'best', 'Interpreter','latex','FontSize', 25)

xlabel('$\triangle r_{bot}$ $(\mu m)$', 'Interpreter','latex')
title('Fraction of Wavelenghts that exceed measurement uncertainty', 'Interpreter','latex')
set(gcf, 'Position', [0 0 1200 700])




%% Plot 






