%% Retreive vertical profiles in a loop using simulated HySICS reflectance measurements
% The for loop changes the simulated HySICS measurement file. Provide
% multiple file names to loop through.

% This script retrieves 4 variables: r_top, r_bot, tau_c, and cwvs



% By Andrew John Buggee



%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


%% Define the HySICS folders for data and storage


folder_paths = define_folderPaths_for_HySICS(0);

which_computer = folder_paths.which_computer;
%% print status updates?

% Would you like each function to print messages to the command terminal on
% the status of the computation?
print_status_updates = false;

% Would you like to print long libRadtran error messages?
print_libRadtran_err = false;



%% -- Start Parallel pool

start_parallel_pool(which_computer)



%%   Delete old files?

% First, delete files in the HySICS INP folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% delete old water vapor profiles
delete([folder_paths.atm_folder_path, '*-aboveCloud.DAT'])

% delete old MIE files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])

%% LOAD SIMULATED HYSICS DATA

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


    % load all filenames in the folder defined above.
    filenames = dir([folder_paths.HySICS_simulated_spectra, '*.mat']);


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % load all filenames in the folder defined above.
    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % load all filenames in the folder defined above that start with a
    % specific block of text
    filenames = dir([folder_paths.HySICS_simulated_spectra,...
        'simulated_spectra_HySICS_reflectance_66bands_0.001%_uncert_rTop_10_rBot_5_tauC_11_tcwv_14_vza_7*.mat']);


end

%% Step through each measurement and perform the retrieval

tic

for nn = 1:size(filenames, 1)




    % Load the simulated measurement
    simulated_measurements = load([folder_paths.HySICS_simulated_spectra,filenames(nn).name]);


    % *** Check to see if these measure have added uncertainty or not ***

    if isfield(simulated_measurements, 'Refl_model_with_noise')==true

        if print_status_updates==true
            disp([newline, 'Using measurements with added uncertianty...', newline])
        end

        % Then we're using measurements with noise and we set this to be the
        % Reflectance measurements
        simulated_measurements.Refl_model = simulated_measurements.Refl_model_with_noise;

    end

    %% Create the name of the file to save all output to

    rev = 1;


    folder_paths.saveOutput_filename = [folder_paths.HySICS_retrievals,'dropletRetrieval_HySICS_',...
        num2str(numel(simulated_measurements.inputs.bands2run)), 'bands_',...
        num2str(100*simulated_measurements.inputs.measurement.uncert), '%_uncert',...
        '_rTop_', num2str(simulated_measurements.changing_variables(1,1)),...
        '_rBot_', num2str(simulated_measurements.changing_variables(1,2)),...
        '_tauC_', num2str(simulated_measurements.changing_variables(1,3)),...
        '_tcwv_', num2str(simulated_measurements.changing_variables(1,4)),...
        '_vza_', num2str(round(simulated_measurements.inputs.RT.vza)),...
        '_vaz_', num2str(round(simulated_measurements.inputs.RT.vaz)),...
        '_sza_', num2str(round(simulated_measurements.inputs.RT.sza)),...
        '_saz_', num2str(round(simulated_measurements.inputs.RT.phi0)),...
        '_sim-ran-on-',char(datetime("today")),'_1.mat'];




    while isfile(folder_paths.saveOutput_filename)
        rev = rev+1;
        if rev<10
            folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-5), num2str(rev),'.mat'];
        elseif rev>10
            folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-6), num2str(rev),'.mat'];
        end
    end


    %% Compute the Two-Band Look-up Table retrieval of effective radius and optical depth

    if print_status_updates==true
        disp([newline, 'Computing the TBLUT retrieval...', newline])
        tic
    end


    tblut_retrieval = TBLUT_for_HySICS_ver2(simulated_measurements, folder_paths,print_status_updates, print_libRadtran_err);


    if print_status_updates==true
        disp([newline, 'TBLUT retrieval completed in ', num2str(toc), ' seconds', newline])
    end


    %% Compute the multispectral estimate of the above cloud column water vapor

    if print_status_updates==true
        disp([newline, 'Computing the ACPW retrieval...', newline])
        tic
    end


    acpw_retrieval = ACPW_retrieval_for_HySICS(simulated_measurements, tblut_retrieval, folder_paths, print_status_updates, print_libRadtran_err);


    if print_status_updates==true
        disp([newline, 'ACPW retrieval completed in ', num2str(toc), ' seconds', newline])
    end


    %% CREATE GAUSS-NEWTON INPUTS

    % Create inputs to retrieve r_top, r_bot, tau_c, cwv
    GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS_ver2(simulated_measurements);

    % if print_status_updates==true
    %     disp('Dont forget to check the inputs and change if needed!!')
    % end
    % 
    % GN_inputs.calc_type = 'forward_model_calcs_forRetrieval';
    % 
    % % what was the assumed above cloud column water vapor path?
    % 
    % %% We're retrieving above cloud column water vapor. Make sure input settings are correct
    % 
    % GN_inputs.RT.modify_total_columnWaterVapor = false;             % don't modify the full column
    % GN_inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % modify the column above the cloud
    % 
    % 
    % 
    % %% CREATE MODEL PRIOR AND COVARIANCE MATRIX AND MEASUREMENT COVARIANCE
    % 
    % % I don't need anything but the covariance matrix and the expected values
    % %inputs = create_model_prior(inputs,data_inputs);
    % 
    % % -------------------------------------------------------
    % % do you want to use your estimates or the MODIS estimate?
    % % -------------------------------------------------------
    % 
    % use_TBLUT_estimates = true;
    % 
    % % Create inputs to retrieve r_top, r_bot, tau_c, cwv
    % GN_inputs = create_model_prior_covariance_HySICS_ver2(GN_inputs, tblut_retrieval, use_TBLUT_estimates);
    % 
    % 
    % GN_inputs = create_HySICS_measurement_covariance(GN_inputs, simulated_measurements);
    % 
    % 
    % %% CALCULATE RETRIEVAL PARAMETERS
    % 
    % tic
    % 
    % % --------------------------------------------------------------
    % % ---------------- Retrieve Vertical Profile! ------------------
    % % --------------------------------------------------------------
    % [GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_HySICS_ver2(GN_inputs, simulated_measurements, folder_paths, print_status_updates);
    % % --------------------------------------------------------------
    % % --------------------------------------------------------------
    % 
    % if print_status_updates==true
    %     disp([newline, 'Hyperspectral retrieval completed in ', num2str(toc), ' seconds', newline])
    % end


    %%
    % ----------------------------------------------
    % ------------ SAVE OUTPUT STRUCTURE -----------
    % ----------------------------------------------

    % Save the version without an measurement uncertainty. Then we can add
    % uncertainty and save the new file


    % If the folder path doesn't exit, create a new directory
    if ~exist(folder_paths.HySICS_retrievals, 'dir')

        mkdir(folder_paths.HySICS_retrievals)

    end

    % if exist(folder_paths.saveOutput_filename, 'file')==true
    %     % append
    %     save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "folder_paths", "tblut_retrieval", '-append');
    % 
    % else
    %     save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "folder_paths", "tblut_retrieval");
    % 
    % end


    if exist(folder_paths.saveOutput_filename, 'file')==true
        % append
        save(folder_paths.saveOutput_filename, "folder_paths", "GN_inputs", '-append');

    else
        save(folder_paths.saveOutput_filename, "folder_paths", "tblut_retrieval", "acpw_retrieval", "GN_inputs");

    end




    %% Clear variables and start again!

    if nn~=size(filenames,1)

        clear simulated_measurements tblut_retrieval acpw_retrieval GN_inputs GN_outputs

    end


end


if print_status_updates==true
    disp([newline, 'Total time to run retrieval on ', num2str(size(filenames,1)), ' files was ', num2str(toc), ' seconds', newline])
end


