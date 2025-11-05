%% Retreive vertical profiles using simulated HySICS reflectance measurements
% This is a function that retrieves a vertical droplet profile and the
% integrated column water vapor amount above cloud

% This script retrieves 4 variables: ln(r_top), ln(r_bot), ln(tau_c), and acpw


% INPUTS

% (1) filename - .mat file that containts the simulated HySICS measurements

% (2) folder_paths - structure containing all directories needed for
% writing and reading libRadtran INP/OUT

% (3) print_status_updates - a true or false variable that tells the
% function to print status messages into the command window, if true

% (4) print_libRadtran_err - a true or false that tells the function to
% write libRadtran error messages, if true


% By Andrew John Buggee


function [tblut_retrieval, acpw_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropletProfile_HySICS_ver4_logState_lowUncertainty(filenames,...
    folder_paths, print_status_updates, print_libRadtran_err)


%% -- Start Parallel pool

start_parallel_pool(folder_paths.which_computer)


for ff = 1:length(filenames)


    if print_status_updates==true
        disp([newline, 'Processing file: ', filenames{ff}, '...', newline])
    end


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



    %% Load the simualted HySICS measurements


    % Load the simulated measurement
    simulated_measurements = load([folder_paths.HySICS_simulated_spectra, filenames{ff}]);



    % *** Check to see if these measure have added uncertainty or not ***

    if isfield(simulated_measurements, 'Refl_model')==true && isfield(simulated_measurements, 'Refl_model_with_noise')==true

        if print_status_updates==true
            disp([newline, 'Using measurements with added uncertianty...', newline])
        end

        % Then we're using measurements with noise and we set this to be the
        % Reflectance measurements
        simulated_measurements.Refl_model = simulated_measurements.Refl_model_with_noise;

    elseif isfield(simulated_measurements, 'Refl_model')==false && isfield(simulated_measurements, 'Refl_model_with_noise')==true

        if print_status_updates==true
            disp([newline, 'Using measurements with added uncertianty...', newline])
        end

        % Then we're using measurements with noise and we set this to be the
        % Reflectance measurements
        simulated_measurements.Refl_model = simulated_measurements.Refl_model_with_noise;


    end



    %% Create the name of the file to save all output to

    rev = 1;


    if size(simulated_measurements.changing_variables,2)>6

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

    elseif size(simulated_measurements.changing_variables,2)<=6

        folder_paths.saveOutput_filename = [folder_paths.HySICS_retrievals,'dropletRetrieval_HySICS_',...
            num2str(numel(simulated_measurements.inputs.bands2run)), 'bands_',...
            num2str(100*simulated_measurements.inputs.measurement.uncert), '%_uncert',...
            '_rTop_', num2str(simulated_measurements.changing_variables(1,1)),...
            '_rBot_', num2str(simulated_measurements.changing_variables(1,2)),...
            '_tauC_', num2str(simulated_measurements.changing_variables(1,3)),...
            '_tcwv_', num2str(simulated_measurements.inputs.RT.waterVapor_column),...
            '_vza_', num2str(round(simulated_measurements.inputs.RT.vza)),...
            '_vaz_', num2str(round(simulated_measurements.inputs.RT.vaz)),...
            '_sza_', num2str(round(simulated_measurements.inputs.RT.sza)),...
            '_saz_', num2str(round(simulated_measurements.inputs.RT.phi0)),...
            '_sim-ran-on-',char(datetime("today")),'_1.mat'];

    end




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


    tblut_retrieval = TBLUT_for_HySICS_ver2(simulated_measurements, folder_paths, print_status_updates, print_libRadtran_err);


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
    GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS_ver2(simulated_measurements, print_libRadtran_err);

    if print_status_updates==true
        disp('Dont forget to check the inputs and change if needed!!')
    end

    GN_inputs.calc_type = 'forward_model_calcs_forRetrieval';

    % what was the assumed above cloud column water vapor path?

    %% We're retrieving above cloud column water vapor. Make sure input settings are correct

    GN_inputs.RT.modify_total_columnWaterVapor = false;             % don't modify the full column
    GN_inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % modify the column above the cloud

    %% override optical depth

    % Do you want to manually set the optical depth?
    GN_inputs.RT.modify_wc_opticalDepth = true;


    %% CREATE MODEL PRIOR AND COVARIANCE MATRIX AND MEASUREMENT COVARIANCE

    % I don't need anything but the covariance matrix and the expected values
    %inputs = create_model_prior(inputs,data_inputs);

    % -------------------------------------------------------
    % do you want to use your estimates or the MODIS estimate?
    % -------------------------------------------------------

    use_TBLUT_estimates = true;

    % Create inputs to retrieve r_top, r_bot, tau_c, acpw
    %     GN_inputs = create_model_prior_covariance_HySICS_ver2(GN_inputs, tblut_retrieval, use_TBLUT_estimates, acpw_retrieval);
    GN_inputs = create_model_prior_covariance_HySICS_ver2_lowUncertainty1(GN_inputs, tblut_retrieval, use_TBLUT_estimates, acpw_retrieval);
    %     GN_inputs = create_model_prior_covariance_HySICS_ver2_highUncertainty1(GN_inputs, tblut_retrieval, use_TBLUT_estimates, acpw_retrieval);


    GN_inputs = create_HySICS_measurement_covariance(GN_inputs, simulated_measurements);


    %% CALCULATE RETRIEVAL PARAMETERS

    tic

    % --------------------------------------------------------------
    % ---------------- Retrieve Vertical Profile! ------------------
    % --------------------------------------------------------------
    [GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_HySICS_ver2(GN_inputs, simulated_measurements, folder_paths, print_status_updates);
    % --------------------------------------------------------------
    % --------------------------------------------------------------

    if print_status_updates==true
        disp([newline, 'Hyperspectral retrieval completed in ', num2str(toc), ' seconds', newline])
    end


    %%
    % ----------------------------------------------
    % ------------ SAVE OUTPUT STRUCTURE -----------
    % ----------------------------------------------

    % Save the version without an measurement uncertainty. Then we can add
    % uncertainty and save the new file


    % If the folder path doesn't exit, create a new directory
    % 7 means a directory exists with the name defined below
    if exist(folder_paths.HySICS_retrievals, 'dir')~=7

        mkdir(folder_paths.HySICS_retrievals)

    end

    % 2 means the file exists with a .mat extension
    if exist(folder_paths.saveOutput_filename, 'file')==2
        % append
        save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "folder_paths", '-append');

    else

        save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "folder_paths", "tblut_retrieval", "acpw_retrieval");

    end


    % if exist(folder_paths.saveOutput_filename, 'file')==2
    %     % append
    %     save(folder_paths.saveOutput_filename, "folder_paths", "GN_inputs", '-append');
    %
    % else
    %     save(folder_paths.saveOutput_filename, "folder_paths", "tblut_retrieval", "acpw_retrieval", "GN_inputs");
    %
    % end








    if print_status_updates==true
        disp([newline, 'Total time to run retrieval on was ', num2str(toc), ' seconds', newline])
    end



    %% Clear variables and start again!


    if length(filenames)>1 && ff~=length(filenames)

        clear simulated_measurements tblut_retrieval acpw_retrieval GN_inputs GN_outputs

    end


end



end