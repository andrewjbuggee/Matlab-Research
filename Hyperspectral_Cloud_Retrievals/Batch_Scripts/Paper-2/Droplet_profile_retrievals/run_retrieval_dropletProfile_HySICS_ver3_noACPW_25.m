%% Retreive vertical profiles using simulated HySICS reflectance measurements
% This is a function that retrieves a vertical droplet profile

% This script retrieves 4 variables: r_top, r_bot, tau_c


% INPUTS

% (1) filename - .mat file that containts the simulated HySICS measurements

% (2) folder_paths - structure containing all directories needed for
% writing and reading libRadtran INP/OUT

% (3) print_status_updates - a true or false variable that tells the
% function to print status messages into the command window, if true

% (4) print_libRadtran_err - a true or false that tells the function to
% write libRadtran error messages, if true


% By Andrew John Buggee


function [tblut_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropletProfile_HySICS_ver3_noACPW_25(filenames, folder_paths, print_status_updates, print_libRadtran_err)


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


    % delete old MIE files
    delete([folder_paths.libRadtran_mie_folder, '*.INP'])
    delete([folder_paths.libRadtran_mie_folder, '*.OUT'])



    %% Load the simualted HySICS measurements


    % Load the simulated measurement
    simulated_measurements = load([folder_paths.HySICS_simulated_spectra, filenames{ff}]);


    %% select 35 wavelengths from paper 1 that do not sample regions with strong water vapor absorption

    % Paper 1 - Figures 7 and 8 - 35 spectral channels that avoid water vapor
    % and other gaseous absorbers
    % inputs.bands2run = [49, 57, 69, 86, 103, 166, 169, 171, 174, 217, 220,...
    %     222, 224, 227, 237, 288, 290, 293, 388, 390, 393,...
    %     426, 434, 436, 570, 574, 577, 579, 582, 613, 616,...
    %     618, 620, 623, 625]';


    if size(simulated_measurements.spec_response.wavelength,1)>35

        paper1_wl_midBand = [498.95, 523.45, 560.25, 612.25, 664.35, 857.35, 866.55,...
            872.656, 881.856, 1013.551,	1022.751, 1028.851, 1035.051, 1044.251, 1074.851,...
            1231.051, 1237.151, 1246.351,	1537.35,	1543.55,	1552.65,	1653.75,...
            1678.25,	1684.45,	2094.85,	2107.15,	2116.35,	2122.45,...
            2131.65,	2226.55,	2235.75,	2241.95,	2248.05,	2257.25,	2263.35];

        idx_35 = zeros(numel(paper1_wl_midBand), 1);

        % step through the measurements and select the 35 bands closest to
        % these midPoints
        % first, compute the mid points of the bands within the simualted
        % measurement
        sim_meas_midBand = mean(simulated_measurements.spec_response.wavelength, 2);

        for nn = 1:length(paper1_wl_midBand)


            [~, idx_35(nn)] = min(abs(sim_meas_midBand - paper1_wl_midBand(nn)));


        end

        % use the above index to grab the subset of measurements applicable
        % for this retrieval
        simulated_measurements.Refl_model = simulated_measurements.Refl_model(idx_35);
        simulated_measurements.Refl_model_uncert = simulated_measurements.Refl_model_uncert(idx_35);
        simulated_measurements.Refl_model_with_noise = simulated_measurements.Refl_model_with_noise(idx_35);
        simulated_measurements.changing_variables = simulated_measurements.changing_variables(idx_35, :);
        simulated_measurements.spec_response.value = simulated_measurements.spec_response.value(idx_35, :);
        simulated_measurements.spec_response.wavelength = simulated_measurements.spec_response.wavelength(idx_35, :);

        simulated_measurements.inputs.bands2run = simulated_measurements.inputs.bands2run(idx_35);
        simulated_measurements.inputs.RT.wavelengths2run = simulated_measurements.inputs.RT.wavelengths2run(idx_35, :);


    end

    %% Check to see if there is uncertainty



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

        folder_paths.saveOutput_filename = [folder_paths.HySICS_retrievals,'dropletRetrieval_noACPW_HySICS_',...
            num2str(numel(simulated_measurements.inputs.bands2run)), 'bands_',...
            num2str(100*simulated_measurements.inputs.measurement.uncert), '%_uncert',...
            '_rTop_', num2str(simulated_measurements.changing_variables(1,1)),...
            '_rBot_', num2str(simulated_measurements.changing_variables(1,2)),...
            '_tauC_', num2str(simulated_measurements.changing_variables(1,3)),...
            '_tcwv_', num2str(simulated_measurements.changing_variables(1,4)),...
            '_tcwv-assumption_25',...
            '_vza_', num2str(round(simulated_measurements.inputs.RT.vza)),...
            '_vaz_', num2str(round(simulated_measurements.inputs.RT.vaz)),...
            '_sza_', num2str(round(simulated_measurements.inputs.RT.sza)),...
            '_saz_', num2str(round(simulated_measurements.inputs.RT.phi0)),...
            '_sim-ran-on-',char(datetime("today")),'_1.mat'];

    elseif size(simulated_measurements.changing_variables,2)<=6

        folder_paths.saveOutput_filename = [folder_paths.HySICS_retrievals,'dropletRetrieval_noACPW_HySICS_',...
            num2str(numel(simulated_measurements.inputs.bands2run)), 'bands_',...
            num2str(100*simulated_measurements.inputs.measurement.uncert), '%_uncert',...
            '_rTop_', num2str(simulated_measurements.changing_variables(1,1)),...
            '_rBot_', num2str(simulated_measurements.changing_variables(1,2)),...
            '_tauC_', num2str(simulated_measurements.changing_variables(1,3)),...
            '_tcwv_', num2str(simulated_measurements.inputs.RT.waterVapor_column),...
            '_tcwv-assumption_25',...
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



    %% CREATE GAUSS-NEWTON INPUTS

    % Create inputs to retrieve r_top, r_bot, tau_c, cwv
    GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS(simulated_measurements, print_libRadtran_err);

    if print_status_updates==true
        disp('Dont forget to check the inputs and change if needed!!')
    end

    GN_inputs.calc_type = 'forward_model_calcs_forRetrieval';

    % what was the assumed above cloud column water vapor path?

    %% This retrieval does NOT retrieve column water vapor. What should the forward model assumption be?

    GN_inputs.RT.modify_total_columnWaterVapor = true;             % don't modify the full column
    GN_inputs.RT.waterVapor_column = 25;   % mm - milimeters of water condensed in a column

    GN_inputs.RT.modify_aboveCloud_columnWaterVapor = false;         % modify the column above the cloud



    %% CREATE MODEL PRIOR AND COVARIANCE MATRIX AND MEASUREMENT COVARIANCE

    % I don't need anything but the covariance matrix and the expected values
    %inputs = create_model_prior(inputs,data_inputs);

    % -------------------------------------------------------
    % do you want to use your estimates or the MODIS estimate?
    % -------------------------------------------------------

    use_TBLUT_estimates = true;

    % Create inputs to retrieve r_top, r_bot, tau_c
    GN_inputs = create_model_prior_covariance_HySICS(GN_inputs, tblut_retrieval, use_TBLUT_estimates);


    GN_inputs = create_HySICS_measurement_covariance(GN_inputs, simulated_measurements);


    %% CALCULATE RETRIEVAL PARAMETERS

    tic

    % --------------------------------------------------------------
    % ---------------- Retrieve Vertical Profile! ------------------
    % --------------------------------------------------------------
    [GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_HySICS(GN_inputs, simulated_measurements, folder_paths, print_status_updates);
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

        save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "folder_paths", "tblut_retrieval");

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

        clear simulated_measurements tblut_retrieval GN_inputs GN_outputs

    end


end



end