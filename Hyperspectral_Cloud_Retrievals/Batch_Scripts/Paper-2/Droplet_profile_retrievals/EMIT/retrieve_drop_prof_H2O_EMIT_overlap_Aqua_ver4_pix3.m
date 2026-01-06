%% Optimal estimation of a vertical droplet profile using EMIT data

% Retrieve a droplet profile from EMIT data
% Use spectral reflectances outside of water vapor bands


% By Andrew John Buggee

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


%% Define EMIT Data locations and LibRadTran paths

folder_paths = define_EMIT_dataPath_and_saveFolders(2);
which_computer = folder_paths.which_computer;


%% Would you like to print status updates and/or the libRadtran error file?

print_status_updates = true;

print_libRadtran_err = true;


%% Define the folder of the coincident data set between EMIT and Aqau

% ---------------------------------------------
% ---------- PICK COINCIDENT DATA SET  --------
% ---------------------------------------------


% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folder_paths.coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/'];

    folder_paths.coincident_dataFolder = '2024-09-12/';

elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/'];

    % folder_paths.coincident_dataFolder = '2024-09-12/';

    % folder_paths.coincident_dataFolder = '2024_05_17-T1835/';

    % EMIT pixels masked out
    % folder_paths.coincident_dataFolder = '2023_9_16_T191106_2/';

    % 11 Pixels with H less than 1.6
    % folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';

    % 2 Pixels with H less than 1.1
    % folder_paths.coincident_dataFolder = '2023_9_16_T191130_1/';

    % 14 Pixels with H less than 1.35
    folder_paths.coincident_dataFolder = '2023_9_16_T191142_1/';

elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];


    % 14 Pixels with H less than 1.35
    folder_paths.coincident_dataFolder = '2023_9_16_T191142_1/';

end






%% Open Aqau data and look for overlapping pixels between EMIT and Aqua that meet certain criteria

criteria.cld_phase = 'water';
criteria.cld_cvr = 1;   % cloud fraction
criteria.cld_tau_min = 3;   % cloud optical depth
criteria.cld_tau_max = 30;   % cloud optical depth
criteria.H = 0.1;         % horizontal inhomogeneity index

% plot flag
plot_data = false;

% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - Then compute: overlap.time_difference_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))

[overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
    criteria, plot_data);

% ** If there aren't any pixels found ... **
% Increase the horizontal inhomogeneity index

while isempty(overlap_pixels.modis.linear_idx) == true

    disp([newline, 'No overlaping pixels that meet defined criteria. Increasing H index....', newline])
    criteria.H = criteria.H + 0.25;         % horizontal inhomogeneity index

    % recompute
    [overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
        criteria, plot_data);

    if isempty(overlap_pixels.modis.linear_idx) == false

        % print the H value used
        disp([newline, 'Horizontal Inhomogeneity Index - H = ', num2str(criteria.H), newline])

    end

end

%% Plot all three swaths

figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.cloud.effRadius17,[],1),'.');
hold on; geoscatter(emit.radiance.geo.lat(:), emit.radiance.geo.long(:), 10, 'r.')
hold on; geoscatter(airs.geo.Latitude(:), airs.geo.Longitude(:), 10, 'c.')
hold on; geoscatter(amsr.geo.Latitude(:), amsr.geo.Longitude(:), 10, 'k.')

%% Plot the pixel footprints on the Earth to see the overlap
% Add an RGB true color image for context
clear options
options.use_radiance = false;
[rgb_img, rgb_lat, rgb_lon] = create_modis_true_color(modis, options);

options.show_rgb = true;
options.rgb_image = rgb_img;
options.rgb_lat = rgb_lat;
options.rgb_lon = rgb_lon;
options.latlim = [-30, -20];  % Only show -30° to -20° latitude
% options.lonlim = [-75, -65];  % Only show -75° to -65° longitude

% ** Plot with RGB Image **
% fig = plot_instrument_footprints(modis, emit, amsr, overlap_pixels, options);
fig1 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);

% ** Plot without RGB Image **
options.show_rgb = false;
fig2 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);


%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);




%% Set libRadtran INP directory

folder_paths.libRadtran_inp = [folder_paths.libRadtran_inp, 'EMIT_', folder_paths.coincident_dataFolder(1:end-1), '_',...
    folder_paths.L1B_fileName_emit{1}(27:30),'/'];


% If the folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_inp, 'dir')

    mkdir(folder_paths.libRadtran_inp)

end


%%   Delete old files?

% First, delete files in the HySICS folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% delete old atm files
delete([folder_paths.atm_folder_path, '*.DAT'])

% delete old mie files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])






%% RUN FOR LOOP OVER ALL PIXELS

for pp = 1:length(overlap_pixels.modis.linear_idx)


    %% Check to see if all chosen EMIT pixels have been masked out

    if all(isnan(emit.radiance.measurements(:, pp)))

        warning([newline, 'EMIT pixels are masked out. Moving to pixel ', num2str(pp), '...', newline])

        continue

    else

        disp([newline, 'Retrieving Profile for pixel ', num2str(pp), '...', newline])

    end



    if pp>1

        clear GN_inputs GN_outputs spec_response

    end


    %% Create an input structure that helps write the INP files

    % this is a built-in function that is defined at the bottom of this script
    GN_inputs = create_gauss_newton_inputs_for_emit_ver4_log(emit, print_libRadtran_err);
    %inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

    % *** Check Inputs ***


    %% Adjust inputs that vary per pixel

    % -----------------------------
    % define the solar zenith angle
    % -----------------------------
    GN_inputs.RT.sza = emit.obs.solar.zenith(pp);           % degree



    % -------------------------------
    % define the solar azimuith angle
    % -------------------------------
    % this is how we map EMIT-defined solar azimuth to the LibRadTran
    % definition.
    % LibRadTran defines the solar azimuth clockwise from South.
    % So at 0deg the Sun is due south, 90deg the Sun is due West,
    % and so on. EMIT defines the solar azimuth clockwise from due North.
    % So we need to add 180deg to the EMIT values, but modulo 360, since it
    % needs to wrap around.
    GN_inputs.RT.phi0 = mod(emit.obs.solar.azimuth(pp) + 180, 360);         % degree





    % --------------------------------
    % define the viewing azimuth angle
    % --------------------------------
    % we need the cosine of the zenith viewing angle
    % positive values solve for upwelling radiance, where the sensor is
    % defined to be looking down towrads the Earth's surface. negative
    % values solve for downwelling radiance, where the sensor is looking
    % upwards towards the sky

    % Define the cosine of the zenith viewing angle
    % ------------------------------------------------
    % define the viewing zenith angle
    GN_inputs.RT.vza = double(emit.obs.sensor.zenith(pp)); % values are in degrees;                        % degree




    % --------------------------------
    % define the viewing azimuth angle
    % --------------------------------
    % LibRadTran defines the viewing azimuth clockwise from North.
    % So at 0deg the sensor is in the north looking south, at 90deg
    % the sensor is in the east looking West, and so on.
    % EMIT defines the sensor azimuth clockwise from due North.
    % So we don't need to change the emit values

    GN_inputs.RT.vaz = emit.obs.sensor.azimuth(pp);



    %% Set the wavelenghts!!

    % --- Indexes using same 35 as above, in addition to 29 water vapor bands ---
    % This set has a total of 64 bands. They are not exactly the same set as
    % the 66 HySICS bands used to retrieve column water vapor because the
    % HySICS channels are more narrow.
    % GN_inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 71, 74, 78, 86, 87, 88, 89, 90,...
    %     94, 97, 99, 101, 105, 115, 116, 117, 139, 141, 142, 145, 147, 148, 149, 151, 156,...
    %     157, 158, 172, 175, 176, 187, 189, 190, 210, 212, 213, 214, 215, 216, 217, 218, 219,...
    %     220, 222, 231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';

    % --- Use all 285 spectral channels -
    GN_inputs.bands2run = (1:285)';


    %% Override input settings with MODIS derived values


    % --------------------------------------------
    % *** Use MODIS Cloud Top Height Retrieval ***
    % --------------------------------------------
    % Override the cloud depth
    GN_inputs.RT.H = 0.3;           % km

    % override the cloud top height
    % ** MODIS cloud top height listed in meters is the geopotential height **
    GN_inputs.RT.z_topBottom = [modis.cloud.topHeight(pp)/1e3,...
        (modis.cloud.topHeight(pp)/1e3 - GN_inputs.RT.H)];    % km

    % Update the height vector based on the MODIS cloud top height
    GN_inputs.RT.z_edges = linspace(GN_inputs.RT.z_topBottom(2),...
        GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers+1);   % km - the edges of each layer

    GN_inputs.RT.z = linspace(GN_inputs.RT.z_topBottom(2),...
        GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers);        % km - altitude above ground vector


    % ------------------------------------------------
    % *** Use MODIS above cloud column water vapor ***
    % ------------------------------------------------
    % ** MODIS NIR above cloud column water vapor listed in cm **
    GN_inputs.RT.waterVapor_column = modis.vapor.col_nir(pp) * 10;    % mm



    % ----------------------------------------------------
    % *** Use AIRS temp/press and water vapor profiles ***
    % ----------------------------------------------------

    % first, write a radiosonde.dat file with airs temperature, pressure and
    % relative humidity
    GN_inputs.RT.use_radiosonde_file = true;
    GN_inputs.RT.radiosonde_num_vars = 3;

    GN_inputs.RT.radiosonde_file_T_P_RH = write_AIRS_radiosonde_DAT(airs, folder_paths, pp, [],...
        GN_inputs.RT.radiosonde_num_vars);
    GN_inputs.RT.radiosonde_file_T_P = write_AIRS_radiosonde_DAT(airs, folder_paths, pp, [],...
        GN_inputs.RT.radiosonde_num_vars-1);






    %% This retrieval does retrieve column water vapor.

    GN_inputs.RT.modify_total_columnWaterVapor = false;             % don't modify the full column

    % *** Retreive Column Water Vapor! ***
    GN_inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % modify the column above the cloud

    %% Set output filename

    rev = 1;


    folder_paths.saveOutput_filename = [folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder,...
        'Droplet_profile_retrievals/',...
        num2str(numel(GN_inputs.bands2run)),...
        'bands_EMIT_dropRetrieval_pixel_', num2str(pp),...
        '_ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];




    while isfile(folder_paths.saveOutput_filename)
        rev = rev+1;
        if rev<10
            folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-5), num2str(rev),'.mat'];
        elseif rev>10
            folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-6), num2str(rev),'.mat'];
        end
    end


    %% Define the spectral response function of EMIT for the desired Bands

    % create the spectral response functions
    [GN_inputs, spec_response] = create_EMIT_specResponse(emit, GN_inputs);


    %% Define the solar source file name and read in the solar source data

    % ********* IMPORTANT *************
    % The source flux is integrated with the EMIT spectral response function

    % define the source file using the input resolution
    GN_inputs = define_source_for_EMIT(GN_inputs, emit);


    %% Convert radiance measurements to TOA reflectance for the desired pixels

    if pp==1

        emit = convert_EMIT_radiance_2_reflectance(emit, GN_inputs);


        %% Compute the radiance measurement uncertainty

        [emit.radiance.uncertainty, emit.radiance.uncertainty_percent_perChannel] = compute_EMIT_radiance_uncertainty(emit);


        %% Compute the reflectance uncertainty

        emit.reflectance.uncertainty = compute_EMIT_reflectance_uncertainty(emit, GN_inputs);


    end



    %%  *** Start parallel pool ***

    % Is parpool running?
    p = gcp('nocreate');
    if isempty(p)==true

        % first read the local number of workers avilabile.
        p = parcluster('local');
        % start the cluster with the number of workers available
        if p.NumWorkers>64
            % Likely the amilan128c partition with 2.1 GB per core
            % Leave some cores for overhead
            parpool(p.NumWorkers - 8);

        elseif p.NumWorkers<=64 && p.NumWorkers>10

            parpool(p.NumWorkers);

        elseif p.NumWorkers<=10

            parpool(p.NumWorkers);

        end

    end




    %% Check the thermodynamic phase of the defined pixels

    %GN_inputs.cloudPhase = determine_cloud_phase_emit(emit, pixels2use);


    %% Compute the TBLUT retrieval estimate

    if print_status_updates==true
        disp([newline, 'Computing the TBLUT retrieval...', newline])
        tic
    end

    use_MODIS_AIRS_data = true;

    tblut_retrieval = TBLUT_forEMIT_perPixel(emit, spec_response, folder_paths, print_libRadtran_err, print_status_updates,...
        GN_inputs, use_MODIS_AIRS_data, pp);

    if print_status_updates==true
        disp([newline, 'TBLUT retrieval completed in ', num2str(toc), ' seconds', newline])
    end




    %% Compute the ACPW retrieval estimate

    if print_status_updates==true
        disp([newline, 'Computing the ACPW retrieval...', newline])
        tic
    end

    use_MODIS_AIRS_data = true;

    acpw_retrieval = ACPW_retrieval_for_EMIT_perPixel(emit, spec_response, tblut_retrieval, folder_paths, use_MODIS_AIRS_data,...
        GN_inputs, print_libRadtran_err, print_status_updates, pp);

    if print_status_updates==true
        disp([newline, 'ACPW retrieval completed in ', num2str(toc), ' seconds', newline])
    end

    %% Create the Model and Measurement prior

    use_TBLUT_estimate = true;

    GN_inputs = create_model_prior_covariance_EMIT_top_bottom_ver4_log(GN_inputs, tblut_retrieval,...
        acpw_retrieval, use_TBLUT_estimate);

    GN_inputs = create_EMIT_measurement_cov_ver4_log_no_FM_uncert_perPixel(GN_inputs, emit, pp);


    %% Create the forward model covariance matrix and transform it into measurement space
    % S_b' = K_b * S_b * K_b' (Maahn et al. 2020 E1516)

    % Just re profile uncertainty
    % GN_inputs = create_EMIT_forMod_cov_ver4_log_reProf(GN_inputs);

    % Including re_profile and cloud top height uncertainty
    GN_inputs = create_EMIT_forMod_cov_ver4_log_reProf_cloudTop(GN_inputs);


    %% Use the tblut retrieval as the initial guess for the hyperspectral retrieval

    % --------------------------------------------------------------
    % ---------------- Retrieve Vertical Profile! ------------------
    % --------------------------------------------------------------

    disp([newline, 'Running Multispectral retrieval... ', newline])

    [GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_EMIT_ver4_log_forMo_uncert_perPixel(GN_inputs, emit, spec_response,...
        folder_paths, print_status_updates, pp);


    % --------------------------------------------------------------
    % --------------------------------------------------------------


    %%
    % ----------------------------------------------
    % ------------ SAVE OUTPUT STRUCTURE -----------
    % ----------------------------------------------

    % Save the version without an measurement uncertainty. Then we can add
    % uncertainty and save the new file



    if exist(folder_paths.saveOutput_filename, 'file')==2
        % append
        save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "overlap_pixels", "folder_paths", '-append');

    else
        save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "overlap_pixels", "folder_paths");

    end


    %% Make plot of the retrieved profile

    % plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)
    % plot_EMIT_retrieved_vertProf_with_MODIS_AIRS_AMSR(GN_outputs, GN_inputs, modis, [], amsr)



end


