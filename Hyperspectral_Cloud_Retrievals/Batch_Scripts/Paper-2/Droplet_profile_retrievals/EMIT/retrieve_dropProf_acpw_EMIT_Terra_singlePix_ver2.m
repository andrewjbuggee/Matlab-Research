
% ** Retrieving 4 variables: log(r_top), log(r_bot), log(tau_c), log(cwv) **

% *** CURRENT FORWARD MODEL UNCERTAINTIES CONSIDERED ***
%   (1) Adiabatic droplet profile assumption
%   (2) Cloud top height assumption
%   (3) Droplet distribution effective variance assumption




function [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = retrieve_dropProf_acpw_EMIT_Terra_singlePix_ver2(emit,...
            modis, era5, overlap_pixels, folder_paths, print_libRadtran_err, print_status_updates, pixel_num)



    %% Create an input structure that helps write the INP files

    % this is a built-in function that is defined at the bottom of this script
    GN_inputs = create_gauss_newton_inputs_for_emit_ver4_log(emit, print_libRadtran_err);
    %inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

    % *** Check Inputs ***

    %% Use the new custom mie tables created on 1 Feb 2026
    GN_inputs.RT.use_custom_mie_calcs = true;


    %% Adjust inputs that vary per pixel

    % -----------------------------
    % define the solar zenith angle
    % -----------------------------
    GN_inputs.RT.sza = emit.obs.solar.zenith(pixel_num);           % degree



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
    GN_inputs.RT.phi0 = mod(emit.obs.solar.azimuth(pixel_num) + 180, 360);         % degree





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
    GN_inputs.RT.vza = double(emit.obs.sensor.zenith(pixel_num)); % values are in degrees;                        % degree




    % --------------------------------
    % define the viewing azimuth angle
    % --------------------------------
    % LibRadTran defines the viewing azimuth clockwise from North.
    % So at 0deg the sensor is in the north looking south, at 90deg
    % the sensor is in the east looking West, and so on.
    % EMIT defines the sensor azimuth clockwise from due North.
    % So we don't need to change the emit values

    GN_inputs.RT.vaz = emit.obs.sensor.azimuth(pixel_num);



    %% Set the wavelenghts!!

    % --- Indexes using same 35 as above, in addition to 29 water vapor bands ---
    % This set has a total of 64 bands. They are not exactly the same set as
    % the 66 HySICS bands used to retrieve column water vapor because the
    % HySICS channels are more narrow.
    % GN_inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 71, 74, 78, 86, 87, 88, 89, 90,...
    %     94, 97, 99, 101, 105, 115, 116, 117, 139, 141, 142, 145, 147, 148, 149, 151, 156,...
    %     157, 158, 172, 175, 176, 187, 189, 190, 210, 212, 213, 214, 215, 216, 217, 218, 219,...
    %     220, 222, 231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';

    % --- Use all 285 spectral channels ---
    % GN_inputs.bands2run = (1:285)';

    % --- Mie tables don't go out far enough!! Can only use first 283 spectral channles ---
    % (2/7/2026) - custom mie tables vary from 300 to 2500 nm, but the last
    % two spectral channels of EMIT have non-zero spectral response
    % functions beyond 2500 nm.
    % GN_inputs.bands2run = (1:283)';

    % *** There are calibration issues below 500 nm and order sorting
    % filter seams to ignore ***
    % index 17 has a center wavelength of 499
    % skip the order sorting filter region between 1245 and 1320
    % skip last two bands due to range of the new custom mie tables (2/9/2026)
    GN_inputs.bands2run = [18:117, 127:283]';

    %% Override input settings with MODIS derived values


    % ** use the MODIS measurement closest to EMIT **
    unique_modis_pix = unique(overlap_pixels.modis.linear_idx);
    unique_modis_pix_idx = zeros(1, length(overlap_pixels.modis.linear_idx));
    for xx = 1:length(unique_modis_pix_idx)

        unique_modis_pix_idx(xx) = find(unique_modis_pix==overlap_pixels.modis.linear_idx(xx));

    end


    % --------------------------------------------
    % *** Use MODIS Cloud Top Height Retrieval ***
    % --------------------------------------------
    % Override the cloud depth
    GN_inputs.RT.cloud_depth = 0.3;           % km

    % override the cloud top height
    % ** MODIS cloud top height listed in meters is the geopotential height **
    GN_inputs.RT.z_topBottom = [modis.cloud.topHeight(unique_modis_pix_idx(pixel_num))/1e3,...
        (modis.cloud.topHeight(unique_modis_pix_idx(pixel_num))/1e3 - GN_inputs.RT.cloud_depth)];    % km

    % Update the height vector based on the MODIS cloud top height
    GN_inputs.RT.z_edges = linspace(GN_inputs.RT.z_topBottom(2),...
        GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers+1);   % km - the edges of each layer

    GN_inputs.RT.z = linspace(GN_inputs.RT.z_topBottom(2),...
        GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers);        % km - altitude above ground vector


    % ------------------------------------------------
    % *** Use MODIS above cloud column water vapor ***
    % ------------------------------------------------
    % ** MODIS NIR above cloud column water vapor listed in cm **
    GN_inputs.RT.waterVapor_column = modis.vapor.col_nir(unique_modis_pix_idx(pixel_num)) * 10;    % mm



    %% Override input settings with ERA5 derived values

    % ----------------------------------------------------
    % *** Use ERA5 temp/press and water vapor profiles ***
    % ----------------------------------------------------

    % first, write a radiosonde.dat file with ERA5 temperature, pressure and
    % relative humidity
    GN_inputs.RT.use_radiosonde_file = true;
    GN_inputs.RT.radiosonde_num_vars = 3;

    [GN_inputs.RT.radiosonde_file_T_P_WV, era5] = write_ERA5_radiosonde_DAT_with_multiPixels(era5,...
        folder_paths, pixel_num, [], GN_inputs.RT.radiosonde_num_vars, overlap_pixels,...
        GN_inputs.RT.atm_file, print_status_updates);

    [GN_inputs.RT.radiosonde_file_T_P, ~] = write_ERA5_radiosonde_DAT_with_multiPixels(era5,...
        folder_paths, pixel_num, [], GN_inputs.RT.radiosonde_num_vars-1, overlap_pixels,...
        GN_inputs.RT.atm_file, print_status_updates);





    %% Override input settings to define the vertical profile of effective varaince

    % ** Computed mean vertical effective variance profile from VOCALS-REx data **

    % load the set of VOCALS-REx in-situ observations
    if strcmp(folder_paths.which_computer, 'anbu8374')==true

        % --- load the effective variance observations ---
        effVar_obs = load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
            'Presentations_and_Papers/paper_2/',...
            'VR_effective_variance_at_normalized_altitudes_20-levels_19-Jan-2026.mat']);

        % --- define the directory for where the custom pre-computed mie tables are ---
        custom_mie_tables_dir = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/custom_mieTables/';


    elseif strcmp(folder_paths.which_computer, 'andrewbuggee')==true

        % --- load the effective variance observations ---
        effVar_obs = load(['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
            'Presentations_and_Papers/paper_2/',...
            'VR_effective_variance_at_normalized_altitudes_20-levels_19-Jan-2026.mat']);

        % --- define the directory for where the custom pre-computed mie tables are ---
        custom_mie_tables_dir = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/wc/custom_mieTables/';
        

    elseif strcmp(folder_paths.which_computer, 'curc')==true

        % --- load the effective variance observations ---
        effVar_obs = load(['/projects/anbu8374/Matlab-Research/Presentations_and_Papers/',...
            'paper_2/VR_effective_variance_at_normalized_altitudes_20-levels_19-Jan-2026.mat']);

        % --- define the directory for where the custom pre-computed mie tables are ---
        custom_mie_tables_dir = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/custom_mieTables/';

    end


    % set the directory for the custom mie tables
     GN_inputs.RT.mie_table_directory = custom_mie_tables_dir;

    % first, read all the filenames in the custom mie tables dir
    mie_table_filenames = dir(fullfile(custom_mie_tables_dir, '*.cdf'));
    
    % step through each filename and extract the alpha value, which is in
    % the filename
    % Extract alpha values from the filenames and store them
    mie_table_alpha_values = zeros( length(mie_table_filenames), 1);

    for ff = 1:length(mie_table_filenames)
        
        alpha_str = extractBetween(mie_table_filenames(ff).name, 'alpha_', '.cdf');

        mie_table_alpha_values(ff) = str2double(alpha_str{1}); % Assuming alpha is in the filename

    end

    % update the effective variance assumption for each cloud layer using
    % the fit using log-normal fit of alpha values
    % the distribution variance fits start at cloud base and move towards
    % cloud top
    GN_inputs.RT.distribution_var_profile = zeros(GN_inputs.RT.n_layers, 1);
    GN_inputs.RT.distribution_var_profile_std = zeros(GN_inputs.RT.n_layers, 1);
    GN_inputs.RT.distribution_var_profile_closest_2file = zeros(GN_inputs.RT.n_layers, 1);
    GN_inputs.RT.mie_table_filename = cell(GN_inputs.RT.n_layers, 1);

    for ll = 1:GN_inputs.RT.n_layers

        % convert the log normal mu and std parameter to the arithmetic
        % mean
        GN_inputs.RT.distribution_var_profile(ll) = exp( effVar_obs.alpha_fit_lognormal(ll).mu +...
            (effVar_obs.alpha_fit_lognormal(ll).sigma^2)/2 );

        % convert the log normal mu and std parameter to the arithmetic
        % standard deviation
        GN_inputs.RT.distribution_var_profile_std(ll) = sqrt( (exp(effVar_obs.alpha_fit_lognormal(ll).sigma^2) -1)...
            * exp(2*effVar_obs.alpha_fit_lognormal(ll).mu + effVar_obs.alpha_fit_lognormal(ll).sigma^2));

        % for each value, find the pre-computed mie table with the closest
        % alpha value. Save this filename
        % find the closest alpha value in the pre-computed mie tables
        [~, idx_min] = min( abs( GN_inputs.RT.distribution_var_profile(ll) - mie_table_alpha_values ) );

        % store the closest alpha value for use in the radiative transfer
        GN_inputs.RT.distribution_var_profile_closest_2file(ll) = mie_table_alpha_values(idx_min);

        GN_inputs.RT.mie_table_filename{ll} = mie_table_filenames(idx_min).name;

        
    end

    % flip the effective variance profile so the first value is cloud top
    % and the last value is cloud bottom
    GN_inputs.RT.distribution_var_profile = flipud( GN_inputs.RT.distribution_var_profile);
    GN_inputs.RT.distribution_var_profile_std = flipud( GN_inputs.RT.distribution_var_profile_std);

    % Lastly, take a mean of the vertical profile of effective variance
    % This is the value that will be used in all calculations, since only a
    % single mie table can be used in the libRadtran uvSpec calculations
    GN_inputs.RT.mean_distribution_var = mean(GN_inputs.RT.distribution_var_profile);
    GN_inputs.RT.mean_distribution_var_std = mean(GN_inputs.RT.distribution_var_profile_std);   % this is the average standard deviation of the alpha value for each cloud layer.

    % find the pre-computed mie table with the closest
    % alpha value to the mean value and save it
    % find the closest alpha value in the pre-computed mie tables
    [~, idx_min] = min( abs( GN_inputs.RT.mean_distribution_var - mie_table_alpha_values ) );

    % store the closest alpha value for use in the radiative transfer
    GN_inputs.RT.mean_distribution_var_closest_2file = mie_table_alpha_values(idx_min);
    % store the filename of the closest alpha value
    GN_inputs.RT.mean_distribution_var_closest_filename = [mie_table_filenames(idx_min).folder, '/', mie_table_filenames(idx_min).name];


    % *** Set the distribution variance that is used to compute mie
    % calculations at the distribution profile ***
    GN_inputs.RT.distribution_var = GN_inputs.RT.distribution_var_profile;




    %% This retrieval does retrieve column water vapor.

    GN_inputs.RT.modify_total_columnWaterVapor = false;             % don't modify the full column

    % *** Retreive Column Water Vapor! ***
    GN_inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % modify the column above the cloud

    %% Set output filename

    rev = 1;

    
    folder_paths.saveOutput_directory = [folder_paths.coincident_dataPath,...
        'Droplet_profile_retrievals/take_1/'];


    folder_paths.saveOutput_filename = [folder_paths.saveOutput_directory,...
        num2str(numel(GN_inputs.bands2run)),...
        'bands_EMIT_dropRetrieval_', folder_paths.coincident_dataFolder(1:end-3),...
        '_pixel_', num2str(pixel_num),...
        '_ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];


    
    % If the directory doesn't exist, create it!
    if ~exist(folder_paths.saveOutput_directory, 'dir')
    
        mkdir(folder_paths.saveOutput_directory)
    
    end




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
    
    
    emit = convert_EMIT_radiance_2_reflectance(emit, GN_inputs);
    
    
    %% Compute the radiance measurement uncertainty
    
    [emit.radiance.uncertainty, emit.radiance.uncertainty_percent_perChannel] = compute_EMIT_radiance_uncertainty(emit);
    
    
    %% Compute the reflectance uncertainty
    
    emit.reflectance.uncertainty = compute_EMIT_reflectance_uncertainty(emit, GN_inputs);



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

    use_MODIS_ERA5_data = true;

    tblut_retrieval = TBLUT_forEMIT_with_MODIS_retrievals_perPixel(emit, spec_response, folder_paths,...
        print_libRadtran_err, print_status_updates,GN_inputs, use_MODIS_ERA5_data, pixel_num);

    if print_status_updates==true
        disp([newline, 'TBLUT retrieval completed in ', num2str(toc), ' seconds', newline])
    end




    %% Compute the ACPW retrieval estimate

    if print_status_updates==true
        disp([newline, 'Computing the ACPW retrieval...', newline])
        tic
    end

    use_MODIS_ERA5_data = true;

    acpw_retrieval = ACPW_retrieval_for_EMIT_perPixel(emit, spec_response, tblut_retrieval, folder_paths, use_MODIS_ERA5_data,...
        GN_inputs, print_libRadtran_err, print_status_updates, pixel_num, era5.datProfiles);

    if print_status_updates==true
        disp([newline, 'ACPW retrieval completed in ', num2str(toc), ' seconds', newline])
    end

    %% Create the Model and Measurement prior

    use_TBLUT_estimate = true;

    GN_inputs = create_model_prior_covariance_EMIT_top_bottom_ver4_log(GN_inputs, tblut_retrieval,...
        acpw_retrieval, use_TBLUT_estimate);

    GN_inputs = create_EMIT_measurement_cov_ver4_log_no_FM_uncert_perPixel(GN_inputs, emit, pixel_num);


    %% Create the forward model covariance matrix and transform it into measurement space
    % S_b' = K_b * S_b * K_b' (Maahn et al. 2020 E1516)

    % Just re profile uncertainty
    % GN_inputs = create_EMIT_forMod_cov_ver4_log_reProf(GN_inputs);

    % Including re_profile and cloud top height uncertainty
    % GN_inputs = create_EMIT_forMod_cov_ver4_log_reProf_cloudTop(GN_inputs);

    % Including re_profile, cloud top height, and effective variance uncertainty
    GN_inputs = create_EMIT_forMod_cov_ver4_log_reProf_CTH_effVar(GN_inputs);




    %% Use the tblut retrieval as the initial guess for the hyperspectral retrieval

    % --------------------------------------------------------------
    % ---------------- Retrieve Vertical Profile! ------------------
    % --------------------------------------------------------------

    disp([newline, 'Running Multispectral retrieval... ', newline])

    [GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_EMIT_ver4_log_forMo_uncert_perPixel(GN_inputs, emit, spec_response,...
        folder_paths, print_status_updates, pixel_num, era5.datProfiles);


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



end