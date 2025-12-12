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

    folder_paths.coincident_dataFolder = '2024_05_17-T1835/';

elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;

end






%% Open Aqau data and look for overlapping pixels between EMIT and Aqua that meet certain criteria

criteria.cld_phase = 'water';
criteria.cld_cvr = 1;   % cloud fraction
criteria.cld_tau_min = 3;   % cloud optical depth
criteria.cld_tau_max = 30;   % cloud optical depth
criteria.H = 2;         % horizontal inhomogeneity index

% plot flag
plot_data = false;

% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - Then compute: overlap.time_difference_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))

[overlap_pixels, emit, modis, airs, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths, criteria, plot_data);


%% Plot all three swaths

figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.cloud.effRadius17,[],1),'.');
hold on; geoscatter(emit.radiance.geo.lat(:), emit.radiance.geo.long(:), 10, 'r.')
hold on; geoscatter(airs.geo.Latitude(:), airs.geo.Longitude(:), 10, 'c.')

%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

%% Set INP filename

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




%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
GN_inputs = create_gauss_newton_inputs_for_emit_ver4_log(emit, print_libRadtran_err);
%inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

% *** Check Inputs ***


%% Set the wavelenghts!!

% --- Indexes using same 35 as above, in addition to 29 water vapor bands ---
% This set has a total of 64 bands. They are not exactly the same set as
% the 66 HySICS bands used to retrieve column water vapor because the
% HySICS channels are more narrow.
GN_inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 71, 74, 78, 86, 87, 88, 89, 90,...
    94, 97, 99, 101, 105, 115, 116, 117, 139, 141, 142, 145, 147, 148, 149, 151, 156,...
    157, 158, 172, 175, 176, 187, 189, 190, 210, 212, 213, 214, 215, 216, 217, 218, 219,...
    220, 222, 231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';


%% Override input settings with MODIS derived values

pixel_2_use = 1;

% --------------------------------------------
% *** Use MODIS Cloud Top Height Retrieval ***
% --------------------------------------------
% Override the cloud depth
GN_inputs.RT.H = 0.3;           % km
% override the cloud top height
GN_inputs.RT.z_topBottom = [modis.cloud.topHeight(pixel_2_use)/1e3,...
    (modis.cloud.topHeight(pixel_2_use)/1e3 - GN_inputs.RT.H)];    % km

% Update the height vector based on the MODIS cloud top height
GN_inputs.RT.z_edges = linspace(GN_inputs.RT.z_topBottom(2), GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers+1);   % km - the edges of each layer
GN_inputs.RT.z = linspace(GN_inputs.RT.z_topBottom(2), GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers);        % km - altitude above ground vector


% ------------------------------------------------
% *** Use MODIS above cloud column water vapor ***
% ------------------------------------------------
GN_inputs.RT.waterVapor_column = modis.vapor.col_nir(pixel_2_use) * 10;    % mm 



% ----------------------------------------------------
% *** Use AIRS temp/press and water vapor profiles ***
% ----------------------------------------------------

% first, write a radiosonde.dat file with airs temperature, pressure and
% relative humidity
GN_inputs.RT.use_radiosonde_file = true;
GN_inputs.RT.radiosonde_file = write_AIRS_radiosonde_DAT(airs, folder_paths, pixel_2_use);




%% This retrieval does retrieve column water vapor. 

GN_inputs.RT.modify_total_columnWaterVapor = false;             % don't modify the full column

% *** Retreive Column Water Vapor! ***
GN_inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % modify the column above the cloud

%% Set output filename

rev = 1;


folder_paths.saveOutput_filename = [folder_paths.coincident_dataPath, folder_paths.coincident_dataFolder,...
    'Droplet_profile_retrievals/',...
    num2str(numel(GN_inputs.bands2run)),...
    'bands_EMIT_dropRetrieval_ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];




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

use_MODIS_AIRS_data = true;

tblut_retrieval = TBLUT_forEMIT(emit, spec_response, folder_paths, print_libRadtran_err,...
    GN_inputs, use_MODIS_AIRS_data);

if print_status_updates==true
    disp([newline, 'TBLUT retrieval completed in ', num2str(toc), ' seconds', newline])
end




%% Compute the ACPW retrieval estimate

if print_status_updates==true
    disp([newline, 'Computing the TBLUT retrieval...', newline])
    tic
end

use_MODIS_AIRS_data = true;

acpw_retrieval = ACPW_retrieval_forEMIT(emit, spec_response, folder_paths, print_libRadtran_err,...
    GN_inputs, use_MODIS_AIRS_data);

if print_status_updates==true
    disp([newline, 'TBLUT retrieval completed in ', num2str(toc), ' seconds', newline])
end

%% Create the Model and Measurement prior


GN_inputs = create_model_prior_covariance_EMIT_top_bottom_ver2(GN_inputs, tblut_retrieval, true);
%inputs = create_model_prior_covariance_EMIT_top_middle(inputs, pixels2use, tblut_retrieval, true);

GN_inputs = create_EMIT_measurement_covariance(GN_inputs, emit);


%% Use the tblut retrieval as the initial guess for the hyperspectral retrieval

tic
% --------------------------------------------------------------
% ---------------- Retrieve Vertical Profile! ------------------
% --------------------------------------------------------------

[GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_4EMIT_top_bottom_ver2(GN_inputs, emit, spec_response, folder_paths);

disp([newline, 'Multispectral retrieval took ', num2str(toc), 'seconds to run', newline])

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
    save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "pixels2use", "folder_paths", '-append');

else
    save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "pixels2use", "folder_paths");

end


%% Make plot of the retrieved profile

plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)


